"""Matched terminal-substitution experiment and truth-free ambiguity reporting.

Run separately from the frozen 270-query baseline. The protocol is locked in
extension.json before evaluating its outcomes; scoring stays in study.json.
"""

import argparse
import hashlib
import json
import platform
import random
from pathlib import Path

import Bio
import numpy as np
import pandas as pd

from experiment import study_targets
from fragment import AMINO, MODES, ROOT, locate, locate_candidates, make_aligner


GROUPS = ("interior", "first", "last", "both")


def keyed_rng(*parts):
    digest = hashlib.sha256(json.dumps(parts).encode()).digest()
    return random.Random(int.from_bytes(digest, "big"))


def validate_protocol(protocol):
    for field in ("seeds", "lengths", "replacement_percent"):
        values = protocol[field]
        if (not values or len(values) != len(set(values))
                or any(type(value) is not int for value in values)):
            raise ValueError(f"{field} must contain distinct integers")
    for field in ("starts_per_protein_per_length", "max_alignments", "bootstrap_resamples"):
        if type(protocol[field]) is not int or protocol[field] < 1:
            raise ValueError(f"{field} must be positive")
    if protocol["terminal_groups"] != list(GROUPS):
        raise ValueError("Expected all four terminal groups in the fixed order")
    for length in protocol["lengths"]:
        for pct in protocol["replacement_percent"]:
            if length < 4 or not 0 < pct <= 100 or not 2 <= length * pct // 100 <= length - 2:
                raise ValueError("Matched design needs 2 <= substitutions <= length-2")


def generate_controlled_queries(targets, protocol):
    """Four queries share a source interval and all unrelocated substitutions.

    k-2 background interior positions stay fixed. Two further interior changes
    may move to the first/last position. Mutation identities are shared wherever
    the position is shared. Starts use the entire valid range for EACH length.
    """
    validate_protocol(protocol)
    for seed in protocol["seeds"]:
        for accession, target in sorted(targets.items()):
            for length in protocol["lengths"]:
                possible = len(target) - length + 1
                count = protocol["starts_per_protein_per_length"]
                if possible < count:
                    raise ValueError(f"{accession}: too short for {count} starts at length {length}")
                rng = keyed_rng(seed, accession, length, "starts")
                for start in sorted(rng.sample(range(possible), count)):
                    original = target[start:start + length]
                    rng = keyed_rng(seed, accession, length, start, "mutations")
                    interior = rng.sample(range(1, length - 1), length - 2)
                    replacements = [rng.choice(AMINO.replace(a, "")) for a in original]
                    for pct in protocol["replacement_percent"]:
                        k = length * pct // 100
                        background = interior[:k - 2]
                        left, right = interior[k - 2:k]
                        positions = {
                            "interior": background + [left, right],
                            "first": background + [0, right],
                            "last": background + [left, length - 1],
                            "both": background + [0, length - 1],
                        }
                        pair_id = f"s{seed}_{accession}_{start:03}_{length:02}_{pct:02}"
                        for group in GROUPS:
                            query = list(original)
                            for pos in positions[group]:
                                query[pos] = replacements[pos]
                            yield dict(query_id=f"{pair_id}_{group}", pair_id=pair_id,
                                       seed=seed, accession=accession, true_start=start,
                                       true_end=start+length, length=length, pct=pct,
                                       terminal_group=group, original=original, query="".join(query),
                                       substitutions=" ".join(map(str, sorted(positions[group]))))


def evaluate(queries, targets, scoring, max_alignments):
    aligners = {mode: make_aligner(mode, scoring) for mode in MODES}
    rows = []
    for q in queries.itertuples():
        truth = q.true_start, q.true_end
        for mode, aligner in aligners.items():
            target = targets[q.accession]
            report = locate_candidates(aligner, target, q.query, max_alignments)
            first = report.pop("first")
            # Independently check the first solution's coordinates using columns.
            _, alignment = locate(aligner, target, q.query)
            if alignment is not None:
                indices = alignment.indices
                paired = indices[:, (indices >= 0).all(axis=0)]
                interval = ((int(paired[0, 0]), int(paired[0, -1])+1)
                            if paired.shape[1] else (None, None))
                if ((first["start"], first["end"]) != interval
                        or first["paired"] != paired.shape[1]
                        or first["coverage"] != paired.shape[1]/q.length):
                    raise ValueError(f"Column audit failed: {q.query_id}, {mode}")
            predicted = first["start"], first["end"]
            intervals = report.pop("intervals")
            rows.append(dict(
                query_id=q.query_id, mode=mode, **first, exact=predicted == truth,
                error=None if first["start"] is None else
                (abs(first["start"]-truth[0]) + abs(first["end"]-truth[1]))/2,
                **report, accepted=report["status"] == "unique",
                truth_in_observed_intervals=any((i["start"], i["end"]) == truth for i in intervals),
                intervals=json.dumps(intervals, separators=(",", ":")),
            ))
    return pd.DataFrame(rows)


def reporting_summary(results, group_columns):
    rows = []
    for keys, group in results.groupby(group_columns, sort=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        accepted = group.loc[group.accepted]
        rows.append(dict(zip(group_columns, keys)) | dict(
            queries=len(group), first_correct=int(group.exact.sum()), first_exact=group.exact.mean(),
            accepted=len(accepted), acceptance_rate=len(accepted)/len(group),
            accepted_wrong=int((~accepted.exact).sum()),
            accepted_error_rate=(~accepted.exact).mean() if len(accepted) else None,
            ambiguous=int(group.status.eq("ambiguous").sum()),
            incomplete=int(group.status.eq("incomplete").sum()),
            no_hit=int(group.status.eq("no_hit").sum()),
            truth_in_observed_intervals=int(group.truth_in_observed_intervals.sum()),
        ))
    return pd.DataFrame(rows)


def paired_effects(joined, protocol):
    keys = ["pair_id", "seed", "accession", "mode", "length", "pct"]
    wide = joined.pivot(index=keys, columns="terminal_group", values="exact").reset_index()
    if wide[list(GROUPS)].isna().any().any():
        raise ValueError("Incomplete matched quartet")
    wide[list(GROUPS)] = wide[list(GROUPS)].astype(int)
    proteins = sorted(wide.accession.unique())
    # Resample proteins with all of their queries, lengths and seeds together.
    indices = np.random.default_rng(protocol["bootstrap_seed"]).integers(
        0, len(proteins), size=(protocol["bootstrap_resamples"], len(proteins)))
    rows, per_protein = [], []
    for mode, group in wide.groupby("mode"):
        for terminal in GROUPS[1:]:
            delta = group[terminal] - group.interior
            changed = group.assign(delta=delta)
            cluster = changed.groupby("accession").delta.mean().reindex(proteins)
            low, high = np.quantile(cluster.to_numpy()[indices].mean(axis=1), [0.025, 0.975])
            seeds = changed.groupby("seed").delta.mean()
            rows.append(dict(mode=mode, terminal_group=terminal, pairs=len(group),
                             gains=int((delta == 1).sum()), losses=int((delta == -1).sum()),
                             delta=cluster.mean(), ci_low=low, ci_high=high,
                             seed_min=seeds.min(), seed_max=seeds.max(),
                             protein_min=cluster.min(), protein_max=cluster.max()))
            for accession, value in cluster.items():
                per_protein.append(dict(mode=mode, terminal_group=terminal,
                                        accession=accession, delta=value))
    return pd.DataFrame(rows), pd.DataFrame(per_protein)


def figures(summary, effects, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"global": "#96451e", "local": "#146782", "semiglobal": "#176b49"}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), layout="constrained")
    for offset, mode in enumerate(MODES):
        group = summary.loc[summary["mode"].eq(mode)].set_index("terminal_group").reindex(GROUPS)
        axes[0].plot(range(4), group.first_exact*100, "o-", label=mode, color=colors[mode])
        effect = effects.loc[effects["mode"].eq(mode)].set_index("terminal_group").reindex(GROUPS[1:])
        x = np.arange(3) + (offset-1)*0.13
        axes[1].errorbar(x, effect.delta*100,
                         yerr=np.array([effect.delta-effect.ci_low, effect.ci_high-effect.delta])*100,
                         fmt="o", capsize=4, color=colors[mode], label=mode)
    axes[0].set(xticks=range(4), xticklabels=["Interior", "First", "Last", "Both"],
                xlabel="Substitution placement", ylabel="First-optimum exact endpoints (%)", ylim=(0, 105))
    axes[1].axhline(0, color="#888888", linewidth=0.7)
    axes[1].set(xticks=range(3), xticklabels=["First", "Last", "Both"],
                xlabel="Placement compared with Interior", ylabel="Paired difference (percentage points)")
    axes[0].legend()
    axes[1].set_title("95% protein-cluster bootstrap intervals", fontsize=10)
    fig.savefig(output / "controlled_endpoints.png", dpi=180)
    plt.close(fig)


def markdown_table(frame, columns):
    lines = ["| " + " | ".join(label for _, label in columns) + " |",
             "| " + " | ".join("---" for _ in columns) + " |"]
    for _, row in frame.iterrows():
        lines.append("| " + " | ".join(str(row[key]) for key, _ in columns) + " |")
    return "\n".join(lines)


def write_report(output, summary, effects, baseline, protocol, protein_count):
    table = summary.copy()
    table["exact_display"] = table.first_exact.map(lambda x: f"{x*100:.2f}%")
    table["wrong_display"] = table.accepted_error_rate.map(lambda x: f"{x*100:.2f}%")
    base = baseline.copy()
    base["acceptance"] = base.acceptance_rate.map(lambda x: f"{x*100:.2f}%")
    base["risk"] = base.accepted_error_rate.map(lambda x: f"{x*100:.2f}%")
    effect = effects.copy()
    effect["display"] = effect.apply(lambda r: f"{r.delta*100:.2f} [{r.ci_low*100:.2f}, {r.ci_high*100:.2f}]", axis=1)
    effect["seed_range"] = effect.apply(lambda r: f"{r.seed_min*100:.2f} 至 {r.seed_max*100:.2f}", axis=1)
    text = "# 扩展实验：端部扰动与定位歧义\n\n"
    text += "本文件由 `extension.py` 从本地真实运行结果生成。原 270 查询基线保持独立。\n\n"
    query_count = int(summary.queries.sum()) // len(MODES)
    text += (f"协议：{len(protocol['seeds'])} 个种子、{protein_count} 条冻结蛋白、每种长度各 "
             f"{protocol['starts_per_protein_per_length']} 个起点、{len(protocol['lengths'])} 种长度、"
             f"{len(protocol['replacement_percent'])} 档替换比例、4 种位置组，共 {query_count} 个查询、"
             f"{query_count*len(MODES)} 次比对。每组 {query_count//4} 个查询。起点按各自长度的完整合法范围抽取。\n\n")
    text += "四组共享母片段、总替换数、k-2 个内部替换及共享位置的替换身份。将另外两个内部替换分别移至首端、末端或两端，比较一次定位成败的配对变化。变化同时包含替换位置和端部残基身份，结论限于此人工干预。\n\n"
    text += "## 原基线上的歧义报告\n\n"
    text += "接受规则仅使用输入序列：枚举结束且观察到唯一端点区间才输出唯一位置。这里的唯一只指固定评分下的最优端点，不保证生物学正确。\n\n"
    text += markdown_table(base, [("mode", "方法"), ("queries", "查询"), ("accepted", "接受"),
                                   ("acceptance", "接受率"), ("accepted_wrong", "接受中错误"),
                                   ("risk", "接受后错误率"), ("ambiguous", "歧义"), ("incomplete", "截断")])
    text += "\n\n## 匹配扰动结果\n\n"
    text += markdown_table(table, [("mode", "方法"), ("terminal_group", "位置组"),
                                    ("queries", "查询"), ("first_correct", "首解正确"),
                                    ("exact_display", "首解准确率"), ("accepted", "接受"),
                                    ("accepted_wrong", "接受中错误"), ("incomplete", "截断")])
    text += "\n\n## 相对内部替换组的配对变化\n\n"
    text += (f"以蛋白为簇等权重采样 {protocol['bootstrap_resamples']} 次；保留该蛋白的全部种子、长度与查询。"
             "区间为百分位 bootstrap 描述区间，依赖便利样本和可交换簇假设，不代表总体保证，也未将相关查询当作独立生物样本。单位为百分点。\n\n")
    text += markdown_table(effect, [("mode", "方法"), ("terminal_group", "位置组"), ("pairs", "配对数"),
                                     ("gains", "由错变对"), ("losses", "由对变错"),
                                     ("display", "变化 [95%区间]"), ("seed_range", "逐种子变化范围")])
    text += "\n\n![端部位置干预与配对差值](controlled_endpoints.png)\n\n"
    text += "推断前固定配置；不按结果挑选种子或评分参数。多个种子是在同一批蛋白上的计算重复，不增加独立蛋白数。未测跨物种、插入缺失、未知母蛋白检索或评分敏感性。\n"
    (output / "report.md").write_text(text, encoding="utf-8", newline="\n")


def run(output):
    protocol = json.loads((ROOT / "extension.json").read_text())
    scoring = json.loads((ROOT / "study.json").read_text())
    targets = study_targets()
    queries = pd.DataFrame(generate_controlled_queries(targets, protocol))
    output.mkdir(parents=True, exist_ok=True)
    results = evaluate(queries, targets, scoring, protocol["max_alignments"])
    joined = results.merge(queries, on="query_id", validate="many_to_one")
    summary = reporting_summary(joined, ["mode", "terminal_group"])
    effects, per_protein_effects = paired_effects(joined, protocol)
    baseline_queries = pd.read_csv(ROOT / "results/queries.csv", keep_default_na=False)
    baseline_results = evaluate(baseline_queries, targets, scoring, protocol["max_alignments"])
    frozen = pd.read_csv(ROOT / "results/results.csv")
    columns = list(frozen.columns)
    pd.testing.assert_frame_equal(frozen, baseline_results[columns], check_dtype=False,
                                  rtol=1e-10, atol=1e-12)
    baseline_summary = reporting_summary(baseline_results, ["mode"])
    tables = dict(controlled_queries=queries, controlled_results=results, controlled_summary=summary,
                  paired_effects=effects, per_protein_effects=per_protein_effects,
                  by_condition=reporting_summary(joined, ["mode", "length", "pct", "terminal_group"]),
                  by_seed=reporting_summary(joined, ["mode", "seed", "terminal_group"]),
                  baseline_candidates=baseline_results, baseline_summary=baseline_summary)
    for name, table in tables.items():
        table.to_csv(output / f"{name}.csv", index=False, lineterminator="\n")
    files = ["extension.json", "study.json", "data/study_manifest.csv", "fragment.py", "extension.py", "experiment.py"]
    metadata = dict(protocol=protocol, scoring=scoring, proteins=len(targets), queries=len(queries),
                    alignments=len(results), python=platform.python_version(), biopython=Bio.__version__,
                    numpy=np.__version__, pandas=pd.__version__, platform=platform.system(),
                    sha256={f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest() for f in files})
    (output / "run_manifest.json").write_text(json.dumps(metadata, indent=2)+"\n", encoding="utf-8", newline="\n")
    figures(summary, effects, output)
    write_report(output, summary, effects, baseline_summary, protocol, len(targets))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "analysis/extension")
    print(run(parser.parse_args().output).to_string(index=False))
