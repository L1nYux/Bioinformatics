"""Audit predictions against the frozen queries and render the paper figures."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from experiment import study_targets, summarize
from fragment import MODES, ROOT, make_aligner


def audit(folder):
    config = json.loads((ROOT / "study.json").read_text())
    targets = study_targets()
    queries = pd.read_csv(folder / "queries.csv", keep_default_na=False)
    results = pd.read_csv(folder / "results.csv")
    if not queries.query_id.is_unique or queries.empty:
        raise ValueError("Expected nonempty, unique query IDs")
    expected = {(qid, mode) for qid in queries.query_id for mode in MODES}
    if len(results) != len(expected) or set(zip(results.query_id, results["mode"])) != expected:
        raise ValueError("Expected one result per query and method")
    aligners = {mode: make_aligner(mode, config) for mode in MODES}
    tie_rows = []
    for row in results.merge(queries, validate="many_to_one", on="query_id").itertuples():
        target = targets[row.accession]
        changes = [i for i, (a, b) in enumerate(zip(row.original, row.query, strict=True)) if a != b]
        if (not 0 <= row.true_start < row.true_end <= len(target)
                or target[row.true_start:row.true_end] != row.original
                or len(row.query) != row.length or row.true_end-row.true_start != row.length
                or changes != list(map(int, row.substitutions.split()))
                or len(changes) != row.length * row.pct // 100):
            raise ValueError(f"Bad truth for {row.query_id}")
        alignments = aligners[row.mode].align(target, row.query)
        bounds = set()
        for i, alignment in enumerate(alignments):
            # The study is small. Prevent unbounded enumeration on new inputs.
            if i >= 1000:
                raise ValueError(f"More than 1000 optimal paths in {row.query_id}: review ties")
            indices = alignment.indices
            pairs = indices[0, (indices >= 0).all(axis=0)]
            if not len(pairs):
                continue
            interval = int(pairs[0]), int(pairs[-1]) + 1
            bounds.add(interval)
            if i == 0:
                error = (abs(interval[0]-row.true_start) + abs(interval[1]-row.true_end)) / 2
                if ((row.start, row.end) != interval or len(pairs) != row.paired
                        or not np.isclose(len(pairs)/row.length, row.coverage)
                        or (interval == (row.true_start,row.true_end)) != row.exact or error != row.error
                        or not np.isclose(row.score, alignment.score)):
                    raise ValueError(f"Metric mismatch for {row.query_id}, {row.mode}")
        if not bounds and not (pd.isna(row.start) and pd.isna(row.end) and row.paired == 0
                               and row.coverage == 0 and not row.exact and pd.isna(row.error)
                               and pd.isna(row.score)):
            raise ValueError(f"Incorrect missing hit for {row.query_id}")
        if len(bounds) > 1:
            tie_rows.append(dict(query_id=row.query_id, mode=row.mode, intervals=sorted(bounds),
                                 exact_first=row.exact, exact_available=(row.true_start,row.true_end) in bounds))
    for name, expected_table in zip(("per_protein", "summary"), summarize(queries, results)):
        try:
            pd.testing.assert_frame_equal(pd.read_csv(folder / f"{name}.csv"), expected_table,
                                          check_dtype=False, rtol=1e-10, atol=1e-12)
        except AssertionError as exc:
            raise ValueError(f"Stale or incorrect {name}.csv") from exc
    pd.DataFrame(tie_rows, columns=["query_id","mode","intervals","exact_first","exact_available"]
                 ).to_csv(folder / "ambiguous.csv", index=False)


def figures(folder):
    summary = pd.read_csv(folder / "summary.csv")
    colors = dict(zip(MODES, ["#96451e", "#146782", "#176b49"]))
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.5), sharey=True, layout="constrained")
    for length, ax in zip(sorted(summary.length.unique()), axes, strict=True):
        for mode in MODES:
            group = summary.loc[(summary.length == length) & (summary["mode"] == mode)].sort_values("pct")
            ax.plot(group.pct, group.exact*100, marker="o", color=colors[mode], label=mode)
        ax.set(title=f"{length} residues", xlabel="Replaced residues (%)", xticks=[0,10,30], ylim=(0,105))
    axes[0].set_ylabel("Exact endpoints (%)")
    axes[-1].legend(loc="lower left")
    fig.savefig(folder / "exact.png", dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(9,3.6), layout="constrained")
    for ax, metric, label, scale in zip(axes, ["coverage","error"],
                                         ["Query coverage (%)","Mean endpoint error (residues)"], [100,1]):
        for mode in MODES:
            group = summary.loc[summary["mode"] == mode].groupby("pct")[metric].mean()
            ax.plot(group.index, group*scale, marker="o", color=colors[mode], label=mode)
        ax.set(xlabel="Replaced residues (%)", ylabel=label, xticks=[0,10,30])
    axes[1].legend()
    fig.savefig(folder / "coverage_error.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results")
    folder = parser.parse_args().output
    audit(folder)
    figures(folder)
