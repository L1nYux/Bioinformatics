"""Check inputs, outputs, dependency versions and local document links.

Use --reproduce to rerun both studies in a disposable project-local directory
and require byte-identical outputs. The frozen result directories are read only.
"""

import argparse
import hashlib
from importlib import metadata
import json
import re
import shutil
import tempfile
from urllib.parse import unquote, urlsplit

import numpy as np
import pandas as pd

from analyze import audit, figures
from experiment import run as run_baseline, study_targets
from extension import generate_controlled_queries, paired_effects, reporting_summary, run as run_extension
from fragment import MODES, ROOT


def check_hashes():
    files = []
    for line in (ROOT / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", 1)
        path = (ROOT / name).resolve()
        if not path.is_relative_to(ROOT.resolve()) or name in files:
            raise ValueError(f"Unsafe or repeated checksum path: {name}")
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Checksum mismatch: {name}")
        files.append(name)
    if not files:
        raise ValueError("Empty checksum manifest")
    return files


def check_links(files):
    count = 0
    for name in files:
        if not name.endswith(".md"):
            continue
        path = ROOT / name
        for link in re.findall(r"\[[^\]]*\]\(([^\n]+?)\)", path.read_text(encoding="utf-8")):
            link = link.strip("<>")
            parsed = urlsplit(link)
            if parsed.scheme or not parsed.path:
                continue
            destination = (path.parent / unquote(parsed.path)).resolve()
            if not destination.is_file():
                raise ValueError(f"Broken local link in {name}: {link}")
            count += 1
    return count


def check_dependencies():
    for line in (ROOT / "requirements-lock.txt").read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        name, version = line.split("==")
        if metadata.version(name) != version:
            raise ValueError(f"Dependency version mismatch: {name}; expected {version}")


def compare_table(actual, expected, label):
    try:
        pd.testing.assert_frame_equal(actual, expected, check_dtype=False, rtol=1e-10, atol=1e-12)
    except AssertionError as exc:
        raise ValueError(f"Incorrect table: {label}") from exc


def check_extended_results(queries, results, cap):
    expected = {(q, m) for q in queries.query_id for m in MODES}
    if (len(results) != len(expected) or set(zip(results.query_id, results["mode"])) != expected):
        raise ValueError("Expected exactly one extended result per query and mode")
    joined = results.merge(queries, on="query_id", validate="many_to_one")
    for r in joined.itertuples():
        intervals = json.loads(r.intervals)
        keys = {(i["start"], i["end"]) for i in intervals}
        truth = r.true_start, r.true_end
        hit = pd.notna(r.start)
        exact = hit and (r.start, r.end) == truth
        error = (abs(r.start-truth[0])+abs(r.end-truth[1]))/2 if hit else np.nan
        status = ("incomplete" if not r.enumeration_complete else "no_hit" if not keys
                  else "unique" if len(keys) == 1 else "ambiguous")
        if (r.exact != exact or not np.isclose(error, r.error, equal_nan=True)
                or not np.isclose(r.coverage, r.paired/r.length)
                or not 0 <= r.coverage <= 1 or r.interval_count_observed != len(keys)
                or len(intervals) != len(keys) or status != r.status
                or r.accepted != (status == "unique")
                or r.truth_in_observed_intervals != (truth in keys)
                or not 0 <= r.alignments_examined <= cap
                or (hit and (r.start, r.end) not in keys)):
            raise ValueError(f"Invalid extended metrics: {r.query_id}, {r.mode}")
    return joined


def verify(reproduce=False):
    files = check_hashes()
    check_dependencies()
    links = check_links(files)
    targets = study_targets()
    folder = ROOT / "analysis/extension"
    manifest = json.loads((folder / "run_manifest.json").read_text())
    for name, digest in manifest["sha256"].items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"Stale extension run: {name} changed")
    protocol = json.loads((ROOT / "extension.json").read_text())
    queries = pd.read_csv(folder / "controlled_queries.csv", keep_default_na=False)
    compare_table(queries, pd.DataFrame(generate_controlled_queries(targets, protocol)), "controlled queries")
    results = pd.read_csv(folder / "controlled_results.csv")
    joined = check_extended_results(queries, results, protocol["max_alignments"])
    effects, per_protein = paired_effects(joined, protocol)
    tables = dict(controlled_summary=reporting_summary(joined, ["mode", "terminal_group"]),
                  paired_effects=effects, per_protein_effects=per_protein,
                  by_condition=reporting_summary(joined, ["mode", "length", "pct", "terminal_group"]),
                  by_seed=reporting_summary(joined, ["mode", "seed", "terminal_group"]))
    for name, expected in tables.items():
        compare_table(pd.read_csv(folder / f"{name}.csv"), expected, name)
    baseline_queries = pd.read_csv(ROOT / "results/queries.csv", keep_default_na=False)
    candidates = pd.read_csv(folder / "baseline_candidates.csv")
    check_extended_results(baseline_queries, candidates, protocol["max_alignments"])
    compare_table(pd.read_csv(folder / "baseline_summary.csv"),
                  reporting_summary(candidates, ["mode"]), "baseline summary")
    cache = ROOT / ".pytest_cache"
    cache.mkdir(exist_ok=True)
    if not cache.resolve().is_relative_to(ROOT.resolve()):
        raise ValueError("Verification directory must stay inside the project")
    reproduced = []
    with tempfile.TemporaryDirectory(prefix="verify-", dir=cache) as name:
        from pathlib import Path
        temporary = Path(name)
        baseline = temporary / "baseline"
        baseline.mkdir()
        if reproduce:
            run_baseline(baseline)
        else:
            for file in ("queries.csv", "results.csv", "per_protein.csv", "summary.csv"):
                shutil.copyfile(ROOT / "results" / file, baseline / file)
        audit(baseline)
        if reproduce:
            figures(baseline)
            for original in sorted((ROOT / "results").iterdir()):
                if original.is_file() and original.read_bytes() != (baseline / original.name).read_bytes():
                    raise ValueError(f"Baseline reproduction differs: {original.name}")
                reproduced.append(f"results/{original.name}")
            extension = temporary / "extension"
            run_extension(extension)
            for original in sorted(folder.iterdir()):
                if original.is_file() and original.read_bytes() != (extension / original.name).read_bytes():
                    raise ValueError(f"Extension reproduction differs: {original.name}")
                reproduced.append(f"analysis/extension/{original.name}")
        elif (baseline / "ambiguous.csv").read_bytes() != (ROOT / "results/ambiguous.csv").read_bytes():
            raise ValueError("Stale baseline ambiguity table")
    return dict(checksum_files=len(files), local_links=links, proteins=len(targets),
                baseline_queries=len(baseline_queries), controlled_queries=len(queries),
                controlled_alignments=len(results), reproduced_identically=reproduced)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reproduce", action="store_true")
    print(json.dumps(verify(parser.parse_args().reproduce), indent=2))
