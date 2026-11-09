"""Generate the frozen study queries, align each one and summarize by protein."""

import argparse
import csv
import hashlib
import json
import random
from pathlib import Path

import pandas as pd

from fragment import AMINO, MODES, ROOT, locate, make_aligner, read_protein


def study_targets():
    targets = {}
    with (ROOT / "data/study_manifest.csv").open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            path = ROOT / "data" / row["raw_file"]
            if hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
                raise ValueError(f"File checksum mismatch: {path}")
            record, seq = read_protein(path)
            key = row["accession"]
            if (record.id.split("|")[1] != key or len(seq) != int(row["length"])
                    or key in targets or seq in targets.values()):
                raise ValueError(f"Invalid or repeated manifest entry: {key}")
            targets[key] = seq
    return targets


def generate_queries(targets, config):
    rng = random.Random(config["seed"])
    for accession, target in targets.items():
        last_start = len(target) - max(config["lengths"])
        starts = sorted(rng.sample(range(last_start + 1), config["starts_per_protein"]))
        for start in starts:
            for length in config["lengths"]:
                original = target[start:start + length]
                for pct in config["replacement_percent"]:
                    query = list(original)
                    count = length * pct // 100
                    positions = sorted(rng.sample(range(length), count))
                    for position in positions:
                        query[position] = rng.choice(AMINO.replace(query[position], ""))
                    yield dict(query_id=f"{accession}_{start:03}_{length:02}_{pct:02}",
                               accession=accession, true_start=start, true_end=start + length,
                               length=length, pct=pct, original=original, query="".join(query),
                               substitutions=" ".join(str(p) for p in positions))


def run(output):
    config = json.loads((ROOT / "study.json").read_text())
    targets = study_targets()
    queries = pd.DataFrame(generate_queries(targets, config))
    aligners = {mode: make_aligner(mode, config) for mode in MODES}
    results = []
    for q in queries.itertuples():
        for mode, aligner in aligners.items():
            result, _ = locate(aligner, targets[q.accession], q.query)
            start, end = result["start"], result["end"]
            error = None if start is None else (abs(start-q.true_start) + abs(end-q.true_end)) / 2
            results.append(dict(query_id=q.query_id, mode=mode, **result,
                                exact=(start, end) == (q.true_start, q.true_end), error=error))
    results = pd.DataFrame(results)
    joined = results.merge(queries, on="query_id", validate="many_to_one")
    groups = ["length", "pct", "mode"]
    per_protein = joined.groupby(["accession", *groups]).agg(
        queries=("exact", "size"), hits=("start", "count"), exact=("exact", "mean"),
        coverage=("coverage", "mean"), error=("error", "mean"),
    ).reset_index()
    summary = per_protein.groupby(groups).agg(
        queries=("queries", "sum"), hits=("hits", "sum"), proteins=("accession", "size"),
        exact=("exact", "mean"), coverage=("coverage", "mean"), error=("error", "mean"),
    ).reset_index()
    output.mkdir(parents=True, exist_ok=True)
    for name, table in [("queries", queries), ("results", results), ("per_protein", per_protein),
                        ("summary", summary)]:
        table.to_csv(output / f"{name}.csv", index=False)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results")
    print(run(parser.parse_args().output).to_string(index=False))
