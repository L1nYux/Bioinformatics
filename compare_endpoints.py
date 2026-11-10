"""Compare global and local first-optimum endpoints in the existing study CSVs."""

import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
OUTCOMES = ["both_exact", "local_only", "global_only", "neither_exact"]
TERMINI = {(False, False): "none", (True, False): "first",
           (False, True): "last", (True, True): "both"}


def paired_queries(directory):
    queries = pd.read_csv(directory / "queries.csv", keep_default_na=False)
    queries = queries.set_index("query_id")
    results = pd.read_csv(directory / "results.csv")
    for mode in ("global", "local"):
        predictions = results.loc[results["mode"] == mode].set_index("query_id")
        if not queries.index.is_unique or not predictions.index.is_unique:
            raise ValueError(f"{mode}: duplicate query ID")
        if set(predictions.index) != set(queries.index):
            raise ValueError(f"{mode}: query IDs differ between queries.csv and results.csv")
        columns = ["start", "end", "paired", "exact"]
        queries = queries.join(predictions[columns].add_prefix(f"{mode}_"))
        exact = (queries[f"{mode}_start"] == queries.true_start) & (queries[f"{mode}_end"] == queries.true_end)
        if not exact.equals(queries[f"{mode}_exact"]):
            raise ValueError(f"{mode}: exact flags disagree with endpoint coordinates")

    termini = []
    for q in queries.itertuples():
        positions = [i for i, (a, b) in enumerate(zip(q.original, q.query, strict=True)) if a != b]
        if len(q.query) != q.length or positions != list(map(int, q.substitutions.split())):
            raise ValueError(f"{q.Index}: substitutions disagree with query sequence")
        termini.append(TERMINI[0 in positions, q.length - 1 in positions])
    queries["terminal_substitution"] = termini
    outcomes = {(True, True): "both_exact", (False, True): "local_only",
                (True, False): "global_only", (False, False): "neither_exact"}
    queries["outcome"] = [outcomes[g, l] for g, l in zip(queries.global_exact, queries.local_exact)]

    ties = pd.read_csv(directory / "ambiguous.csv")
    alternative = ties.loc[(ties["mode"] == "global") & ties.exact_available, "query_id"]
    queries["global_truth_optimal"] = queries.global_exact | queries.index.isin(alternative)
    return queries.sort_index()


def counts(group, outcome):
    table = pd.crosstab(group, outcome).reindex(columns=OUTCOMES, fill_value=0)
    table["queries"] = table.sum(axis=1)
    return table


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "results")
    parser.add_argument("--details", type=Path, help="write the discordant queries to a separate CSV")
    args = parser.parse_args()
    paired = paired_queries(args.input)
    print(counts(paired.accession, paired.outcome).to_string())
    print(counts(paired.terminal_substitution, paired.outcome).to_string())
    print(counts(paired.terminal_substitution.ne("none").rename("any_terminal_substitution"),
                 paired.outcome).to_string())
    if args.details:
        columns = ["accession", "length", "pct", "substitutions", "terminal_substitution",
                   "true_start", "true_end", "global_start", "global_end", "global_paired",
                   "local_start", "local_end", "local_paired", "outcome", "global_truth_optimal"]
        discordant = paired.loc[paired.global_exact != paired.local_exact, columns]
        args.details.parent.mkdir(parents=True, exist_ok=True)
        discordant.to_csv(args.details, lineterminator="\n")
