import json
import shutil
import subprocess
import sys

import pandas as pd
import pytest
from Bio import Align

import analyze
from extension import GROUPS, generate_controlled_queries, paired_effects, reporting_summary
from fragment import ROOT, locate_candidates, make_aligner


CONFIG = json.loads((ROOT / "study.json").read_text())
PROTOCOL = json.loads((ROOT / "extension.json").read_text())


def test_report_distinguishes_complete_unique_ambiguous_and_capped():
    aligner = make_aligner("semiglobal", CONFIG)
    unique = locate_candidates(aligner, "ACDEFGHIKLMN", "EFGH", 1)
    assert unique["status"] == "unique" and unique["enumeration_complete"]
    repeated = locate_candidates(aligner, "EFGHAEFGH", "EFGH")
    assert repeated["status"] == "ambiguous"
    assert {(i["start"], i["end"]) for i in repeated["intervals"]} == {(0, 4), (5, 9)}
    capped = locate_candidates(aligner, "EFGHAEFGH", "EFGH", 1)
    assert capped["status"] == "incomplete" and not capped["enumeration_complete"]
    assert capped["interval_count_observed"] == capped["alignments_examined"] == 1


def test_paths_with_same_endpoints_are_one_interval():
    aligner = Align.PairwiseAligner(mode="global", match_score=2, mismatch_score=-10,
                                   open_gap_score=-1, extend_gap_score=-1)
    report = locate_candidates(aligner, "AAACAAA", "AAAGAAA")
    assert report["status"] == "unique"
    assert report["alignments_examined"] > 1
    assert report["interval_count_observed"] == 1


def test_no_hit_and_invalid_inputs():
    aligner = make_aligner("local", CONFIG)
    report = locate_candidates(aligner, "DDDD", "WWWW")
    assert report["status"] == "no_hit" and report["first"]["coverage"] == 0
    for cap in (0, -1, True, 1.5):
        with pytest.raises(ValueError, match="positive integer"):
            locate_candidates(aligner, "ACD", "ACD", cap)
    for query in ("", "ABZ", "acd"):
        with pytest.raises(ValueError, match="amino-acid"):
            locate_candidates(aligner, "ACD", query)


def test_controlled_quartets_have_equal_counts_and_shared_mutations():
    protocol = PROTOCOL | {"seeds": [11, 12], "lengths": [20, 40],
                           "starts_per_protein_per_length": 2}
    targets = {"one": "ACDEFGHIKLMNPQRSTVWY" * 4, "two": "WVTSRQPNMLKIHGFEDCAW" * 4}
    rows = list(generate_controlled_queries(targets, protocol))
    assert len(rows) == 2*2*2*2*2*4
    assert len({r["query_id"] for r in rows}) == len(rows)
    assert rows == list(generate_controlled_queries(dict(reversed(list(targets.items()))), protocol))
    for begin in range(0, len(rows), 4):
        group = rows[begin:begin+4]
        assert [r["terminal_group"] for r in group] == list(GROUPS)
        assert len({r["pair_id"] for r in group}) == 1
        positions = []
        for row in group:
            changed = {i for i, (a, b) in enumerate(zip(row["original"], row["query"])) if a != b}
            assert changed == set(map(int, row["substitutions"].split()))
            assert len(changed) == row["length"] * row["pct"] // 100
            positions.append(changed)
        end = group[0]["length"] - 1
        assert [p & {0, end} for p in positions] == [set(), {0}, {end}, {0, end}]
        assert len(set.intersection(*positions)) == len(positions[0])-2
        for a in group:
            for b in group:
                shared = set(map(int, a["substitutions"].split())) & set(map(int, b["substitutions"].split()))
                assert all(a["query"][i] == b["query"][i] for i in shared)


def test_full_length_source_and_impossible_protocol():
    protocol = PROTOCOL | {"seeds": [1], "lengths": [20], "starts_per_protein_per_length": 1}
    rows = list(generate_controlled_queries({"whole": "ACDEFGHIKLMNPQRSTVWY"}, protocol))
    assert {r["true_start"] for r in rows} == {0}
    assert {r["true_end"] for r in rows} == {20}
    with pytest.raises(ValueError, match="too short"):
        list(generate_controlled_queries({"short": "ACD"}, protocol))
    with pytest.raises(ValueError, match="Matched design"):
        list(generate_controlled_queries({"whole": "ACDEFGHIKLMNPQRSTVWY"}, protocol | {"replacement_percent": [0]}))


def test_abstention_keeps_denominators_and_zero_accepted_risk_missing():
    frame = pd.DataFrame(dict(mode=["local", "local"], accepted=[False, False], exact=[True, False],
                              status=["ambiguous", "incomplete"], truth_in_observed_intervals=[True, False]))
    row = reporting_summary(frame, ["mode"]).iloc[0]
    assert row.queries == 2 and row.first_exact == 0.5 and row.accepted == 0
    assert pd.isna(row.accepted_error_rate) and row.incomplete == row.ambiguous == 1


def test_cluster_effect_uses_proteins_and_preserves_direction():
    rows = []
    for accession, n, interior in [("a", 1, True), ("b", 3, False)]:
        for i in range(n):
            for terminal in GROUPS:
                rows.append(dict(pair_id=f"{accession}{i}", seed=1, accession=accession,
                                 mode="local", length=20, pct=10, terminal_group=terminal,
                                 exact=interior if terminal == "interior" else not interior))
    effects, _ = paired_effects(pd.DataFrame(rows), PROTOCOL)
    # Per-query mean is +0.5; equal-protein mean is zero.
    assert effects.delta.tolist() == [0, 0, 0]
    assert effects.gains.tolist() == [3, 3, 3] and effects.losses.tolist() == [1, 1, 1]


@pytest.mark.parametrize("filename,column", [("results.csv", "score"), ("summary.csv", "exact")])
def test_audit_rejects_tampered_scores_and_stale_aggregates(tmp_path, filename, column):
    for name in ("queries.csv", "results.csv", "per_protein.csv", "summary.csv"):
        shutil.copyfile(ROOT / "results" / name, tmp_path / name)
    table = pd.read_csv(tmp_path / filename)
    table.loc[0, column] += 1
    table.to_csv(tmp_path / filename, index=False)
    with pytest.raises(ValueError, match="Metric mismatch|Stale or incorrect"):
        analyze.audit(tmp_path)


def test_cli_json_is_parseable_and_input_error_is_actionable():
    command = [sys.executable, str(ROOT / "fragment.py"), str(ROOT / "data/raw/P69905.fasta"),
               str(ROOT / "examples/query.fasta"), "--json"]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    report = json.loads(result.stdout)
    assert report["status"] == "unique" and report["first"]["start"] == 20
    bad = subprocess.run(command + ["--max-alignments", "0"], capture_output=True, text=True)
    assert bad.returncode == 2 and "must be positive" in bad.stderr
