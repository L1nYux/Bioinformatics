import json

import pytest

from experiment import generate_queries, study_targets
from fragment import MODES, ROOT, locate, make_aligner


CONFIG = json.loads((ROOT / "study.json").read_text())


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("query,start,end", [("ACDE", 0, 4), ("EFGH", 3, 7), ("KLMN", 8, 12)])
def test_unique_substring(mode, query, start, end):
    result, _ = locate(make_aligner(mode, CONFIG), "ACDEFGHIKLMN", query)
    assert (result["start"], result["end"], result["coverage"]) == (start, end, 1)


def test_free_target_overhangs_only():
    aligner = make_aligner("semiglobal", CONFIG)
    base = aligner.score("EFGH", "EFGH")
    assert aligner.score("ACDEFGHIKLMN", "EFGH") == base
    assert aligner.score("EFGH", "ACDEFGHIKLMN") < base


def test_no_positive_local_and_tied_positions():
    aligner = make_aligner("local", CONFIG)
    result, _ = locate(aligner, "DDDD", "WWWW")
    assert result["start"] is None and result["coverage"] == 0
    result, _ = locate(aligner, "ACDEFGHKEFGH", "EFGH")
    assert (result["start"], result["end"]) in {(3, 7), (8, 12)}


def test_paired_columns_exclude_internal_gap():
    result, _ = locate(make_aligner("global", CONFIG), "ACDEFG", "ACDWEFG")
    assert result["paired"] == 6 and result["coverage"] == pytest.approx(6/7)


def test_study_queries():
    targets = study_targets()
    rows = list(generate_queries(targets, CONFIG))
    assert len({row["query_id"] for row in rows}) == len(rows) == 270
    assert rows == list(generate_queries(targets, CONFIG))
    for row in rows:
        original = targets[row["accession"]][row["true_start"]:row["true_end"]]
        assert original == row["original"]
        changes = [i for i, (a,b) in enumerate(zip(original,row["query"], strict=True)) if a != b]
        assert changes == list(map(int, row["substitutions"].split()))
        assert len(changes) == row["length"] * row["pct"] // 100
