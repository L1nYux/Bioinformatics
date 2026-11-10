"""Locate a protein fragment within a target FASTA sequence."""

import argparse
import json
from itertools import islice
from pathlib import Path

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices


ROOT = Path(__file__).resolve().parent
MODES = ("global", "local", "semiglobal")
AMINO = "ACDEFGHIKLMNPQRSTVWY"


def validate_sequences(target, query):
    for label, sequence in (("target", target), ("query", query)):
        if not sequence or set(sequence) - set(AMINO):
            raise ValueError(f"{label}: expected nonempty standard amino-acid sequence")


def read_protein(path):
    record = SeqIO.read(path, "fasta")
    sequence = str(record.seq).upper()
    if not sequence or set(sequence) - set(AMINO):
        raise ValueError(f"{path}: expected one protein with the 20 standard amino acids")
    return record, sequence


def make_aligner(mode, config):
    if mode not in MODES:
        raise ValueError(f"Unknown alignment mode: {mode}")
    aligner = Align.PairwiseAligner(
        mode="local" if mode == "local" else "global",
        substitution_matrix=substitution_matrices.load(config["matrix"]),
        open_gap_score=config["open_gap"], extend_gap_score=config["extend_gap"],
    )
    if mode == "semiglobal":
        # Deletions in the query let the target extend past either fragment end.
        aligner.open_left_deletion_score = aligner.extend_left_deletion_score = 0
        aligner.open_right_deletion_score = aligner.extend_right_deletion_score = 0
    return aligner


def prediction_from_alignment(alignment, query_length):
    if alignment is None or not len(alignment.aligned[0]):
        return dict(start=None, end=None, paired=0, coverage=0.0, score=None)
    blocks = alignment.aligned[0]
    paired = int(sum(end - start for start, end in blocks))
    return dict(start=int(blocks[0, 0]), end=int(blocks[-1, 1]), paired=paired,
                coverage=paired / query_length, score=alignment.score)


def locate(aligner, target, query):
    validate_sequences(target, query)
    alignment = next(iter(aligner.align(target, query)), None)
    return prediction_from_alignment(alignment, len(query)), alignment


def locate_candidates(aligner, target, query, max_alignments=1000):
    """Report observed optimal intervals without using truth or counting all paths.

    One extra path is requested to distinguish exhaustion from truncation. A
    unique interval is asserted only after complete enumeration. Different gap
    paths with the same target endpoints count as one interval.
    """
    validate_sequences(target, query)
    if type(max_alignments) is not int or max_alignments < 1:
        raise ValueError("max_alignments must be a positive integer")
    first = prediction_from_alignment(None, len(query))
    intervals = {}
    examined, complete = 0, True
    iterator = iter(aligner.align(target, query))
    for index, alignment in enumerate(islice(iterator, max_alignments + 1)):
        if index == max_alignments:
            complete = False
            break
        examined += 1
        result = prediction_from_alignment(alignment, len(query))
        if index == 0:
            first = result
        if result["start"] is None:
            continue
        key = result["start"], result["end"]
        item = intervals.setdefault(key, dict(start=key[0], end=key[1],
                                             coverage_min=result["coverage"],
                                             coverage_max=result["coverage"]))
        item["coverage_min"] = min(item["coverage_min"], result["coverage"])
        item["coverage_max"] = max(item["coverage_max"], result["coverage"])
    status = ("incomplete" if not complete else "no_hit" if not intervals
              else "unique" if len(intervals) == 1 else "ambiguous")
    return dict(first=first, status=status, enumeration_complete=complete,
                alignments_examined=examined, interval_count_observed=len(intervals),
                intervals=[intervals[key] for key in sorted(intervals)])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", type=Path)
    parser.add_argument("query", type=Path)
    parser.add_argument("--mode", choices=MODES, default="semiglobal")
    parser.add_argument("--json", action="store_true", help="emit a machine-readable candidate report")
    parser.add_argument("--max-alignments", type=int, default=1000,
                        help="maximum optimal paths to examine (default: 1000)")
    args = parser.parse_args()
    try:
        _, target = read_protein(args.target)
        _, query = read_protein(args.query)
        if args.max_alignments < 1:
            raise ValueError("--max-alignments must be positive")
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    config = json.loads((ROOT / "study.json").read_text())
    aligner = make_aligner(args.mode, config)
    report = locate_candidates(aligner, target, query, args.max_alignments)
    if args.json:
        print(json.dumps(dict(schema_version=1, mode=args.mode, scoring=config, **report)))
        return
    prediction, alignment = locate(aligner, target, query)
    print(f"Mode: {args.mode}; matrix: {config['matrix']}; algorithm: {aligner.algorithm}")
    print(aligner)
    print(json.dumps(prediction))
    print(f"Location status: {report['status']}; observed intervals: {report['interval_count_observed']}; "
          f"enumeration complete: {report['enumeration_complete']}")
    print(json.dumps(report["intervals"]))
    if prediction["start"] is not None:
        print(f"Target residues (1-based, inclusive): {prediction['start']+1}..{prediction['end']}")
        print(alignment)


if __name__ == "__main__":
    main()
