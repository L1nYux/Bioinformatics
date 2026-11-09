"""Locate a protein fragment within a target FASTA sequence."""

import argparse
import json
from pathlib import Path

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices


ROOT = Path(__file__).resolve().parent
MODES = ("global", "local", "semiglobal")
AMINO = "ACDEFGHIKLMNPQRSTVWY"


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


def locate(aligner, target, query):
    alignment = next(iter(aligner.align(target, query)), None)
    if alignment is None or not len(alignment.aligned[0]):
        return dict(start=None, end=None, paired=0, coverage=0.0, score=None), alignment
    blocks = alignment.aligned[0]
    paired = int(sum(end - start for start, end in blocks))
    return dict(start=int(blocks[0, 0]), end=int(blocks[-1, 1]), paired=paired,
                coverage=paired / len(query), score=alignment.score), alignment


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", type=Path)
    parser.add_argument("query", type=Path)
    parser.add_argument("--mode", choices=MODES, default="semiglobal")
    args = parser.parse_args()
    try:
        _, target = read_protein(args.target)
        _, query = read_protein(args.query)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    config = json.loads((ROOT / "study.json").read_text())
    aligner = make_aligner(args.mode, config)
    prediction, alignment = locate(aligner, target, query)
    print(f"Mode: {args.mode}; matrix: {config['matrix']}; algorithm: {aligner.algorithm}")
    print(aligner)
    print(json.dumps(prediction))
    if prediction["start"] is not None:
        print(f"Target residues (1-based, inclusive): {prediction['start']+1}..{prediction['end']}")
        print(alignment)


if __name__ == "__main__":
    main()
