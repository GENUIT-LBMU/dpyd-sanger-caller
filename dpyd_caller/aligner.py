from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq


COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


@dataclass
class AlignmentResult:
    score: float
    target_start: int
    target_end: int
    query_start: int
    query_end: int
    aligned_blocks: list[tuple[tuple[int, int], tuple[int, int]]]
    oriented_query: str


def _make_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    return aligner


def align_read(read_seq: str, amplicon_sense: str, is_reverse: bool) -> AlignmentResult:
    oriented = reverse_complement(read_seq) if is_reverse else read_seq
    aligner = _make_aligner()
    alignments = aligner.align(amplicon_sense, oriented)
    best = alignments[0]

    target_blocks, query_blocks = best.aligned
    blocks = [((int(t0), int(t1)), (int(q0), int(q1)))
              for (t0, t1), (q0, q1) in zip(target_blocks, query_blocks)]

    return AlignmentResult(
        score=float(best.score),
        target_start=int(target_blocks[0][0]),
        target_end=int(target_blocks[-1][1]),
        query_start=int(query_blocks[0][0]),
        query_end=int(query_blocks[-1][1]),
        aligned_blocks=blocks,
        oriented_query=oriented,
    )


def map_amplicon_to_query_index(alignment: AlignmentResult, amplicon_pos: int) -> Optional[int]:
    for (t_start, t_end), (q_start, q_end) in alignment.aligned_blocks:
        if t_start <= amplicon_pos < t_end:
            return q_start + (amplicon_pos - t_start)
    return None


def query_index_to_original_read_index(q_index: int, original_length: int, is_reverse: bool) -> int:
    if is_reverse:
        return original_length - 1 - q_index
    return q_index
