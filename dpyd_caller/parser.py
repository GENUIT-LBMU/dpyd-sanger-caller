from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Union

from Bio import SeqIO


@dataclass
class SangerRead:
    sample_name: str
    sequence: str
    traces: dict[str, list[int]]
    peak_locations: list[int]
    quality: list[int]
    base_order: str


def _decode(value) -> str:
    return value.decode() if isinstance(value, (bytes, bytearray)) else str(value)


def parse_fsa(path: Union[str, Path]) -> SangerRead:
    record = SeqIO.read(str(path), "abi")
    raw = record.annotations["abif_raw"]

    base_order = _decode(raw.get("FWO_1", "GATC"))
    if len(base_order) != 4:
        base_order = "GATC"

    trace_keys = ("DATA9", "DATA10", "DATA11", "DATA12")
    traces = {base_order[i]: list(raw[trace_keys[i]]) for i in range(4)}

    peak_locations = list(raw.get("PLOC2") or raw.get("PLOC1") or [])
    quality = record.letter_annotations.get("phred_quality", [])

    return SangerRead(
        sample_name=record.name,
        sequence=str(record.seq),
        traces=traces,
        peak_locations=peak_locations,
        quality=list(quality),
        base_order=base_order,
    )
