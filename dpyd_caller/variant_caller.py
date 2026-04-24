from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Optional

from .parser import SangerRead
from .aligner import (
    COMPLEMENT,
    AlignmentResult,
    align_read,
    map_amplicon_to_query_index,
    query_index_to_original_read_index,
)


HETERO_MIN_FRAC = 0.25
DOMINANT_FRAC = 0.70
MIN_ALIGNMENT_SCORE = 50.0


@dataclass
class ReadCall:
    direction: str
    call: str
    heights: dict[str, int]
    ref_channel: str
    alt_channel: str
    ref_height: int
    alt_height: int
    ref_frac: float
    alt_frac: float
    peak_location: int
    read_index: int
    alignment_score: float
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class VariantAnalysis:
    variant_id: str
    ref_fwd: Optional[ReadCall]
    ref_rev: Optional[ReadCall]
    sample_fwd: Optional[ReadCall]
    sample_rev: Optional[ReadCall]
    reference_qc_pass: bool
    sample_concordant: bool
    sample_genotype: str
    analysis_mode: str = "bidirectional"   # bidirectional | fwd_only | rev_only | missing
    errors: list[str] = field(default_factory=list)


def _call_single_read(
    read: SangerRead,
    is_reverse: bool,
    amplicon_sense: str,
    variant_pos: int,
    ref_base_sense: str,
    alt_base_sense: str,
    direction_label: str,
) -> ReadCall:
    notes: list[str] = []
    alignment = align_read(read.sequence, amplicon_sense, is_reverse)

    if alignment.score < MIN_ALIGNMENT_SCORE:
        notes.append(f"Alineamiento débil (score={alignment.score:.0f})")

    q_index = map_amplicon_to_query_index(alignment, variant_pos)
    if q_index is None:
        raise ValueError(
            f"La posición de la variante (offset {variant_pos}) no está cubierta por la lectura {direction_label}."
        )

    read_index = query_index_to_original_read_index(q_index, len(read.sequence), is_reverse)
    if read_index < 0 or read_index >= len(read.peak_locations):
        raise ValueError(
            f"Índice de base fuera de rango en la lectura {direction_label} (read_index={read_index})."
        )

    peak_loc = read.peak_locations[read_index]

    if is_reverse:
        ref_channel = COMPLEMENT[ref_base_sense]
        alt_channel = COMPLEMENT[alt_base_sense]
    else:
        ref_channel = ref_base_sense
        alt_channel = alt_base_sense

    heights = {b: int(read.traces[b][peak_loc]) for b in read.traces}
    total = sum(heights.values()) or 1
    ref_h = heights[ref_channel]
    alt_h = heights[alt_channel]
    ref_frac = ref_h / total
    alt_frac = alt_h / total

    if ref_frac >= DOMINANT_FRAC:
        call = "hom_ref"
    elif alt_frac >= DOMINANT_FRAC:
        call = "hom_alt"
    elif ref_frac >= HETERO_MIN_FRAC and alt_frac >= HETERO_MIN_FRAC:
        call = "het"
    else:
        call = "unclear"
        notes.append("Señal ambigua: picos ref/alt por debajo del umbral de heterocigota")

    return ReadCall(
        direction=direction_label,
        call=call,
        heights=heights,
        ref_channel=ref_channel,
        alt_channel=alt_channel,
        ref_height=ref_h,
        alt_height=alt_h,
        ref_frac=ref_frac,
        alt_frac=alt_frac,
        peak_location=peak_loc,
        read_index=read_index,
        alignment_score=alignment.score,
        notes=notes,
    )


def analyze_variant(
    variant_id: str,
    variant: dict,
    ref_fwd: Optional[SangerRead] = None,
    ref_rev: Optional[SangerRead] = None,
    sample_fwd: Optional[SangerRead] = None,
    sample_rev: Optional[SangerRead] = None,
) -> VariantAnalysis:
    """Analyze a variant with any combination of reads present.

    All four reads are optional. If both sample_fwd and sample_rev are None,
    the variant is reported as "missing". Single-direction analysis is supported
    (analysis_mode will be fwd_only or rev_only).
    """
    amplicon_sense = variant["amplicon_sense"]
    variant_pos = variant["variant_offset_in_amplicon"]
    ref_base = variant["ref_allele_sense"].upper()
    alt_base = variant["alt_allele_sense"].upper()

    errors: list[str] = []
    calls: dict[str, Optional[ReadCall]] = {
        "ref_fwd": None, "ref_rev": None,
        "sample_fwd": None, "sample_rev": None,
    }
    inputs = [
        ("ref_fwd", ref_fwd, False),
        ("ref_rev", ref_rev, True),
        ("sample_fwd", sample_fwd, False),
        ("sample_rev", sample_rev, True),
    ]
    for key, read, is_rev in inputs:
        if read is None:
            continue
        try:
            calls[key] = _call_single_read(
                read, is_rev, amplicon_sense, variant_pos, ref_base, alt_base, key,
            )
        except Exception as exc:
            errors.append(f"{key}: {exc}")

    # Reference QC: pass if every available reference read shows hom_ref.
    ref_calls_present = [c for c in (calls["ref_fwd"], calls["ref_rev"]) if c is not None]
    if ref_calls_present:
        ref_qc = all(c.call == "hom_ref" for c in ref_calls_present)
    else:
        ref_qc = False

    sample_fwd_call = calls["sample_fwd"]
    sample_rev_call = calls["sample_rev"]

    if sample_fwd_call is None and sample_rev_call is None:
        analysis_mode = "missing"
        genotype = "missing"
        sample_concordant = False
    elif sample_fwd_call is not None and sample_rev_call is not None:
        analysis_mode = "bidirectional"
        fc, rc = sample_fwd_call.call, sample_rev_call.call
        if fc == rc and fc != "unclear":
            genotype = fc
            sample_concordant = True
        elif fc == "unclear" and rc == "unclear":
            genotype = "unclear"
            sample_concordant = False
        elif fc == "unclear":
            genotype = rc
            sample_concordant = False
            errors.append("fwd ambigua — llamado basado sólo en rev")
        elif rc == "unclear":
            genotype = fc
            sample_concordant = False
            errors.append("rev ambigua — llamado basado sólo en fwd")
        else:
            genotype = "discordant"
            sample_concordant = False
    elif sample_fwd_call is not None:
        analysis_mode = "fwd_only"
        genotype = sample_fwd_call.call
        sample_concordant = False
        errors.append("Sólo forward disponible — sin confirmación bidireccional")
    else:
        analysis_mode = "rev_only"
        genotype = sample_rev_call.call
        sample_concordant = False
        errors.append("Sólo reverse disponible — sin confirmación bidireccional")

    return VariantAnalysis(
        variant_id=variant_id,
        ref_fwd=calls["ref_fwd"],
        ref_rev=calls["ref_rev"],
        sample_fwd=sample_fwd_call,
        sample_rev=sample_rev_call,
        reference_qc_pass=ref_qc,
        sample_concordant=sample_concordant,
        sample_genotype=genotype,
        analysis_mode=analysis_mode,
        errors=errors,
    )
