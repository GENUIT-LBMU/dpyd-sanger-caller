"""Verify that an uploaded .fsa/.ab1 matches a specific variant's amplicon."""
from __future__ import annotations

import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Union

from Bio import SeqIO

from .aligner import align_read


MIN_VERIFY_SCORE = 200.0   # below this: file does NOT correspond to this amplicon
STRONG_VERIFY_SCORE = 350.0  # high confidence match


@dataclass
class FileVerification:
    ok: bool
    score: float
    detected_orientation: str        # "forward" | "reverse" | "unknown"
    expected_orientation: str        # from the slot ("fwd" | "rev")
    orientation_matches: bool
    sample_name: str
    abif_sample_name: str
    message: str
    level: str                        # "ok" | "warn" | "error"


def _parse_abif(path_or_bytes: Union[str, Path, bytes]):
    if isinstance(path_or_bytes, (bytes, bytearray)):
        with tempfile.NamedTemporaryFile(suffix=".fsa", delete=False) as tmp:
            tmp.write(path_or_bytes)
            tmp_path = tmp.name
        try:
            record = SeqIO.read(tmp_path, "abi")
        finally:
            Path(tmp_path).unlink(missing_ok=True)
    else:
        record = SeqIO.read(str(path_or_bytes), "abi")
    return record


def _decode(value) -> str:
    return value.decode() if isinstance(value, (bytes, bytearray)) else str(value) if value else ""


def verify_file_matches_variant(
    file_bytes: bytes,
    variant: dict,
    expected_direction: str,
) -> FileVerification:
    """Verify that a Sanger file corresponds to the given variant's amplicon.

    Tries both forward and reverse-complement orientations and picks the better one.
    Returns ok=True if the best alignment score exceeds MIN_VERIFY_SCORE.
    Flags orientation mismatch if the best orientation doesn't match the slot.
    """
    assert expected_direction in ("fwd", "rev")

    try:
        record = _parse_abif(file_bytes)
    except Exception as exc:
        return FileVerification(
            ok=False, score=0.0,
            detected_orientation="unknown",
            expected_orientation=expected_direction,
            orientation_matches=False,
            sample_name="", abif_sample_name="",
            message=f"No pude parsear el archivo: {exc}",
            level="error",
        )

    abif_sample_name = _decode((record.annotations.get("abif_raw") or {}).get("SMPL1", "")).strip()
    sample_name = record.name or ""
    amplicon = variant.get("amplicon_sense")
    if not amplicon:
        return FileVerification(
            ok=False, score=0.0,
            detected_orientation="unknown",
            expected_orientation=expected_direction,
            orientation_matches=False,
            sample_name=sample_name, abif_sample_name=abif_sample_name,
            message=f"La variante {variant.get('name')} no tiene amplicón cargado.",
            level="error",
        )

    seq = str(record.seq)
    if len(seq) < 50:
        return FileVerification(
            ok=False, score=0.0,
            detected_orientation="unknown",
            expected_orientation=expected_direction,
            orientation_matches=False,
            sample_name=sample_name, abif_sample_name=abif_sample_name,
            message=f"Archivo con secuencia muy corta ({len(seq)} bases).",
            level="error",
        )

    try:
        fwd_score = align_read(seq, amplicon, is_reverse=False).score
        rev_score = align_read(seq, amplicon, is_reverse=True).score
    except Exception as exc:
        return FileVerification(
            ok=False, score=0.0,
            detected_orientation="unknown",
            expected_orientation=expected_direction,
            orientation_matches=False,
            sample_name=sample_name, abif_sample_name=abif_sample_name,
            message=f"Error en alineamiento: {exc}",
            level="error",
        )

    best_score = max(fwd_score, rev_score)
    if fwd_score >= rev_score:
        detected_orientation = "forward"
    else:
        detected_orientation = "reverse"

    if best_score < MIN_VERIFY_SCORE:
        return FileVerification(
            ok=False, score=best_score,
            detected_orientation="unknown",
            expected_orientation=expected_direction,
            orientation_matches=False,
            sample_name=sample_name, abif_sample_name=abif_sample_name,
            message=(
                f"El archivo NO corresponde al amplicón de {variant['name']} "
                f"(score máx {best_score:.0f} en ambas orientaciones, umbral {MIN_VERIFY_SCORE:.0f}). "
                f"Probable archivo mal asignado."
            ),
            level="error",
        )

    expected_dir_word = "forward" if expected_direction == "fwd" else "reverse"
    orientation_matches = (detected_orientation == expected_dir_word)

    if not orientation_matches:
        return FileVerification(
            ok=False, score=best_score,
            detected_orientation=detected_orientation,
            expected_orientation=expected_direction,
            orientation_matches=False,
            sample_name=sample_name, abif_sample_name=abif_sample_name,
            message=(
                f"El archivo alinea a {variant['name']} pero en orientación **{detected_orientation}**, "
                f"no **{expected_dir_word}** como indica el slot. Probable archivo en slot equivocado."
            ),
            level="warn",
        )

    confidence = "alta" if best_score >= STRONG_VERIFY_SCORE else "media"
    return FileVerification(
        ok=True, score=best_score,
        detected_orientation=detected_orientation,
        expected_orientation=expected_direction,
        orientation_matches=True,
        sample_name=sample_name, abif_sample_name=abif_sample_name,
        message=(
            f"Corresponde a {variant['name']} ({detected_orientation}, score {best_score:.0f}, confianza {confidence})."
            + (f" Sample name: {abif_sample_name}." if abif_sample_name else "")
        ),
        level="ok",
    )
