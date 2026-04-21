from __future__ import annotations

from pathlib import Path
from typing import Optional

REFERENCES_DIR = Path(__file__).resolve().parent.parent / "data" / "references"


def find_reference(variant_id: str, direction: str) -> Optional[Path]:
    """Return the reference .fsa/.ab1 file for a variant+direction if present."""
    assert direction in ("fwd", "rev")
    for ext in (".fsa", ".ab1"):
        p = REFERENCES_DIR / f"{variant_id}_ref_{direction}{ext}"
        if p.exists():
            return p
    return None


def references_status(variant_ids: list[str]) -> dict[str, dict[str, Optional[Path]]]:
    """Return {variant_id: {'fwd': Path|None, 'rev': Path|None}} for each variant."""
    return {
        vid: {"fwd": find_reference(vid, "fwd"), "rev": find_reference(vid, "rev")}
        for vid in variant_ids
    }


def all_references_present(variant_ids: list[str]) -> bool:
    status = references_status(variant_ids)
    return all(paths["fwd"] and paths["rev"] for paths in status.values())
