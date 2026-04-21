from __future__ import annotations

import shutil
from pathlib import Path
from typing import Optional

PATIENT_FILES_ROOT = Path(__file__).resolve().parent.parent / "data" / "patient_files"


def patient_dir(patient_id: str) -> Path:
    return PATIENT_FILES_ROOT / patient_id


def patient_file_path(patient_id: str, variant_id: str, direction: str, ext: str = "fsa") -> Path:
    assert direction in ("fwd", "rev")
    return patient_dir(patient_id) / f"{variant_id}_{direction}.{ext}"


def save_patient_file(
    patient_id: str,
    variant_id: str,
    direction: str,
    content: bytes,
    original_name: str = "",
) -> Path:
    d = patient_dir(patient_id)
    d.mkdir(parents=True, exist_ok=True)
    ext = "ab1" if original_name.lower().endswith(".ab1") else "fsa"
    # Remove any competing file with the other extension so find_patient_file stays deterministic
    for e in ("fsa", "ab1"):
        alt = patient_file_path(patient_id, variant_id, direction, e)
        if alt.exists() and e != ext:
            alt.unlink()
    p = patient_file_path(patient_id, variant_id, direction, ext)
    p.write_bytes(content)
    return p


def find_patient_file(patient_id: str, variant_id: str, direction: str) -> Optional[Path]:
    for ext in ("fsa", "ab1"):
        p = patient_file_path(patient_id, variant_id, direction, ext)
        if p.exists():
            return p
    return None


def patient_files_status(patient_id: str, variant_ids: list[str]) -> dict[str, dict[str, Optional[Path]]]:
    return {
        vid: {
            "fwd": find_patient_file(patient_id, vid, "fwd"),
            "rev": find_patient_file(patient_id, vid, "rev"),
        }
        for vid in variant_ids
    }


def delete_patient_files(patient_id: str) -> None:
    d = patient_dir(patient_id)
    if d.exists():
        shutil.rmtree(d)
