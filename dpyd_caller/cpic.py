from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

CPIC_TABLE_PATH = Path(__file__).resolve().parent.parent / "data" / "cpic_table.json"


def load_cpic_table(path: Optional[Path] = None) -> dict:
    return json.loads((path or CPIC_TABLE_PATH).read_text())


def calculate_activity_score(variant_calls: dict[str, dict]) -> dict:
    """Compute a diplotype activity score assuming trans-phasing for compound hets.

    variant_calls: {variant_id: {"genotype": str, "activity_value": float}}
    genotype ∈ {"hom_ref", "het", "hom_alt", "unclear", "discordant", "failed"}

    Returns {"score": float|None, "allele1": float, "allele2": float, "reliable": bool,
             "notes": [str]}.
    """
    allele1_values = [1.0]
    allele2_values = [1.0]
    notes = []
    reliable = True

    for vid, call in variant_calls.items():
        gt = call["genotype"]
        av = call["activity_value"]
        if gt in ("unclear", "discordant", "failed"):
            reliable = False
            notes.append(f"{vid}: genotipo no confiable ({gt}) — no computado en el score")
            continue
        if gt == "hom_alt":
            allele1_values.append(av)
            allele2_values.append(av)
        elif gt == "het":
            # trans assumption: asignar a la alelo "más sano" que tenga
            if min(allele1_values) >= min(allele2_values):
                allele1_values.append(av)
            else:
                allele2_values.append(av)

    allele1 = min(allele1_values)
    allele2 = min(allele2_values)
    score = round(allele1 + allele2, 2) if reliable else None
    return {"score": score, "allele1": allele1, "allele2": allele2, "reliable": reliable, "notes": notes}


def phenotype_from_activity_score(activity_score: float, table: Optional[dict] = None) -> Optional[dict]:
    if activity_score is None:
        return None
    table = table or load_cpic_table()
    for row in table["activity_score_to_phenotype"]:
        if row["min"] <= activity_score <= row["max"]:
            return row
    return None


def dosing_for(activity_score: Optional[float], phenotype_code: Optional[str], table: Optional[dict] = None) -> dict:
    if phenotype_code is None or activity_score is None:
        return {}
    table = table or load_cpic_table()
    dosing = table["dosing_recommendations"]
    if phenotype_code == "IM":
        return dosing["IM_1.5"] if activity_score == 1.5 else dosing["IM_1.0"]
    return dosing.get(phenotype_code, {})
