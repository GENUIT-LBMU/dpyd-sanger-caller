from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

CPIC_TABLE_PATH = Path(__file__).resolve().parent.parent / "data" / "cpic_table.json"


def load_cpic_table(path: Optional[Path] = None) -> dict:
    return json.loads((path or CPIC_TABLE_PATH).read_text())


def calculate_activity_score(
    variant_calls: dict[str, dict],
    all_variant_ids: Optional[list[str]] = None,
) -> dict:
    """Compute a diplotype activity score assuming trans-phasing for compound hets.

    variant_calls: {variant_id: {"genotype": str, "activity_value": float}}
    genotype ∈ {"hom_ref", "het", "hom_alt", "unclear", "discordant", "failed", "missing"}

    - "missing" = variant not tested for this patient; assumed wild-type for the score
      calculation (most likely scenario given low CPIC variant frequencies), but flagged
      in `untested` and `notes` so the operator knows the interpretation is incomplete.
    - "unclear" / "discordant" / "failed" = tested but inconclusive → reliable=False.

    Returns dict with score, allele activities, reliability flag, notes, tested count,
    and list of untested variants.
    """
    allele1_values = [1.0]
    allele2_values = [1.0]
    notes = []
    reliable = True

    tested_reliably = []
    for vid, call in variant_calls.items():
        gt = call["genotype"]
        av = call["activity_value"]
        if gt == "missing":
            continue
        if gt in ("unclear", "discordant", "failed"):
            reliable = False
            notes.append(f"{vid}: genotipo no confiable ({gt}) — no computado en el score")
            continue
        tested_reliably.append(vid)
        if gt == "hom_alt":
            allele1_values.append(av)
            allele2_values.append(av)
        elif gt == "het":
            if min(allele1_values) >= min(allele2_values):
                allele1_values.append(av)
            else:
                allele2_values.append(av)

    allele1 = min(allele1_values)
    allele2 = min(allele2_values)
    score = round(allele1 + allele2, 2) if reliable else None

    untested: list[str] = []
    if all_variant_ids:
        untested = [
            vid for vid in all_variant_ids
            if vid not in variant_calls or variant_calls[vid]["genotype"] == "missing"
        ]
        if untested:
            notes.append(
                f"Variantes no testeadas (asumidas wild-type para el cálculo): {', '.join(untested)}. "
                f"Interpretación CPIC **incompleta**; un Poor Metabolizer podría no detectarse."
            )

    return {
        "score": score,
        "allele1": allele1,
        "allele2": allele2,
        "reliable": reliable,
        "notes": notes,
        "tested_count": len(tested_reliably),
        "total_expected": len(all_variant_ids) if all_variant_ids else len(variant_calls),
        "untested": untested,
    }


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
