from __future__ import annotations

import json
from dataclasses import asdict
from typing import Any, Optional

from .schema import get_connection, init_db


PATIENT_FIELDS = [
    "patient_id", "full_name", "sex", "birth_date", "dni",
    "diagnosis", "tumor_stage", "referring_physician", "referring_institution",
    "sample_date", "notes",
]


def _row_to_dict(row) -> dict:
    return dict(row) if row is not None else None


def upsert_patient(data: dict) -> str:
    """Insert or update a patient. Returns patient_id."""
    init_db()
    fields = {k: data.get(k) for k in PATIENT_FIELDS}
    if not fields["patient_id"]:
        raise ValueError("patient_id es requerido")

    with get_connection() as conn:
        existing = conn.execute(
            "SELECT patient_id FROM patients WHERE patient_id = ?",
            (fields["patient_id"],),
        ).fetchone()
        if existing:
            placeholders = ", ".join(f"{k} = ?" for k in PATIENT_FIELDS if k != "patient_id")
            values = [fields[k] for k in PATIENT_FIELDS if k != "patient_id"]
            conn.execute(
                f"UPDATE patients SET {placeholders}, updated_at = DATETIME('now') "
                f"WHERE patient_id = ?",
                (*values, fields["patient_id"]),
            )
        else:
            cols = ", ".join(PATIENT_FIELDS)
            qmarks = ", ".join("?" for _ in PATIENT_FIELDS)
            values = [fields[k] for k in PATIENT_FIELDS]
            conn.execute(
                f"INSERT INTO patients ({cols}) VALUES ({qmarks})",
                values,
            )
        conn.commit()
    return fields["patient_id"]


def get_patient(patient_id: str) -> Optional[dict]:
    init_db()
    with get_connection() as conn:
        row = conn.execute(
            "SELECT * FROM patients WHERE patient_id = ?",
            (patient_id,),
        ).fetchone()
    return _row_to_dict(row)


def list_patients(search: str = "", status: Optional[str] = None, limit: int = 500) -> list[dict]:
    init_db()
    query = "SELECT * FROM patients WHERE 1=1"
    params: list[Any] = []
    if search:
        query += " AND (patient_id LIKE ? OR full_name LIKE ? OR dni LIKE ?)"
        like = f"%{search}%"
        params.extend([like, like, like])
    if status:
        query += " AND status = ?"
        params.append(status)
    query += " ORDER BY entry_date DESC LIMIT ?"
    params.append(limit)
    with get_connection() as conn:
        rows = conn.execute(query, params).fetchall()
    return [dict(r) for r in rows]


def delete_patient(patient_id: str) -> None:
    init_db()
    with get_connection() as conn:
        conn.execute("DELETE FROM patients WHERE patient_id = ?", (patient_id,))
        conn.commit()
    try:
        from dpyd_caller.patient_files import delete_patient_files
        delete_patient_files(patient_id)
    except Exception:
        pass


def save_variant_call(patient_id: str, variant_id: str, analysis) -> None:
    """Save a VariantAnalysis object to the DB."""
    init_db()
    details = {
        "ref_fwd": asdict(analysis.ref_fwd) if analysis.ref_fwd else None,
        "ref_rev": asdict(analysis.ref_rev) if analysis.ref_rev else None,
        "sample_fwd": asdict(analysis.sample_fwd) if analysis.sample_fwd else None,
        "sample_rev": asdict(analysis.sample_rev) if analysis.sample_rev else None,
    }
    fwd_call = analysis.sample_fwd.call if analysis.sample_fwd else None
    rev_call = analysis.sample_rev.call if analysis.sample_rev else None
    with get_connection() as conn:
        conn.execute(
            """
            INSERT INTO variant_calls
              (patient_id, variant_id, genotype, ref_qc_pass, sample_concordant,
               fwd_call, rev_call, details_json, errors_json, analyzed_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, DATETIME('now'))
            ON CONFLICT(patient_id, variant_id) DO UPDATE SET
              genotype = excluded.genotype,
              ref_qc_pass = excluded.ref_qc_pass,
              sample_concordant = excluded.sample_concordant,
              fwd_call = excluded.fwd_call,
              rev_call = excluded.rev_call,
              details_json = excluded.details_json,
              errors_json = excluded.errors_json,
              analyzed_at = DATETIME('now')
            """,
            (
                patient_id,
                variant_id,
                analysis.sample_genotype,
                int(bool(analysis.reference_qc_pass)),
                int(bool(analysis.sample_concordant)),
                fwd_call,
                rev_call,
                json.dumps(details),
                json.dumps(analysis.errors),
            ),
        )
        conn.commit()


def get_variant_calls(patient_id: str) -> dict[str, dict]:
    init_db()
    with get_connection() as conn:
        rows = conn.execute(
            "SELECT * FROM variant_calls WHERE patient_id = ?",
            (patient_id,),
        ).fetchall()
    result = {}
    for r in rows:
        d = dict(r)
        d["details"] = json.loads(d.pop("details_json")) if d.get("details_json") else {}
        d["errors"] = json.loads(d.pop("errors_json")) if d.get("errors_json") else []
        result[d["variant_id"]] = d
    return result


def save_cpic_report(patient_id: str, cpic_result: dict, phenotype: dict, dosing: dict) -> None:
    init_db()
    with get_connection() as conn:
        conn.execute(
            """
            INSERT INTO cpic_reports
              (patient_id, activity_score, allele1_av, allele2_av, phenotype, phenotype_code,
               dosing_label, dosing_recommendation, dosing_strength, reliable, notes_json, reported_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, DATETIME('now'))
            ON CONFLICT(patient_id) DO UPDATE SET
              activity_score = excluded.activity_score,
              allele1_av = excluded.allele1_av,
              allele2_av = excluded.allele2_av,
              phenotype = excluded.phenotype,
              phenotype_code = excluded.phenotype_code,
              dosing_label = excluded.dosing_label,
              dosing_recommendation = excluded.dosing_recommendation,
              dosing_strength = excluded.dosing_strength,
              reliable = excluded.reliable,
              notes_json = excluded.notes_json,
              reported_at = DATETIME('now')
            """,
            (
                patient_id,
                cpic_result.get("score"),
                cpic_result.get("allele1"),
                cpic_result.get("allele2"),
                (phenotype or {}).get("phenotype"),
                (phenotype or {}).get("code"),
                (dosing or {}).get("label"),
                (dosing or {}).get("recommendation"),
                (dosing or {}).get("strength"),
                int(bool(cpic_result.get("reliable"))),
                json.dumps(cpic_result.get("notes", [])),
            ),
        )
        conn.commit()


def get_cpic_report(patient_id: str) -> Optional[dict]:
    init_db()
    with get_connection() as conn:
        row = conn.execute(
            "SELECT * FROM cpic_reports WHERE patient_id = ?",
            (patient_id,),
        ).fetchone()
    if row is None:
        return None
    d = dict(row)
    d["notes"] = json.loads(d.pop("notes_json")) if d.get("notes_json") else []
    return d


def mark_reported(patient_id: str) -> None:
    init_db()
    with get_connection() as conn:
        conn.execute(
            "UPDATE patients SET status = 'reported', report_date = DATETIME('now'), "
            "updated_at = DATETIME('now') WHERE patient_id = ?",
            (patient_id,),
        )
        conn.commit()


def mark_pending(patient_id: str) -> None:
    init_db()
    with get_connection() as conn:
        conn.execute(
            "UPDATE patients SET status = 'pending', report_date = NULL, "
            "updated_at = DATETIME('now') WHERE patient_id = ?",
            (patient_id,),
        )
        conn.commit()


def dashboard_stats() -> dict:
    init_db()
    with get_connection() as conn:
        total = conn.execute("SELECT COUNT(*) AS n FROM patients").fetchone()["n"]
        pending = conn.execute(
            "SELECT COUNT(*) AS n FROM patients WHERE status = 'pending'"
        ).fetchone()["n"]
        reported = conn.execute(
            "SELECT COUNT(*) AS n FROM patients WHERE status = 'reported'"
        ).fetchone()["n"]
        pheno_rows = conn.execute(
            "SELECT phenotype_code, COUNT(*) AS n FROM cpic_reports "
            "WHERE phenotype_code IS NOT NULL GROUP BY phenotype_code"
        ).fetchall()
        phenotypes = {r["phenotype_code"]: r["n"] for r in pheno_rows}
    return {
        "total_patients": total,
        "pending": pending,
        "reported": reported,
        "phenotypes": phenotypes,
    }


def variant_frequencies(variant_ids: list[str]) -> list[dict]:
    """For each variant, count genotype calls across all patients."""
    init_db()
    results = []
    with get_connection() as conn:
        for vid in variant_ids:
            rows = conn.execute(
                "SELECT genotype, COUNT(*) AS n FROM variant_calls "
                "WHERE variant_id = ? GROUP BY genotype",
                (vid,),
            ).fetchall()
            counts = {r["genotype"]: r["n"] for r in rows}
            hom_ref = counts.get("hom_ref", 0)
            het = counts.get("het", 0)
            hom_alt = counts.get("hom_alt", 0)
            total = hom_ref + het + hom_alt
            carrier = het + hom_alt
            alt_alleles = het + 2 * hom_alt
            total_alleles = 2 * total
            results.append({
                "variant_id": vid,
                "hom_ref": hom_ref,
                "het": het,
                "hom_alt": hom_alt,
                "unclear": counts.get("unclear", 0),
                "discordant": counts.get("discordant", 0),
                "failed": counts.get("failed", 0),
                "total_confident": total,
                "carrier_count": carrier,
                "carrier_pct": round(100 * carrier / total, 2) if total else None,
                "maf_pct": round(100 * alt_alleles / total_alleles, 2) if total_alleles else None,
            })
    return results


def recent_patients(limit: int = 10) -> list[dict]:
    init_db()
    with get_connection() as conn:
        rows = conn.execute(
            """
            SELECT p.patient_id, p.full_name, p.status, p.entry_date, p.report_date,
                   r.phenotype, r.phenotype_code, r.activity_score
            FROM patients p
            LEFT JOIN cpic_reports r ON r.patient_id = p.patient_id
            ORDER BY p.entry_date DESC
            LIMIT ?
            """,
            (limit,),
        ).fetchall()
    return [dict(r) for r in rows]
