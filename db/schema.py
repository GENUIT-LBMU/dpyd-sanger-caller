from __future__ import annotations

import sqlite3
from pathlib import Path

DB_PATH = Path(__file__).resolve().parent.parent / "data" / "dpyd_sanger.db"


SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS patients (
    patient_id TEXT PRIMARY KEY,
    full_name TEXT NOT NULL,
    sex TEXT NOT NULL,
    birth_date TEXT NOT NULL,
    dni TEXT,
    diagnosis TEXT NOT NULL,
    tumor_stage TEXT,
    referring_physician TEXT,
    referring_institution TEXT,
    sample_date TEXT NOT NULL,
    notes TEXT,
    status TEXT NOT NULL DEFAULT 'pending' CHECK (status IN ('pending', 'reported')),
    entry_date TEXT NOT NULL DEFAULT (DATETIME('now')),
    report_date TEXT,
    created_at TEXT NOT NULL DEFAULT (DATETIME('now')),
    updated_at TEXT NOT NULL DEFAULT (DATETIME('now'))
);

CREATE TABLE IF NOT EXISTS variant_calls (
    patient_id TEXT NOT NULL,
    variant_id TEXT NOT NULL,
    genotype TEXT NOT NULL,
    ref_qc_pass INTEGER NOT NULL,
    sample_concordant INTEGER NOT NULL,
    fwd_call TEXT,
    rev_call TEXT,
    details_json TEXT NOT NULL,
    errors_json TEXT,
    analyzed_at TEXT NOT NULL DEFAULT (DATETIME('now')),
    PRIMARY KEY (patient_id, variant_id),
    FOREIGN KEY (patient_id) REFERENCES patients(patient_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS cpic_reports (
    patient_id TEXT PRIMARY KEY,
    activity_score REAL,
    allele1_av REAL,
    allele2_av REAL,
    phenotype TEXT,
    phenotype_code TEXT,
    dosing_label TEXT,
    dosing_recommendation TEXT,
    dosing_strength TEXT,
    reliable INTEGER NOT NULL,
    notes_json TEXT,
    reported_at TEXT NOT NULL DEFAULT (DATETIME('now')),
    FOREIGN KEY (patient_id) REFERENCES patients(patient_id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_patients_status ON patients(status);
CREATE INDEX IF NOT EXISTS idx_patients_entry_date ON patients(entry_date);
CREATE INDEX IF NOT EXISTS idx_patients_full_name ON patients(full_name);
CREATE INDEX IF NOT EXISTS idx_variant_calls_variant ON variant_calls(variant_id);
CREATE INDEX IF NOT EXISTS idx_cpic_phenotype ON cpic_reports(phenotype_code);
"""


def get_connection() -> sqlite3.Connection:
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def init_db() -> None:
    with get_connection() as conn:
        conn.executescript(SCHEMA_SQL)
        conn.commit()
