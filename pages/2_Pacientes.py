"""Lista buscable de pacientes con genotipos por variante."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import streamlit as st

from dpyd_caller.auth import require_auth, logout_button
from db import (
    get_cpic_report,
    get_variant_calls,
    init_db,
    list_patients,
)

VARIANTS_PATH = Path(__file__).parent.parent / "data" / "variants.json"

CALL_EMOJI = {
    "hom_ref": "🟢", "het": "🟡", "hom_alt": "🔴",
    "unclear": "⚪", "discordant": "⚠️", "failed": "❌",
}
CALL_LABELS = {
    "hom_ref": "WT/WT", "het": "WT/var", "hom_alt": "var/var",
    "unclear": "ambig", "discordant": "disc", "failed": "fail",
}


@st.cache_data(ttl=30)
def load_variants() -> dict:
    return json.loads(VARIANTS_PATH.read_text())


def build_table(patients: list[dict], variants: dict) -> pd.DataFrame:
    rows = []
    for p in patients:
        calls = get_variant_calls(p["patient_id"])
        cpic = get_cpic_report(p["patient_id"])
        row = {
            "ID": p["patient_id"],
            "Nombre": p["full_name"],
            "Sexo": p["sex"],
            "Dx": p["diagnosis"],
            "Status": "🟡 Pendiente" if p["status"] == "pending" else "✅ Reportado",
            "Ingreso": p["entry_date"][:10] if p["entry_date"] else "—",
            "Reporte": p["report_date"][:10] if p["report_date"] else "—",
        }
        for vid, v in variants.items():
            call = calls.get(vid)
            if call:
                g = call["genotype"]
                row[v["name"]] = f"{CALL_EMOJI.get(g,'')} {CALL_LABELS.get(g, g)}"
            else:
                row[v["name"]] = "—"
        row["Fenotipo"] = (cpic or {}).get("phenotype", "—") or "—"
        row["AS"] = (cpic or {}).get("activity_score") if cpic else None
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    st.set_page_config(page_title="Pacientes", layout="wide", page_icon="🧑‍⚕️")
    require_auth()
    init_db()
    variants = load_variants()

    with st.sidebar:
        st.markdown("## 🧬 DPYD Sanger")
        st.divider()
        logout_button()

    st.title("🧑‍⚕️ Pacientes")

    top_cols = st.columns([3, 1, 1])
    with top_cols[0]:
        search = st.text_input("Buscar por ID, nombre o DNI", placeholder="Empezá a tipear...")
    with top_cols[1]:
        status_filter = st.selectbox("Estado", ["Todos", "Pendiente", "Reportado"])
    with top_cols[2]:
        st.write("")
        st.write("")
        if st.button("➕ Nuevo paciente", use_container_width=True):
            st.switch_page("pages/1_Nuevo_paciente.py")

    status_map = {"Todos": None, "Pendiente": "pending", "Reportado": "reported"}
    patients = list_patients(search=search.strip(), status=status_map[status_filter])

    if not patients:
        st.info("No hay pacientes que coincidan con el criterio.")
        return

    st.caption(f"{len(patients)} paciente(s) listado(s). Seleccioná una fila y presioná **Abrir detalle**.")
    df = build_table(patients, variants)

    event = st.dataframe(
        df,
        use_container_width=True,
        hide_index=True,
        on_select="rerun",
        selection_mode="single-row",
    )

    if event.selection.rows:
        idx = event.selection.rows[0]
        pid = df.iloc[idx]["ID"]
        if st.button(f"📋 Abrir detalle de {pid}", type="primary", use_container_width=True):
            st.session_state["viewing_patient_id"] = pid
            st.switch_page("pages/3_Detalle_paciente.py")


if __name__ == "__main__":
    main()
