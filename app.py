"""DPYD Sanger Variant Caller — Dashboard (home)."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import streamlit as st

from dpyd_caller.auth import require_auth, logout_button
from dpyd_caller.references import references_status
from db import (
    dashboard_stats,
    init_db,
    recent_patients,
    truncate_all_data,
    variant_frequencies,
)

VARIANTS_PATH = Path(__file__).parent / "data" / "variants.json"

CALL_EMOJI = {
    "hom_ref": "🟢", "het": "🟡", "hom_alt": "🔴",
    "unclear": "⚪", "discordant": "⚠️", "failed": "❌",
    "missing": "⏸️",
}
CALL_LABELS_SHORT = {
    "hom_ref": "WT/WT", "het": "WT/var", "hom_alt": "var/var",
    "unclear": "ambig", "discordant": "disc", "failed": "fail",
    "missing": "no test",
}


@st.cache_data(ttl=30)
def load_variants() -> dict:
    return json.loads(VARIANTS_PATH.read_text())


def render_sidebar(variants: dict):
    with st.sidebar:
        st.markdown("## 🧬 DPYD Sanger")
        st.divider()
        st.markdown("### Referencias")
        status = references_status(list(variants.keys()))
        ok = sum(1 for p in status.values() if p["fwd"] and p["rev"])
        st.metric("Refs completas", f"{ok}/{len(variants)}")
        for vid, paths in status.items():
            fwd_ok = "✅" if paths["fwd"] else "❌"
            rev_ok = "✅" if paths["rev"] else "❌"
            st.caption(f"{variants[vid]['name']}: fwd {fwd_ok} · rev {rev_ok}")
        st.divider()
        logout_button()


def render_header():
    st.title("📊 Dashboard")
    st.caption(
        "DPYD Sanger Variant Caller — panel de control. "
        "Total de pacientes, fenotipos CPIC y frecuencias de polimorfismos."
    )


def render_cards(stats: dict):
    total = stats["total_patients"]
    pending = stats["pending"]
    reported = stats["reported"]
    pheno = stats["phenotypes"]

    row1 = st.columns(3)
    row1[0].metric("Pacientes totales", total)
    row1[1].metric("Pendientes de análisis", pending)
    row1[2].metric("Reportados", reported)

    st.markdown("### Distribución de fenotipos CPIC")
    row2 = st.columns(3)
    row2[0].metric("🟢 Normal Metabolizer", pheno.get("NM", 0))
    row2[1].metric("🟡 Intermediate Metabolizer", pheno.get("IM", 0))
    row2[2].metric("🔴 Poor Metabolizer", pheno.get("PM", 0))


def render_variant_frequencies(variants: dict):
    st.markdown("### Frecuencia de polimorfismos en la base")
    freqs = variant_frequencies(list(variants.keys()))
    if not any(f["total_confident"] for f in freqs):
        st.info("Todavía no hay análisis cargados. Creá un paciente y corré el panel.")
        return
    rows = []
    for f in freqs:
        v = variants[f["variant_id"]]
        rows.append({
            "Variante": v["name"],
            "HGVS": v["hgvs_coding"],
            "rsID": v["rsid"],
            "🟢 hom_ref": f["hom_ref"],
            "🟡 het": f["het"],
            "🔴 hom_alt": f["hom_alt"],
            "⚠️ inválidos": f["unclear"] + f["discordant"] + f["failed"],
            "n (válidos)": f["total_confident"],
            "Portadores": f"{f['carrier_count']} ({f['carrier_pct']}%)" if f["carrier_pct"] is not None else "—",
            "MAF": f"{f['maf_pct']}%" if f["maf_pct"] is not None else "—",
        })
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)


def render_recent_patients():
    st.markdown("### Últimos 10 pacientes ingresados")
    recent = recent_patients(10)
    if not recent:
        st.info("Sin pacientes ingresados aún.")
        return
    rows = []
    for r in recent:
        status_emoji = "🟡 Pendiente" if r["status"] == "pending" else "✅ Reportado"
        rows.append({
            "ID": r["patient_id"],
            "Nombre": r["full_name"],
            "Status": status_emoji,
            "Ingreso": r["entry_date"],
            "Reporte": r["report_date"] or "—",
            "Fenotipo": r["phenotype"] or "—",
            "AS": r["activity_score"] if r["activity_score"] is not None else "—",
        })
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)


def main():
    st.set_page_config(page_title="DPYD Dashboard", layout="wide", page_icon="🧬")
    require_auth()
    init_db()
    variants = load_variants()

    render_sidebar(variants)
    render_header()

    stats = dashboard_stats()
    render_cards(stats)
    st.divider()
    render_variant_frequencies(variants)
    st.divider()
    render_recent_patients()
    st.divider()
    render_admin_panel(stats)


def render_admin_panel(stats: dict):
    with st.expander("⚙️ Administración"):
        st.caption(
            "Acciones críticas. El reset borra **todos los pacientes, análisis, reportes CPIC "
            "y archivos `.fsa`** subidos. Las referencias y los amplicones **no** se tocan."
        )
        st.markdown(f"**Estado actual**: {stats['total_patients']} paciente(s) — "
                    f"{stats['pending']} pendientes, {stats['reported']} reportados")
        confirm1 = st.checkbox(
            "Confirmo que quiero eliminar TODOS los datos de pacientes",
            key="admin_reset_confirm_1",
        )
        confirm2 = False
        if confirm1:
            confirm2 = st.checkbox(
                "Entiendo que esta acción es **irreversible**",
                key="admin_reset_confirm_2",
            )
        if confirm1 and confirm2:
            if st.button("🗑️ Eliminar todos los datos ahora", type="secondary"):
                counts = truncate_all_data()
                st.session_state.pop("admin_reset_confirm_1", None)
                st.session_state.pop("admin_reset_confirm_2", None)
                st.success(
                    f"✅ Base reseteada. Eliminados: {counts['patients']} pacientes, "
                    f"{counts['variant_calls']} llamados de variante, "
                    f"{counts['cpic_reports']} reportes CPIC, "
                    f"{counts['patient_files']} archivos `.fsa`."
                )
                st.rerun()


if __name__ == "__main__":
    main()
