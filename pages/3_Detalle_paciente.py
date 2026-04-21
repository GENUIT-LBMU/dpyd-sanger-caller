"""Detalle individual de paciente: datos, reporte CPIC, análisis por variante, visor de electroferogramas, acción re-analizar."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from dpyd_caller.auth import require_auth, logout_button
from dpyd_caller.cpic import (
    calculate_activity_score,
    dosing_for,
    load_cpic_table,
    phenotype_from_activity_score,
)
from dpyd_caller.parser import SangerRead, parse_fsa
from dpyd_caller.patient_files import (
    find_patient_file,
    patient_files_status,
    save_patient_file,
)
from dpyd_caller.references import find_reference, references_status
from dpyd_caller.variant_caller import analyze_variant
from db import (
    delete_patient,
    get_cpic_report,
    get_patient,
    get_variant_calls,
    init_db,
    mark_pending,
    mark_reported,
    save_cpic_report,
    save_variant_call,
)

VARIANTS_PATH = Path(__file__).parent.parent / "data" / "variants.json"

CALL_EMOJI = {
    "hom_ref": "🟢", "het": "🟡", "hom_alt": "🔴",
    "unclear": "⚪", "discordant": "⚠️", "failed": "❌",
}
CALL_LABELS = {
    "hom_ref": "Homo wild-type",
    "het": "Heterocigota",
    "hom_alt": "Homo variante",
    "unclear": "Ambigua",
    "discordant": "Discordante fwd/rev",
    "failed": "Falló",
}
CHANNEL_COLORS = {"A": "#2ca02c", "C": "#1f77b4", "G": "#111111", "T": "#d62728"}


@st.cache_data(ttl=30)
def load_variants() -> dict:
    return json.loads(VARIANTS_PATH.read_text())


@st.cache_data(show_spinner=False)
def parse_fsa_cached(path_str: str, mtime: float) -> SangerRead:
    """Cache keyed by (path, mtime) so re-analyses bust the cache automatically."""
    return parse_fsa(path_str)


def read_from_disk(path: Path | None) -> SangerRead | None:
    if path is None or not path.exists():
        return None
    return parse_fsa_cached(str(path), path.stat().st_mtime)


def plot_variant_region(
    read: SangerRead,
    peak_location: int,
    ref_channel: str,
    alt_channel: str,
    window_bases: int = 10,
    title: str = "",
):
    half_win = window_bases * 12
    max_len = len(next(iter(read.traces.values())))
    start = max(0, peak_location - half_win)
    end = min(max_len, peak_location + half_win)
    fig, ax = plt.subplots(figsize=(10, 2.8))
    x = list(range(start, end))
    for base, color in CHANNEL_COLORS.items():
        if base in read.traces:
            is_focus = base in (ref_channel, alt_channel)
            ax.plot(
                x,
                read.traces[base][start:end],
                color=color,
                label=f"{base}{' ◀' if is_focus else ''}",
                linewidth=2.0 if is_focus else 1.0,
                alpha=1.0 if is_focus else 0.45,
            )
    ax.axvspan(peak_location - 3, peak_location + 3, color="#f1c40f", alpha=0.25, label="_nolegend_")
    ax.axvline(peak_location, color="#8e44ad", linestyle="--", alpha=0.9, linewidth=1.5)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Posición en traza")
    ax.set_ylabel("Intensidad")
    ax.legend(loc="upper right", fontsize=8, ncols=4)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    return fig


def render_demographics(patient: dict):
    st.subheader("Datos demográficos")
    c1, c2, c3 = st.columns(3)
    c1.markdown(f"**Nombre**  \n{patient['full_name']}")
    c1.markdown(f"**Sexo**  \n{patient['sex']}")
    c1.markdown(f"**DNI**  \n{patient.get('dni') or '—'}")
    c2.markdown(f"**Fecha de nacimiento**  \n{patient['birth_date']}")
    c2.markdown(f"**Fecha de muestra**  \n{patient['sample_date']}")
    c2.markdown(f"**Diagnóstico**  \n{patient['diagnosis']}")
    c3.markdown(f"**Estadío**  \n{patient.get('tumor_stage') or '—'}")
    c3.markdown(f"**Médico solicitante**  \n{patient.get('referring_physician') or '—'}")
    c3.markdown(f"**Institución**  \n{patient.get('referring_institution') or '—'}")
    if patient.get("notes"):
        st.markdown(f"**Notas**  \n{patient['notes']}")


def render_cpic_card(cpic: dict | None):
    st.subheader("Reporte CPIC")
    if not cpic:
        st.info("Aún no se generó reporte CPIC para este paciente.")
        return
    c1, c2, c3 = st.columns(3)
    c1.metric("Activity Score", cpic["activity_score"] if cpic["activity_score"] is not None else "—")
    c2.metric("Fenotipo", cpic["phenotype"] or "—")
    c3.metric("Alelos (AV)", f"{cpic['allele1_av']} + {cpic['allele2_av']}" if cpic["allele1_av"] is not None else "—")
    if cpic["dosing_recommendation"]:
        st.success(f"**{cpic['dosing_label'] or ''}**")
        st.markdown(f"💊 **Recomendación**: {cpic['dosing_recommendation']}")
        st.caption(f"Fuerza de evidencia CPIC: {cpic.get('dosing_strength') or '—'}")
    if not cpic.get("reliable"):
        st.warning("⚠️ Reporte no confiable — hay genotipos ambiguos/discordantes/fallidos.")
    for note in cpic.get("notes", []):
        st.caption(f"⚠️ {note}")
    st.caption(f"Reportado: {cpic.get('reported_at', '—')}")


def render_variant_calls_table(variants: dict, calls: dict[str, dict]):
    st.subheader("Variantes analizadas")
    if not calls:
        st.info("Aún no hay análisis de variantes para este paciente.")
        return
    rows = []
    for vid, v in variants.items():
        call = calls.get(vid)
        if call is None:
            rows.append({
                "Variante": v["name"],
                "HGVS": v["hgvs_coding"],
                "Genotipo": "—",
                "QC ref": "—",
                "Concord.": "—",
                "Analizado": "—",
            })
            continue
        g = call["genotype"]
        rows.append({
            "Variante": v["name"],
            "HGVS": v["hgvs_coding"],
            "Genotipo": f"{CALL_EMOJI.get(g,'')} {CALL_LABELS.get(g, g)}",
            "QC ref": "✅" if call["ref_qc_pass"] else "❌",
            "Concord.": "✅" if call["sample_concordant"] else "⚠️",
            "Analizado": call.get("analyzed_at", "—"),
        })
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)


def render_electropherogram_viewer(patient_id: str, variants: dict, calls: dict[str, dict]):
    st.subheader("🔬 Visor de electroferogramas")

    analyzed_vids = [vid for vid in variants if vid in calls]
    if not analyzed_vids:
        st.info("No hay análisis aún. Corré el análisis para ver los electroferogramas.")
        return

    files_status = patient_files_status(patient_id, analyzed_vids)

    variant_id = st.selectbox(
        "Seleccioná variante",
        options=analyzed_vids,
        format_func=lambda vid: f"{variants[vid]['name']} — {CALL_EMOJI.get(calls[vid]['genotype'],'')} {CALL_LABELS.get(calls[vid]['genotype'], calls[vid]['genotype'])}",
        key="viewer_variant",
    )

    call = calls[variant_id]
    v = variants[variant_id]
    details = call.get("details", {})
    ref_channel = (details.get("sample_fwd") or {}).get("ref_channel") or v["ref_allele_sense"]
    alt_channel = (details.get("sample_fwd") or {}).get("alt_channel") or v["alt_allele_sense"]

    header_cols = st.columns(5)
    header_cols[0].metric("Genotipo", CALL_LABELS.get(call["genotype"], call["genotype"]))
    header_cols[1].metric("HGVS", v["hgvs_coding"])
    header_cols[2].metric("Ubicación", v["location"])
    header_cols[3].metric("Ref / Alt", f"{v['ref_allele_sense']} → {v['alt_allele_sense']}")
    header_cols[4].metric("Función", v["function"])

    sample_fwd_path = files_status[variant_id]["fwd"]
    sample_rev_path = files_status[variant_id]["rev"]
    ref_fwd_path = find_reference(variant_id, "fwd")
    ref_rev_path = find_reference(variant_id, "rev")

    if not sample_fwd_path or not sample_rev_path:
        st.warning(
            "⚠️ Los archivos `.fsa` del paciente no están en disco. "
            "Esto puede pasar con pacientes analizados antes de esta versión. "
            "Usá **Re-analizar** más abajo para recuperar los electroferogramas."
        )
        return

    window = st.slider(
        "Ventana alrededor del pico (bases)", min_value=5, max_value=40, value=12, step=1,
        help="Cuántas bases mostrar a cada lado del pico variante.",
    )

    tabs = st.tabs(["Muestra forward", "Muestra reverse", "Referencia forward", "Referencia reverse"])
    plot_configs = [
        ("sample_fwd", sample_fwd_path, tabs[0], "Muestra — Forward"),
        ("sample_rev", sample_rev_path, tabs[1], "Muestra — Reverse"),
        ("ref_fwd", ref_fwd_path, tabs[2], "Referencia — Forward"),
        ("ref_rev", ref_rev_path, tabs[3], "Referencia — Reverse"),
    ]

    for read_key, path, tab, label in plot_configs:
        with tab:
            if path is None:
                st.info("Archivo no disponible.")
                continue
            read = read_from_disk(path)
            if read is None:
                st.error(f"No pude parsear {path.name}")
                continue
            d = details.get(read_key)
            if not d:
                st.warning("No hay detalle guardado para esta lectura.")
                continue
            peak_location = d["peak_location"]
            per_read_call = d.get("call", "")
            cols = st.columns(4)
            cols[0].metric("Llamado", CALL_LABELS.get(per_read_call, per_read_call))
            cols[1].metric(f"Pico ref ({d.get('ref_channel')})", d.get("ref_height"))
            cols[2].metric(f"Pico alt ({d.get('alt_channel')})", d.get("alt_height"))
            cols[3].metric("Alt fraction", f"{d.get('alt_frac', 0):.2f}")
            fig = plot_variant_region(
                read,
                peak_location=peak_location,
                ref_channel=d.get("ref_channel") or ref_channel,
                alt_channel=d.get("alt_channel") or alt_channel,
                window_bases=window,
                title=f"{label} — pico en {peak_location}  (call: {per_read_call})",
            )
            st.pyplot(fig)
            plt.close(fig)


def render_analysis_uploader(patient_id: str, variants: dict, mode_label: str):
    st.subheader(mode_label)
    ref_status = references_status(list(variants.keys()))
    missing_refs = [vid for vid, p in ref_status.items() if not (p["fwd"] and p["rev"])]
    if missing_refs:
        st.error(f"Faltan referencias para: {', '.join(missing_refs)}. Completar en `data/references/`.")
        return

    st.caption("Subí fwd + rev de las 4 variantes. Las referencias se toman de `data/references/`.")
    uploads = {}
    for vid, v in variants.items():
        with st.expander(f"{v['name']} — `{v['hgvs_coding']}`", expanded=False):
            c1, c2 = st.columns(2)
            uploads[(vid, "fwd")] = c1.file_uploader(f"Forward — {v['name']}", key=f"reanal_{vid}_fwd", type=["fsa", "ab1"])
            uploads[(vid, "rev")] = c2.file_uploader(f"Reverse — {v['name']}", key=f"reanal_{vid}_rev", type=["fsa", "ab1"])

    all_uploaded = all(u is not None for u in uploads.values())
    missing_count = sum(1 for u in uploads.values() if u is None)
    st.caption(f"{8 - missing_count}/8 archivos subidos")

    if not st.button("Correr análisis", type="primary", disabled=not all_uploaded, key="run_analysis_btn"):
        return

    variant_calls_for_cpic = {}
    progress = st.progress(0.0, text="Iniciando...")
    for i, (vid, v) in enumerate(variants.items()):
        progress.progress((i + 0.1) / len(variants), text=f"Analizando {v['name']}...")
        fwd_upload = uploads[(vid, "fwd")]
        rev_upload = uploads[(vid, "rev")]
        try:
            sample_fwd_path = save_patient_file(
                patient_id, vid, "fwd", fwd_upload.read(), original_name=fwd_upload.name,
            )
            sample_rev_path = save_patient_file(
                patient_id, vid, "rev", rev_upload.read(), original_name=rev_upload.name,
            )
            reads = {
                "ref_fwd": parse_fsa(ref_status[vid]["fwd"]),
                "ref_rev": parse_fsa(ref_status[vid]["rev"]),
                "sample_fwd": parse_fsa(sample_fwd_path),
                "sample_rev": parse_fsa(sample_rev_path),
            }
            analysis = analyze_variant(
                vid, v,
                reads["ref_fwd"], reads["ref_rev"],
                reads["sample_fwd"], reads["sample_rev"],
            )
            save_variant_call(patient_id, vid, analysis)
            variant_calls_for_cpic[vid] = {
                "genotype": analysis.sample_genotype,
                "activity_value": v["activity_value"],
            }
        except Exception as exc:
            st.error(f"Error analizando {vid}: {exc}")
            variant_calls_for_cpic[vid] = {
                "genotype": "failed",
                "activity_value": v["activity_value"],
            }
        progress.progress((i + 1) / len(variants), text=f"{v['name']} listo")

    cpic_result = calculate_activity_score(variant_calls_for_cpic)
    phenotype = phenotype_from_activity_score(cpic_result["score"]) if cpic_result["reliable"] else None
    dosing = dosing_for(cpic_result["score"], phenotype["code"] if phenotype else None) if phenotype else None
    save_cpic_report(patient_id, cpic_result, phenotype or {}, dosing or {})
    mark_reported(patient_id)
    parse_fsa_cached.clear()  # bust cache so viewer loads new traces
    st.session_state["_flash_detail"] = "✅ Análisis actualizado correctamente."
    st.rerun()


def resolve_patient_id() -> str | None:
    pid = st.session_state.get("viewing_patient_id")
    if not pid:
        pid = st.query_params.get("patient_id")
    if pid:
        st.session_state["viewing_patient_id"] = pid
    return pid


def main():
    st.set_page_config(page_title="Detalle paciente", layout="wide", page_icon="📋")
    require_auth()
    init_db()
    variants = load_variants()

    with st.sidebar:
        st.markdown("## 🧬 DPYD Sanger")
        st.divider()
        logout_button()

    st.title("📋 Detalle de paciente")

    flash = st.session_state.pop("_flash_detail", None)
    if flash:
        st.success(flash)

    patient_id = resolve_patient_id()
    if not patient_id:
        st.error("No hay paciente seleccionado. Elegí uno desde la lista.")
        if st.button("🧑‍⚕️ Ir a Pacientes"):
            st.switch_page("pages/2_Pacientes.py")
        return

    patient = get_patient(patient_id)
    if not patient:
        st.error(f"No existe un paciente con ID `{patient_id}`.")
        if st.button("🧑‍⚕️ Ir a Pacientes"):
            st.switch_page("pages/2_Pacientes.py")
        return

    status_badge = "🟡 Pendiente" if patient["status"] == "pending" else "✅ Reportado"
    hdr1, hdr2 = st.columns([4, 1])
    hdr1.markdown(f"### `{patient['patient_id']}` — {patient['full_name']} ({status_badge})")
    with hdr2:
        if patient["status"] == "reported":
            if st.button("↩️ Marcar pendiente"):
                mark_pending(patient_id)
                st.rerun()
        else:
            if st.button("✅ Marcar reportado"):
                mark_reported(patient_id)
                st.rerun()

    calls = get_variant_calls(patient_id)
    cpic = get_cpic_report(patient_id)

    render_demographics(patient)
    st.divider()
    render_cpic_card(cpic)
    st.divider()
    render_variant_calls_table(variants, calls)
    st.divider()
    render_electropherogram_viewer(patient_id, variants, calls)
    st.divider()

    mode = "Re-analizar muestras" if calls else "Analizar muestras (paciente pendiente)"
    render_analysis_uploader(patient_id, variants, mode)

    st.divider()
    with st.expander("⚠️ Acciones administrativas"):
        confirm = st.checkbox(f"Confirmo que quiero eliminar al paciente `{patient_id}` y todos sus análisis")
        if confirm and st.button("🗑️ Eliminar paciente", type="secondary"):
            delete_patient(patient_id)
            st.session_state.pop("viewing_patient_id", None)
            st.success(f"Paciente `{patient_id}` eliminado.")
            st.switch_page("pages/2_Pacientes.py")


if __name__ == "__main__":
    main()
