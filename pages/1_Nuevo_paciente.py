"""Alta de paciente + (opcional) análisis inmediato."""
from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import streamlit as st

from dpyd_caller.auth import require_auth, logout_button
from dpyd_caller.cpic import (
    calculate_activity_score,
    dosing_for,
    load_cpic_table,
    phenotype_from_activity_score,
)
from dpyd_caller.parser import parse_fsa
from dpyd_caller.patient_files import save_patient_file
from dpyd_caller.references import references_status
from dpyd_caller.variant_caller import analyze_variant
from db import (
    get_patient,
    init_db,
    mark_reported,
    save_cpic_report,
    save_variant_call,
    upsert_patient,
)

VARIANTS_PATH = Path(__file__).parent.parent / "data" / "variants.json"


@st.cache_data(ttl=30)
def load_variants() -> dict:
    return json.loads(VARIANTS_PATH.read_text())


def patient_form():
    st.subheader("Datos demográficos")
    col1, col2 = st.columns(2)
    with col1:
        patient_id = st.text_input("ID paciente *", placeholder="PAC-2026-0001")
        full_name = st.text_input("Nombre completo *")
        sex = st.selectbox("Sexo *", ["F", "M", "Otro"])
        birth_date = st.date_input(
            "Fecha de nacimiento *",
            value=None,
            min_value=date(1900, 1, 1),
            max_value=date.today(),
        )
        dni = st.text_input("DNI")
        sample_date = st.date_input(
            "Fecha de muestra *",
            value=date.today(),
            max_value=date.today(),
        )
    with col2:
        diagnosis = st.text_input("Diagnóstico oncológico *", placeholder="Ca colorrectal estadío III")
        tumor_stage = st.text_input("Estadío tumoral")
        referring_physician = st.text_input("Médico solicitante")
        referring_institution = st.text_input("Institución solicitante")
        notes = st.text_area("Notas", height=150)

    return {
        "patient_id": patient_id.strip() if patient_id else "",
        "full_name": full_name.strip() if full_name else "",
        "sex": sex,
        "birth_date": birth_date.isoformat() if birth_date else "",
        "dni": dni.strip() if dni else None,
        "diagnosis": diagnosis.strip() if diagnosis else "",
        "tumor_stage": tumor_stage.strip() if tumor_stage else None,
        "referring_physician": referring_physician.strip() if referring_physician else None,
        "referring_institution": referring_institution.strip() if referring_institution else None,
        "sample_date": sample_date.isoformat() if sample_date else "",
        "notes": notes.strip() if notes else None,
    }


def validate(data: dict) -> list[str]:
    errors = []
    required = ["patient_id", "full_name", "birth_date", "diagnosis", "sample_date"]
    for field in required:
        if not data.get(field):
            errors.append(f"Campo obligatorio: {field}")
    return errors


def analysis_uploader(variants: dict):
    st.subheader("Archivos del paciente (8 en total)")
    st.caption("Subí fwd + rev para cada una de las 4 variantes. Las referencias se toman de `data/references/`.")
    uploads = {}
    for vid, v in variants.items():
        with st.expander(f"{v['name']} — `{v['hgvs_coding']}`", expanded=False):
            c1, c2 = st.columns(2)
            uploads[(vid, "fwd")] = c1.file_uploader(f"Forward — {v['name']}", key=f"{vid}_fwd", type=["fsa", "ab1"])
            uploads[(vid, "rev")] = c2.file_uploader(f"Reverse — {v['name']}", key=f"{vid}_rev", type=["fsa", "ab1"])
    return uploads


def run_panel_analysis(patient_id: str, variants: dict, uploads: dict, ref_status: dict):
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
    progress.progress(1.0, text="Calculando CPIC...")

    cpic_result = calculate_activity_score(variant_calls_for_cpic)
    phenotype = phenotype_from_activity_score(cpic_result["score"]) if cpic_result["reliable"] else None
    dosing = dosing_for(cpic_result["score"], phenotype["code"] if phenotype else None) if phenotype else None
    save_cpic_report(patient_id, cpic_result, phenotype or {}, dosing or {})
    mark_reported(patient_id)


def main():
    st.set_page_config(page_title="Nuevo paciente", layout="wide", page_icon="➕")
    require_auth()
    init_db()
    variants = load_variants()

    with st.sidebar:
        st.markdown("## 🧬 DPYD Sanger")
        st.divider()
        logout_button()

    st.title("➕ Nuevo paciente")

    data = patient_form()
    st.divider()

    analyze_now = st.checkbox("Analizar secuencias ahora (si ya tenés los 8 archivos del paciente)")

    uploads = {}
    ref_status = references_status(list(variants.keys()))
    missing_refs = [vid for vid, p in ref_status.items() if not (p["fwd"] and p["rev"])]

    if analyze_now:
        if missing_refs:
            st.error(
                f"Faltan referencias para: {', '.join(missing_refs)}. "
                f"No se puede analizar hasta completarlas en `data/references/`."
            )
            analyze_now = False
        else:
            uploads = analysis_uploader(variants)

    all_uploaded = analyze_now and all(u is not None for u in uploads.values()) if uploads else False

    st.divider()
    submit_label = "Registrar + analizar" if (analyze_now and all_uploaded) else "Registrar como pendiente"
    submitted = st.button(submit_label, type="primary")

    if not submitted:
        if analyze_now and not all_uploaded:
            missing_count = sum(1 for u in uploads.values() if u is None)
            st.caption(f"{8 - missing_count}/8 archivos subidos")
        return

    errors = validate(data)
    if errors:
        for e in errors:
            st.error(e)
        return

    if get_patient(data["patient_id"]):
        st.error(f"Ya existe un paciente con ID `{data['patient_id']}`.")
        return

    upsert_patient(data)
    st.success(f"✅ Paciente `{data['patient_id']}` registrado.")

    if analyze_now and all_uploaded:
        with st.spinner("Analizando panel DPYD..."):
            run_panel_analysis(data["patient_id"], variants, uploads, ref_status)
        st.session_state["_flash_detail"] = (
            f"✅ Paciente `{data['patient_id']}` registrado y análisis completado."
        )
    else:
        st.session_state["_flash_detail"] = (
            f"🟡 Paciente `{data['patient_id']}` registrado como **pendiente**. "
            f"Corré el análisis cuando tengas los `.fsa`."
        )

    st.session_state["viewing_patient_id"] = data["patient_id"]
    st.switch_page("pages/3_Detalle_paciente.py")


if __name__ == "__main__":
    main()
