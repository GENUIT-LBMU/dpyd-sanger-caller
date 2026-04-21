from __future__ import annotations

import os

import streamlit as st


def _get_expected_password() -> str | None:
    try:
        if "app_password" in st.secrets:
            return st.secrets["app_password"]
    except Exception:
        pass
    return os.environ.get("APP_PASSWORD")


def require_auth() -> None:
    """Block page rendering until the user provides the correct password.

    Reads the expected password from `.streamlit/secrets.toml` (key: `app_password`)
    or from the env var `APP_PASSWORD`. If neither is configured, shows a setup error.
    """
    if st.session_state.get("authenticated"):
        return

    expected = _get_expected_password()
    if not expected:
        st.error(
            "⚠️ No hay contraseña configurada. Definir `app_password` en "
            "`.streamlit/secrets.toml` o la variable de entorno `APP_PASSWORD`."
        )
        st.stop()

    st.title("🔒 DPYD Sanger Variant Caller")
    st.caption("Acceso restringido")
    with st.form("auth_form"):
        pwd = st.text_input("Contraseña", type="password")
        submitted = st.form_submit_button("Ingresar")
    if not submitted:
        st.stop()

    if pwd == expected:
        st.session_state["authenticated"] = True
        st.rerun()
    else:
        st.error("Contraseña incorrecta.")
        st.stop()


def logout_button():
    if st.sidebar.button("Cerrar sesión", use_container_width=True):
        st.session_state["authenticated"] = False
        st.rerun()
