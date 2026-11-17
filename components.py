import streamlit as st


def render_navbar():
    """Renders the sticky navigation bar at the top of the app."""

    st.markdown(
        """
    <style>
        /* Make the first two rows sticky at the top */
        div.stMainBlockContainer > div.stVerticalBlock > div.stHorizontalBlock:nth-of-type(1),
        div.stMainBlockContainer > div.stVerticalBlock > div.stHorizontalBlock:nth-of-type(2) {
            position: sticky !important;
            top: 0px !important;
            z-index: 99999 !important;
            background-color: var(--background-color, #090D16) !important;
            padding-top: 8px !important;
            padding-bottom: 8px !important;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.2);
        }
        
        /* The 2nd row sits right below the 1st */
        div.stMainBlockContainer > div.stVerticalBlock > div.stHorizontalBlock:nth-of-type(2) {
            top: 48px !important; 
        }
    </style>
    """,
        unsafe_allow_html=True,
    )

    # --- Row 1 (Home + Tools 1-3) ---
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.page_link("app.py", label="🏠 Home Page")
    with col2:
        st.page_link("pages/1_lipinski.py", label="💊 Lipinski Ro5")
    with col3:
        st.page_link("pages/2_veber.py", label="📐 Veber Rules")
    with col4:
        st.page_link("pages/3_admet.py", label="🛡️ ADMET & Toxicity")

    # --- Row 2 (Tools 4-6) ---
    col5, col6, col7, col8 = st.columns(4)
    with col5:
        st.page_link("pages/4_ghose.py", label="🧪 Ghose Filter")
    with col6:
        st.page_link("pages/5_qed.py", label="📈 QED Evaluator")
    with col7:
        st.page_link("pages/6_names_to_smiles.py", label="🔤 Name ⇄ SMILES")
    with col8:
        st.empty()

    # Small spacer after navbar before page content begins
    st.markdown("<div style='margin-bottom: 25px;'></div>", unsafe_allow_html=True)


def render_footer():
    """Renders a clean, professional footer with contact info across all pages."""
    st.markdown(
        """
        <hr style="margin-top: 50px; margin-bottom: 20px; border: 0; border-top: 1px solid rgba(128, 128, 128, 0.2);">
        <div style="text-align: center; font-size: 14px; padding-bottom: 30px;">
            <p>🧬 <b>biochemtools.tech</b> | Free Drug Molecules Screening Suite</p>
            <p>📧 Contact & Support: <a href="mailto:info@biochemtools.tech" style="color: #4CAF50; text-decoration: none;">info@biochemtools.tech</a></p>
            <p style="font-size: 12px; margin-top: 5px; opacity: 0.6;">© 2026 BiochemTools. All rights reserved.</p>
        </div>
        """,
        unsafe_allow_html=True,
    )
