import streamlit as st
from components import render_navbar
from components import render_footer

st.set_page_config(
    page_title="biochemtools.tech | Enterprise Cheminformatics",
    page_icon="🧬",
    layout="wide",
)

render_navbar()

# --- HERO SECTION ---
st.title("🧬 biochemtools.tech")
st.markdown(
    "### A modern, fast, and free tool suite for drug discovery"
)
st.markdown("---")


# --- GRID NAVIGATION WITH THE 6 TOOLS ---
st.subheader("🛠️ Available Suite Tools")

c1, c2 = st.columns(2)

with c1:
    st.markdown("#### 1. Lipinski's Rule of 5")
    st.markdown("Evaluate drug-likeness and core physicochemical properties.")
    st.page_link("pages/1_lipinski.py", label="Open Lipinski", icon="💊")

    st.markdown("#### 3. Simple ADMET")
    st.markdown(
        "Simple ADMET screening, ESOL, PAINS, CYP450, and toxicity alerts."
    )
    st.page_link("pages/3_admet.py", label="Open ADMET", icon="🛡️")

    st.markdown("#### 5. QED Evaluator")
    st.markdown(
        "Calculate quantitative estimate of drug-likeness (0 to 1 score)."
    )
    st.page_link("pages/5_qed.py", label="Open QED", icon="📈")

with c2:
    st.markdown("#### 2. Veber's Rules")
    st.markdown(
        "Analyze molecular flexibility and polar surface area (TPSA) for"
        " bioavailability."
    )
    st.page_link("pages/2_veber.py", label="Open Veber", icon="📐")

    st.markdown("#### 4. Ghose Filter")
    st.markdown(
        "Filter based on MW, LogP, heavy atoms, and Molar Refractivity."
    )
    st.page_link("pages/4_ghose.py", label="Open Ghose", icon="🧪")

    st.markdown("#### 6. Name ⇄ SMILES")
    st.markdown(
        "Bidirectional batch conversion between compound names (PubChem) and"
        " SMILES."
    )
    st.page_link(
        "pages/6_names_to_smiles.py", label="Open Name ⇄ SMILES", icon="🔤"
    )

st.markdown("---")

render_footer()
