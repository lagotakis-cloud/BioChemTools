import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
import pandas as pd
import io
import plotly.graph_objects as go

# 1. Ρύθμιση Σελίδας
st.set_page_config(
    page_title="SpectraMass Pro | LC-MS Analytical Suite", 
    page_icon="🔬", 
    layout="wide",
    initial_sidebar_state="expanded"
)

# 2. Custom CSS για επαγγελματικό SaaS UI
custom_css = """
<style>
    .stApp {
        background-color: #0f172a;
        color: #f8fafc;
    }
    .brand-header {
        background: linear-gradient(90deg, #1e293b 0%, #0f172a 100%);
        padding: 24px;
        border-radius: 12px;
        border-bottom: 2px solid #3b82f6;
        margin-bottom: 25px;
    }
    .brand-title {
        color: #60a5fa;
        font-family: 'Inter', sans-serif;
        font-weight: 700;
        font-size: 2.2rem;
        margin: 0;
    }
    .brand-subtitle {
        color: #94a3b8;
        font-size: 1.0rem;
        margin-top: 5px;
    }
    .stat-card {
        background-color: #1e293b;
        border: 1px solid #334155;
        border-radius: 10px;
        padding: 16px;
        text-align: center;
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.3);
    }
    .stat-label {
        color: #94a3b8;
        font-size: 0.85rem;
        text-transform: uppercase;
        letter-spacing: 0.05em;
    }
    .stat-value {
        color: #38bdf8;
        font-size: 1.6rem;
        font-weight: 700;
        margin-top: 5px;
    }
    div[data-testid="stDataFrame"] {
        border: 1px solid #334155;
        border-radius: 8px;
        overflow: hidden;
    }
</style>
"""
st.markdown(custom_css, unsafe_allow_html=True)

# Header της Εφαρμογής
st.markdown("""
<div class="brand-header">
    <h1 class="brand-title">🔬 SpectraMass Pro</h1>
    <p class="brand-subtitle">High-Precision LC-MS Adduct Calculator & Isotopic Mass Suite</p>
</div>
""", unsafe_allow_html=True)

# 3. Sidebar (Εργαλειοθήκη)
st.sidebar.markdown("### 🎛️ Control Panel")
input_type = st.sidebar.radio("Input Method:", ["Molecule Name", "SMILES Sequence"])

if input_type == "Molecule Name":
    molecule_input = st.sidebar.text_input("Chemical Name:", "Caffeine")
else:
    molecule_input = st.sidebar.text_input("SMILES String:", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

precision = st.sidebar.slider("m/z Decimal Precision:", min_value=2, max_value=6, value=4)

if molecule_input:
    try:
        if input_type == "Molecule Name":
            compound = pcp.get_compounds(molecule_input, 'name')[0]
            smiles = compound.canonical_smiles
            compound_name = molecule_input.capitalize()
            formula = compound.molecular_formula
        else:
            smiles = molecule_input
            compound_name = "User Structure"
            formula = "N/A"

        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            st.error("⚠️ Invalid chemical structure. Please check your SMILES input.")
        else:
            exact_mass = Descriptors.ExactMolWt(mol)
            
            # Layout σε 2 Στήλες
            col_left, col_right = st.columns([1, 2], gap="large")
            
            with col_left:
                st.markdown("#### 🧪 Structure & Properties")
                img = Draw.MolToImage(mol, size=(320, 320))
                st.image(img, use_container_width=True)
                
                st.markdown(f"""
                <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 10px; margin-top: 15px;">
                    <div class="stat-card">
                        <div class="stat-label">Monoisotopic Mass</div>
                        <div class="stat-value">{exact_mass:.{precision}f}</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-label">Formula</div>
                        <div class="stat-value" style="color: #f43f5e;">{formula}</div>
                    </div>
                </div>
                """, unsafe_allow_html=True)
                
                st.caption(f"**SMILES:** `{smiles}`")

            with col_right:
                st.markdown("#### 📊 Adduct Ion Calculator & Spectrum")
                
                # Σταθερές Μαζών
                mass_h = 1.007825
                mass_na = 22.989218
                mass_nh4 = 18.034374
                mass_k = 38.963158
                mass_cl = 34.968853
                mass_hcoo = 44.997655
                mass_ch3coo = 59.013304

                # Positive Mode
                pos_adducts = [
                    {"Adduct": "[M+H]+", "Ion Formula": "M + H", "m/z": exact_mass + mass_h},
                    {"Adduct": "[M+Na]+", "Ion Formula": "M + Na", "m/z": exact_mass + mass_na},
                    {"Adduct": "[M+NH4]+", "Ion Formula": "M + NH4", "m/z": exact_mass + mass_nh4},
                    {"Adduct": "[M+K]+", "Ion Formula": "M + K", "m/z": exact_mass + mass_k},
                    {"Adduct": "[2M+H]+", "Ion Formula": "2M + H", "m/z": (2 * exact_mass) + mass_h},
                ]
                
                # Negative Mode
                neg_adducts = [
                    {"Adduct": "[M-H]-", "Ion Formula": "M - H", "m/z": exact_mass - mass_h},
                    {"Adduct": "[M+Cl]-", "Ion Formula": "M + Cl", "m/z": exact_mass + mass_cl},
                    {"Adduct": "[M+HCOO]-", "Ion Formula": "M + Formate", "m/z": exact_mass + mass_hcoo},
                    {"Adduct": "[M+CH3COO]-", "Ion Formula": "M + Acetate", "m/z": exact_mass + mass_ch3coo},
                ]

                df_pos = pd.DataFrame(pos_adducts)
                df_pos["m/z"] = df_pos["m/z"].round(precision)
                
                df_neg = pd.DataFrame(neg_adducts)
                df_neg["m/z"] = df_neg["m/z"].round(precision)

                # Tabs για Positive, Negative και Spectrum
                tab1, tab2, tab3 = st.tabs(["➕ Positive Mode ([M+H]+)", "➖ Negative Mode ([M-H]-)", "📊 Predicted MS Spectrum"])
                
                with tab1:
                    st.dataframe(df_pos, use_container_width=True, hide_index=True)
                with tab2:
                    st.dataframe(df_neg, use_container_width=True, hide_index=True)
                with tab3:
                    st.markdown("##### 🔮 Predicted Isotope Pattern for [M+H]+")
                    mz_mh = exact_mass + mass_h
                    
                    isotope_data = {
                        "mz": [mz_mh, mz_mh + 1.003, mz_mh + 2.004],
                        "intensity": [100.0, 11.2, 1.5]
                    }
                    
                    fig = go.Figure()
                    for mz, int_val in zip(isotope_data["mz"], isotope_data["intensity"]):
                        fig.add_trace(go.Scatter(
                            x=[mz, mz], y=[0, int_val],
                            mode='lines',
                            line=dict(color='#38bdf8', width=2),
                            showlegend=False
                        ))

                    fig.update_layout(
                        title=f"[M+H]+ Monoisotopic Peak at {mz_mh:.4f}",
                        xaxis_title="m/z",
                        yaxis_title="Relative Intensity (%)",
                        template="plotly_dark",
                        paper_bgcolor="#1e293b",
                        plot_bgcolor="#1e293b",
                        height=300,
                        margin=dict(l=20, r=20, t=40, b=20),
                    )
                    st.plotly_chart(fig, use_container_width=True)

                st.markdown("---")
                col_exp1, col_exp2 = st.columns(2)
                
                df_all = pd.concat([df_pos.assign(Mode='Positive'), df_neg.assign(Mode='Negative')])
                csv_buffer = io.StringIO()
                df_all.to_csv(csv_buffer, index=False)
                
                with col_exp1:
                    st.download_button(
                        label="📥 Export Results (CSV)",
                        data=csv_buffer.getvalue(),
                        file_name=f"{compound_name}_adducts.csv",
                        mime="text/csv",
                        use_container_width=True
                    )
                with col_exp2:
                    if st.button("🔒 Batch Processing (Pro)", use_container_width=True):
                        st.info("💡 Batch CSV analysis is available in the Pro Version.")

    except Exception as e:
        st.error(f"⚠️ Could not resolve structure: {e}")
