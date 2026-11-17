import io
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
import streamlit as st
from components import render_navbar

st.set_page_config(
    page_title="ADMET & Toxicity Screener | BiochemTools",
    page_icon="assets/organic.png",
    layout="wide",
)

render_navbar()

# --- TITLE & DESCRIPTION ---
st.title("ADMET & Toxicity Screener (Enterprise Edition)")
st.markdown("---")

# --- INITIALIZE SESSION STATE ---
if "admet_results" not in st.session_state:
    st.session_state.admet_results = None


# --- FUNCTIONS ---
def get_smiles_from_pubchem(compound_name):
    try:
        compounds = pcp.get_compounds(compound_name.strip(), "name")
        if compounds:
            return compounds[0].canonical_smiles
        return None
    except Exception:
        return None


def calculate_esol(mol):
    """Calculate aqueous solubility LogS (Delaney equation)."""
    logp = Descriptors.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    ar_atoms = len([a for a in mol.GetAtoms() if a.GetIsAromatic()])
    ar_prop = ar_atoms / heavy_atoms if heavy_atoms > 0 else 0

    return 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rotb - 0.74 * ar_prop


def classify_solubility(log_s):
    """Classify aqueous solubility based on LogS."""
    if log_s > 0:
        return "Highly Soluble"
    elif log_s > -2:
        return "Very Soluble"
    elif log_s > -4:
        return "Soluble"
    elif log_s > -6:
        return "Moderately Soluble"
    elif log_s > -10:
        return "Poorly Soluble"
    else:
        return "Insoluble"


def check_structural_alerts(mol):
    """Detect known toxic / reactive substructures via SMARTS."""
    alerts = {
        "Epoxide (Alkylating)": "[OX1]1C2C1C2",
        "Michael Acceptor": "[C,c]=C-[C,c]=[O,N,S]",
        "Nitroaromatic (Mutagenic)": "c[N+](=O)[O-]",
        "Aromatic Amine": "c[NH2]",
        "Acyl Halide": "[CX3](=[OX1])[F,Cl,Br,I]",
        "Alkyl Halide": "[CX4][Cl,Br,I]",
        "Aldehyde (Reactive)": "[CX3H1](=O)[#6]",
        "Hydrazine": "[NX3][NX3]",
    }
    found_alerts = []
    for name, smarts in alerts.items():
        pat = Chem.MolFromSmarts(smarts)
        if pat and mol.HasSubstructMatch(pat):
            found_alerts.append(name)
    return found_alerts


def check_cyp450_inhibition(mol):
    """Estimate potential hepatic CYP450 enzyme inhibition via structural motifs."""
    cyp_patterns = {
        "CYP3A4 Inhibitor Motif": "c1ccccc1[NX3](C)C",
        "CYP2D6 Inhibitor Motif": "[NX3](C)(C)CCC",
        "CYP2C9 Inhibitor Motif": "c1ccc(cc1)C(=O)O",
    }
    inhibited = []
    for name, smarts in cyp_patterns.items():
        pat = Chem.MolFromSmarts(smarts)
        if pat and mol.HasSubstructMatch(pat):
            inhibited.append(name.split()[0])
    return inhibited if inhibited else ["Low Risk"]


def highlight_admet(row):
    styles = [""] * len(row)
    name_idx = row.index.get_loc("Name / SMILES")
    score_idx = row.index.get_loc("ADMET Score (%)")

    if row["ADMET Score (%)"] >= 80:
        score_style = (
            "background-color: #d4edda; color: #155724; font-weight: bold;"
        )
    elif row["ADMET Score (%)"] >= 50:
        score_style = (
            "background-color: #fff3cd; color: #856404; font-weight: bold;"
        )
    else:
        score_style = (
            "background-color: #f8d7da; color: #721c24; font-weight: bold;"
        )

    styles[name_idx] = score_style
    styles[score_idx] = score_style

    if "Poorly" in str(row["Solubility"]) or "Insoluble" in str(
        row["Solubility"]
    ):
        styles[row.index.get_loc("Solubility")] = (
            "background-color: #fff3cd; color: #856404;"
        )

    if row["PAINS Match"] != "None detected":
        styles[row.index.get_loc("PAINS Match")] = (
            "background-color: #f8d7da; color: #721c24;"
        )

    if row["Reactive Alerts"] != "None detected":
        styles[row.index.get_loc("Reactive Alerts")] = (
            "background-color: #f8d7da; color: #721c24;"
        )

    if row["CYP450 Profile"] != "Low Risk":
        styles[row.index.get_loc("CYP450 Profile")] = (
            "background-color: #fff3cd; color: #856404;"
        )

    return styles


# ==========================================
# 1. DATA INPUT
# ==========================================
if st.session_state.admet_results is None:
    input_option = st.radio(
        "Choose data input:",
        (
            "Test with pre-existing names",
            "Type Compound Names",
            "Type SMILES Strings",
            "Upload Compound Names File",
            "Upload SMILES Strings File",
        ),
    )

    input_data_list = []
    is_smiles_mode = False

    if input_option == "Test with pre-existing names":
        # Χρησιμοποιούμε ασφαλή default ζευγάρια (Όνομα, SMILES) για εγγυημένη εκτέλεση
        input_data_list = [
            ("Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O"),
            ("Paracetamol", "CC(=O)NC1=CC=C(C=C1)O"),
            ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
            ("Ibuprofen", "CC(C)CC1=CC=CC(=C1)C(C)C(=O)O"),
        ]
        st.write("Test List: Aspirin, Paracetamol, Caffeine, Ibuprofen")
        is_preset_tuples = True

    elif input_option == "Type Compound Names":
        is_preset_tuples = False
        user_input = st.text_area(
            "Type Compound Names (one per line):",
            value="Aspirin\nParacetamol\nCaffeine\nIbuprofen",
        )
        input_data_list = [
            name.strip() for name in user_input.split("\n") if name.strip()
        ]

    elif input_option == "Type SMILES Strings":
        is_smiles_mode = True
        is_preset_tuples = False
        user_input = st.text_area(
            "Type SMILES strings (one per line):",
            value=(
                "CC(=O)OC1=CC=CC=C1C(=O)O\n"
                "CC(=O)NC1=CC=C(C=C1)O\n"
                "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
            ),
        )
        input_data_list = [
            smi.strip() for smi in user_input.split("\n") if smi.strip()
        ]

    elif input_option == "Upload Compound Names File":
        is_preset_tuples = False
        uploaded_file = st.file_uploader(
            "Upload File (CSV with 'Name' column or TXT with one name per"
            " line):",
            type=["csv", "txt"],
        )
        if uploaded_file is not None:
            try:
                if uploaded_file.name.endswith(".csv"):
                    df_upload = pd.read_csv(uploaded_file)
                    col_name = (
                        "Name"
                        if "Name" in df_upload.columns
                        else df_upload.columns[0]
                    )
                    input_data_list = (
                        df_upload[col_name].dropna().astype(str).tolist()
                    )
                else:
                    stringdata = uploaded_file.getvalue().decode("utf-8")
                    input_data_list = [
                        line.strip()
                        for line in stringdata.splitlines()
                        if line.strip()
                    ]
                st.success(
                    f"Successfully uploaded {len(input_data_list)} compound"
                    " names!"
                )
            except Exception as e:
                st.error(f"Error reading file: {e}")

    else:  # Upload SMILES Strings File
        is_smiles_mode = True
        is_preset_tuples = False
        uploaded_file = st.file_uploader(
            "Upload File (CSV with 'SMILES' column or TXT with one SMILES per"
            " line):",
            type=["csv", "txt"],
        )
        if uploaded_file is not None:
            try:
                if uploaded_file.name.endswith(".csv"):
                    df_upload = pd.read_csv(uploaded_file)
                    col_name = (
                        "SMILES"
                        if "SMILES" in df_upload.columns
                        else df_upload.columns[0]
                    )
                    input_data_list = (
                        df_upload[col_name].dropna().astype(str).tolist()
                    )
                else:
                    stringdata = uploaded_file.getvalue().decode("utf-8")
                    input_data_list = [
                        line.strip()
                        for line in stringdata.splitlines()
                        if line.strip()
                    ]
                st.success(
                    f"Successfully uploaded {len(input_data_list)} SMILES"
                    " strings!"
                )
            except Exception as e:
                st.error(f"Error reading file: {e}")

    if st.button("Search & Calculate ADMET", type="primary"):
        if not input_data_list:
            st.warning("Please enter or upload data.")
        else:
            results = []
            progress_bar = st.progress(0)

            try:
                params = FilterCatalogParams()
                params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
                catalog = FilterCatalog(params)
            except Exception:
                catalog = None

            for i, item in enumerate(input_data_list):
                if input_option == "Test with pre-existing names":
                    display_name, smi = item
                elif is_smiles_mode:
                    smi = item
                    display_name = item
                else:
                    display_name = item
                    smi = get_smiles_from_pubchem(item)
                    if smi is None:
                        # Fallback to direct name search simulation if PubChem fails
                        smi = get_smiles_from_pubchem(display_name)

                if smi is None:
                    progress_bar.progress((i + 1) / len(input_data_list))
                    continue

                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    progress_bar.progress((i + 1) / len(input_data_list))
                    continue

                esol = calculate_esol(mol)
                sol_class = classify_solubility(esol)

                pains_descriptions = []
                if catalog is not None:
                    try:
                        pains_matches = catalog.GetMatches(mol)
                        pains_descriptions = [
                            m.GetDescription() for m in pains_matches
                        ]
                    except Exception:
                        pass

                alerts = check_structural_alerts(mol)
                cyp_inhibitions = check_cyp450_inhibition(mol)

                score = 100
                if pains_descriptions:
                    score -= 40
                score -= len(alerts) * 25
                if "Insoluble" in sol_class or "Poorly" in sol_class:
                    score -= 15
                score = max(0, score)

                failure_reasons = []
                if pains_descriptions:
                    failure_reasons.append(
                        f"PAINS: {', '.join(pains_descriptions)}"
                    )
                if alerts:
                    failure_reasons.append(f"Alerts: {', '.join(alerts)}")

                toxicity_status = (
                    " | ".join(failure_reasons)
                    if failure_reasons
                    else "Safe (No Alerts)"
                )

                results.append({
                    "Name / SMILES": display_name,
                    "SMILES": smi,
                    "ADMET Score (%)": score,
                    "LogS (ESOL)": round(esol, 2),
                    "Solubility": sol_class,
                    "CYP450 Profile": ", ".join(cyp_inhibitions),
                    "PAINS Match": (
                        ", ".join(pains_descriptions)
                        if pains_descriptions
                        else "None detected"
                    ),
                    "Reactive Alerts": (
                        ", ".join(alerts) if alerts else "None detected"
                    ),
                    "Overall Status": toxicity_status,
                })

                progress_bar.progress((i + 1) / len(input_data_list))

            if results:
                st.session_state.admet_results = results
                st.rerun()
            else:
                st.error(
                    "Could not process compounds. Please try entering SMILES"
                    " directly."
                )


# ==========================================
# 2. RESULTS & DOWNLOAD
# ==========================================
else:
    if st.button("🔄 New Search / Change Data"):
        st.session_state.admet_results = None
        st.rerun()

    df_results = pd.DataFrame(st.session_state.admet_results)

    # --- KPI METRICS ---
    total_mols = len(df_results)
    safe_mols = len(
        df_results[df_results["Overall Status"].str.contains("Safe")]
    )
    alert_mols = total_mols - safe_mols
    avg_score = round(df_results["ADMET Score (%)"].mean(), 1)

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Total Molecules", total_mols)
    col2.metric("Safe (No Alerts)", safe_mols)
    col3.metric("Alerts / PAINS", alert_mols)
    col4.metric("Avg ADMET Score", f"{avg_score}%")

    st.markdown("---")
    st.subheader("📊 Detailed ADMET Results Table")

    df_display = df_results.drop(columns=["SMILES"])

    styled_df = df_display.style.apply(highlight_admet, axis=1).format(
        {"LogS (ESOL)": "{:.2f}"}
    )

    st.dataframe(styled_df, use_container_width=True)

    # --- DOWNLOAD OPTIONS ---
    st.markdown("### 📥 Download Results")
    file_format = st.radio(
        "Choose file format for download:",
        options=["CSV (.csv)", "Text (.txt)", "Excel (.xlsx)"],
        horizontal=True,
        key="download_format_radio_admet",
    )

    if file_format == "CSV (.csv)":
        file_data = df_display.to_csv(index=False).encode("utf-8")
        file_name = "admet_results.csv"
        mime_type = "text/csv"

    elif file_format == "Text (.txt)":
        file_data = df_display.to_csv(sep="\t", index=False).encode("utf-8")
        file_name = "admet_results.txt"
        mime_type = "text/plain"

    elif file_format == "Excel (.xlsx)":
        buffer = io.BytesIO()
        with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
            df_display.to_excel(writer, index=False, sheet_name="ADMET_Results")
        file_data = buffer.getvalue()
        file_name = "admet_results.xlsx"
        mime_type = (
            "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )

    st.download_button(
        label=f"⬇️ Download file ({file_name})",
        data=file_data,
        file_name=file_name,
        mime=mime_type,
        type="primary",
    )

st.markdown("""
### What is ADMET & Toxicity Screening?

Evaluating **ADMET** (Absorption, Distribution, Metabolism, Excretion, and Toxicity) properties early in drug discovery is critical to reducing late-stage clinical attrition.
""")
