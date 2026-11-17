import io
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
import streamlit as st
from components import render_navbar

st.set_page_config(
    page_title="Name ⇄ SMILES Converter | BiochemTools",
    page_icon="assets/organic.png",
    layout="wide",
)

render_navbar()

st.title("Name ⇄ SMILES Converter")
st.markdown(
    "Convert chemical names to SMILES structures or resolve SMILES strings"
    " back to compound names via PubChem."
)
st.markdown("---")

# --- INITIALIZE SESSION STATE ---
if "converter_results" not in st.session_state:
    st.session_state.converter_results = None


# --- FUNCTIONS ---
def get_smiles_from_pubchem(compound_name):
    try:
        compounds = pcp.get_compounds(compound_name.strip(), "name")
        if compounds:
            return compounds[0].canonical_smiles
        return None
    except Exception:
        return None


def get_name_from_pubchem(smiles_string):
    try:
        compounds = pcp.get_compounds(smiles_string.strip(), "smiles")
        if compounds:
            # Get IUPAC name or synonym if available
            synonyms = compounds[0].synonyms
            if synonyms:
                return synonyms[0]
        return "Unknown Compound"
    except Exception:
        return "Unknown Compound"


# ==========================================
# 1. DATA INPUT
# ==========================================
if st.session_state.converter_results is None:
    input_option = st.radio(
        "Choose conversion mode:",
        (
            "Type Names to SMILES",
            "Type SMILES to Names",
            "Upload Names.txt to SMILES",
            "Upload SMILES.txt to Names",
        ),
    )

    input_data_list = []
    is_reverse_mode = "to Names" in input_option
    is_upload_mode = "Upload" in input_option

    user_input = ""
    uploaded_file = None

    if not is_upload_mode:
        if not is_reverse_mode:
            user_input = st.text_area(
                "Type Compound Names (one per line):",
                value="Aspirin\nParacetamol\nIbuprofen\nCaffeine",
            )
        else:
            user_input = st.text_area(
                "Type SMILES strings (one per line):",
                value=(
                    "CC(=O)OC1=CC=CC=C1C(=O)O\n"
                    "CC(=O)NC1=CC=C(C=C1)O\n"
                    "CC(C)CC1=CC=CC(=C1)C(C)C(=O)O"
                ),
            )
        input_data_list = [
            line.strip() for line in user_input.split("\n") if line.strip()
        ]
    else:
        uploaded_file = st.file_uploader(
            "Upload text file (one entry per line):", type=["txt", "csv"]
        )
        if uploaded_file is not None:
            try:
                stringdata = uploaded_file.getvalue().decode("utf-8")
                input_data_list = [
                    line.strip()
                    for line in stringdata.splitlines()
                    if line.strip()
                ]
                st.success(
                    f"Successfully loaded {len(input_data_list)} entries from"
                    " file!"
                )
            except Exception as e:
                st.error(f"Error reading file: {e}")

    if st.button("Execute Conversion", type="primary"):
        if not input_data_list:
            st.warning("Please enter or upload data before executing.")
        else:
            results = []
            progress_bar = st.progress(0)

            for i, item in enumerate(input_data_list):
                if not is_reverse_mode:
                    name = item
                    smi = get_smiles_from_pubchem(item)
                    if smi is None:
                        smi = "Not Found in PubChem"
                else:
                    smi = item
                    # Validate SMILES via RDKit first
                    mol = Chem.MolFromSmiles(smi)
                    if mol is None:
                        name = "Invalid SMILES"
                    else:
                        name = get_name_from_pubchem(smi)

                results.append({
                    "Compound Name": name,
                    "SMILES String": smi,
                })

                progress_bar.progress((i + 1) / len(input_data_list))

            if results:
                st.session_state.converter_results = results
                st.rerun()

# ==========================================
# 2. RESULTS & DOWNLOAD
# ==========================================
else:
    if st.button("🔄 New Conversion / Change Data"):
        st.session_state.converter_results = None
        st.rerun()

    df_results = pd.DataFrame(st.session_state.converter_results)

    st.subheader("📊 Conversion Results")
    st.dataframe(df_results, use_container_width=True)

    # --- DOWNLOAD OPTIONS ---
    st.markdown("### 📥 Download Results")
    file_format = st.radio(
        "Choose file format for download:",
        options=["CSV (.csv)", "Text (.txt)", "Excel (.xlsx)"],
        horizontal=True,
        key="download_format_radio_converter",
    )

    if file_format == "CSV (.csv)":
        file_data = df_results.to_csv(index=False).encode("utf-8")
        file_name = "name_smiles_conversion.csv"
        mime_type = "text/csv"

    elif file_format == "Text (.txt)":
        file_data = df_results.to_csv(sep="\t", index=False).encode("utf-8")
        file_name = "name_smiles_conversion.txt"
        mime_type = "text/plain"

    elif file_format == "Excel (.xlsx)":
        buffer = io.BytesIO()
        with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
            df_results.to_excel(writer, index=False, sheet_name="Conversion")
        file_data = buffer.getvalue()
        file_name = "name_smiles_conversion.xlsx"
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
### What is Name ⇄ SMILES Conversion?

In cheminformatics and drug discovery, chemical structures are represented in standardized machine-readable formats—most notably **SMILES (Simplified Molecular Input Line Entry System)**—whereas literature and databases frequently reference compounds by common names, trade names, or IUPAC nomenclature.

Bidirectional translation between chemical names and SMILES strings is essential for data cleaning, library merging, and automated screening workflows.

### How Does This Tool Work?

The **Name ⇄ SMILES Converter** is a high-throughput batch translation utility powered by PubChemPy and RDKit integration, supporting four flexible conversion modes:

1. **Type Names to SMILES:** Resolves typed chemical names into canonical SMILES strings via PubChem.
2. **Type SMILES to Names:** Resolves input SMILES strings into standardized compound names or synonyms.
3. **Upload Names.txt to SMILES:** Batch processes text files containing compound names.
4. **Upload SMILES.txt to Names:** Batch processes text files containing SMILES strings.

### What Does the Tool Calculate and Resolve?

For each entry in the input set, the tool:

- Queries the PubChem database for robust name resolution and synonym matching.
- Validates chemical structural integrity locally using RDKit parsers.
- Outputs clean, tabular results containing paired compound names and SMILES strings without visual clutter, ready for immediate export in `.csv`, `.txt`, or `.xlsx` formats.

### Limitations

Name resolution relies on external database indexing (PubChem). Ambiguous trade names, proprietary identifiers, or highly novel synthesized compounds without PubChem entries may require direct SMILES input.

### Scientific Reference

Kim, S., Thiessen, P. A., Bolton, E. E., Chen, J., Fu, G., Gindulyte, A., ... & Wang, Y. (2016). PubChem Substance and Compound databases. *Nucleic Acids Research, 44*(D1), D1202–D1213.
""")
