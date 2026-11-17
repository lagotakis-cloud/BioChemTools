import io
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors
import streamlit as st
from components import render_navbar

st.set_page_config(
    page_title="Veber's Rules Checker | BiochemTools",
    page_icon="assets/organic.png",
    layout="wide",
)

render_navbar()

# --- TITLE & DESCRIPTION ---
st.title("🧬 Veber's Rules Checker (TPSA & Rotatable Bonds)")
st.write(
    "Evaluate molecular flexibility and membrane / Blood-Brain Barrier (BBB)"
    " permeability based on Veber's rules."
)
st.markdown("---")

# --- TPSA THRESHOLD SELECTION ---
tpsa_option = st.selectbox(
    "Select permeability target (TPSA threshold):",
    (
        "Cell Membranes (TPSA ≤ 140 Å²)",
        "Blood-Brain Barrier / BBB (TPSA ≤ 90 Å²)",
    ),
)
tpsa_limit = 140 if "140" in tpsa_option else 90

# --- INITIALIZE SESSION STATE ---
if "veber_results" not in st.session_state:
    st.session_state.veber_results = None


# --- FUNCTIONS ---
def get_smiles_from_pubchem(compound_name):
    try:
        compounds = pcp.get_compounds(compound_name.strip(), "name")
        if compounds:
            return compounds[0].canonical_smiles
        else:
            return None
    except Exception:
        return None


def highlight_veber(row):
    styles = [""] * len(row)
    name_idx = row.index.get_loc("Name / SMILES")

    if "Passed" in str(row["Status"]):
        styles[name_idx] = (
            "background-color: #d4edda; color: #155724; font-weight: bold;"
        )
    else:
        styles[name_idx] = (
            "background-color: #f8d7da; color: #721c24; font-weight: bold;"
        )

    if row["TPSA (Å²)"] > tpsa_limit:
        styles[row.index.get_loc("TPSA (Å²)")] = (
            "background-color: #f8d7da; color: #721c24;"
        )
    if row["Rotatable Bonds (RB)"] > 10:
        styles[row.index.get_loc("Rotatable Bonds (RB)")] = (
            "background-color: #f8d7da; color: #721c24;"
        )

    return styles


# ==========================================
# 1. DATA INPUT
# ==========================================
if st.session_state.veber_results is None:
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
        input_data_list = [
            "Aspirin",
            "Cyclosporin A",
            "Ibuprofen",
            "Caffeine",
            "Diazepam",
        ]
        st.write("Test List:", ", ".join(input_data_list))

    elif input_option == "Type Compound Names":
        user_input = st.text_area(
            "Type Compound Names (one per line):",
            value="Aspirin\nDiazepam\nVancomycin",
        )
        input_data_list = [
            name.strip() for name in user_input.split("\n") if name.strip()
        ]

    elif input_option == "Type SMILES Strings":
        is_smiles_mode = True
        user_input = st.text_area(
            "Type SMILES strings (one per line):",
            value=(
                "CC(=O)OC1=CC=CC=C1C(=O)O\nCN1C=NC2=C1C(=O)N(C(=O)N2C)C\nCC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=CC(=C3)OC"
            ),
        )
        input_data_list = [
            smi.strip() for smi in user_input.split("\n") if smi.strip()
        ]

    elif input_option == "Upload Compound Names File":
        uploaded_file = st.file_uploader(
            "Upload File (CSV with 'Name' column or TXT with one name per line):",
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

    if st.button("Search & Calculate Veber", type="primary"):
        if not input_data_list:
            st.warning("Please enter or upload data.")
        else:
            results = []
            progress_bar = st.progress(0)

            for i, item in enumerate(input_data_list):
                if is_smiles_mode:
                    smi = item
                    display_name = item
                else:
                    display_name = item
                    smi = get_smiles_from_pubchem(item)

                if smi is None:
                    st.warning(
                        "Compound Name / SMILES:"
                        f" **{display_name}** not found in PubChem. Compound"
                        " excluded."
                    )
                    progress_bar.progress((i + 1) / len(input_data_list))
                    continue

                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    st.warning(
                        f"Invalid SMILES string: **{display_name}**. Compound"
                        " excluded."
                    )
                    progress_bar.progress((i + 1) / len(input_data_list))
                    continue

                tpsa = Descriptors.TPSA(mol)
                rot_bonds = Descriptors.NumRotatableBonds(mol)

                violations = 0
                if tpsa > tpsa_limit:
                    violations += 1
                if rot_bonds > 10:
                    violations += 1

                status = (
                    "Passed (Veber-like)" if violations == 0 else "Failed"
                )

                results.append({
                    "Name / SMILES": display_name,
                    "SMILES": smi,
                    "TPSA (Å²)": round(tpsa, 2),
                    "Rotatable Bonds (RB)": rot_bonds,
                    "TPSA Limit": tpsa_limit,
                    "Violations": violations,
                    "Status": status,
                })

                progress_bar.progress((i + 1) / len(input_data_list))

            if results:
                st.session_state.veber_results = results
                st.rerun()


# ==========================================
# 2. RESULTS & DOWNLOAD
# ==========================================
else:
    if st.button("🔄 New Search / Change Data"):
        st.session_state.veber_results = None
        st.rerun()

    df_results = pd.DataFrame(st.session_state.veber_results)

    st.subheader("📊 Veber Analysis Results")
    st.info(
        "Active TPSA criterion: "
        f"{'≤ 140 Å² (Cell Membranes)' if tpsa_limit == 140 else '≤ 90 Å² (BBB)'}"
    )

    df_display = df_results.drop(columns=["SMILES", "TPSA Limit"])

    styled_df = df_display.style.apply(highlight_veber, axis=1).format(
        {"TPSA (Å²)": "{:.2f}"}
    )

    st.dataframe(styled_df, use_container_width=True)

    # --- DOWNLOAD OPTIONS ---
    st.markdown("### 📥 Download Veber Results")
    file_format = st.radio(
        "Choose file format for download:",
        options=["CSV (.csv)", "Text (.txt)", "Excel (.xlsx)"],
        horizontal=True,
        key="download_format_radio_veber",
    )

    if file_format == "CSV (.csv)":
        file_data = df_display.to_csv(index=False).encode("utf-8")
        file_name = "veber_results.csv"
        mime_type = "text/csv"

    elif file_format == "Text (.txt)":
        file_data = df_display.to_csv(sep="\t", index=False).encode("utf-8")
        file_name = "veber_results.txt"
        mime_type = "text/plain"

    elif file_format == "Excel (.xlsx)":
        buffer = io.BytesIO()
        with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
            df_display.to_excel(writer, index=False, sheet_name="Veber_Rules")
        file_data = buffer.getvalue()
        file_name = "veber_results.xlsx"
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
### What are Veber's Rules?

Veber's rules are key guidelines formulated by Valerie J. Veber and colleagues at GlaxoSmithKline (2002) to evaluate oral bioavailability and membrane permeability in drug discovery. 

While Lipinski's Rule of Five focuses primarily on molecular weight, lipophilicity, and hydrogen bonding, Veber's work demonstrated that molecular flexibility and polar surface area are exceptionally powerful predictors of oral bioavailability in rats.

According to Veber's criteria, compounds are more likely to exhibit good oral bioavailability if they satisfy:

- **Rotatable Bonds (NRotB):** $\le 10$ (measuring molecular flexibility)
- **Topological Polar Surface Area (TPSA):** $\le 140 \text{ \AA}^2$ (measuring polar atoms contribution, often $\le 60 \text{ \AA}^2$ for Blood-Brain Barrier penetration)

### How Does This Tool Work?

The **Veber's Rules & TPSA Screener** allows users to evaluate chemical compounds through multiple input methods:

**1. Compound Name Search (PubChem)**
Users enter compound names, and the tool queries the PubChem database to retrieve their canonical SMILES strings.

**2. Direct SMILES Input**
Users can manually input SMILES representations for instant, local analysis.

**3. File Upload**
Users can upload `.csv` or `.txt` files containing compound names or SMILES strings for high-throughput batch screening.

### What Does the Tool Calculate?

For each successfully processed compound, the tool evaluates:

- **Rotatable Bond Count:** Counts single bonds that allow free rotation (excluding terminal bonds).
- **Topological Polar Surface Area (TPSA):** Computes the surface sum over all polar atoms (primarily oxygen and nitrogen, including attached hydrogens).
- **Blood-Brain Barrier (BBB) Screening:** Evaluates whether the TPSA and lipophilicity profiles suggest potential penetration into the central nervous system.
- **Compliance Status:** Determines whether the compound satisfies Veber's criteria for oral bioavailability.

### Scientific Interpretation and Limitations

Veber's rules serve as an invaluable early-stage filter for molecular optimization. However, they represent statistical guidelines rather than absolute biological laws. Active transporter mechanisms, protein binding, and metabolic stability can allow certain molecules to bypass traditional topological constraints.

### Scientific Reference

Veber, D. F., Johnson, S. R., Cheng, H. Y., Smith, B. R., Ward, K. W., & Kopple, K. D. (2002). Molecular properties that influence oral bioavailability and drug likeness of candidates. *Journal of Medicinal Chemistry, 45*(12), 2615–2623.
""")


