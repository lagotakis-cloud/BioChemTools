import io
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors
import streamlit as st
from components import render_navbar

render_navbar()

# --- TITLE & DESCRIPRTION ---
st.set_page_config(
    page_title="Lipinski Rule of 5 | BiochemTools",
    page_icon="assets/organic.png",
    layout="wide",
)
st.title("Lipinski Rule of 5 Checker")

st.markdown("---")

# --- INITIALIZE SESSION STATE ---
if "lipinski_results" not in st.session_state:
  st.session_state.lipinski_results = None


# --- FUNCTIONS ---
def get_smiles_from_pubchem(compound_name):
  try:
    compounds = pcp.get_compounds(compound_name.strip(), "name")
    if compounds:
      return compounds[0].canonical_smiles
    else:
      return None
  except Exception as e:
    return None


def highlight_lipinski(row):
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

  if row["Molecular Weight (MW)"] > 500:
    styles[row.index.get_loc("Molecular Weight (MW)")] = (
        "background-color: #f8d7da; color: #721c24;"
    )
  if row["LogP"] > 5:
    styles[row.index.get_loc("LogP")] = "background-color: #f8d7da; color: #721c24;"
  if row["H Bond Donors (HBD)"] > 5:
    styles[row.index.get_loc("H Bond Donors (HBD)")] = (
        "background-color: #f8d7da; color: #721c24;"
    )
  if row["H Bond Acceptors (HBA)"] > 10:
    styles[row.index.get_loc("H Bond Acceptors (HBA)")] = (
        "background-color: #f8d7da; color: #721c24;"
    )

  return styles


# ==========================================
# 1. DATA INPUT
# ==========================================
if st.session_state.lipinski_results is None:
  input_option = st.radio(
      "Choose data input:",
      (
          "Test with pre-existing names",
          "Type PubChem Names",
          "Type SMILES Strings",
          "Upload PubChem Names File",
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
        "Vancomycin",
        "Omeprazole",
    ]
    st.write("Test List:", ", ".join(input_data_list))

  elif input_option == "Type Names":
    user_input = st.text_area(
        "Type Compound Names (one in each line):",
        value="Aspirin\nSertraline\nVancomycin",
    )
    input_data_list = [
        name.strip() for name in user_input.split("\n") if name.strip()
    ]

  elif input_option == "Type SMILES":
    is_smiles_mode = True
    user_input = st.text_area(
        "Type SMILES strings (one in each line):",
        value=(
            "CC(=O)OC1=CC=CC=C1C(=O)O\nCC(C)CC1=CC=CC(=C1)C(C)C(=O)O\nCN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        ),
    )
    input_data_list = [
        smi.strip() for smi in user_input.split("\n") if smi.strip()
    ]

  elif input_option == "Upload Names File":
    uploaded_file = st.file_uploader(
        "Upload File (CSV with 'Name' column ή TXT with one name in each line):",
        type=["csv", "txt"],
    )
    if uploaded_file is not None:
      try:
        if uploaded_file.name.endswith(".csv"):
          df_upload = pd.read_csv(uploaded_file)
          col_name = (
              "Name" if "Name" in df_upload.columns else df_upload.columns[0]
          )
          input_data_list = (
              df_upload[col_name].dropna().astype(str).tolist()
          )
        else:
          stringdata = uploaded_file.getvalue().decode("utf-8")
          input_data_list = [
              line.strip() for line in stringdata.splitlines() if line.strip()
          ]

        st.success(
            f"Successful {len(input_data_list)} compound names uploaded!"
        )
      except Exception as e:
        st.error(f"Error in uploading: {e}")

  else:  # Ανέβασμα αρχείου SMILES
    is_smiles_mode = True
    uploaded_file = st.file_uploader(
        "Upload File (CSV with 'SMILES' column ή TXT with one SMILES string in each line):",
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
              line.strip() for line in stringdata.splitlines() if line.strip()
          ]

        st.success(
            f"Successful {len(input_data_list)} SMILES strings uploaded!"
        )
      except Exception as e:
        st.error(f"Error in uploading: {e}")

  if st.button("Search & Calculate", type="primary"):
    if not input_data_list:
      st.warning("Please type or upload data.")
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
              f"Compound Name / SMILES: **{display_name}** not found in PubChem."
              "Compound excluded."
          )
          continue

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
          st.warning(
              f"Invalid SMILES string: **{display_name}**. Compound excluded."
          )
          continue

        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)

        violations = 0
        if mw > 500:
          violations += 1
        if logp > 5:
          violations += 1
        if hbd > 5:
          violations += 1
        if hba > 10:
          violations += 1

        status = "Passed (Drug-like)" if violations <= 1 else "Failed"

        results.append({
            "Name/ SMILES": display_name,
            "SMILES": smi,
            "Molecular Weight(MW)": round(mw, 2),
            "LogP": round(logp, 2),
            "H Bond Donors (HBD)": hbd,
            "HBA_val": hba,
            "H Bond Acceptors (HBA)": hba,
            "Violations": violations,
            "Status": status,
        })

        progress_bar.progress((i + 1) / len(input_data_list))

      if results:
        st.session_state.lipinski_results = results
        st.rerun()

# ==========================================
# 2. RESULTS & DOWNLOAD
# ==========================================
else:
  if st.button("🔄 New Search / Change Data"):
    st.session_state.lipinski_results = None
    st.rerun()

  df_results = pd.DataFrame(st.session_state.lipinski_results)

  st.subheader("📊 Analysis Results")

  df_display = df_results.drop(columns=["SMILES", "HBA_val"])

  styled_df = (
      df_display.style.apply(highlight_lipinski, axis=1)
      .format({"Μοριακό Βάρος (MW)": "{:.2f}", "LogP": "{:.2f}"})
  )

  st.dataframe(styled_df, use_container_width=True)

  # --- DOWNLOAD OPTIONS ---
  st.markdown("### 📥 Download Results")
  file_format = st.radio(
      "Choose file type for download:",
      options=["CSV (.csv)", "Text (.txt)", "Excel (.xlsx)"],
      horizontal=True,
      key="download_format_radio",
  )

  if file_format == "CSV (.csv)":
    file_data = df_display.to_csv(index=False).encode("utf-8")
    file_name = "lipinski_results.csv"
    mime_type = "text/csv"

  elif file_format == "Text (.txt)":
    file_data = df_display.to_csv(sep="\t", index=False).encode("utf-8")
    file_name = "lipinski_results.txt"
    mime_type = "text/plain"

  elif file_format == "Excel (.xlsx)":
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
      df_display.to_excel(writer, index=False, sheet_name="Lipinski_Ro5")
    file_data = buffer.getvalue()
    file_name = "lipinski_results.xlsx"
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
### What is Lipinski's Rule of Five?

Lipinski's Rule of Five (Ro5) is a widely used guideline in drug discovery and development for evaluating the drug-likeness of small molecules, particularly in relation to their potential oral bioavailability.

The rule was introduced by Christopher A. Lipinski and colleagues and is based on the observation that compounds with certain physicochemical properties are more likely to exhibit favorable absorption and permeability characteristics.

According to Lipinski's Rule of Five, a compound is more likely to have favorable oral absorption when it satisfies the following criteria:

- **Molecular Weight (MW):** ≤ 500 Da
- **Lipophilicity (LogP):** ≤ 5
- **Hydrogen Bond Donors (HBD):** ≤ 5
- **Hydrogen Bond Acceptors (HBA):** ≤ 10

The rule is commonly used during the early stages of drug discovery to identify compounds with physicochemical properties that may be compatible with oral drug development.

### How Does This Tool Work?

The **Lipinski's Rule of Five Calculator** allows users to evaluate one or multiple chemical compounds through three different input methods:

**1. Compound Name Search (PubChem)**

Users can enter the name of a chemical compound, and the tool performs a search through the PubChem database to retrieve the corresponding molecular information.

**2. Direct SMILES Input**

Users can manually enter the SMILES (Simplified Molecular Input Line Entry System) representation of a compound. This allows direct analysis without requiring a database search.

**3. File Upload**

Users can upload a `.csv` or a `.txt` file containing SMILES strings, enabling the analysis of multiple compounds in a single operation.

### What Does the Tool Calculate?

For each successfully processed compound, the calculator evaluates the relevant physicochemical properties and checks them against Lipinski's Rule of Five criteria.

The results allow users to:

- Review the calculated molecular properties.
- Identify which Lipinski criteria are satisfied or violated.
- Determine the total number of rule violations for each compound.
- Perform an initial assessment of the compound's drug-likeness.

### Scientific Interpretation and Limitations

Lipinski's Rule of Five is a screening guideline rather than a definitive predictor of drug efficacy, safety, or oral bioavailability.

A compound that satisfies all four criteria is not necessarily orally bioavailable, while compounds that violate one or more criteria should not automatically be considered unsuitable for drug development.

The results should therefore be interpreted alongside additional physicochemical properties, ADMET predictions, experimental data, and other relevant drug discovery considerations.

### Scientific Reference

Lipinski, C. A., Lombardo, F., Dominy, B. W., & Feeney, P. J. (2001). Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. *Advanced Drug Delivery Reviews, 46*(1–3), 3–26.
""")
