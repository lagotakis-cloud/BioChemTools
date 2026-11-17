import io
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, QED
import streamlit as st

st.set_page_config(
    page_title="QED Drug-Likeness Evaluator | BiochemTools",
    page_icon="assets/organic.png",
    layout="wide",
)

from components import render_navbar

render_navbar()

st.title("QED Drug-Likeness Evaluator")
st.markdown("---")

# --- INITIALIZE SESSION STATE ---
if "qed_results" not in st.session_state:
  st.session_state.qed_results = None


# --- FUNCTIONS ---
def get_smiles_from_pubchem(compound_name):
  try:
    compounds = pcp.get_compounds(compound_name.strip(), "name")
    if compounds:
      return compounds[0].canonical_smiles
    return None
  except Exception:
    return None


def highlight_qed(row):
  styles = [""] * len(row)
  name_idx = row.index.get_loc("Name / SMILES")
  qed_idx = row.index.get_loc("QED Score")

  if row["QED Score"] >= 0.67:
    color = "background-color: #d4edda; color: #155724; font-weight: bold;"
  elif row["QED Score"] >= 0.35:
    color = "background-color: #fff3cd; color: #856404; font-weight: bold;"
  else:
    color = "background-color: #f8d7da; color: #721c24; font-weight: bold;"

  styles[name_idx] = color
  styles[qed_idx] = color
  return styles


# ==========================================
# 1. DATA INPUT
# ==========================================
if st.session_state.qed_results is None:
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
        "Vancomycin",
        "Omeprazole",
    ]
    st.write("Test List:", ", ".join(input_data_list))

  elif input_option == "Type Compound Names":
    user_input = st.text_area(
        "Type Compound Names (one per line):",
        value="Aspirin\nSertraline\nVancomycin",
    )
    input_data_list = [
        name.strip() for name in user_input.split("\n") if name.strip()
    ]

  elif input_option == "Type SMILES Strings":
    is_smiles_mode = True
    user_input = st.text_area(
        "Type SMILES strings (one per line):",
        value=(
            "CC(=O)OC1=CC=CC=C1C(=O)O\nCC(C)CC1=CC=CC(=C1)C(C)C(=O)O\nCN1C=NC2=C1C(=O)N(C(=O)N2C)C"
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
            f"Successfully uploaded {len(input_data_list)} compound names!"
        )
      except Exception as e:
        st.error(f"Error reading file: {e}")

  else:  # Upload SMILES Strings File
    is_smiles_mode = True
    uploaded_file = st.file_uploader(
        "Upload File (CSV with 'SMILES' column or TXT with one SMILES per line):",
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
            f"Successfully uploaded {len(input_data_list)} SMILES strings!"
        )
      except Exception as e:
        st.error(f"Error reading file: {e}")

  if st.button("Search & Calculate", type="primary"):
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
              f"Compound Name / SMILES: **{display_name}** not found in PubChem."
              " Compound excluded."
          )
          progress_bar.progress((i + 1) / len(input_data_list))
          continue

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
          st.warning(
              f"Invalid SMILES string: **{display_name}**. Compound excluded."
          )
          progress_bar.progress((i + 1) / len(input_data_list))
          continue

        qed_score = round(QED.qed(mol), 3)
        mw = round(Descriptors.MolWt(mol), 2)
        logp = round(Crippen.MolLogP(mol), 2)
        tpsa = round(Descriptors.TPSA(mol), 2)

        results.append({
            "Name / SMILES": display_name,
            "SMILES": smi,
            "Molecular Weight (MW)": mw,
            "LogP": logp,
            "TPSA (Å²)": tpsa,
            "QED Score": qed_score,
        })

        progress_bar.progress((i + 1) / len(input_data_list))

      if results:
        st.session_state.qed_results = results
        st.rerun()

# ==========================================
# 2. RESULTS & DOWNLOAD
# ==========================================
else:
  if st.button("🔄 New Search / Change Data"):
    st.session_state.qed_results = None
    st.rerun()

  df_results = pd.DataFrame(st.session_state.qed_results)

  st.subheader("📊 Analysis Results")

  df_display = df_results.drop(columns=["SMILES"])

  styled_df = df_display.style.apply(highlight_qed, axis=1).format({
      "Molecular Weight (MW)": "{:.2f}",
      "LogP": "{:.2f}",
      "TPSA (Å²)": "{:.2f}",
      "QED Score": "{:.3f}",
  })

  st.dataframe(styled_df, use_container_width=True)

  # --- DOWNLOAD OPTIONS ---
  st.markdown("### 📥 Download Results")
  file_format = st.radio(
      "Choose file format for download:",
      options=["CSV (.csv)", "Text (.txt)", "Excel (.xlsx)"],
      horizontal=True,
      key="download_format_radio_qed",
  )

  if file_format == "CSV (.csv)":
    file_data = df_display.to_csv(index=False).encode("utf-8")
    file_name = "qed_results.csv"
    mime_type = "text/csv"

  elif file_format == "Text (.txt)":
    file_data = df_display.to_csv(sep="\t", index=False).encode("utf-8")
    file_name = "qed_results.txt"
    mime_type = "text/plain"

  elif file_format == "Excel (.xlsx)":
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
      df_display.to_excel(writer, index=False, sheet_name="QED_Score")
    file_data = buffer.getvalue()
    file_name = "qed_results.xlsx"
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
### What is QED (Quantitative Estimate of Drug-Likeness)?

**QED** (Quantitative Estimate of Drug-Likeness) is a widely recognized metric in modern drug discovery and cheminformatics used to evaluate the overall "drug-likeness" or "chemical beauty" of small molecules.

Unlike traditional rule-based filters (such as Lipinski's Rule of Five) which rely on binary thresholds (pass/fail based on strict cutoffs), QED is a continuous numerical score ranging from **0** (least drug-like) to **1** (most drug-like).

The concept was introduced by G. Richard Bickerton and colleagues and is based on the underlying principles of desirability functions. It combines eight fundamental physicochemical descriptors into a unified score:

- **Molecular Weight (MW)**
- **Lipophilicity (ALogP)**
- **Hydrogen Bond Donors (HBD)**
- **Hydrogen Bond Acceptors (HBA)**
- **Topological Polar Surface Area (TPSA)**
- **Rotatable Bonds (ROTB)**
- **Aromatic Rings (AROM)**
- **Structural Alerts (ALERTS)**

### How Does This Tool Work?

The **QED Evaluator** allows researchers to assess one or multiple chemical compounds through three flexible input methods:

**1. Compound Name Search (PubChem)**

Users can enter the commercial or generic name of a chemical compound. The tool queries the PubChem database to automatically retrieve its canonical SMILES string and molecular structure.

**2. Direct SMILES Input**

Users can directly input SMILES strings to calculate and evaluate proprietary or novel chemical structures without requiring database lookups.

**3. File Upload**

Users can upload a `.csv` or `.txt` file containing lists of compound names or SMILES strings, enabling high-throughput in silico screening and batch analysis in a single operation.

### What Does the Tool Calculate?

For each successfully processed compound, the application computes key molecular properties alongside the final QED score.

The results empower researchers to:

- Quantify the overall drug-likeness of candidate molecules on a continuous scale.
- Compare multiple compounds within a unified tabular format.
- Quickly identify high-potential candidates (typically indicated by a QED score ≥ 0.67).
- Export processed reports securely in multiple formats (`.csv`, `.txt`, `.xlsx`).

### Scientific Interpretation and Limitations

While QED provides an excellent statistical estimation of how closely a molecule resembles established oral drugs, it is intended as a guiding metric rather than an absolute rule.

A high QED score suggests favorable physicochemical properties associated with clinical candidate success, but factors such as target-specific binding affinity, metabolic stability, cellular toxicity, and pharmacokinetics must also be evaluated in parallel.

### Scientific Reference

Bickerton, G. R., Paolini, G. V., Besnard, J., Muresan, S., & Hopkins, A. L. (2012). Quantifying the chemical beauty of drugs. *Nature Chemistry, 4*(2), 90–98.
""")
