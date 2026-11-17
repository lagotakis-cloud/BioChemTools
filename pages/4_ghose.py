import io
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors
import streamlit as st

st.set_page_config(
    page_title="Ghose Filter Checker | BiochemTools",
    page_icon="assets/organic.png",
    layout="wide",
)

from components import render_navbar

render_navbar()

st.title("Ghose Filter Checker")
st.markdown("---")

# --- INITIALIZE SESSION STATE ---
if "ghose_results" not in st.session_state:
  st.session_state.ghose_results = None


# --- ΣΥΝΑΡΤΗΣΕΙΣ ---
def get_smiles_from_pubchem(compound_name):
  try:
    compounds = pcp.get_compounds(compound_name.strip(), "name")
    if compounds:
      return compounds[0].canonical_smiles
    return None
  except Exception:
    return None


def highlight_ghose(row):
  styles = [""] * len(row)
  name_idx = row.index.get_loc("Όνομα / SMILES")

  if "Passed" in str(row["Κατάσταση"]):
    styles[name_idx] = (
        "background-color: #d4edda; color: #155724; font-weight: bold;"
    )
  else:
    styles[name_idx] = (
        "background-color: #f8d7da; color: #721c24; font-weight: bold;"
    )

  if not (160 <= row["Μοριακό Βάρος (MW)"] <= 480):
    styles[row.index.get_loc("Μοριακό Βάρος (MW)")] = (
        "background-color: #f8d7da; color: #721c24;"
    )
  if not (-0.4 <= row["LogP"] <= 5.6):
    styles[row.index.get_loc("LogP")] = (
        "background-color: #f8d7da; color: #721c24;"
    )
  if not (20 <= row["Βαριά Άτομα"] <= 70):
    styles[row.index.get_loc("Βαριά Άτομα")] = (
        "background-color: #f8d7da; color: #721c24;"
    )
  if not (40 <= row["Molar Refractivity"] <= 130):
    styles[row.index.get_loc("Molar Refractivity")] = (
        "background-color: #f8d7da; color: #721c24;"
    )

  return styles


# ==========================================
# ΚΑΤΑΣΤΑΣΗ 1: ΕΙΣΑΓΩΓΗ ΔΕΔΟΜΕΝΩΝ
# ==========================================
if st.session_state.ghose_results is None:
  input_option = st.radio(
      "Επιλέξτε τρόπο εισαγωγής μορίων:",
      (
          "Δοκιμή με έτοιμα ονόματα",
          "Πληκτρολόγηση ονομάτων",
          "Πληκτρολόγηση SMILES",
          "Ανέβασμα αρχείου ονομάτων",
          "Ανέβασμα αρχείου SMILES",
      ),
  )

  input_data_list = []
  is_smiles_mode = False

  if input_option == "Δοκιμή με έτοιμα ονόματα":
    input_data_list = [
        "Aspirin",
        "Cyclosporin A",
        "Ibuprofen",
        "Atorvastatin",
        "Omeprazole",
    ]
    st.write("Δοκιμαστική λίστα:", ", ".join(input_data_list))

  elif input_option == "Πληκτρολόγηση ονομάτων":
    user_input = st.text_area(
        "Γράψτε τα ονόματα (ένα ανά γραμμή):",
        value="Aspirin\nSertraline\nAtorvastatin",
    )
    input_data_list = [
        name.strip() for name in user_input.split("\n") if name.strip()
    ]

  elif input_option == "Πληκτρολόγηση SMILES":
    is_smiles_mode = True
    user_input = st.text_area(
        "Γράψτε τα SMILES strings (ένα ανά γραμμή):",
        value=(
            "CC(=O)OC1=CC=CC=C1C(=O)O\nCC(C)CC1=CC=CC(=C1)C(C)C(=O)O\nCN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        ),
    )
    input_data_list = [
        smi.strip() for smi in user_input.split("\n") if smi.strip()
    ]

  elif input_option == "Ανέβασμα αρχείου ονομάτων":
    uploaded_file = st.file_uploader(
        "Ανεβάστε αρχείο (CSV με στήλη 'Name' ή TXT με ένα όνομα ανά γραμμή):",
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
            f"Επιτυχής φόρτωση {len(input_data_list)} εγγραφών από το αρχείο!"
        )
      except Exception as e:
        st.error(f"Σφάλμα κατά την ανάγνωση του αρχείου: {e}")

  else:  # Ανέβασμα αρχείου SMILES
    is_smiles_mode = True
    uploaded_file = st.file_uploader(
        "Ανεβάστε αρχείο (CSV με στήλη 'SMILES' ή TXT με ένα SMILES ανά γραμμή):",
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
            f"Επιτυχής φόρτωση {len(input_data_list)} SMILES από το αρχείο!"
        )
      except Exception as e:
        st.error(f"Σφάλμα κατά την ανάγνωση του αρχείου: {e}")

  if st.button("Αναζήτηση & Υπολογισμός", type="primary"):
    if not input_data_list:
      st.warning("Παρακαλώ εισάγετε ή ανεβάστε δεδομένα.")
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
              f"Δεν βρέθηκε η ένωση / SMILES: **{display_name}** στο PubChem."
              " Παράλειψη."
          )
          progress_bar.progress((i + 1) / len(input_data_list))
          continue

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
          st.warning(
              f"Μη έγκυρο SMILES για την ένωση: **{display_name}**. Παράλειψη."
          )
          progress_bar.progress((i + 1) / len(input_data_list))
          continue

        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        heavy_atoms = mol.GetNumHeavyAtoms()
        mr = Crippen.MolMR(mol)

        violations = 0
        if not (160 <= mw <= 480):
          violations += 1
        if not (-0.4 <= logp <= 5.6):
          violations += 1
        if not (20 <= heavy_atoms <= 70):
          violations += 1
        if not (40 <= mr <= 130):
          violations += 1

        status = "Passed (Drug-like)" if violations == 0 else "Failed"

        results.append({
            "Όνομα / SMILES": display_name,
            "SMILES": smi,
            "Μοριακό Βάρος (MW)": round(mw, 2),
            "LogP": round(logp, 2),
            "Βαριά Άτομα": heavy_atoms,
            "Molar Refractivity": round(mr, 2),
            "Παραβάσεις": violations,
            "Κατάσταση": status,
        })

        progress_bar.progress((i + 1) / len(input_data_list))

      if results:
        st.session_state.ghose_results = results
        st.rerun()

# ==========================================
# ΚΑΤΑΣΤΑΣΗ 2: ΑΠΟΤΕΛΕΣΜΑΤΑ & ΛΗΨΗ
# ==========================================
else:
  if st.button("🔄 Νέα Αναζήτηση / Αλλαγή Δεδομένων"):
    st.session_state.ghose_results = None
    st.rerun()

  df_results = pd.DataFrame(st.session_state.ghose_results)

  st.subheader("📊 Αποτελέσματα Ανάλυσης")

  df_display = df_results.drop(columns=["SMILES"])

  styled_df = df_display.style.apply(highlight_ghose, axis=1).format({
      "Μοριακό Βάρος (MW)": "{:.2f}",
      "LogP": "{:.2f}",
      "Molar Refractivity": "{:.2f}",
  })

  st.dataframe(styled_df, use_container_width=True)

  # --- ΠΟΛΛΑΠΛΕΣ ΕΠΙΛΟΓΕΣ ΛΗΨΗΣ ΑΡΧΕΙΟΥ ---
  st.markdown("### 📥 Λήψη Αποτελεσμάτων")
  file_format = st.radio(
      "Επιλέξτε μορφή αρχείου για λήψη:",
      options=["CSV (.csv)", "Text (.txt)", "Excel (.xlsx)"],
      horizontal=True,
      key="download_format_radio_ghose",
  )

  if file_format == "CSV (.csv)":
    file_data = df_display.to_csv(index=False).encode("utf-8")
    file_name = "ghose_results.csv"
    mime_type = "text/csv"

  elif file_format == "Text (.txt)":
    file_data = df_display.to_csv(sep="\t", index=False).encode("utf-8")
    file_name = "ghose_results.txt"
    mime_type = "text/plain"

  elif file_format == "Excel (.xlsx)":
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
      df_display.to_excel(writer, index=False, sheet_name="Ghose_Filter")
    file_data = buffer.getvalue()
    file_name = "ghose_results.xlsx"
    mime_type = (
        "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    )

  st.download_button(
      label=f"⬇️ Λήψη αρχείου ({file_name})",
      data=file_data,
      file_name=file_name,
      mime=mime_type,
      type="primary",
  )

st.markdown("""
### What is the Ghose Filter?

The **Ghose Filter** is a foundational knowledge-based physicochemical filter introduced by Arup K. Ghose and colleagues to characterize the drug-like chemical space of small molecules. 

While rules like Lipinski's focus primarily on molecular weight, lipophilicity, and hydrogen bonding, the Ghose filter broadens the evaluation framework by examining additional structural descriptors derived from comprehensive statistical analyses of commercial drug databases.

According to the standard Ghose criteria, a compound is considered to fall within a favorable chemical space when its key physicochemical properties satisfy the following explicit boundaries:

- **Molecular Weight (MW):** 160 – 480 Da
- **Lipophilicity (LogP):** -0.4 – 5.6
- **Heavy Atom Count:** 20 – 70 atoms
- **Molar Refractivity (MR):** 40 – 130

### How Does This Tool Work?

The **Ghose Filter Checker** enables researchers and medicinal chemists to evaluate single or multiple chemical compounds through three flexible input mechanisms:

**1. Compound Name Search (PubChem)**

Users can input commercial or generic chemical names. The tool performs an automated query against the PubChem database to fetch the exact canonical SMILES and molecular representations.

**2. Direct SMILES Input**

Users can directly paste custom or novel SMILES strings (one per line) for rapid in silico evaluation without requiring external database lookups.

**3. File Upload**

Users can upload structured `.csv` files or raw `.txt` files containing lists of compound names or SMILES strings, facilitating high-throughput screening and batch analysis in seconds.

### What Does the Tool Calculate?

For each successfully parsed compound, the application computes the required molecular descriptors and evaluates them against the formal Ghose criteria.

The resulting reports allow researchers to:

- Inspect calculated physicochemical values (MW, LogP, Heavy Atoms, Molar Refractivity).
- Instantly identify specific parameter violations.
- Determine overall compliance through an integrated pass/fail classification.
- Export structured reports securely in multiple standard formats (`.csv`, `.txt`, `.xlsx`).

### Scientific Interpretation and Limitations

The Ghose filter serves as an essential preliminary screening tool in early-stage drug discovery to filter out outlier molecules that deviate significantly from the typical physicochemical profiles of successful oral drugs. 

However, like other heuristic filters, it should not be treated as an absolute predictor of biological activity or clinical success. Compounds failing one or more criteria may still possess therapeutic potential depending on specific target classes and mechanisms of action.

### Scientific Reference

Ghose, A. K., Viswanadhan, V. N., & Wendoloski, J. J. (1999). A knowledge-based approach in designing combinatorial or medicinal chemistry libraries for drug discovery. 1. A qualitative and quantitative characterization of known drug databases. *Journal of Combinatorial Chemistry, 1*(1), 55–68.
""")
