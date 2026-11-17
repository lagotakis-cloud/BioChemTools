import io
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors
import streamlit as st

st.set_page_config(
    page_title="Ghose Filter Screener", page_icon="🧪", layout="wide"
)

from components import render_navbar

render_navbar()

st.title("🧪 Ghose Filter Drug-Likeness Screener")
st.markdown(
    "Έλεγχος συμμόρφωσης μορίων με το κλασικό φίλτρο Ghose (MW, LogP, βαριά"
    " άτομα, Molar Refractivity)."
)
st.markdown("---")

# --- INITIALIZE SESSION STATE ---
if "ghose_results" not in st.session_state:
  st.session_state.ghose_results = None

# ==========================================
# ΚΑΤΑΣΤΑΣΗ 1: ΕΙΣΑΓΩΓΗ ΔΕΔΟΜΕΝΩΝ (ΜΕ ΦΟΡΜΑ)
# ==========================================
if st.session_state.ghose_results is None:
  with st.form("ghose_input_form"):
    input_option = st.radio(
        "Επιλέξτε τρόπο εισαγωγής δεδομένων:",
        (
            "Δοκιμή με έτοιμη λίστα",
            "Πληκτρολόγηση Ονομάτων (ένα ανά γραμμή)",
            "Πληκτρολόγηση SMILES (ένα ανά γραμμή)",
            "Ανέβασμα αρχείου Ονομάτων (.csv / .txt)",
            "Ανέβασμα αρχείου SMILES (.csv / .txt)",
        ),
    )

    default_names = "Aspirin\nParacetamol\nIbuprofen\nCaffeine\nEnalapril"
    default_smiles = (
        "CC(=O)Oc1ccccc1C(=O)O\nCC(=O)Nc1ccc(O)cc1\nCC(C)Cc1ccc(cc1)C(C)C(=O)O"
    )

    text_input_data = ""
    uploaded_file = None

    if input_option == "Δοκιμή με έτοιμη λίστα":
      st.markdown(
          "**Προεπιλεγμένα μόρια:** Aspirin, Paracetamol, Ibuprofen, Caffeine,"
          " Enalapril"
      )
    elif input_option == "Πληκτρολόγηση Ονομάτων (ένα ανά γραμμή)":
      text_input_data = st.text_area(
          "Εισάγετε ονόματα (ένα ανά γραμμή):", value=default_names, height=120
      )
    elif input_option == "Πληκτρολόγηση SMILES (ένα ανά γραμμή)":
      text_input_data = st.text_area(
          "Εισάγετε SMILES strings (ένα ανά γραμμή):",
          value=default_smiles,
          height=120,
      )
    elif input_option == "Ανέβασμα αρχείου Ονομάτων (.csv / .txt)":
      uploaded_file = st.file_uploader(
          "Ανεβάστε αρχείο με στήλη ονομάτων:", type=["csv", "txt"]
      )
    else:
      uploaded_file = st.file_uploader(
          "Ανεβάστε αρχείο με στήλη SMILES:", type=["csv", "txt"]
      )

    submitted = st.form_submit_button(
        "Εκτέλεση Έλεγχου Ghose Filter", type="primary"
    )

  if submitted:
    raw_items = []

    if input_option == "Δοκιμή με έτοιμη λίστα":
      raw_items = ["Aspirin", "Paracetamol", "Ibuprofen", "Caffeine", "Enalapril"]
      mode = "name"
    elif "Πληκτρολόγηση Ονομάτων" in input_option:
      raw_items = [
          line.strip() for line in text_input_data.split("\n") if line.strip()
      ]
      mode = "name"
    elif "Πληκτρολόγηση SMILES" in input_option:
      raw_items = [
          line.strip() for line in text_input_data.split("\n") if line.strip()
      ]
      mode = "smiles"
    elif "Ανέβασμα αρχείου Ονομάτων" in input_option:
      mode = "name"
      if uploaded_file is not None:
        try:
          if uploaded_file.name.endswith(".csv"):
            df_up = pd.read_csv(uploaded_file)
            raw_items = df_up.iloc[:, 0].dropna().astype(str).tolist()
          else:
            stringdata = uploaded_file.getvalue().decode("utf-8")
            raw_items = [
                l.strip() for l in stringdata.splitlines() if l.strip()
            ]
        except Exception as e:
          st.error(f"Σφάλμα ανάγνωσης αρχείου: {e}")
    else:
      mode = "smiles"
      if uploaded_file is not None:
        try:
          if uploaded_file.name.endswith(".csv"):
            df_up = pd.read_csv(uploaded_file)
            raw_items = df_up.iloc[:, 0].dropna().astype(str).tolist()
          else:
            stringdata = uploaded_file.getvalue().decode("utf-8")
            raw_items = [
                l.strip() for l in stringdata.splitlines() if l.strip()
            ]
        except Exception as e:
          st.error(f"Σφάλμα ανάγνωσης αρχείου: {e}")

    results = []

    for item in raw_items:
      mol = None
      mol_name = item
      smi = ""

      if mode == "name":
        try:
          comps = pcp.get_compounds(item, "name")
          if comps:
            smi = comps[0].canonical_smiles
            mol = Chem.MolFromSmiles(smi)
        except:
          pass
      else:
        smi = item
        mol = Chem.MolFromSmiles(smi)
        mol_name = "Custom_SMILES"

      if mol is not None:
        mw = round(Descriptors.MolWt(mol), 2)
        logp = round(Crippen.MolLogP(mol), 2)
        num_atoms = mol.GetNumHeavyAtoms()
        mr = round(Crippen.MolMR(mol), 2)

        passed = (
            (160 <= mw <= 480)
            and (-0.4 <= logp <= 5.6)
            and (20 <= num_atoms <= 70)
            and (40 <= mr <= 130)
        )

        results.append({
            "Molecule": mol_name,
            "SMILES": smi,
            "MW": mw,
            "LogP": logp,
            "Heavy Atoms": num_atoms,
            "Molar Refractivity": mr,
            "Ghose Status": "Passed ✅" if passed else "Failed ❌",
        })

    if not results:
      st.warning(
          "Δεν βρέθηκαν έγκυρα μόρια προς επεξεργασία. Ελέγξτε την είσοδο σας."
      )
    else:
      df_results = pd.DataFrame(results)
      st.session_state.ghose_results = df_results
      st.rerun()

# ==========================================
# ΚΑΤΑΣΤΑΣΗ 2: ΑΠΟΤΕΛΕΣΜΑΤΑ & KPI METRICS
# ==========================================
else:
  if st.button("🔄 Νέα Εισαγωγή / Αλλαγή Δεδομένων"):
    st.session_state.ghose_results = None
    st.rerun()

  df_res = st.session_state.ghose_results
  total_mols = len(df_res)
  passed_mols = len(df_res[df_res["Ghose Status"] == "Passed ✅"])
  success_rate = (
      round((passed_mols / total_mols) * 100, 1) if total_mols > 0 else 0
  )

  # --- KPI METRICS ---
  col1, col2, col3 = st.columns(3)
  col1.metric("Συνολικά Μόρια", total_mols)
  col2.metric("Επιτυχίες (Passed)", f"{passed_mols} ({success_rate}%)")
  col3.metric("Αποτυχίες (Failed)", total_mols - passed_mols)

  st.markdown("---")
  st.subheader("📊 Πίνακας Αποτελεσμάτων Ghose Filter")
  st.dataframe(df_res, use_container_width=True)

  # --- ΛΗΨΗ ΑΡΧΕΙΟΥ (CSV, TXT, EXCEL) ---
  st.markdown("### 📥 Λήψη Αναφοράς")
  file_format = st.radio(
      "Επιλέξτε μορφή αρχείου για λήψη:",
      options=["CSV (.csv)", "Text (.txt)", "Excel (.xlsx)"],
      horizontal=True,
      key="download_format_radio_ghose",
  )

  if file_format == "CSV (.csv)":
    file_data = df_res.to_csv(index=False).encode("utf-8")
    file_name = "ghose_filter_report.csv"
    mime_type = "text/csv"
  elif file_format == "Text (.txt)":
    file_data = df_res.to_csv(sep="\t", index=False).encode("utf-8")
    file_name = "ghose_filter_report.txt"
    mime_type = "text/plain"
  else:
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
      df_res.to_excel(writer, index=False, sheet_name="Ghose_Report")
    file_data = buffer.getvalue()
    file_name = "ghose_filter_report.xlsx"
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
