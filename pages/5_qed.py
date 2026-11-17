import io
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors
from rdkit.Chem import QED
import streamlit as st

st.set_page_config(
    page_title="QED Drug-Likeness Screener", page_icon="📈", layout="wide"
)

from components import render_navbar

render_navbar()

st.title("📈 QED (Quantitative Estimate of Drug-Likeness)")
st.markdown(
    "Υπολογισμός του σύγχρονου ποσοτικού δείκτη QED (από 0 έως 1) για την"
    " αξιολόγηση της φαρμακευτικής ελκυστικότητας μορίων."
)
st.markdown("---")

# --- INITIALIZE SESSION STATE ---
if "qed_results" not in st.session_state:
  st.session_state.qed_results = None

# ==========================================
# ΚΑΤΑΣΤΑΣΗ 1: ΕΙΣΑΓΩΓΗ ΔΕΔΟΜΕΝΩΝ (ΜΕ ΦΟΡΜΑ)
# ==========================================
if st.session_state.qed_results is None:
  with st.form("qed_input_form"):
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

    submitted = st.form_submit_button("Υπολογισμός QED Score", type="primary")

  if submitted:
    raw_items = []

    if input_option == "Δοκιμή με έτοιμη λίστα":
      raw_items = [
          "Aspirin",
          "Paracetamol",
          "Ibuprofen",
          "Caffeine",
          "Enalapril",
      ]
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
        except Exception:
          pass
      else:
        smi = item
        mol = Chem.MolFromSmiles(smi)
        mol_name = "Custom_SMILES"

      if mol is not None:
        qed_score = round(QED.qed(mol), 3)
        mw = round(Descriptors.MolWt(mol), 2)
        logp = round(Crippen.MolLogP(mol), 2)

        if qed_score >= 0.67:
          classification = "High (Drug-like) ✅"
        elif qed_score >= 0.35:
          classification = "Moderate ℹ️"
        else:
          classification = "Low ⚠️"

        results.append({
            "Molecule": mol_name,
            "SMILES": smi,
            "MW": mw,
            "LogP": logp,
            "QED Score": qed_score,
            "Drug-Likeness": classification,
        })

    if not results:
      st.warning(
          "Δεν βρέθηκαν έγκυρα μόρια προς επεξεργασία. Ελέγξτε την είσοδό σας."
      )
    else:
      df_results = pd.DataFrame(results)
      st.session_state.qed_results = df_results
      st.rerun()

# ==========================================
# ΚΑΤΑΣΤΑΣΗ 2: ΑΠΟΤΕΛΕΣΜΑΤΑ & KPI METRICS
# ==========================================
else:
  if st.button("🔄 Νέα Εισαγωγή / Αλλαγή Δεδομένων"):
    st.session_state.qed_results = None
    st.rerun()

  df_res = st.session_state.qed_results
  total_mols = len(df_res)
  avg_qed = round(df_res["QED Score"].mean(), 3) if total_mols > 0 else 0.0
  high_qed_count = len(df_res[df_res["QED Score"] >= 0.67])

  # --- KPI METRICS ---
  col1, col2, col3 = st.columns(3)
  col1.metric("Συνολικά Μόρια", total_mols)
  col2.metric("Μέσο QED Score", avg_qed)
  col3.metric("Υψηλό QED (≥ 0.67)", high_qed_count)

  st.markdown("---")
  st.subheader("📊 Πίνακας Αποτελεσμάτων QED")
  st.dataframe(df_res, use_container_width=True)

  # --- ΛΗΨΗ ΑΡΧΕΙΟΥ (CSV, TXT, EXCEL) ---
  st.markdown("### 📥 Λήψη Αναφοράς")
  file_format = st.radio(
      "Επιλέξτε μορφή αρχείου για λήψη:",
      options=["CSV (.csv)", "Text (.txt)", "Excel (.xlsx)"],
      horizontal=True,
      key="download_format_radio_qed",
  )

  if file_format == "CSV (.csv)":
    file_data = df_res.to_csv(index=False).encode("utf-8")
    file_name = "qed_evaluation_report.csv"
    mime_type = "text/csv"
  elif file_format == "Text (.txt)":
    file_data = df_res.to_csv(sep="\t", index=False).encode("utf-8")
    file_name = "qed_evaluation_report.txt"
    mime_type = "text/plain"
  else:
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
      df_res.to_excel(writer, index=False, sheet_name="QED_Report")
    file_data = buffer.getvalue()
    file_name = "qed_evaluation_report.xlsx"
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
