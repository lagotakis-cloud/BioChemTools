import io
import numpy as np
import pandas as pd
import plotly.express as px
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import streamlit as st

st.set_page_config(
    page_title="Multi-Molecule Similarity Matrix", page_icon="📈", layout="wide"
)

from components import render_navbar

render_navbar()

st.title("📈 Multi-Molecule Tanimoto Similarity Matrix")
st.markdown(
    "Υπολογισμός μήτρας ομοιότητας Tanimoto για πολλαπλά μόρια ταυτόχρονα και"
    " οπτικοποίηση με Heatmap."
)
st.markdown("---")

# --- INITIALIZE SESSION STATE ---
if "tanimoto_results" not in st.session_state:
  st.session_state.tanimoto_results = None

# ==========================================
# ΚΑΤΑΣΤΑΣΗ 1: ΕΙΣΑΓΩΓΗ ΔΕΔΟΜΕΝΩΝ (ΜΕ ΦΟΡΜΑ)
# ==========================================
if st.session_state.tanimoto_results is None:
  with st.form("tanimoto_input_form"):
    input_option = st.radio(
        "Επιλέξτε τρόπο εισαγωγής μορίων:",
        (
            "Δοκιμή με έτοιμη λίστα",
            "Πληκτρολόγηση λίστας (SMILES, Όνομα)",
            "Ανέβασμα αρχείου (CSV / TXT)",
        ),
    )

    default_smiles = (
        "CC(=O)Oc1ccccc1C(=O)O, Aspirin\nCC(=O)Nc1ccc(O)cc1,"
        " Paracetamol\nCC(C)Cc1ccc(cc1)C(C)C(=O)O, Ibuprofen\nCN1C=NC2=C1C(=O)N(C(=O)N2C)C,"
        " Caffeine\nCCO, Ethanol"
    )

    user_text_input = ""
    uploaded_file = None

    if input_option == "Δοκιμή με έτοιμη λίστα":
      st.markdown(
          "**Προεπιλεγμένα μόρια:** Aspirin, Paracetamol, Ibuprofen, Caffeine,"
          " Ethanol"
      )

    elif input_option == "Πληκτρολόγηση λίστας (SMILES, Όνομα)":
      user_text_input = st.text_area(
          "Εισάγετε λίστα μορίων (Μορφή ανά γραμμή: `SMILES, Όνομα`):",
          value=default_smiles,
          height=150,
      )

    else:  # Ανέβασμα αρχείου
      uploaded_file = st.file_uploader(
          "Ανεβάστε αρχείο (CSV ή TXT με στήλες SMILES, Name):",
          type=["csv", "txt"],
      )

    submitted = st.form_submit_button(
        "Υπολογισμός Μήτρας Ομοιότητας", type="primary"
    )

  # --- ΕΠΕΞΕΡΓΑΣΙΑ ΜΕΤΑ ΤΗΝ ΥΠΟΒΟΛΗ ---
  if submitted:
    lines = []

    if input_option == "Δοκιμή με έτοιμη λίστα":
      lines = default_smiles.strip().split("\n")
    elif input_option == "Πληκτρολόγηση λίστας (SMILES, Όνομα)":
      lines = [
          line.strip() for line in user_text_input.split("\n") if line.strip()
      ]
    elif input_option == "Ανέβασμα αρχείου (CSV / TXT)":
      if uploaded_file is not None:
        try:
          if uploaded_file.name.endswith(".csv"):
            df_upload = pd.read_csv(uploaded_file)
            # Υποθέτουμε πρώτη στήλη SMILES, δεύτερη Name (ή προσαρμογή)
            smi_col = df_upload.columns[0]
            name_col = (
                df_upload.columns[1] if len(df_upload.columns) > 1 else smi_col
            )
            lines = [
                f"{row[smi_col]},{row[name_col]}"
                for _, row in df_upload.iterrows()
            ]
          else:
            stringdata = uploaded_file.getvalue().decode("utf-8")
            lines = [
                line.strip() for line in stringdata.splitlines() if line.strip()
            ]
        except Exception as e:
          st.error(f"Σφάλμα κατά την ανάγνωση του αρχείου: {e}")

    molecules = []
    names = []

    for line in lines:
      parts = line.split(",")
      if len(parts) >= 1:
        smi = parts[0].strip()
        name = parts[1].strip() if len(parts) > 1 else smi
        mol = Chem.MolFromSmiles(smi)
        if mol:
          molecules.append(mol)
          names.append(name)

    if len(molecules) < 2:
      st.warning(
          "Παρακαλώ εισάγετε τουλάχιστον 2 έγκυρα μόρια για τη δημιουργία της"
          " μήτρας ομοιότητας."
      )
    else:
      # Υπολογισμός Morgan/RDKit Fingerprints και Matrix Tanimoto
      fps = [FingerprintMols.FingerprintMol(m) for m in molecules]
      n = len(molecules)
      matrix = [[0.0 for _ in range(n)] for _ in range(n)]

      for i in range(n):
        for j in range(n):
          matrix[i][j] = round(
              DataStructs.TanimotoSimilarity(fps[i], fps[j]), 3
          )

      df_matrix = pd.DataFrame(matrix, index=names, columns=names)

      # Αποθήκευση στο session_state
      st.session_state.tanimoto_results = {
          "df_matrix": df_matrix,
          "n": n,
          "matrix": matrix,
      }
      st.rerun()

# ==========================================
# ΚΑΤΑΣΤΑΣΗ 2: ΑΠΟΤΕΛΕΣΜΑΤΑ & ΛΗΨΗ
# ==========================================
else:
  if st.button("🔄 Νέα Εισαγωγή / Αλλαγή Δεδομένων"):
    st.session_state.tanimoto_results = None
    st.rerun()

  data = st.session_state.tanimoto_results
  df_matrix = data["df_matrix"]
  n = data["n"]
  matrix = data["matrix"]

  # --- KPI METRICS ---
  arr = np.array(matrix)
  mask = ~np.eye(arr.shape[0], dtype=bool)
  avg_similarity = round(arr[mask].mean(), 3) if arr[mask].size > 0 else 0.0

  col1, col2, col3 = st.columns(3)
  col1.metric("Συνολικά Μόρια", n)
  col2.metric("Μέση Ομοιότητα (Tanimoto)", avg_similarity)
  col3.metric("Διαστάσεις Μήτρας", f"{n} x {n}")

  st.markdown("---")
  st.subheader("📊 Πίνακας Ομοιότητας Tanimoto")
  st.dataframe(df_matrix, use_container_width=True)

  # --- ΛΗΨΗ ΑΡΧΕΙΟΥ ---
  st.markdown("### 📥 Λήψη Μήτρας Ομοιότητας")
  file_format = st.radio(
      "Επιλέξτε μορφή αρχείου για λήψη:",
      options=["CSV (.csv)", "Text (.txt)", "Excel (.xlsx)"],
      horizontal=True,
      key="download_format_radio_tanimoto",
  )

  if file_format == "CSV (.csv)":
    file_data = df_matrix.to_csv().encode("utf-8")
    file_name = "tanimoto_matrix.csv"
    mime_type = "text/csv"
  elif file_format == "Text (.txt)":
    file_data = df_matrix.to_csv(sep="\t").encode("utf-8")
    file_name = "tanimoto_matrix.txt"
    mime_type = "text/plain"
  else:
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
      df_matrix.to_excel(writer, index=True, sheet_name="Tanimoto_Matrix")
    file_data = buffer.getvalue()
    file_name = "tanimoto_matrix.xlsx"
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

  st.markdown("---")

  # Plotly Heatmap
  st.markdown("### 🌡️ Heatmap Οπτικοποίηση")
  fig = px.imshow(
      df_matrix,
      text_auto=True,
      color_continuous_scale="Blues",
      aspect="auto",
      labels=dict(color="Tanimoto Similarity"),
  )
  fig.update_layout(
      title="Σύγκριση Ομοιότητας Ζεύγους Μορίων",
      xaxis_title="",
      yaxis_title="",
      plot_bgcolor="rgba(0,0,0,0)",
      paper_bgcolor="rgba(0,0,0,0)",
  )
  st.plotly_chart(fig, use_container_width=True)
