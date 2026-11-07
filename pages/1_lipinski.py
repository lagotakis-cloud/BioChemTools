import io
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors
from rdkit.Chem import Draw
import streamlit as st
from components import render_navbar

render_navbar()

# --- ΕΠΙΚΕΦΑΛΙΔΑ & ΠΕΡΙΓΡΑΦΗ ---
st.title("🧪 Lipinski Rule of 5 Checker")
st.write(
    "Ελέγξτε τη φαρμακολογική καταλληλότητα μορίων αυτόματα μέσω PubChem, SMILES &"
    " RDKit."
)
st.markdown("---")

# --- INITIALIZE SESSION STATE ---
if "lipinski_results" not in st.session_state:
  st.session_state.lipinski_results = None


# --- ΣΥΝΑΡΤΗΣΕΙΣ ---
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
  name_idx = row.index.get_loc("Όνομα / SMILES")

  if "Passed" in str(row["Κατάσταση"]):
    styles[name_idx] = (
        "background-color: #d4edda; color: #155724; font-weight: bold;"
    )
  else:
    styles[name_idx] = (
        "background-color: #f8d7da; color: #721c24; font-weight: bold;"
    )

  if row["Μοριακό Βάρος (MW)"] > 500:
    styles[row.index.get_loc("Μοριακό Βάρος (MW)")] = (
        "background-color: #f8d7da; color: #721c24;"
    )
  if row["LogP"] > 5:
    styles[row.index.get_loc("LogP")] = "background-color: #f8d7da; color: #721c24;"
  if row["Δότες H (HBD)"] > 5:
    styles[row.index.get_loc("Δότες H (HBD)")] = (
        "background-color: #f8d7da; color: #721c24;"
    )
  if row["Δέκτες H (HBA)"] > 10:
    styles[row.index.get_loc("Δέκτες H (HBA)")] = (
        "background-color: #f8d7da; color: #721c24;"
    )

  return styles


# ==========================================
# ΚΑΤΑΣΤΑΣΗ 1: ΕΙΣΑΓΩΓΗ ΔΕΔΟΜΕΝΩΝ
# ==========================================
if st.session_state.lipinski_results is None:
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
        "Vancomycin",
        "Omeprazole",
    ]
    st.write("Δοκιμαστική λίστα:", ", ".join(input_data_list))

  elif input_option == "Πληκτρολόγηση ονομάτων":
    user_input = st.text_area(
        "Γράψτε τα ονόματα (ένα ανά γραμμή):",
        value="Aspirin\nSertraline\nVancomycin",
    )
    input_data_list = [
        name.strip() for name in user_input.split("\n") if name.strip()
    ]

  elif input_option == "Πληκτρολόγηση SMILES":
    is_smiles_mode = True
    user_input = st.text_area(
        "Γράψτε τα SMILES strings (ένα ανά γραμμή):",
        value=(
            "CC(=O)OC1=CC=CC=C1C(=O)O\nCC(C)CC1=CC=C(C=C1)C(C)C(=O)O\nCN1C=NC2=C1C(=O)N(C(=O)N2C)C"
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
          display_name = item  # Ως όνομα εμφανίζουμε το SMILES string
        else:
          display_name = item
          smi = get_smiles_from_pubchem(item)

        if smi is None:
          st.warning(
              f"Δεν βρέθηκε η ένωση / SMILES: **{display_name}** στο PubChem."
              " Παράλειψη."
          )
          continue

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
          st.warning(
              f"Μη έγκυρο SMILES για την ένωση: **{display_name}**. Παράλειψη."
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
            "Όνομα / SMILES": display_name,
            "SMILES": smi,
            "Μοριακό Βάρος (MW)": round(mw, 2),
            "LogP": round(logp, 2),
            "Δότες H (HBD)": hbd,
            "HBA_val": hba,
            "Δέκτες H (HBA)": hba,
            "Παραβάσεις": violations,
            "Κατάσταση": status,
            "Molecule": mol,
        })

        progress_bar.progress((i + 1) / len(input_data_list))

      if results:
        st.session_state.lipinski_results = results
        st.rerun()

# ==========================================
# ΚΑΤΑΣΤΑΣΗ 2: ΑΠΟΤΕΛΕΣΜΑΤΑ & ΛΗΨΗ
# ==========================================
else:
  if st.button("🔄 Νέα Αναζήτηση / Αλλαγή Δεδομένων"):
    st.session_state.lipinski_results = None
    st.rerun()

  df_results = pd.DataFrame(st.session_state.lipinski_results)

  st.subheader("📊 Αποτελέσματα Ανάλυσης")

  # Αφαίρεση βοηθητικών στηλών
  df_display = df_results.drop(columns=["SMILES", "Molecule", "HBA_val"])

  # Εφαρμογή styling και μορφοποίηση σε 2 δεκαδικά ψηφία
  styled_df = (
      df_display.style.apply(highlight_lipinski, axis=1)
      .format({"Μοριακό Βάρος (MW)": "{:.2f}", "LogP": "{:.2f}"})
  )

  st.dataframe(styled_df, use_container_width=True)

  # --- ΠΟΛΛΑΠΛΕΣ ΕΠΙΛΟΓΕΣ ΛΗΨΗΣ ΑΡΧΕΙΟΥ ---
  st.markdown("### 📥 Λήψη Αποτελεσμάτων")
  file_format = st.radio(
      "Επιλέξτε μορφή αρχείου για λήψη:",
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
      label=f"⬇️ Λήψη αρχείου ({file_name})",
      data=file_data,
      file_name=file_name,
      mime=mime_type,
      type="primary",
  )

  st.markdown("---")
  st.subheader("🔬 2D Δομές Μορίων")
  for index, row in df_results.iterrows():
    col1, col2 = st.columns([1, 2])
    with col1:
      st.write(f"**{row['Όνομα / SMILES']}**")
      st.write(f"Κατάσταση: {row['Κατάσταση']}")
    with col2:
      img = Draw.MolToImage(row["Molecule"], size=(300, 150))
      st.image(img)
