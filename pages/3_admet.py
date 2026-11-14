import io
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
import streamlit as st
from components import render_navbar

render_navbar()

# --- ΕΠΙΚΕΦΑΛΙΔΑ & ΠΕΡΙΓΡΑΦΗ ---
st.title("🧪 ADMET & Toxicity Screener")
st.write(
    "Αξιολόγηση υδατικής διαλυτότητας (ESOL), σάρωση PAINS filters και ανίχνευση"
    " δραστικών/τοξικών χημικών ομάδων (Toxicophores)."
)
st.markdown("---")

# --- INITIALIZE SESSION STATE ---
if "admet_results" not in st.session_state:
  st.session_state.admet_results = None


# --- ΣΥΝΑΡΤΗΣΕΙΣ ---
def get_smiles_from_pubchem(compound_name):
  try:
    compounds = pcp.get_compounds(compound_name.strip(), "name")
    if compounds:
      return compounds[0].canonical_smiles
    return None
  except Exception:
    return None


def calculate_esol(mol):
  """Υπολογισμός υδατικής διαλυτότητας LogS (Delaney equation)."""
  logp = Descriptors.MolLogP(mol)
  mw = Descriptors.MolWt(mol)
  rotb = Descriptors.NumRotatableBonds(mol)
  heavy_atoms = mol.GetNumHeavyAtoms()
  ar_atoms = len([a for a in mol.GetAtoms() if a.GetIsAromatic()])
  ar_prop = ar_atoms / heavy_atoms if heavy_atoms > 0 else 0

  log_s = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rotb - 0.74 * ar_prop
  return log_s


def classify_solubility(log_s):
  """Κατηγοριοποίηση υδατικής διαλυτότητας βάσει LogS."""
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
  """Ανίχνευση γνωστών τοξικών / δραστικών υποδομών μέσω SMARTS."""
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


def highlight_admet(row):
  styles = [""] * len(row)
  status_idx = row.index.get_loc("Συνολική Κατάσταση")
  name_idx = row.index.get_loc("Όνομα / SMILES")

  if "Safe" in str(row["Συνολική Κατάσταση"]):
    color_style = (
        "background-color: #d4edda; color: #155724; font-weight: bold;"
    )
  else:
    color_style = (
        "background-color: #f8d7da; color: #721c24; font-weight: bold;"
    )

  styles[name_idx] = color_style
  styles[status_idx] = color_style
  return styles


# ==========================================
# ΚΑΤΑΣΤΑΣΗ 1: ΕΙΣΑΓΩΓΗ ΔΕΔΟΜΕΝΩΝ
# ==========================================
if st.session_state.admet_results is None:
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
        "Paracetamol",
        "Caffeine",
        "Glycidol",
        "Emodin",
        "Nitrobenzene",
    ]
    st.write("Δοκιμαστική λίστα:", ", ".join(input_data_list))

  elif input_option == "Πληκτρολόγηση ονομάτων":
    user_input = st.text_area(
        "Γράψτε τα ονόματα (ένα ανά γραμμή):",
        value="Aspirin\nParacetamol\nCaffeine\nGlycidol",
    )
    input_data_list = [
        name.strip() for name in user_input.split("\n") if name.strip()
    ]

  elif input_option == "Πληκτρολόγηση SMILES":
    is_smiles_mode = True
    user_input = st.text_area(
        "Γράψτε τα SMILES strings (ένα ανά γραμμή):",
        value=(
            "CC(=O)OC1=CC=CC=C1C(=O)O\n"
            "CC(=O)NC1=CC=C(C=C1)O\n"
            "C1CO1CO\n"
            "CC1=CC2=C(C(=C1)O)C(=O)C3=C(C2=O)C=C(C=C3O)O"
        ),
    )
    input_data_list = [
        smi.strip() for smi in user_input.split("\n") if smi.strip()
    ]

  elif input_option == "Ανέβασμα αρχείου ονομάτων":
    uploaded_file = st.file_uploader(
        "Ανεβάστε αρχείο (CSV ή TXT):", type=["csv", "txt"]
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
        "Ανεβάστε αρχείο SMILES (CSV ή TXT):", type=["csv", "txt"]
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

  if st.button("Εκτέλεση ADMET & Toxicity Screening", type="primary"):
    if not input_data_list:
      st.warning("Παρακαλώ εισάγετε ή ανεβάστε δεδομένα.")
    else:
      results = []
      not_found = []
      progress_bar = st.progress(0)

      # Φόρτωση καταλόγου PAINS (A, B, C)
      params = FilterCatalogParams()
      params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
      catalog = FilterCatalog(params)

      for i, item in enumerate(input_data_list):
        if is_smiles_mode:
          smi = item
          display_name = item
        else:
          display_name = item
          smi = get_smiles_from_pubchem(item)

        if smi is None:
          not_found.append(display_name)
          progress_bar.progress((i + 1) / len(input_data_list))
          continue

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
          not_found.append(display_name)
          progress_bar.progress((i + 1) / len(input_data_list))
          continue

        # Υπολογισμοί ADMET
        esol = calculate_esol(mol)
        sol_class = classify_solubility(esol)
        has_pains = catalog.HasMatch(mol)
        alerts = check_structural_alerts(mol)

        # Συνολική Κατάσταση
        if has_pains or len(alerts) > 0:
          toxicity_status = "Alert / PAINS Match"
        else:
          toxicity_status = "Safe (No Alerts)"

        results.append({
            "Όνομα / SMILES": display_name,
            "SMILES": smi,
            "LogS (ESOL)": round(esol, 2),
            "Κατηγορία Διαλυτότητας": sol_class,
            "PAINS Match": "Yes" if has_pains else "No",
            "Reactive Alerts": (
                ", ".join(alerts) if alerts else "None detected"
            ),
            "Συνολική Κατάσταση": toxicity_status,
        })

        progress_bar.progress((i + 1) / len(input_data_list))

      if not_found:
        st.warning(
            f"Δεν βρέθηκαν ή δεν ήταν έγκυρα τα εξής: {', '.join(not_found)}"
        )

      if results:
        st.session_state.admet_results = results
        st.rerun()

# ==========================================
# ΚΑΤΑΣΤΑΣΗ 2: ΑΠΟΤΕΛΕΣΜΑΤΑ & ΛΗΨΗ
# ==========================================
else:
  if st.button("🔄 Νέα Αναζήτηση / Αλλαγή Δεδομένων"):
    st.session_state.admet_results = None
    st.rerun()

  df_results = pd.DataFrame(st.session_state.admet_results)

  # --- KPI METRICS ---
  total_mols = len(df_results)
  safe_mols = len(
      df_results[df_results["Συνολική Κατάσταση"].str.contains("Safe")]
  )
  alert_mols = total_mols - safe_mols

  col1, col2, col3 = st.columns(3)
  col1.metric("Συνολικά Μόρια", total_mols)
  col2.metric("Ασφαλή (No Alerts)", safe_mols)
  col3.metric("Alerts / PAINS", alert_mols)

  st.markdown("---")

  st.subheader("📊 Αποτελέσματα ADMET & Toxicity Screening")
  st.info(
      "**Οδηγός LogS (ESOL):** Τιμές > 0 (Highly Soluble), -2 έως 0 (Very"
      " Soluble), -4 έως -2 (Soluble), -6 έως -4 (Moderately Soluble), < -6"
      " (Poorly Soluble / Insoluble)."
  )

  df_display = df_results.drop(columns=["SMILES"])

  styled_df = df_display.style.apply(highlight_admet, axis=1).format(
      {"LogS (ESOL)": "{:.2f}"}
  )

  st.dataframe(styled_df, use_container_width=True)

  # --- ΛΗΨΗ ΑΡΧΕΙΟΥ ---
  st.markdown("### 📥 Λήψη Αποτελεσμάτων ADMET")
  file_format = st.radio(
      "Επιλέξτε μορφή αρχείου για λήψη:",
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
      label=f"⬇️ Λήψη αρχείου ({file_name})",
      data=file_data,
      file_name=file_name,
      mime=mime_type,
      type="primary",
  )
