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
    "Αξιολόγηση υδατικής διαλυτότητας (ESOL), σάρωση PAINS filters με λεπτομερείς"
    " περιγραφές και ανίχνευση δραστικών/τοξικών χημικών ομάδων (Toxicophores)."
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
        value="Aspirin\nParacetamol\nCaffeine\nGlycidol\nEmodin",
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
