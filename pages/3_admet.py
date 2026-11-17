import io
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
import streamlit as st

st.set_page_config(
    page_title="ADMET & Toxicity Screener | BiochemTools",
    page_icon="assets/organic.png",
    layout="wide",
)

from components import render_navbar

render_navbar()

# --- ΕΠΙΚΕΦΑΛΙΔΑ & ΠΕΡΙΓΡΑΦΗ ---
st.title("ADMET & Toxicity Screener (Enterprise Edition)")
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

  return 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rotb - 0.74 * ar_prop


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


def check_cyp450_inhibition(mol):
  """Εκτίμηση πιθανής αναστολής ηπατικών ενζύμων CYP450 μέσω υποδομών."""
  cyp_patterns = {
      "CYP3A4 Inhibitor Motif": "c1ccccc1[NX3](C)C",
      "CYP2D6 Inhibitor Motif": "[NX3](C)(C)CCC",
      "CYP2C9 Inhibitor Motif": "c1ccc(cc1)C(=O)O",
  }
  inhibited = []
  for name, smarts in cyp_patterns.items():
    pat = Chem.MolFromSmarts(smarts)
    if pat and mol.HasSubstructMatch(pat):
      inhibited.append(name.split()[0])
  return inhibited if inhibited else ["Low Risk"]


def highlight_admet(row):
  styles = [""] * len(row)
  name_idx = row.index.get_loc("Όνομα / SMILES")
  score_idx = row.index.get_loc("ADMET Score (%)")

  if row["ADMET Score (%)"] >= 80:
    score_style = (
        "background-color: #d4edda; color: #155724; font-weight: bold;"
    )
  elif row["ADMET Score (%)"] >= 50:
    score_style = (
        "background-color: #fff3cd; color: #856404; font-weight: bold;"
    )
  else:
    score_style = (
        "background-color: #f8d7da; color: #721c24; font-weight: bold;"
    )

  styles[name_idx] = score_style
  styles[score_idx] = score_style

  if (
      "Poorly" in str(row["Διαλυτότητα"])
      or "Insoluble" in str(row["Διαλυτότητα"])
  ):
    styles[row.index.get_loc("Διαλυτότητα")] = (
        "background-color: #fff3cd; color: #856404;"
    )

  if row["PAINS Match"] != "None detected":
    styles[row.index.get_loc("PAINS Match")] = (
        "background-color: #f8d7da; color: #721c24;"
    )

  if row["Reactive Alerts"] != "None detected":
    styles[row.index.get_loc("Reactive Alerts")] = (
        "background-color: #f8d7da; color: #721c24;"
    )

  if row["CYP450 Profile"] != "Low Risk":
    styles[row.index.get_loc("CYP450 Profile")] = (
        "background-color: #fff3cd; color: #856404;"
    )

  return styles


# ==========================================
# ΚΑΤΑΣΤΑΣΗ 1: ΕΙΣΑΓΩΓΗ ΔΕΔΟΜ
