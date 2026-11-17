import io
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import QED
import streamlit as st

st.set_page_config(
    page_title="QED Drug-Likeness Score", page_icon="📈", layout="wide"
)

from components import render_navbar

render_navbar()

st.title("📈 QED (Quantitative Estimate of Drug-Likeness)")
st.markdown(
    "Υπολογισμός του σύγχρονου αριθμητικού δείκτη QED (από 0 έως 1) για την"
    " εκτίμηση της φαρμακευτικής ελκυστικότητας."
)
st.markdown("---")

if "qed_results" not in st.session_state:
  st.session_state.qed_results = None

if st.session_state.qed_results is None:
  with st.form("qed_input_form"):
    input_option = st.radio(
        "Επιλέξτε τρόπο εισαγωγής:",
        ("SMILES", "Όνομα Μορίου (από PubChem)"),
    )
    smiles_input = ""
    name_input = ""

    if input_option == "SMILES":
      smiles_input = st.text_input("Εισάγετε SMILES:", "CC(=O)Oc1ccccc1C(=O)O")
    else:
      name_input = st.text_input("Όνομα μορίου:", "Aspirin")

    submitted = st.form_submit_button("Υπολογισμός QED Score", type="primary")

  if submitted:
    mol = None
    resolved_smiles = ""
    molecule_name = "molecule"

    if input_option == "SMILES":
      resolved_smiles = smiles_input.strip()
      mol = Chem.MolFromSmiles(resolved_smiles)
    else:
      try:
        comps = pcp.get_compounds(name_input.strip(), "name")
        if comps:
          resolved_smiles = comps[0].canonical_smiles
          mol = Chem.MolFromSmiles(resolved_smiles)
          molecule_name = name_input.strip().lower().replace(" ", "_")
      except:
        mol = None

    if mol is None:
      st.error("Μη έγκυρο μόριο ή δεν βρέθηκε.")
    else:
      qed_score = round(QED.qed(mol), 3)
      mw = round(Descriptors.MolWt(mol), 2)
      logp = round(Descriptors.MolLogP(mol), 2)

      st.session_state.qed_results = {
          "molecule_name": molecule_name,
          "qed_score": qed_score,
          "mw": mw,
          "logp": logp,
          "mol": mol,
      }
      st.rerun()
else:
  if st.button("🔄 Νέος Υπολογισμός"):
    st.session_state.qed_results = None
    st.rerun()

  res = st.session_state.qed_results
  score = res["qed_score"]
  mol = res["mol"]

  col1, col2, col3 = st.columns(3)
  col1.metric("QED Score (0 - 1)", score)
  col2.metric("Μοριακό Βάρος", res["mw"])
  col3.metric("LogP", res["logp"])

  st.markdown("---")
  c1, c2 = st.columns([1, 2])
  with c1:
    from rdkit.Chem import Draw

    img = Draw.MolToImage(mol, size=(300, 300))
    st.image(img, use_container_width=True)
  with c2:
    st.subheader("📊 Ερμηνεία Δείκτη QED")
    if score >= 0.67:
      st.success("✨ **Υψηλή φαρμακευτική ελκυστικότητα** (Drug-like).")
    elif score >= 0.35:
      st.info("ℹ️ **Μέτρια φαρμακευτική ελκυστικότητα**.")
    else:
      st.warning("⚠️ **Χαμηλή φαρμακευτική ελκυστικότητα**.")

    df_out = pd.DataFrame([{
        "Molecule": res["molecule_name"],
        "QED_Score": score,
        "MW": res["mw"],
        "LogP": res["logp"],
    }])
    st.download_button(
        "⬇️ Λήψη αποτελέσματος (CSV)",
        data=df_out.to_csv(index=False).encode("utf-8"),
        file_name="qed_result.csv",
        mime="text/csv",
        type="primary",
    )
