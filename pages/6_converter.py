import io
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
import streamlit as st

st.set_page_config(
    page_title="Chemical Format Converter", page_icon="🔄", layout="wide"
)

from components import render_navbar

render_navbar()

st.title("🔄 Chemical Format Converter (SMILES ⇄ InChI ⇄ SDF)")
st.markdown(
    "Μετατροπή χημικών δομών μεταξύ διαφόρων μορφών (SMILES, InChI, InChIKey,"
    " SDF) και δυνατότητα λήψης αρχείων."
)
st.markdown("---")

# --- INITIALIZE SESSION STATE ---
if "converter_results" not in st.session_state:
  st.session_state.converter_results = None

# ==========================================
# ΚΑΤΑΣΤΑΣΗ 1: ΕΙΣΑΓΩΓΗ ΔΕΔΟΜΕΝΩΝ (ΜΕ ΦΟΡΜΑ)
# ==========================================
if st.session_state.converter_results is None:
  with st.form("converter_input_form"):
    input_option = st.radio(
        "Επιλέξτε τρόπο εισαγωγής:",
        ("SMILES", "Όνομα Μορίου (από PubChem)"),
    )

    smiles_input = ""
    name_input = ""

    if input_option == "SMILES":
      smiles_input = st.text_input("Εισάγετε SMILES:", "CCO")  # Ethanol default
    else:
      name_input = st.text_input("Όνομα μορίου:", "Aspirin")

    submitted = st.form_submit_button(
        "Μετατροπή & Ανάλυση Μορίου", type="primary"
    )

  if submitted:
    mol = None
    molecule_name = "molecule"
    resolved_smiles = ""

    if input_option == "SMILES":
      resolved_smiles = smiles_input.strip()
      mol = Chem.MolFromSmiles(resolved_smiles)
      molecule_name = "custom_molecule"
    else:
      try:
        compounds = pcp.get_compounds(name_input.strip(), "name")
        if compounds:
          resolved_smiles = compounds[0].canonical_smiles
          mol = Chem.MolFromSmiles(resolved_smiles)
          molecule_name = name_input.strip().lower().replace(" ", "_")
      except Exception:
        mol = None

    if mol is None:
      st.error(
          "Μη έγκυρο μόριο ή δεν βρέθηκε. Ελέγξτε την είσοδο και δοκιμάστε"
          " ξανά."
      )
    else:
      # Υπολογισμός μορφών
      canonical_smiles = Chem.MolToSmiles(mol)
      inchi = Chem.MolInchi.MolToInchi(mol) if hasattr(Chem, "MolInchi") else Chem.MolToInchi(mol)
      inchikey = Chem.MolToInchiKey(mol)
      mol_block = Chem.MolToMolBlock(mol)

      st.session_state.converter_results = {
          "molecule_name": molecule_name,
          "canonical_smiles": canonical_smiles,
          "inchi": inchi,
          "inchikey": inchikey,
          "mol_block": mol_block,
          "mol": mol,
      }
      st.rerun()

# ==========================================
# ΚΑΤΑΣΤΑΣΗ 2: ΑΠΟΤΕΛΕΣΜΑΤΑ & ΛΗΨΗ
# ==========================================
else:
  if st.button("🔄 Νέα Μετατροπή / Αλλαγή Μορίου"):
    st.session_state.converter_results = None
    st.rerun()

  res = st.session_state.converter_results
  molecule_name = res["molecule_name"]
  canonical_smiles = res["canonical_smiles"]
  inchi = res["inchi"]
  inchikey = res["inchikey"]
  mol_block = res["mol_block"]
  mol = res["mol"]

  # --- KPI METRICS ---
  col_m1, col_m2, col_m3 = st.columns(3)
  col_m1.metric("Βαριά Άτομα", mol.GetNumHeavyAtoms())
  col_m2.metric("Μοριακό Βάρος", round(Chem.Descriptors.MolWt(mol), 2))
  col_m3.metric("Αριθμός Δακτυλίων", Chem.CalcNumRings(mol))

  st.markdown("---")
  col1, col2 = st.columns([1, 2])

  with col1:
    st.subheader("🖼️ Δομή")
    img = Draw.MolToImage(mol, size=(300, 300))
    st.image(img, use_container_width=True)

  with col2:
    st.subheader("📄 Μετατροπές Μορφών")

    st.text_input("Canonical SMILES", value=canonical_smiles)
    st.text_input("InChI", value=inchi)
    st.text_input("InChIKey", value=inchikey)

    st.markdown("### 📥 Λήψη Αρχείων (Export)")
    col_d1, col_d2, col_d3 = st.columns(3)

    with col_d1:
      st.download_button(
          label="📥 Λήψη .SDF",
          data=mol_block,
          file_name=f"{molecule_name}.sdf",
          mime="chemical/x-mdl-sdfile",
          type="primary",
      )

    with col_d2:
      st.download_button(
          label="📥 Λήψη .smi",
          data=canonical_smiles,
          file_name=f"{molecule_name}.smi",
          mime="text/plain",
      )

    with col_d3:
      st.download_button(
          label="📥 Λήψη InChI",
          data=inchi,
          file_name=f"{molecule_name}_inchi.txt",
          mime="text/plain",
      )
