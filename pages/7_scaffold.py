import streamlit as st
from components import render_navbar
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Scaffolds import MurckoScaffold
import pubchempy as pcp

render_navbar()

st.set_page_config(page_title="Molecular Scaffold Splitter", page_icon="🧬", layout="wide")

st.title("🧬 Molecular Scaffold Splitter (Bemis-Murcko)")
st.markdown("Απομόνωση και ανάλυση του κεντρικού δομικού σκελετού (Scaffold / Core) οργανικών μορίων.")

st.markdown("---")

input_type = st.radio("Επιλέξτε τρόπο εισαγωγής:", ["SMILES", "Όνομα Μορίου (από PubChem)"])

mol = None
molecule_name = "Molecule"

if input_type == "SMILES":
    smiles_input = st.text_input("Εισάγετε SMILES:", "CC(=O)Oc1ccccc1C(=O)O") # Aspirin default
    if smiles_input:
        mol = Chem.MolFromSmiles(smiles_input)
        molecule_name = "Custom Molecule"
        if not mol:
            st.error("Μη έγκυρο SMILES string.")
else:
    name_input = st.text_input("Όνομα μορίου:", "Ibuprofen")
    if name_input:
        try:
            compounds = pcp.get_compounds(name_input, 'name')
            if compounds:
                smiles_input = compounds[0].canonical_smiles
                mol = Chem.MolFromSmiles(smiles_input)
                molecule_name = name_input
                st.success(f"Βρέθηκε: {name_input} (`{smiles_input}`)")
            else:
                st.error("Δεν βρέθηκε το μόριο στην PubChem.")
        except Exception as e:
            st.error(f"Σφάλμα: {e}")

if mol:
    st.markdown("---")
    
    # Υπολογισμός Murcko Scaffold
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    scaffold_smiles = Chem.MolToSmiles(scaffold)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("🖼️ Πλήρες Μόριο")
        img_mol = Draw.MolToImage(mol, size=(350, 350))
        st.image(img_mol, use_container_width=True)
        st.caption(f"SMILES: `{Chem.MolToSmiles(mol)}`")
        
    with col2:
        st.subheader("🦴 Bemis-Murcko Scaffold (Πυρήνας)")
        if scaffold and scaffold.GetNumAtoms() > 0:
            img_scaffold = Draw.MolToImage(scaffold, size=(350, 350))
            st.image(img_scaffold, use_container_width=True)
            st.code(scaffold_smiles, language="text")
        else:
            st.info("Το μόριο δεν περιέχει αρωματικούς ή αλεικυκλικούς δακτυλίους για τη δημιουργία Murcko Scaffold (είναι αλυσίδα).")

    st.markdown("---")
    st.subheader("📊 Στατιστικά Δομής")
    c1, c2, c3 = st.columns(3)
    c1.metric("Άτομα Πλήρους Μορίου", mol.GetNumAtoms())
    c2.metric("Άτομα Scaffold", scaffold.GetNumAtoms() if scaffold else 0)
    c3.metric("Ποσοστό Πυρήνα/Μορίου", f"{(scaffold.GetNumAtoms()/mol.GetNumAtoms())*100:.1f}%" if mol.GetNumAtoms() > 0 else "0%")
