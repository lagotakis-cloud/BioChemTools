import streamlit as st
from components import render_navbar
from rdkit import Chem
from rdkit.Chem import Draw
import pubchempy as pcp

render_navbar()

st.set_page_config(page_title="Chemical Format Converter", page_icon="🔄", layout="wide")

st.title("🔄 Chemical Format Converter (SMILES ⇄ InChI ⇄ SDF)")
st.markdown("Μετατροπή χημικών δομών μεταξύ διαφόρων μορφών (SMILES, InChI, InChIKey, SDF) και δυνατότητα λήψης αρχείων.")

st.markdown("---")

input_type = st.radio("Επιλέξτε τρόπο εισαγωγής:", ["SMILES", "Όνομα Μορίου (από PubChem)"])

mol = None
molecule_name = "molecule"

if input_type == "SMILES":
    smiles_input = st.text_input("Εισάγετε SMILES:", "CCO") # Ethanol default
    if smiles_input:
        mol = Chem.MolFromSmiles(smiles_input)
        molecule_name = "custom_molecule"
        if not mol:
            st.error("Μη έγκυρο SMILES string.")
else:
    name_input = st.text_input("Όνομα μορίου:", "Aspirin")
    if name_input:
        try:
            compounds = pcp.get_compounds(name_input, 'name')
            if compounds:
                smiles_input = compounds[0].canonical_smiles
                mol = Chem.MolFromSmiles(smiles_input)
                molecule_name = name_input.lower().replace(" ", "_")
                st.success(f"Βρέθηκε: {name_input} (`{smiles_input}`)")
            else:
                st.error("Δεν βρέθηκε το μόριο στην PubChem.")
        except Exception as e:
            st.error(f"Σφάλμα: {e}")

if mol:
    st.markdown("---")
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.subheader("🖼️ Δομή")
        img = Draw.MolToImage(mol, size=(300, 300))
        st.image(img, use_container_width=True)
        
    with col2:
        st.subheader("📄 Μετατροπές Μορφών")
        
        # Υπολογισμός μορφών
        canonical_smiles = Chem.MolToSmiles(mol)
        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.MolToInchiKey(mol)
        mol_block = Chem.MolToMolBlock(mol)
        
        st.text_input("Canonical SMILES", value=canonical_smiles)
        st.text_input("InChI", value=inchi)
        st.text_input("InChIKey", value=inchikey)
        
        st.markdown("### 📥 Λήψη Αρχείων (Export)")
        col_d1, col_d2 = st.columns(2)
        
        with col_d1:
            st.download_button(
                label="📥 Λήψη ως .SDF",
                data=mol_block,
                file_name=f"{molecule_name}.sdf",
                mime="chemical/x-mdl-sdfile"
            )
            
        with col_d2:
            st.download_button(
                label="📥 Λήψη ως .smi",
                data=canonical_smiles,
                file_name=f"{molecule_name}.smi",
                mime="text/plain"
            )
