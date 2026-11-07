import streamlit as st
from components import render_navbar
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors
import pubchempy as pcp

render_navbar()

st.set_page_config(page_title="Functional Groups Analyzer", page_icon="🧪", layout="wide")

st.title("🧪 Functional Groups & Substructure Analyzer")
st.markdown("Αυτόματη ανίχνευση και καταμέτρηση βασικών λειτουργικών ομάδων σε οργανικά μόρια.")

st.markdown("---")

input_type = st.radio("Επιλέξτε τρόπο εισαγωγής:", ["SMILES", "Όνομα Μορίου (από PubChem)"])

mol = None

if input_type == "SMILES":
    smiles_input = st.text_input("Εισάγετε SMILES:", "CC(=O)Oc1ccccc1C(=O)O") # Aspirin
    if smiles_input:
        mol = Chem.MolFromSmiles(smiles_input)
        if not mol:
            st.error("Μη έγκυρο SMILES.")
else:
    name_input = st.text_input("Όνομα μορίου:", "Paracetamol")
    if name_input:
        try:
            compounds = pcp.get_compounds(name_input, 'name')
            if compounds:
                smiles_input = compounds[0].canonical_smiles
                mol = Chem.MolFromSmiles(smiles_input)
                st.success(f"Βρέθηκε: {name_input} (`{smiles_input}`)")
            else:
                st.error("Δεν βρέθηκε το μόριο.")
        except Exception as e:
            st.error(f"Σφάλμα: {e}")

if mol:
    st.markdown("---")
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("🖼️ Μοριακή Δομή")
        img = Draw.MolToImage(mol, size=(400, 400))
        st.image(img, use_container_width=True)
        
    with col2:
        st.subheader("🔍 Ανίχνευση Ομάδων")
        
        # Ορισμός βασικών SMARTS patterns για λειτουργικές ομάδες
        patterns = {
            "Υδροξυλική Ομάδα (-OH αλκοόλης)": "[OX2H][CX4]",
            "Καρβοξυλικό Οξύ (-COOH)": "C(=O)[OX2H1]",
            "Εστέρας (-COO-)": "C(=O)O[CX4]",
            "Πρωτοταγής Αμίνη (-NH2)": "[NX3H2]",
            "Αρωματικός Δακτύλιος (Benzene-like)": "c1ccccc1",
            "Καρβονύλιο (Κετόνη/Αλδεΰδη C=O)": "[CX3](=O)[#6]",
            "Αλογόνα (F, Cl, Br, I)": "[F,Cl,Br,I]"
        }
        
        results = []
        for name, smarts in patterns.items():
            patt = Chem.MolFromSmarts(smarts)
            if patt:
                matches = mol.GetSubstructMatches(patt)
                count = len(matches)
                if count > 0:
                    results.append({"Λειτουργική Ομάδα": name, "Εμφανίσεις": count})
                    
        if results:
            import pandas as pd
            df = pd.DataFrame(results)
            st.dataframe(df, use_container_width=True)
        else:
            st.info("Δεν εντοπίστηκαν οι βασικές ομάδες της λίστας αναφοράς.")
