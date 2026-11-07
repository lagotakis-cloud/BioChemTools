import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, Draw
import pubchempy as pcp

st.set_page_config(page_title="Fragment-Based Drug Design (Ro3)", page_icon="🔬", layout="wide")

st.title("🔬 FBDD & Rule of 3 (Ro3) Checker")
st.markdown("Αξιολόγηση μορίων και θραυσμάτων (Fragments) βάσει του κανόνα **Rule of 3** (Congreve et al.) για Fragment-Based Drug Discovery.")

st.markdown("---")

input_type = st.radio("Επιλέξτε τρόπο εισαγωγής:", ["SMILES", "Όνομα Μορίου (από PubChem)"])

mol = None
molecule_name = ""

if input_type == "SMILES":
    smiles_input = st.text_input("Εισάγετε SMILES θραύσματος:", "c1ccccc1O") # Phenol
    if smiles_input:
        mol = Chem.MolFromSmiles(smiles_input)
        molecule_name = "Custom Fragment"
        if not mol:
            st.error("Μη έγκυρο SMILES.")
else:
    name_input = st.text_input("Όνομα μορίου/θραύσματος:", "Phenol")
    if name_input:
        try:
            compounds = pcp.get_compounds(name_input, 'name')
            if compounds:
                smiles_input = compounds[0].canonical_smiles
                mol = Chem.MolFromSmiles(smiles_input)
                molecule_name = name_input
                st.success(f"Βρέθηκε: {molecule_name} (`{smiles_input}`)")
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
        st.subheader("📊 Υπολογισμός Κανόνων Rule of 3")
        
        # Υπολογισμός ιδιοτήτων
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        h_donors = Lipinski.NumHDonors(mol)
        h_acceptors = Lipinski.NumHAcceptors(mol)
        rot_bonds = Descriptors.NumRotatableBonds(mol)
        
        # Έλεγχος Ro3
        checks = [
            ("Molecular Weight ≤ 300 Da", mw, mw <= 300, f"{mw:.1f} Da"),
            ("LogP ≤ 3", logp, logp <= 3, f"{logp:.2f}"),
            ("H-Bond Donors ≤ 3", h_donors, h_donors <= 3, str(h_donors)),
            ("H-Bond Acceptors ≤ 3", h_acceptors, h_acceptors <= 3, str(h_acceptors)),
            ("Rotatable Bonds ≤ 3", rot_bonds, rot_bonds <= 3, str(rot_bonds))
        ]
        
        passed_count = 0
        for desc, val, passed, val_str in checks:
            status = "✅ Pass" if passed else "⚠️ Fail"
            if passed:
                passed_count += 1
            st.write(f"**{desc}**: `{val_str}` — {status}")
            
        st.markdown("---")
        if passed_count == 5:
            st.success("🎯 **Εξαιρετικό Fragment!** Πληροί απόλυτα όλες τις προϋποθέσεις του Rule of 3 (Ro3).")
        else:
            st.warning(f"⚠️ Το μόριο πληροί {passed_count}/5 κριτήρια του Rule of 3. Ενδέχεται να είναι μεγαλύτερο από το τυπικό 'fragment'.")
