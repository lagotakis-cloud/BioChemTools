
import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Draw

st.title("Έλεγχος Κανόνων Lipinski")
st.write("Εισάγετε το **όνομα της ένωσης** (π.χ. *Aspirin*, *Ibuprofen*, *Caffeine*, *Paracetamol*) για να δείτε αυτόματα τη δομή και τις ιδιότητές της.")

# Πεδίο κειμένου για το όνομα
compound_name = st.text_input("Όνομα Ένωσης:", "Aspirin")

if compound_name:
    with st.spinner("Αναζήτηση στη βάση δεδομένων PubChem..."):
        try:
            # Αναζήτηση SMILES βάσει ονόματος μέσω PubChem
            compounds = pcp.get_compounds(compound_name, 'name')
            if compounds:
                smiles_input = compounds[0].canonical_smiles
            else:
                smiles_input = None
        except Exception as e:
            smiles_input = None

    if not smiles_input:
        st.error(f"❌ Δεν βρέθηκε η ένωση '{compound_name}' στο PubChem. Δοκιμάστε άλλο όνομα (στα αγγλικά).")
    else:
        st.success(f"✅ Βρέθηκε! (SMILES: `{smiles_input}`)")
        
        mol = Chem.MolFromSmiles(smiles_input)
        
        if mol is None:
            st.error("❌ Σφάλμα κατά τη δημιουργία του μορίου.")
        else:
            # Εμφάνιση εικόνας μορίου
            st.subheader("Μοριακή Δομή")
            img = Draw.MolToImage(mol, size=(300, 300))
            st.image(img)
            
            # Υπολογισμοί
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            
            # Εμφάνιση ιδιοτήτων
            st.subheader("Φυσικοχημικές Ιδιότητες")
            c1, c2 = st.columns(2)
            c1.metric("Μοριακό Βάρος (MW)", f"{mw:.2f} g/mol")
            c1.metric("LogP", f"{logp:.2f}")
            c2.metric("Δότες δεσμών H (HBD)", hbd)
            c2.metric("Δέκτες δεσμών H (HBA)", hba)
            
            # Αξιολόγηση Lipinski
            st.subheader("Αξιολόγηση Lipinski (Rule of 5)")
            if mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10:
                st.markdown("🟢 **Το μόριο ικανοποιεί πλήρως τους κανόνες του Lipinski!**")
            else:
                st.markdown("🔴 **Το μόριο παραβιάζει κάποιον/ους από τους κανόνες του Lipinski.**")
