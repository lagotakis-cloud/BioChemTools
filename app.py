
import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Draw

st.title("🧪 Υπολογιστής & Έλεγχος Κανόνων Lipinski")
st.write("Αναζητήστε ένα μόριο με αυτόματη συμπλήρωση (autocomplete) και ελέγξτε τη φαρμακευτική του συμμόρφωση.")

# Μια πλούσια λίστα με δημοφιλή μόρια/φάρμακα για το autocomplete
molecule_database = [
    "Caffeine",
    "Cocaine",
    "Aspirin",
    "Paracetamol",
    "Ibuprofen",
    "Cyclosporine A",
    "Morphine",
    "Atorvastatin",
    "Metformin",
    "Omeprazole",
    "Amoxicillin",
    "Diazepam"
]

# Το st.selectbox λειτουργεί αυτόματα ως autocomplete καθώς ο χρήστης πληκτρολογεί
selected_molecule = st.selectbox(
    "Πληκτρολογήστε ή επιλέξτε μόριο (π.χ. c, ca...):",
    options=[""] + sorted(molecule_database),
    index=0
)

# Αν ο χρήστης επιλέξει ή γράψει κάτι:
if selected_molecule:
    st.info(f"Αναζήτηση για: **{selected_molecule}**")
    
    try:
        # Αναζήτηση στο PubChem
        compound = pcp.get_compounds(selected_molecule, 'name')[0]
        smiles = compound.canonical_smiles
        st.success(f"Βρέθηκε SMILES: `{smiles}`")
        
        # RDKit Υπολογισμοί
        mol = Chem.MolFromSmiles(smiles)
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        h_donors = Descriptors.NumHDonors(mol)
        h_acceptors = Descriptors.NumHAcceptors(mol)
        
        # Εμφάνιση ιδιοτήτων
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Μοριακό Βάρος (MW)", f"{mw:.2f} g/mol")
            st.metric("H-Donors", h_donors)
        with col2:
            st.metric("LogP", f"{logp:.2f}")
            st.metric("H-Acceptors", h_acceptors)
            
        # Έλεγχος Lipinski (Rule of 5)
        violations = 0
        if mw > 500: violations += 1
        if logp > 5: violations += 1
        if h_donors > 5: violations += 1
        if h_acceptors > 10: violations += 1
        
        st.subheader("📋 Αποτελέσματα Lipinski Rule of 5")
        if violations == 0:
            st.success("✅ Εξαιρετική συμμόρφωση (0 παραβιάσεις). Κατάλληλο για από του στόματος χορήγηση!")
        else:
            st.warning(f"⚠️ Εντοπίστηκαν **{violations}** παραβιάσεις των κανόνων του Lipinski.")
            
        # Σχεδίαση μορίου
        st.subheader("Molecule Structure")
        img = Draw.MolToImage(mol, size=(300, 300))
        st.image(img)
        
    except Exception as e:
        st.error("⚠️ Δεν βρέθηκε το μόριο στη βάση δεδομένων. Δοκιμάστε άλλη ονομασία.")
