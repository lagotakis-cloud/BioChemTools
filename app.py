import streamlit as st
import requests
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Draw

st.title("🧪 Παγκόσμιος Υπολογιστής & Έλεγχος Lipinski")
st.write("Αναζητήστε **οποιοδήποτε** μόριο παγκοσμίως με δυναμική αναζήτηση στη βάση δεδομένων του PubChem.")

# Πεδίο ελεύθερης πληκτρολόγησης
search_query = st.text_input("Πληκτρολογήστε όνομα μορίου (π.χ. caffeine, aspirin, paracetamol...):", "")

# Αν ο χρήστης έχει πληκτρολογήσει τουλάχιστον 2 γράμματα, ζητάμε προτάσεις από το PubChem
suggestions = []
if search_query and len(search_query) >= 2:
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/{search_query}/json"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            # Το PubChem επιστρέφει τα αποτελέσματα στο κλειδί 'dictionary_terms' -> 'compound'
            suggestions = data.get("dictionary_terms", {}).get("compound", [])
    except Exception:
        suggestions = []

# Αν βρέθηκαν προτάσεις, τις δείχνουμε σε έναν επιλογέα
selected_molecule = search_query
if suggestions:
    st.caption(f"Βρέθηκαν {len(suggestions)} πιθανές προτάσεις από το PubChem:")
    selected_molecule = st.selectbox("Επιλέξτε από τις προτάσεις:", options=suggestions, index=0)

# Κουμπί ή αυτόματη εκτέλεση αν υπάρχει επιλεγμένο μόριο
if selected_molecule:
    st.info(f"Ανάλυση για: **{selected_molecule}**")
    
    try:
        # Αναζήτηση και φόρτωση μορίου από PubChem
        compound = pcp.get_compounds(selected_molecule, 'name')[0]
        smiles = compound.canonical_smiles
        st.success(f"Βρέθηκε SMILES: `{smiles}`")
        
        # RDKit Υπολογισμοί
        mol = Chem.MolFromSmiles(smiles)
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        h_donors = Descriptors.NumHDonors(mol)
        h_acceptors = Descriptors.NumHAcceptors(mol)
        
        # Εμφάνιση ιδιοτήτων σε στήλες
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
        st.error("⚠️ Δεν ήταν δυνατή η φόρτωση του μορίου. Ελέγξτε την ονομασία.")
