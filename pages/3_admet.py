import streamlit as st
from components import render_navbar
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski
import pubchempy as pcp

render_navbar()

st.set_page_config(page_title="ADMET & Physicochemical Profiler", page_icon="💊", layout="wide")

st.title("💊 ADMET & Physicochemical Profiler")
st.markdown("Υπολογισμός βασικών φυσικοχημικών ιδιοτήτων και προφίλ ADMET για την εκτίμηση της φαρμακευτικότητας ενός μορίου.")

st.markdown("---")

# Επιλογή τρόπου εισαγωγής
input_type = st.radio("Επιλέξτε τρόπο εισαγωγής μορίου:", ["SMILES", "Όνομα Μορίου (από PubChem)"])

mol = None
molecule_name = ""

if input_type == "SMILES":
    smiles_input = st.text_input("Εισάγετε SMILES:", "CC(=O)OC1=CC=CC=C1C(=O)O") # Default Aspirin
    if smiles_input:
        try:
            mol = Chem.MolFromSmiles(smiles_input)
            molecule_name = "Custom SMILES"
            if mol is None:
                st.error("Μη έγκυρο SMILES string.")
        except Exception as e:
            st.error(f"Σφάλμα: {e}")

else:
    name_input = st.text_input("Εισάγετε όνομα μορίου (π.χ. Aspirin, Caffeine, Paracetamol):", "Caffeine")
    if name_input:
        try:
            compounds = pcp.get_compounds(name_input, 'name')
            if compounds:
                smiles_input = compounds[0].canonical_smiles
                mol = Chem.MolFromSmiles(smiles_input)
                molecule_name = name_input
                st.success(f"Βρέθηκε: {molecule_name} (SMILES: `{smiles_input}`)")
            else:
                st.error("Το μόριο δεν βρέθηκε στην PubChem.")
        except Exception as e:
            st.error(f"Σφάλμα επικοινωνίας με PubChem: {e}")

# Αν έχουμε έγκυρο μόριο, υπολογίζουμε τις ιδιότητες
if mol:
    st.markdown("---")
    st.subheader(f"📊 Αποτελέσματα Ανάλυσης για: {molecule_name}")
    
    # Υπολογισμός Δεικτών
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    h_donors = Lipinski.NumHDonors(mol)
    h_acceptors = Lipinski.NumHAcceptors(mol)
    rot_bonds = Descriptors.NumRotatableBonds(mol)
    
    # Εμφάνιση σε Metrics
    col1, col2, col3 = st.columns(3)
    col1.metric("Molecular Weight (Da)", f"{mw:.2f}")
    col2.metric("LogP (Lipophilicity)", f"{logp:.2f}")
    col3.metric("TPSA (Å²)", f"{tpsa:.2f}")

    col4, col5, col6 = st.columns(3)
    col4.metric("H-Bond Donors", h_donors)
    col5.metric("H-Bond Acceptors", h_acceptors)
    col6.metric("Rotatable Bonds", rot_bonds)
    
    st.markdown("### 🧪 ADMET / Drug-likeness Preliminary Check")
    
    # Απλός έλεγχος ADMET / Veber / Egan rules
    alerts = []
    if tpsa > 140:
        alerts.append("⚠️ Υψηλή TPSA (> 140 Å²): Πιθανή χαμηλή διαπερατότητα από τον αιματοεγκεφαλικό φραγμό ή χαμηλή απορρόφηση από το έντερο.")
    if logp > 5:
        alerts.append("⚠️ Υψηλό LogP (> 5): Πιθανότητα κακής υδατοδιαλυτότητας και αυξημένης λιποφιλικότητας.")
    if rot_bonds > 10:
        alerts.append("⚠️ Αυξημένος αριθμός περιστρεφόμενων δεσμών (> 10): Μειωμένη ευκαμψία σύνδεσης με υποδοχείς (oral bioavailability risk).")
        
    if not alerts:
        st.success("✅ Το μόριο παρουσιάζει εξαιρετικό προφίλ φαρμακευτικότητας βάσει των κλασικών κανόνων ADMET!")
    else:
        for alert in alerts:
            st.warning(alert)
