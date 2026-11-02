import streamlit as st
import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
from rdkit.Chem import Draw

# --- ΣΥΝΑΡΤΗΣΗ: Αναζήτηση SMILES από όνομα μέσω PubChem API ---
def get_smiles_from_pubchem(compound_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name.strip()}/property/CanonicalSMILES/JSON"
    try:
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            data = response.json()
            smiles = data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
            return smiles
        else:
            return None
    except Exception as e:
        return None

# --- 1. ΤΙΤΛΟΣ & ΠΕΡΙΓΡΑΦΗ ---
st.title("🧪 Lipinski Rule of 5 Checker")
st.write("Πληκτρολόγησε απλά τα ονόματα των φαρμάκων στα αγγλικά και το εργαλείο θα βρει μόνο του τη δομή και τους κανόνες!")

# --- 2. ΕΙΣΑΓΩΓΗ ΔΕΔΟΜΕΝΩΝ ΑΠΟ ΧΡΗΣΤΗ ---
input_option = st.radio("Επιλογή:", ("Δοκιμή με έτοιμα ονόματα", "Πληκτρολόγηση ονομάτων"))

names_list = []

if input_option == "Δοκιμή με έτοιμα ονόματα":
    names_list = ["Aspirin", "Paracetamol", "Caffeine", "Ibuprofen", "Omeprazole"]
    st.write("Δοκιμαστική λίστα:", ", ".join(names_list))
else:
    user_input = st.text_area("Γράψτε τα ονόματα (ένα ανά γραμμή):", 
                              value="Aspirin\nParacetamol\nCaffeine")
    names_list = [name.strip() for name in user_input.split("\n") if name.strip()]

# --- 3. ΥΠΟΛΟΓΙΣΜΟΣ ΚΑΝΟΝΩΝ ---
if st.button("Αναζήτηση & Υπολογισμός"):
    if not names_list:
        st.warning("Παρακαλώ εισάγετε τουλάχιστον ένα όνομα.")
    else:
        results = []
        progress_bar = st.progress(0)
        
        for i, name in enumerate(names_list):
            # 1. Βρίσκουμε το SMILES αυτόματα από το PubChem
            smi = get_smiles_from_pubchem(name)
            
            if smi is None:
                st.error(f"Δεν βρέθηκε η ένωση: **{name}** στο PubChem. Ελέγξτε την ορθογραφία.")
                continue
            
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            
            # 2. Υπολογισμοί Lipinski
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            
            violations = 0
            if mw > 500: violations += 1
            if logp > 5: violations += 1
            if hbd > 5: violations += 1
            if hba > 10: violations += 1
            
            status = "Passed (Drug-like)" if violations <= 1 else "Failed"
            
            results.append({
                "Name": name,
                "SMILES": smi,
                "MW": round(mw, 2),
                "LogP": round(logp, 2),
                "HBD": hbd,
                "HBA": hba,
                "Violations": violations,
                "Status": status,
                "Molecule": mol
            })
            
            progress_bar.progress((i + 1) / len(names_list))
        
        # --- 4. ΠΑΡΟΥΣΙΑΣΗ ΑΠΟΤΕΛΕΣΜΑΤΩΝ ---
        if results:
            df_results = pd.DataFrame(results)
            
            st.subheader("📊 Αποτελέσματα Ανάλυσης")
            st.dataframe(df_results.drop(columns=["Molecule"]))
            
            st.subheader("🔬 2D Δομές Μορίων")
            for index, row in df_results.iterrows():
                col1, col2 = st.columns([1, 2])
                with col1:
                    st.write(f"**{row['Name']}**")
                    st.write(f"Status: {row['Status']}")
                with col2:
                    img = Draw.MolToImage(row["Molecule"], size=(300, 150))
                    st.image(img)
