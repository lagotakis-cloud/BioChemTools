import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
from rdkit.Chem import Draw

# --- 1. ΤΙΤΛΟΣ & ΠΕΡΙΓΡΑΦΗ ---
st.title("🧪 Lipinski Rule of 5 Checker")
st.write("Αυτό το εργαλείο ελέγχει τη φαρμακολογική καταλληλότητα μορίων με βάση τους κανόνες Lipinski.")

# --- 2. ΕΙΣΑΓΩΓΗ ΔΕΔΟΜΕΝΩΝ ---
input_option = st.radio("Επιλέξτε τρόπο εισαγωγής μορίων:", ("Δοκιμή με έτοιμα φάρμακα", "Εισαγωγή SMILES χειροκίνητα"))

smiles_list = []

if input_option == "Δοκιμή με έτοιμα φάρμακα":
    # Μερικά γνωστά φάρμακα ως παράδειγμα
    sample_data = {
        "Name": ["Aspirin", "Paracetamol", "Caffeine", "Ibuprofen"],
        "SMILES": [
            "CC(=O)Oc1ccccc1C(=O)O",
            "CC(=O)Nc1ccc(O)cc1",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
        ]
    }
    df_input = pd.DataFrame(sample_data)
    st.write("Χρησιμοποιούμε τα εξής δείγματα:")
    st.dataframe(df_input)
    smiles_list = list(zip(df_input["Name"], df_input["SMILES"]))

else:
    user_input = st.text_area("Πληκτρολογήστε SMILES (ένα ανά γραμμή, με μορφή: Όνομα,SMILES):", 
                              value="Aspirin,CC(=O)Oc1ccccc1C(=O)O\nParacetamol,CC(=O)Nc1ccc(O)cc1")
    for line in user_input.split("\n"):
        if "," in line:
            name, smi = line.split(",", 1)
            smiles_list.append((name.strip(), smi.strip()))

# --- 3. ΥΠΟΛΟΓΙΣΜΟΣ ΚΑΝΟΝΩΝ ---
if st.button("Υπολογισμός Lipinski Rules"):
    if not smiles_list:
        st.warning("Παρακαλώ εισάγετε τουλάχιστον ένα μόριο.")
    else:
        results = []
        
        for name, smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            
            # Έλεγχος κανόνων Lipinski (επιτρέπεται το πολύ 1 παράβαση)
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
        
        # --- 4. ΠΑΡΟΥΣΙΑΣΗ ΑΠΟΤΕΛΕΣΜΑΤΩΝ ---
        df_results = pd.DataFrame(results)
        
        st.subheader("📊 Αποτελέσματα Ανάλυσης")
        
        # Εμφάνιση πίνακα (χωρίς το αντικείμενο Mol για να μην χαλάσει το dataframe view)
        st.dataframe(df_results.drop(columns=["Molecule"]))
        
        # Εμφάνιση μορίων με εικόνες
        st.subheader("🔬 2D Δομές Μορίων")
        for index, row in df_results.iterrows():
            col1, col2 = st.columns([1, 2])
            with col1:
                st.write(f"**{row['Name']}**")
                st.write(f"Status: {row['Status']}")
            with col2:
                img = Draw.MolToImage(row["Molecule"], size=(300, 150))
                st.image(img)
