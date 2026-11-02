import streamlit as st
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
from rdkit.Chem import Draw

# --- ΣΥΝΑΡΤΗΣΗ: Αναζήτηση SMILES από όνομα μέσω PubChemPy ---
def get_smiles_from_pubchem(compound_name):
    try:
        compounds = pcp.get_compounds(compound_name.strip(), 'name')
        if compounds:
            return compounds[0].canonical_smiles
        else:
            return None
    except Exception as e:
        return None

# --- ΣΥΝΑΡΤΗΣΗ: Χρωματισμός κελιών (Styling) ---
def highlight_lipinski(row):
    styles = [''] * len(row)
    
    # 1. Χρωματισμός του ονόματος (Name) ανάλογα με το Status
    name_idx = row.index.get_loc('Name')
    if "Passed" in str(row['Status']):
        styles[name_idx] = 'background-color: #d4edda; color: #155724; font-weight: bold;'
    else:
        styles[name_idx] = 'background-color: #f8d7da; color: #721c24; font-weight: bold;'
        
    # 2. Χρωματισμός συγκεκριμένων παραμέτρων αν παραβιάζονται
    if row['MW'] > 500:
        idx = row.index.get_loc('MW')
        styles[idx] = 'background-color: #f8d7da; color: #721c24;'
        
    if row['LogP'] > 5:
        idx = row.index.get_loc('LogP')
        styles[idx] = 'background-color: #f8d7da; color: #721c24;'
        
    if row['HBD'] > 5:
        idx = row.index.get_loc('HBD')
        styles[idx] = 'background-color: #f8d7da; color: #721c24;'
        
    if row['HBA'] > 10:
        idx = row.index.get_loc('HBA')
        styles[idx] = 'background-color: #f8d7da; color: #721c24;'
        
    return styles

# --- 1. ΤΙΤΛΟΣ & ΠΕΡΙΓΡΑΦΗ ---
st.title("🧪 Lipinski Rule of 5 Checker")
st.write("Ελέγξτε τη φαρμακολογική καταλληλότητα μορίων αυτόματα μέσω PubChem & RDKit.")

# --- 2. ΕΙΣΑΓΩΓΗ ΔΕΔΟΜΕΝΩΝ ΑΠΟ ΧΡΗΣΤΗ (3 ΕΠΙΛΟΓΕΣ) ---
input_option = st.radio(
    "Επιλέξτε τρόπο εισαγωγής μορίων:", 
    ("Δοκιμή με έτοιμα ονόματα", "Πληκτρολόγηση ονομάτων", "Ανέβασμα αρχείου")
)

names_list = []

if input_option == "Δοκιμή με έτοιμα ονόματα":
    names_list = ["Aspirin", "Cyclosporin A", "Ibuprofen", "Vancomycin", "Omeprazole"]
    st.write("Δοκιμαστική λίστα:", ", ".join(names_list))

elif input_option == "Πληκτρολόγηση ονομάτων":
    user_input = st.text_area("Γράψτε τα ονόματα (ένα ανά γραμμή):", 
                              value="Aspirin\nParacetamol\nCaffeine")
    names_list = [name.strip() for name in user_input.split("\n") if name.strip()]

else:  # Ανέβασμα αρχείου
    uploaded_file = st.file_uploader("Ανεβάστε αρχείο (CSV με στήλη 'Name' ή TXT με ένα όνομα ανά γραμμή):", type=["csv", "txt"])
    
    if uploaded_file is not None:
        try:
            if uploaded_file.name.endswith('.csv'):
                df_upload = pd.read_csv(uploaded_file)
                # Αν υπάρχει στήλη 'Name' την παίρνουμε, αλλιώς παίρνουμε την πρώτη στήλη
                col_name = 'Name' if 'Name' in df_upload.columns else df_upload.columns[0]
                names_list = df_upload[col_name].dropna().astype(str).tolist()
            else:  # Αρχείο TXT
                stringdata = uploaded_file.getvalue().decode("utf-8")
                names_list = [line.strip() for line in stringdata.splitlines() if line.strip()]
            
            st.success(f"Επιτυχής φόρτωση {len(names_list)} ονομάτων από το αρχείο!")
        except Exception as e:
            st.error(f"Σφάλμα κατά την ανάγνωση του αρχείου: {e}")

# --- 3. ΥΠΟΛΟΓΙΣΜΟΣ ΚΑΝΟΝΩΝ ---
if st.button("Αναζήτηση & Υπολογισμός"):
    if not names_list:
        st.warning("Παρακαλώ εισάγετε ή ανεβάστε τουλάχιστον ένα όνομα.")
    else:
        results = []
        progress_bar = st.progress(0)
        
        for i, name in enumerate(names_list):
            smi = get_smiles_from_pubchem(name)
            
            if smi is None:
                st.warning(f"Δεν βρέθηκε η ένωση: **{name}** στο PubChem. Παράλειψη.")
                continue
            
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            
            # Υπολογισμοί Lipinski
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
            
            df_display = df_results.drop(columns=["SMILES", "Molecule"])
            styled_df = df_display.style.apply(highlight_lipinski, axis=1)
            
            st.dataframe(styled_df, use_container_width=True)
            
            # Δυνατότητα λήψης των αποτελεσμάτων σε CSV
            csv_data = df_display.to_csv(index=False).encode('utf-8')
            st.download_button(
                label="📥 Λήψη αποτελεσμάτων ως CSV",
                data=csv_data,
                file_name='lipinski_results.csv',
                mime='text/csv',
            )
            
            st.subheader("🔬 2D Δομές Μορίων")
            for index, row in df_results.iterrows():
                col1, col2 = st.columns([1, 2])
                with col1:
                    st.write(f"**{row['Name']}**")
                    st.write(f"Status: {row['Status']}")
                with col2:
                    img = Draw.MolToImage(row["Molecule"], size=(300, 150))
                    st.image(img)

