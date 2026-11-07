import streamlit as st
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
from rdkit.Chem import Draw

# --- 1. ΔΙΑΧΕΙΡΙΣΗ STATE & ΓΛΩΣΣΑΣ (Default: EN) ---
if 'lang' not in st.session_state:
    st.session_state.lang = 'EN'

# --- 2. ΛΕΞΙΚΟ ΜΕΤΑΦΡΑΣΕΩΝ ---
LANG = {
    "EN": {
        "title": "🧪 Lipinski Rule of 5 Checker",
        "description": "Assess molecular drug-likeness automatically via PubChem & RDKit.",
        "input_label": "Select molecule input method:",
        "opt_preset": "Try sample compounds",
        "opt_manual": "Type compound names",
        "opt_upload": "Upload file",
        "preset_info": "Sample list:",
        "manual_prompt": "Enter compound names (one per line):",
        "upload_prompt": "Upload file (CSV with 'Name' column or TXT with one name per line):",
        "upload_success": "Successfully loaded {} compounds from file!",
        "upload_error": "Error reading file: {}",
        "btn_calc": "Search & Calculate",
        "warn_no_input": "Please enter or upload at least one compound name.",
        "warn_not_found": "Compound **{}** was not found in PubChem. Skipping.",
        "res_header": "📊 Analysis Results",
        "btn_download": "📥 Download results as CSV",
        "struct_header": "🔬 2D Molecular Structures",
        "col_name": "Name",
        "col_mw": "Molecular Weight (MW)",
        "col_logp": "LogP",
        "col_hbd": "H-Bond Donors (HBD)",
        "col_hba": "H-Bond Acceptors (HBA)",
        "col_violations": "Violations",
        "col_status": "Status",
        "status_pass": "Passed (Drug-like)",
        "status_fail": "Failed"
    },
    "EL": {
        "title": "🧪 Lipinski Rule of 5 Checker",
        "description": "Ελέγξτε τη φαρμακολογική καταλληλότητα μορίων αυτόματα μέσω PubChem & RDKit.",
        "input_label": "Επιλέξτε τρόπο εισαγωγής μορίων:",
        "opt_preset": "Δοκιμή με έτοιμα ονόματα",
        "opt_manual": "Πληκτρολόγηση ονομάτων",
        "opt_upload": "Ανέβασμα αρχείου",
        "preset_info": "Δοκιμαστική λίστα:",
        "manual_prompt": "Γράψτε τα ονόματα (ένα ανά γραμμή):",
        "upload_prompt": "Ανεβάστε αρχείο (CSV με στήλη 'Name' ή TXT με ένα όνομα ανά γραμμή):",
        "upload_success": "Επιτυχής φόρτωση {} ονομάτων από το αρχείο!",
        "upload_error": "Σφάλμα κατά την ανάγνωση του αρχείου: {}",
        "btn_calc": "Αναζήτηση & Υπολογισμός",
        "warn_no_input": "Παρακαλώ εισάγετε ή ανεβάστε τουλάχιστον ένα όνομα.",
        "warn_not_found": "Δεν βρέθηκε η ένωση: **{}** στο PubChem. Παράλειψη.",
        "res_header": "📊 Αποτελέσματα Ανάλυσης",
        "btn_download": "📥 Λήψη αποτελεσμάτων ως CSV",
        "struct_header": "🔬 2D Δομές Μορίων",
        "col_name": "Όνομα",
        "col_mw": "Μοριακό Βάρος (MW)",
        "col_logp": "LogP",
        "col_hbd": "Δότες H (HBD)",
        "col_hba": "Δέκτες H (HBA)",
        "col_violations": "Παραβάσεις",
        "col_status": "Κατάσταση",
        "status_pass": "Passed (Drug-like)",
        "status_fail": "Failed"
    }
}

t = LANG[st.session_state.lang]

# --- 3. ΕΠΙΚΕΦΑΛΙΔΑ & ΚΟΥΜΠΙΑ ΣΗΜΑΙΩΝ ΜΕ ΕΙΚΟΝΕΣ ---
header_col1, header_col2 = st.columns([4, 2])

with header_col1:
    st.title(t["title"])

with header_col2:
    st.write("") # Μικρό κενό για ευθυγράμμιση
    # Δημιουργούμε δύο στήλες για τις σημαίες
    flag_col1, flag_col2 = st.columns(2)
    
    with flag_col1:
        # Εμφάνιση μικρής αγγλικής σημαίας και κουμπί από κάτω
        st.image("https://flagcdn.com/w40/gb.png", width=30)
        if st.button("EN", key="btn_en"):
            st.session_state.lang = 'EN'
            st.rerun()
            
    with flag_col2:
        # Εμφάνιση μικρής ελληνικής σημαίας και κουμπί από κάτω
        st.image("https://flagcdn.com/w40/gr.png", width=30)
        if st.button("EL", key="btn_el"):
            st.session_state.lang = 'EL'
            st.rerun()

st.write(t["description"])
st.markdown("---")

# --- ΣΥΝΑΡΤΗΣΕΙΣ ---
def get_smiles_from_pubchem(compound_name):
    try:
        compounds = pcp.get_compounds(compound_name.strip(), 'name')
        if compounds:
            return compounds[0].canonical_smiles
        else:
            return None
    except Exception as e:
        return None

def highlight_lipinski(row):
    styles = [''] * len(row)
    name_idx = row.index.get_loc(t["col_name"])
    
    if t["status_pass"] in str(row[t["col_status"]]):
        styles[name_idx] = 'background-color: #d4edda; color: #155724; font-weight: bold;'
    else:
        styles[name_idx] = 'background-color: #f8d7da; color: #721c24; font-weight: bold;'
        
    if row[t["col_mw"]] > 500:
        idx = row.index.get_loc(t["col_mw"])
        styles[idx] = 'background-color: #f8d7da; color: #721c24;'
        
    if row[t["col_logp"]] > 5:
        idx = row.index.get_loc(t["col_logp"])
        styles[idx] = 'background-color: #f8d7da; color: #721c24;'
        
    if row[t["col_hbd"]] > 5:
        idx = row.index.get_loc(t["col_hbd"])
        styles[idx] = 'background-color: #f8d7da; color: #721c24;'
        
    if row[t["col_hba"]] > 10:
        idx = row.index.get_loc(t["col_hba"])
        styles[idx] = 'background-color: #f8d7da; color: #721c24;'
        
    return styles

# --- 4. ΕΙΣΑΓΩΓΗ ΔΕΔΟΜΕΝΩΝ ---
input_option = st.radio(
    t["input_label"], 
    (t["opt_preset"], t["opt_manual"], t["opt_upload"])
)

names_list = []

if input_option == t["opt_preset"]:
    names_list = ["Aspirin", "Cyclosporin A", "Ibuprofen", "Vancomycin", "Omeprazole"]
    st.write(t["preset_info"], ", ".join(names_list))

elif input_option == t["opt_manual"]:
    user_input = st.text_area(t["manual_prompt"], value="Aspirin\nSertraline\nVancomycin")
    names_list = [name.strip() for name in user_input.split("\n") if name.strip()]

else:
    uploaded_file = st.file_uploader(t["upload_prompt"], type=["csv", "txt"])
    if uploaded_file is not None:
        try:
            if uploaded_file.name.endswith('.csv'):
                df_upload = pd.read_csv(uploaded_file)
                col_name = 'Name' if 'Name' in df_upload.columns else df_upload.columns[0]
                names_list = df_upload[col_name].dropna().astype(str).tolist()
            else:
                stringdata = uploaded_file.getvalue().decode("utf-8")
                names_list = [line.strip() for line in stringdata.splitlines() if line.strip()]
            
            st.success(t["upload_success"].format(len(names_list)))
        except Exception as e:
            st.error(t["upload_error"].format(e))

# --- 5. ΥΠΟΛΟΓΙΣΜΟΙ & ΠΡΟΒΟΛΗ ---
if st.button(t["btn_calc"]):
    if not names_list:
        st.warning(t["warn_no_input"])
    else:
        results = []
        progress_bar = st.progress(0)
        
        for i, name in enumerate(names_list):
            smi = get_smiles_from_pubchem(name)
            
            if smi is None:
                st.warning(t["warn_not_found"].format(name))
                continue
            
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            
            violations = 0
            if mw > 500: violations += 1
            if logp > 5: violations += 1
            if hbd > 5: violations += 1
            if hba > 10: violations += 1
            
            status = t["status_pass"] if violations <= 1 else t["status_fail"]
            
            results.append({
                t["col_name"]: name,
                "SMILES": smi,
                t["col_mw"]: round(mw, 2),
                t["col_logp"]: round(logp, 2),
                t["col_hbd"]: hbd,
                "HBA_val": hba,
                t["col_hba"]: hba,
                t["col_violations"]: violations,
                t["col_status"]: status,
                "Molecule": mol
            })
            
            progress_bar.progress((i + 1) / len(names_list))
        
        if results:
            df_results = pd.DataFrame(results)
            
            st.subheader(t["res_header"])
            
            df_display = df_results.drop(columns=["SMILES", "Molecule", "HBA_val"])
            styled_df = df_display.style.apply(highlight_lipinski, axis=1)
            
            st.dataframe(styled_df, use_container_width=True)
            
            csv_data = df_display.to_csv(index=False).encode('utf-8')
            st.download_button(
                label=t["btn_download"],
                data=csv_data,
                file_name='lipinski_results.csv',
                mime='text/csv',
            )
            
            st.subheader(t["struct_header"])
            for index, row in df_results.iterrows():
                col1, col2 = st.columns([1, 2])
                with col1:
                    st.write(f"**{row[t['col_name']]}**")
                    st.write(f"{t['col_status']}: {row[t['col_status']]}")
                with col2:
                    img = Draw.MolToImage(row["Molecule"], size=(300, 150))
                    st.image(img)
