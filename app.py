import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
import pandas as pd

# Ρύθμιση σελίδας
st.set_page_config(page_title="LC-MS Adduct & Mass Pro", page_icon="🧪", layout="wide")

st.title("🧪 LC-MS Adduct & Mass Pro")
st.markdown("Επαγγελματικό εργαλείο ταυτοποίησης φασματοσκοπίας μάζας για αναλυτικούς χημικούς.")

# Πλαϊνή μπάρα (Sidebar) για επιλογές
st.sidebar.header("⚙️ Ρυθμίσεις Αναζήτησης")
input_type = st.sidebar.radio("Τρόπος Εισόδου:", ["Όνομα Μορίου (π.χ. Caffeine)", "SMILES String"])

molecule_input = ""
if input_type == "Όνομα Μορίου (π.χ. Caffeine)":
    molecule_input = st.sidebar.text_input("Εισάγετε όνομα:", "Caffeine")
else:
    molecule_input = st.sidebar.text_input("Εισάγετε SMILES:", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

if molecule_input:
    try:
        # Εύρεση SMILES και δεδομένων από PubChem αν δόθηκε όνομα
        if input_type == "Όνομα Μορίου (π.χ. Caffeine)":
            compound = pcp.get_compounds(molecule_input, 'name')[0]
            smiles = compound.canonical_smiles
            compound_name = molecule_input
        else:
            smiles = molecule_input
            compound_name = "Custom Molecule"

        # Δημιουργία μορίου RDKit
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            st.error("⚠️ Μη έγκυρη δομή μορίου. Ελέγξτε την είσοδο.")
        else:
            # Υπολογισμός Μονοϊσοτοπικής Μάζας (Exact Mass / Monoisotopic Mass)
            exact_mass = Descriptors.ExactMolWt(mol)
            
            # Κεντρικό Layout
            col1, col2 = st.columns([1, 2])
            
            with col1:
                st.subheader("🧬 Δομή & Βασικές Ιδιότητες")
                st.image(Draw.MolToImage(mol, size=(250, 250)), caption=compound_name)
                st.metric("Exact Monoisotopic Mass", f"{exact_mass:.4f} Da")
                st.metric("Canonical SMILES", smiles)

            with col2:
                st.subheader("📊 Πίνακας LC-MS Adducts")
                st.markdown("Αυτόματος υπολογισμός λόγου μάζας προς φορτίο ($m/z$) για δημοφιλή adducts:")
                
                # Σταθερές μαζών ιόντων / adducts
                mass_electron = 0.00054858
                mass_h = 1.007825
                mass_na = 22.989218
                mass_nh4 = 18.034374
                mass_k = 38.963158
                mass_cl = 34.968853
                mass_hcoo = 44.997655
                mass_ch3coo = 59.013304

                # Υπολογισμός Positive Adducts
                pos_adducts = [
                    {"Adduct": "[M+H]+", "Formula": "M + H", "m/z": exact_mass + mass_h},
                    {"Adduct": "[M+Na]+", "Formula": "M + Na", "m/z": exact_mass + mass_na},
                    {"Adduct": "[M+NH4]+", "Formula": "M + NH4", "m/z": exact_mass + mass_nh4},
                    {"Adduct": "[M+K]+", "Formula": "M + K", "m/z": exact_mass + mass_k},
                    {"Adduct": "[2M+H]+", "Formula": "2M + H", "m/z": (2 * exact_mass) + mass_h},
                ]
                
                # Υπολογισμός Negative Adducts
                neg_adducts = [
                    {"Adduct": "[M-H]-", "Formula": "M - H", "m/z": exact_mass - mass_h},
                    {"Adduct": "[M+Cl]-", "Formula": "M + Cl", "m/z": exact_mass + mass_cl},
                    {"Adduct": "[M+HCOO]-", "Formula": "M + Formate", "m/z": exact_mass + mass_hcoo},
                    {"Adduct": "[M+CH3COO]-", "Formula": "M + Acetate", "m/z": exact_mass + mass_ch3coo},
                ]

                st.markdown("##### ➕ Positive Ionization Mode")
                df_pos = pd.DataFrame(pos_adducts)
                df_pos["m/z"] = df_pos["m/z"].round(4)
                st.dataframe(df_pos, use_container_width=True, hide_index=True)

                st.markdown("##### ➖ Negative Ionization Mode")
                df_neg = pd.DataFrame(neg_adducts)
                df_neg["m/z"] = df_neg["m/z"].round(4)
                st.dataframe(df_neg, use_container_width=True, hide_index=True)

    except Exception as e:
        st.error(f"⚠️ Δεν βρέθηκε το μόριο ή προέκυψε σφάλμα: {e}")
