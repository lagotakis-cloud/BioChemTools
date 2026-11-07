import streamlit as st
from components import render_navbar
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import pandas as pd
import plotly.express as px

render_navbar()

st.set_page_config(page_title="Multi-Molecule Similarity Matrix", page_icon="📈", layout="wide")

st.title("📈 Multi-Molecule Tanimoto Similarity Matrix")
st.markdown("Υπολογισμός μήτρας ομοιότητας Tanimoto για πολλαπλά μόρια ταυτόχρονα και οπτικοποίηση με Heatmap.")

st.markdown("---")

# Προεπιλεγμένα παραδείγματα (Aspirin, Paracetamol, Ibuprofen, Caffeine, Ethanol)
default_smiles = """CC(=O)Oc1ccccc1C(=O)O, Aspirin
CC(=O)Nc1ccc(O)cc1, Paracetamol
CC(C)Cc1ccc(cc1)C(C)C(=O)O, Ibuprofen
CN1C=NC2=C1C(=O)N(C(=O)N2C)C, Caffeine
CCO, Ethanol"""

st.markdown("### Εισάγετε λίστα μορίων (Μορφή ανά γραμμή: `SMILES, Όνομα`)")
input_data = st.text_area("Λίστα:", value=default_smiles, height=150)

if input_data:
    lines = input_data.strip().split("\n")
    molecules = []
    names = []
    
    for line in lines:
        parts = line.split(",")
        if len(parts) >= 1:
            smi = parts[0].strip()
            name = parts[1].strip() if len(parts) > 1 else smi
            mol = Chem.MolFromSmiles(smi)
            if mol:
                molecules.append(mol)
                names.append(name)
                
    if len(molecules) < 2:
        st.warning("Παρακαλώ εισάγετε τουλάχιστον 2 έγκυρα μόρια για τη δημιουργία της μήτρας ομοιότητας.")
    else:
        st.markdown(f"Εντοπίστηκαν **{len(molecules)}** έγκυρα μόρια.")
        
        # Υπολογισμός Morgan/RDKit Fingerprints και Matrix Tanimoto
        fps = [FingerprintMols.FingerprintMol(m) for m in molecules]
        n = len(molecules)
        matrix = [[0.0 for _ in range(n)] for _ in range(n)]
        
        for i in range(n):
            for j in range(n):
                matrix[i][j] = round(DataStructs.TanimotoSimilarity(fps[i], fps[j]), 3)
                
        # Δημιουργία Pandas DataFrame για το Heatmap
        df_matrix = pd.DataFrame(matrix, index=names, columns=names)
        
        st.markdown("---")
        st.subheader("📊 Πίνακας Ομοιότητας Tanimoto")
        st.dataframe(df_matrix, use_container_width=True)
        
        # Plotly Heatmap (Διόρθωση σε color_continuous_scale)
        st.markdown("### 🌡️ Heatmap Οπτικοποίηση")
        fig = px.imshow(
            df_matrix,
            text_auto=True,
            color_continuous_scale="Blues",
            aspect="auto",
            labels=dict(color="Tanimoto Similarity")
        )
        fig.update_layout(title="Σύγκριση Ομοιότητας Ζεύγους Μορίων", xaxis_title="", yaxis_title="")
        st.plotly_chart(fig, use_container_width=True)
    
