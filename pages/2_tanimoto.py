import streamlit as st
from components import render_navbar
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Draw
import pandas as pd
from PIL import Image
import io

render_navbar()

st.set_page_config(page_title="Tanimoto Similarity Search | BiochemTools", page_icon="🎯", layout="wide")

st.title("🎯 Tanimoto Molecular Similarity Search")
st.write("Συγκρίνετε ένα μόριο αναφοράς με μια ομάδα μορίων χρησιμοποιώντας Morgan Fingerprints (ECFP4).")

st.markdown("---")

# 1. Sidebar - Input Controls
st.sidebar.header("⚙️ Παράμετροι Εισαγωγής")

query_input = st.sidebar.text_input(
    "1. Μόριο Αναφοράς (Query SMILES):", 
    value="CC(=O)Oc1ccccc1C(=O)O"  # Aspirin default
)

st.sidebar.markdown("---")
st.sidebar.subheader("2. Λίστα Μορίων προς Σύγκριση")
targets_input = st.sidebar.text_area(
    "Εισάγετε SMILES (ένα ανά γραμμή):",
    value="CC(=O)Oc1ccccc1C(=O)O\nCC(=O)NC1=CC=C(O)C=C1\nCC(C)CC1=CC=C(C=C1)C(C)C(=O)O\nCN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    height=150
)
st.sidebar.caption("Προεπιλογή: Ασπιρίνη, Παρακεταμόλη, Ιβουπροφαίνη, Καφεΐνη")

# 2. Main Processing Logic
if st.button("🚀 Υπολογισμός Ομοιότητας", type="primary"):
    # Validate Query Molecule
    query_mol = Chem.MolFromSmiles(query_input.strip())
    
    if query_mol is None:
        st.error("❌ Το SMILES του μορίου αναφοράς δεν είναι έγκυρο. Παρακαλώ ελέγξτε την εισαγωγή σας.")
    else:
        # Calculate Query Fingerprint
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048)
        
        # Display Query Molecule
        col_q1, col_q2 = st.columns([1, 2])
        with col_q1:
            st.subheader("📌 Μόριο Αναφοράς")
            img = Draw.MolToImage(query_mol, size=(250, 250))
            st.image(img, caption=f"SMILES: {query_input}", use_container_width=True)
            
        with col_q2:
            st.subheader("📊 Αποτελέσματα Σύγκρισης")
            
            target_smiles_list = [s.strip() for s in targets_input.split('\n') if s.strip()]
            results = []

            for idx, smiles in enumerate(target_smiles_list):
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                    similarity = DataStructs.TanimotoSimilarity(query_fp, fp)
                    results.append({
                        "ID": f"Mol_{idx+1}",
                        "SMILES": smiles,
                        "Tanimoto Score": round(similarity, 4),
                        "Molecule Object": mol
                    })
                else:
                    results.append({
                        "ID": f"Mol_{idx+1}",
                        "SMILES": f"INVALID: {smiles}",
                        "Tanimoto Score": 0.0,
                        "Molecule Object": None
                    })

            # Sort by highest similarity
            results = sorted(results, key=lambda x: x["Tanimoto Score"], reverse=True)
            
            # DataFrame display
            df_results = pd.DataFrame(results)[["ID", "SMILES", "Tanimoto Score"]]
            st.dataframe(df_results, use_container_width=True)

        # 3. Grid Visualization of Top Matches
        st.markdown("---")
        st.subheader("🖼️ 2D Απεικόνιση & Σκορ Ομοιότητας")
        
        cols = st.columns(3)
        for i, item in enumerate(results):
            with cols[i % 3]:
                score = item["Tanimoto Score"]
                st.markdown(f"**{item['ID']}** — Score: `{score}`")
                
                # Visual Indicator
                st.progress(score)
                
                if item["Molecule Object"]:
                    img_target = Draw.MolToImage(item["Molecule Object"], size=(200, 200))
                    st.image(img_target, use_container_width=True)
                else:
                    st.warning("Μη έγκυρο μόριο")
