import io
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
import streamlit as st
from components import render_navbar

# st.set_page_config πρέπει να είναι πάντα πρώτη (αν δεν υπάρχει στο κεντρικό app.py για τα pages)
# st.set_page_config(page_title="Functional Groups Analyzer", page_icon="🧪", layout="wide")

render_navbar()

# --- ΕΠΙΚΕΦΑΛΙΔΑ & ΠΕΡΙΓΡΑΦΗ ---
st.title("🧪 Functional Groups & Substructure Analyzer")
st.write(
    "Αυτόματη ανίχνευση και καταμέτρηση βασικών λειτουργικών ομάδων σε οργανικά"
    " μόρια μέσω SMARTS patterns."
)
st.markdown("---")

# --- INPUT METHOD ---
input_type = st.radio(
    "Επιλέξτε τρόπο εισαγωγής:", ("SMILES", "Όνομα Μορίου (από PubChem)")
)

mol = None
smiles_input = ""

if input_type == "SMILES":
  smiles_input = st.text_input(
      "Εισάγετε SMILES:", "CC(=O)Oc1ccccc1C(=O)O"
  )  # Aspirin
  if smiles_input:
    mol = Chem.MolFromSmiles(smiles_input)
    if not mol:
      st.error("Μη έγκυρο SMILES.")
else:
  name_input = st.text_input("Όνομα μορίου:", "Paracetamol")
  if name_input:
    try:
      compounds = pcp.get_compounds(name_input, "name")
      if compounds:
        smiles_input = compounds[0].canonical_smiles
        mol = Chem.MolFromSmiles(smiles_input)
        st.success(f"Βρέθηκε: {name_input} (`{smiles_input}`)")
      else:
        st.error("Δεν βρέθηκε το μόριο.")
    except Exception as e:
      st.error(f"Σφάλμα: {e}")

if mol:
  st.markdown("---")
  col1, col2 = st.columns([1, 1])

  with col1:
    st.subheader("🖼️ Μοριακή Δομή")
    img = Draw.MolToImage(mol, size=(400, 400))
    st.image(img, use_container_width=True)

  with col2:
    st.subheader("🔍 Ανίχνευση Ομάδων")

    # Ορισμός βασικών SMARTS patterns για λειτουργικές ομάδες
    patterns = {
        "Υδροξυλική Ομάδα (-OH αλκοόλης)": "[OX2H][CX4]",
        "Φαινολικό Υδροξύλιο (-OH φαινόλης)": "[OX2H]c",
        "Καρβοξυλικό Οξύ (-COOH)": "C(=O)[OX2H1]",
        "Εστέρας (-COO-)": "C(=O)O[CX4]",
        "Αμίδιο (-C(=O)N-)": "C(=O)[NX3]",
        "Πρωτοταγής Αμίνη (-NH2)": "[NX3H2]",
        "Αρωματικός Δακτύλιος (Benzene-like)": "c1ccccc1",
        "Καρβονύλιο (Κετόνη/Αλδεΰδη C=O)": "[CX3](=O)[#6]",
        "Αλογόνα (F, Cl, Br, I)": "[F,Cl,Br,I]",
        "Νιτροομάδα (-NO2)": "[N+](=O)[O-]",
    }

    results = []
    for name, smarts in patterns.items():
      patt = Chem.MolFromSmarts(smarts)
      if patt:
        matches = mol.GetSubstructMatches(patt)
        count = len(matches)
        if count > 0:
          results.append({"Λειτουργική Ομάδα": name, "Εμφανίσεις": count})

    if results:
      df = pd.DataFrame(results)

      # --- KPI METRICS ---
      total_groups = df["Εμφανίσεις"].sum()
      unique_types = len(df)

      mcol1, mcol2 = st.columns(2)
      mcol1.metric("Μοναδικοί Τύποι Ομάδων", unique_types)
      mcol2.metric("Συνολικές Εμφανίσεις", total_groups)

      st.markdown("---")
      st.dataframe(df, use_container_width=True)

      # --- DOWNLOAD BUTTON (για SaaS αισθητική) ---
      csv_data = df.to_csv(index=False).encode("utf-8")
      st.download_button(
          label="⬇️ Λήψη αποτελεσμάτων (CSV)",
          data=csv_data,
          file_name="functional_groups.csv",
          mime="text/csv",
      )
    else:
      st.info("Δεν εντοπίστηκαν οι βασικές ομάδες της λίστας αναφοράς.")
