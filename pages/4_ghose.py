import io
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors
import streamlit as st

st.set_page_config(
    page_title="Ghose Filter Evaluator", page_icon="🧪", layout="wide"
)

from components import render_navbar

render_navbar()

st.title("🧪 Ghose Filter Drug-Likeness Screener")
st.markdown(
    "Έλεγχος συμμόρφωσης ενός μορίου με το κλασικό φίλτρο Ghose (MW, LogP,"
    " βαριά άτομα, Molar Refractivity)."
)
st.markdown("---")

if "ghose_results" not in st.session_state:
  st.session_state.ghose_results = None

if st.session_state.ghose_results is None:
  with st.form("ghose_input_form"):
    input_option = st.radio(
        "Επιλέξτε τρόπο εισαγωγής:",
        ("SMILES", "Όνομα Μορίου (από PubChem)"),
    )
    smiles_input = ""
    name_input = ""

    if input_option == "SMILES":
      smiles_input = st.text_input("Εισάγετε SMILES:", "CC(=O)Oc1ccccc1C(=O)O")
    else:
      name_input = st.text_input("Όνομα μορίου:", "Aspirin")

    submitted = st.form_submit_button("Έλεγχος Ghose Filter", type="primary")

  if submitted:
    mol = None
    resolved_smiles = ""
    molecule_name = "molecule"

    if input_option == "SMILES":
      resolved_smiles = smiles_input.strip()
      mol = Chem.MolFromSmiles(resolved_smiles)
    else:
      try:
        comps = pcp.get_compounds(name_input.strip(), "name")
        if comps:
          resolved_smiles = comps[0].canonical_smiles
          mol = Chem.MolFromSmiles(resolved_smiles)
          molecule_name = name_input.strip().lower().replace(" ", "_")
      except:
        mol = None

    if mol is None:
      st.error("Μη έγκυρο μόριο ή δεν βρέθηκε.")
    else:
      mw = Descriptors.MolWt(mol)
      logp = Crippen.MolLogP(mol)
      num_atoms = mol.GetNumHeavyAtoms()
      mr = Crippen.MolMR(mol)

      passed = (
          (160 <= mw <= 480)
          and (-0.4 <= logp <= 5.6)
          and (20 <= num_atoms <= 70)
          and (40 <= mr <= 130)
      )

      st.session_state.ghose_results = {
          "molecule_name": molecule_name,
          "mw": round(mw, 2),
          "logp": round(logp, 2),
          "num_atoms": num_atoms,
          "mr": round(mr, 2),
          "passed": passed,
          "mol": mol,
      }
      st.rerun()
else:
  if st.button("🔄 Νέος Έλεγχος"):
    st.session_state.ghose_results = None
    st.rerun()

  res = st.session_state.ghose_results
  passed = res["passed"]
  mol = res["mol"]

  col1, col2, col3 = st.columns(3)
  col1.metric("Ghose Status", "Passed ✅" if passed else "Failed ❌")
  col2.metric("Μοριακό Βάρος", res["mw"])
  col3.metric("Βαριά Άτομα", res["num_atoms"])

  st.markdown("---")
  c1, c2 = st.columns([1, 2])
  with c1:
    from rdkit.Chem import Draw

    img = Draw.MolToImage(mol, size=(300, 300))
    st.image(img, use_container_width=True)
  with c2:
    st.subheader("📊 Αναλυτικός Πίνακας Κριτηρίων Ghose")
    df_details = pd.DataFrame({
        "Παράμετρος": [
            "Μοριακό Βάρος (MW)",
            "LogP (Λιποφιλία)",
            "Βαριά Άτομα",
            "Molar Refractivity (MR)",
        ],
        "Τιμή Μορίου": [
            res["mw"],
            res["logp"],
            res["num_atoms"],
            res["mr"],
        ],
        "Αποδεκτό Εύρος": [
            "160 – 480",
            "-0.4 – 5.6",
            "20 – 70",
            "40 – 130",
        ],
    })
    st.dataframe(df_details, use_container_width=True)

    if passed:
      st.success(
          "✨ Το μόριο πληροί όλες τις προϋποθέσεις του **Ghose Filter**."
      )
    else:
      st.warning(
          "⚠️ Το μόριο αποκλίνει από τα τυπικά όρια του **Ghose Filter**."
      )

    st.download_button(
        "⬇️ Λήψη αναφοράς (CSV)",
        data=df_details.to_csv(index=False).encode("utf-8"),
        file_name="ghose_filter_report.csv",
        mime="text/csv",
        type="primary",
    )
