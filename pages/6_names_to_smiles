import io
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import streamlit as st

st.set_page_config(
    page_title="Name ⇄ SMILES Converter", page_icon="🔤", layout="wide"
)

from components import render_navbar

render_navbar()

st.title("🔤 Name ⇄ SMILES Bidirectional & Batch Converter")
st.markdown(
    "Αμφίδρομη μετατροπή μεταξύ Ονομάτων Χημικών Ενώσεων και SMILES strings,"
    " με δυνατότητα μαζικής επεξεργασίας (batch) μέσω αρχείου ή λίστας."
)
st.markdown("---")

if "name_smiles_results" not in st.session_state:
  st.session_state.name_smiles_results = None

if st.session_state.name_smiles_results is None:
  with st.form("name_smiles_form"):
    conversion_mode = st.radio(
        "Κατεύθυνση Μετατροπής:",
        (
            "Όνομα Μορίου ➔ SMILES (PubChem)",
            "SMILES ➔ Canonical SMILES & Έλεγχος",
            "Μαζική Επεξεργασία (Batch από Αρχείο CSV/TXT)",
        ),
    )

    single_input = st.text_input(
        "Εισάγετε τιμή (για μεμονωμένη μετατροπή):", "Paracetamol"
    )
    uploaded_file = st.file_uploader(
        "Ή ανεβάστε αρχείο για batch μετατροπή (.csv ή .txt με 1 στήλη):",
        type=["csv", "txt"],
    )

    submitted = st.form_submit_button("Εκτέλεση Μετατροπής", type="primary")

  if submitted:
    results_list = []

    # Περίπτωση 1: Batch Upload
    if uploaded_file is not None:
      try:
        if uploaded_file.name.endswith(".csv"):
          df_up = pd.read_csv(uploaded_file)
          items = df_up.iloc[:, 0].dropna().astype(str).tolist()
        else:
          stringdata = uploaded_file.getvalue().decode("utf-8")
          items = [
              line.strip() for line in stringdata.splitlines() if line.strip()
          ]

        for item in items:
          if "Όνομα" in conversion_mode or "Όνομα Μορίου" in conversion_mode:
            try:
              comps = pcp.get_compounds(item, "name")
              smi = comps[0].canonical_smiles if comps else "Not Found"
              results_list.append({"Input": item, "Converted_SMILES": smi})
            except:
              results_list.append({"Input": item, "Converted_SMILES": "Error"})
          else:
            mol = Chem.MolFromSmiles(item)
            smi = Chem.MolToSmiles(mol) if mol else "Invalid SMILES"
            results_list.append({"Input": item, "Canonical_SMILES": smi})

        df_res = pd.DataFrame(results_list)
        st.session_state.name_smiles_results = {
            "type": "batch",
            "df": df_res,
        }
        st.rerun()
      except Exception as e:
        st.error(f"Σφάλμα ανάγνωσης αρχείου: {e}")

    # Περίπτωση 2: Single Input
    else:
      if not single_input.strip():
        st.warning("Παρακαλώ εισάγετε μια τιμή.")
      else:
        if "Όνομα" in conversion_mode:
          try:
            comps = pcp.get_compounds(single_input.strip(), "name")
            if comps:
              smi = comps[0].canonical_smiles
              mol = Chem.MolFromSmiles(smi)
              st.session_state.name_smiles_results = {
                  "type": "single",
                  "input": single_input,
                  "output": smi,
                  "mol": mol,
                  "mode": "name_to_smi",
              }
              st.rerun()
            else:
              st.error("Δεν βρέθηκε το μόριο στην PubChem.")
          except Exception as e:
            st.error(f"Σφάλμα σύνδεσης: {e}")
        else:
          mol = Chem.MolFromSmiles(single_input.strip())
          if mol:
            smi = Chem.MolToSmiles(mol)
            st.session_state.name_smiles_results = {
                "type": "single",
                "input": single_input,
                "output": smi,
                "mol": mol,
                "mode": "smi_to_smi",
            }
            st.rerun()
          else:
            st.error("Μη έγκυρο SMILES string.")

else:
  if st.button("🔄 Νέα Μετατροπή"):
    st.session_state.name_smiles_results = None
    st.rerun()

  res = st.session_state.name_smiles_results

  if res["type"] == "single":
    col1, col2 = st.columns([1, 2])
    with col1:
      st.subheader("🖼️ Δομή Μορίου")
      if res["mol"]:
        img = Draw.MolToImage(res["mol"], size=(300, 300))
        st.image(img, use_container_width=True)
    with col2:
      st.subheader("📊 Αποτελέσματα Μετατροπής")
      st.text_input("Αρχική Είσοδος", value=res["input"])
      st.text_input("Παραγόμενο SMILES", value=res["output"])

      csv_data = pd.DataFrame(
          [{"Input": res["input"], "Output": res["output"]}]
      ).to_csv(index=False)
      st.download_button(
          "⬇️ Λήψη αποτελέσματος (CSV)",
          data=csv_data,
          file_name="conversion_result.csv",
          mime="text/csv",
          type="primary",
      )
  else:
    st.subheader("📊 Αποτελέσματα Μαζικής Μετατροπής (Batch)")
    df = res["df"]
    st.dataframe(df, use_container_width=True)

    csv_data = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        "⬇️ Λήψη ολόκληρης λίστας (CSV)",
        data=csv_data,
        file_name="batch_conversion_results.csv",
        mime="text/csv",
        type="primary",
    )
