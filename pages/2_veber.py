import io
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors
import streamlit as st
from components import render_navbar

render_navbar()

# --- ΕΠΙΚΕΦΑΛΙΔΑ & ΠΕΡΙΓΡΑΦΗ ---
st.title("🧬 Veber's Rules Checker (TPSA & Rotatable Bonds)")
st.write(
    "Αξιολογήστε την ευκαμψία και τη διαπερατότητα μεμβρανών/Αιματοεγκεφαλικού"
    " Φραγμού (BBB) βάσει των κανόνων Veber."
)
st.markdown("---")

# --- ΕΠΙΛΟΓΗ ΟΡΙΟΥ TPSA ---
tpsa_option = st.selectbox(
    "Επιλέξτε στόχο διαπερατότητας (Όριο TPSA):",
    (
        "Cell Membranes (TPSA ≤ 140 Å²)",
        "Blood-Brain Barrier / BBB (TPSA ≤ 90 Å²)",
    ),
)
tpsa_limit = 140 if "140" in tpsa_option else 90

# --- INITIALIZE SESSION STATE ---
if "veber_results" not in st.session_state:
  st.session_state.veber_results = None


# --- ΣΥΝΑΡΤΗΣΕΙΣ ---
def get_smiles_from_pubchem(compound_name):
  try:
    compounds = pcp.get_compounds(compound_name.strip(), "name")
    if compounds:
      return compounds[0].canonical_smiles
    else:
      return None
  except Exception as e:
    return None


def highlight_veber(row):
  styles = [""] * len(row)
  name_idx = row.index.get_loc("Όνομα / SMILES")

  if "Passed" in str(row["Κατάσταση"]):
    styles[name_idx] = (
        "background-color: #d4edda; color: #155724; font-weight: bold;"
    )
  else:
    styles[name_idx] = (
        "background-color: #f8d7da; color: #721c24; font-weight: bold;"
    )

  if row["TPSA (Å²)"] > tpsa_limit:
    styles[row.index.get_loc("TPSA (Å²)")] = (
        "background-color: #f8d7da; color: #721c24;"
    )
  if row["Περιστρ. Δεσμοί (RB)"] > 10:
    styles[row.index.get_loc("Περιστρ. Δεσμοί (RB)")] = (
        "background-color: #f8d7da; color: #721c24;"
    )

  return styles


# ==========================================
# ΚΑΤΑΣΤΑΣΗ 1: ΕΙΣΑΓΩΓΗ ΔΕΔΟΜΕΝΩΝ
# ==========================================
if st.session_state.veber_results is None:
  input_option = st.radio(
      "Επιλέξτε τρόπο εισαγωγής μορίων:",
      (
          "Δοκιμή με έτοιμα ονόματα",
          "Πληκτρολόγηση ονομάτων",
          "Πληκτρολόγηση SMILES",
          "Ανέβασμα αρχείου ονομάτων",
          "Ανέβασμα αρχείου SMILES",
      ),
  )

  input_data_list = []
  is_smiles_mode = False

  if input_option == "Δοκιμή με έτοιμα ονόματα":
    input_data_list = [
        "Aspirin",
        "Cyclosporin A",
        "Ibuprofen",
        "Caffeine",
        "Diazepam",
    ]
    st.write("Δοκιμαστική λίστα:", ", ".join(input_data_list))

  elif input_option == "Πληκτρολόγηση ονομάτων":
    user_input = st.text_area(
        "Γράψτε τα ονόματα (ένα ανά γραμμή):",
        value="Aspirin\nDiazepam\nVancomycin",
    )
    input_data_list = [
        name.strip() for name in user_input.split("\n") if name.strip()
    ]

  elif input_option == "Πληκτρολόγηση SMILES":
    is_smiles_mode = True
    user_input = st.text_area(
        "Γράψτε τα SMILES strings (ένα ανά γραμμή):",
        value=(
            "CC(=O)OC1=CC=CC=C1C(=O)O\nCN1C=NC2=C1C(=O)N(C(=O)N2C)C\nCC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=CC(=C3)OC"
        ),
    )
    input_data_list = [
        smi.strip() for smi in user_input.split("\n") if smi.strip()
    ]

  elif input_option == "Ανέβασμα αρχείου ονομάτων":
    uploaded_file = st.file_uploader(
        "Ανεβάστε αρχείο (CSV με στήλη 'Name' ή TXT με ένα όνομα ανά γραμμή):",
        type=["csv", "txt"],
    )
    if uploaded_file is not None:
      try:
        if uploaded_file.name.endswith(".csv"):
          df_upload = pd.read_csv(uploaded_file)
          col_name = (
              "Name" if "Name" in df_upload.columns else df_upload.columns[0]
          )
          input_data_list = (
              df_upload[col_name].dropna().astype(str).tolist()
          )
        else:
          stringdata = uploaded_file.getvalue().decode("utf-8")
          input_data_list = [
              line.strip() for line in stringdata.splitlines() if line.strip()
          ]
        st.success(
            f"Επιτυχής φόρτωση {len(input_data_list)} εγγραφών από το αρχείο!"
        )
      except Exception as e:
        st.error(f"Σφάλμα κατά την ανάγνωση του αρχείου: {e}")

  else:  # Ανέβασμα αρχείου SMILES
    is_smiles_mode = True
    uploaded_file = st.file_uploader(
        "Ανεβάστε αρχείο (CSV με στήλη 'SMILES' ή TXT με ένα SMILES ανά γραμμή):",
        type=["csv", "txt"],
    )
    if uploaded_file is not None:
      try:
        if uploaded_file.name.endswith(".csv"):
          df_upload = pd.read_csv(uploaded_file)
          col_name = (
              "SMILES"
              if "SMILES" in df_upload.columns
              else df_upload.columns[0]
          )
          input_data_list = (
              df_upload[col_name].dropna().astype(str).tolist()
          )
        else:
          stringdata = uploaded_file.getvalue().decode("utf-8")
          input_data_list = [
              line.strip() for line in stringdata.splitlines() if line.strip()
          ]
        st.success(
            f"Επιτυχής φόρτωση {len(input_data_list)} SMILES από το αρχείο!"
        )
      except Exception as e:
        st.error(f"Σφάλμα κατά την ανάγνωση του αρχείου: {e}")

  if st.button("Αναζήτηση & Υπολογισμός Veber", type="primary"):
    if not input_data_list:
      st.warning("Παρακαλώ εισάγετε ή ανεβάστε δεδομένα.")
    else:
      results = []
      progress_bar = st.progress(0)

      for i, item in enumerate(input_data_list):
        if is_smiles_mode:
          smi = item
          display_name = item
        else:
          display_name = item
          smi = get_smiles_from_pubchem(item)

        if smi is None:
          continue

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
          continue

        tpsa = Descriptors.TPSA(mol)
        rot_bonds = Descriptors.NumRotatableBonds(mol)

        violations = 0
        if tpsa > tpsa_limit:
          violations += 1
        if rot_bonds > 10:
          violations += 1

        status = "Passed (Veber-like)" if violations == 0 else "Failed"

        results.append({
            "Όνομα / SMILES": display_name,
            "SMILES": smi,
            "TPSA (Å²)": round(tpsa, 2),
            "Περιστρ. Δεσμοί (RB)": rot_bonds,
            "Όριο TPSA": tpsa_limit,
            "Παραβάσεις": violations,
            "Κατάσταση": status,
        })

        progress_bar.progress((i + 1) / len(input_data_list))

      if results:
        st.session_state.veber_results = results
        st.rerun()


# ==========================================
# ΚΑΤΑΣΤΑΣΗ 2: ΑΠΟΤΕΛΕΣΜΑΤΑ & ΛΗΨΗ
# ==========================================
else:
  if st.button("🔄 Νέα Αναζήτηση / Αλλαγή Δεδομένων"):
    st.session_state.veber_results = None
    st.rerun()

  df_results = pd.DataFrame(st.session_state.veber_results)

  st.subheader("📊 Αποτελέσματα Ανάλυσης Veber")
  st.info(
      f"Ενεργό κριτήριο TPSA: **{'≤ 140 Å² (Cell Membranes)' if tpsa_limit == 140 else '≤ 90 Å² (BBB)'}**"
  )

  # Αφαίρεση βοηθητικών στηλών
  df_display = df_results.drop(columns=["SMILES", "Όριο TPSA"])

  # Εφαρμογή styling και μορφοποίηση
  styled_df = df_display.style.apply(highlight_veber, axis=1).format(
      {"TPSA (Å²)": "{:.2f}"}
  )

  st.dataframe(styled_df, use_container_width=True)

  # --- ΠΟΛΛΑΠΛΕΣ ΕΠΙΛΟΓΕΣ ΛΗΨΗΣ ΑΡΧΕΙΟΥ ---
  st.markdown("### 📥 Λήψη Αποτελεσμάτων Veber")
  file_format = st.radio(
      "Επιλέξτε μορφή αρχείου για λήψη:",
      options=["CSV (.csv)", "Text (.txt)", "Excel (.xlsx)"],
      horizontal=True,
      key="download_format_radio_veber",
  )

  if file_format == "CSV (.csv)":
    file_data = df_display.to_csv(index=False).encode("utf-8")
    file_name = "veber_results.csv"
    mime_type = "text/csv"

  elif file_format == "Text (.txt)":
    file_data = df_display.to_csv(sep="\t", index=False).encode("utf-8")
    file_name = "veber_results.txt"
    mime_type = "text/plain"

  elif file_format == "Excel (.xlsx)":
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
      df_display.to_excel(writer, index=False, sheet_name="Veber_Rules")
    file_data = buffer.getvalue()
    file_name = "veber_results.xlsx"
    mime_type = (
        "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    )

  st.download_button(
      label=f"⬇️ Λήψη αρχείου ({file_name})",
      data=file_data,
      file_name=file_name,
      mime=mime_type,
      type="primary",
  )


st.markdown("""
## Veber Rules Calculator

### What are Veber's Rules?

Veber's Rules are widely used guidelines in drug discovery for evaluating the physicochemical properties of small molecules that may influence their oral bioavailability and membrane permeability.

The rules were proposed by Daniel F. Veber and colleagues, based on studies investigating the relationship between molecular properties and oral bioavailability.

According to the original Veber criteria, compounds are more likely to exhibit favorable oral bioavailability when they meet the following conditions:

- **Rotatable Bonds (RB):** ≤ 10
- **Topological Polar Surface Area (TPSA):** ≤ 140 Å²

Alternatively, the original study also identified a relationship between oral bioavailability and the sum of hydrogen bond donors and acceptors, with a total of ≤ 12.

These parameters help assess molecular flexibility and polarity, two important factors that can influence a compound's ability to cross biological membranes.

### How Does This Tool Work?

The **Veber Rules Calculator** allows users to evaluate chemical compounds based on their Topological Polar Surface Area (TPSA), using two different assessment modes:

**1. Cell Membranes (TPSA ≤ 140 Å²)**

This option evaluates compounds against a TPSA threshold of 140 Å².

The threshold is based on Veber's original findings, which associated lower molecular polarity with improved oral bioavailability and the ability to permeate biological membranes.

Compounds with TPSA values ≤ 140 Å² meet this particular screening criterion.

**2. Blood-Brain Barrier (BBB) (TPSA ≤ 90 Å²)**

This option evaluates compounds using a more restrictive TPSA threshold of 90 Å².

Lower TPSA values are generally associated with a greater potential for passive penetration across the blood-brain barrier (BBB), which separates the bloodstream from the central nervous system.

Compounds with TPSA values ≤ 90 Å² meet this particular screening criterion.

### What Does the Tool Calculate?

For each successfully processed compound, the calculator determines its TPSA and compares the calculated value against the threshold corresponding to the selected assessment mode.

The results allow users to:

- Review the calculated TPSA value for each compound.
- Assess whether the compound meets the selected TPSA threshold.
- Compare compounds based on their molecular polarity.
- Perform an initial screening for membrane permeability or potential BBB penetration.

### Scientific Interpretation and Limitations

TPSA is an important molecular descriptor associated with membrane permeability and drug absorption. However, it should not be considered a standalone predictor of biological membrane penetration.

Meeting the TPSA ≤ 140 Å² criterion does not guarantee favorable oral bioavailability or membrane permeability.

Similarly, meeting the TPSA ≤ 90 Å² criterion does not guarantee blood-brain barrier penetration. BBB permeability is influenced by several additional factors, including molecular size, lipophilicity, hydrogen bonding, ionization state, transporter interactions, and active efflux mechanisms.

The two thresholds implemented in this tool should therefore be interpreted as preliminary screening criteria rather than definitive predictors of biological activity or pharmacokinetic behavior.

### Scientific Reference

Veber, D. F., Johnson, S. R., Cheng, H.-Y., Smith, B. R., Ward, K. W., & Kopple, K. D. (2002). Molecular properties that influence the oral bioavailability of drug candidates. *Journal of Medicinal Chemistry, 45*(12), 2615–2623.
""")
