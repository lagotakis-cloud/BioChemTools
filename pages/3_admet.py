import io
import pandas as pd
import plotly.express as px
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
import streamlit as st
from components import render_navbar

render_navbar()

# --- ΕΠΙΚΕΦΑΛΙΔΑ & ΠΕΡΙΓΡΑΦΗ ---
st.title("🧪 ADMET & Toxicity Screener (Enterprise Edition)")
st.write(
    "Προηγμένη αξιολόγηση υδατικής διαλυτότητας (ESOL), μεταβολισμού CYP450,"
    " σάρωση PAINS, τοξικών alerts και αυτόματη 2D απεικόνιση δομών."
)
st.markdown("---")

# --- INITIALIZE SESSION STATE ---
if "admet_results" not in st.session_state:
  st.session_state.admet_results = None


# --- ΣΥΝΑΡΤΗΣΕΙΣ ---
def get_smiles_from_pubchem(compound_name):
  try:
    compounds = pcp.get_compounds(compound_name.strip(), "name")
    if compounds:
      return compounds[0].canonical_smiles
    return None
  except Exception:
    return None


def calculate_esol(mol):
  """Υπολογισμός υδατικής διαλυτότητας LogS (Delaney equation)."""
  logp = Descriptors.MolLogP(mol)
  mw = Descriptors.MolWt(mol)
  rotb = Descriptors.NumRotatableBonds(mol)
  heavy_atoms = mol.GetNumHeavyAtoms()
  ar_atoms = len([a for a in mol.GetAtoms() if a.GetIsAromatic()])
  ar_prop = ar_atoms / heavy_atoms if heavy_atoms > 0 else 0

  log_s = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rotb - 0.74 * ar_prop
  return log_s


def classify_solubility(log_s):
  """Κατηγοριοποίηση υδατικής διαλυτότητας βάσει LogS."""
  if log_s > 0:
    return "Highly Soluble"
  elif log_s > -2:
    return "Very Soluble"
  elif log_s > -4:
    return "Soluble"
  elif log_s > -6:
    return "Moderately Soluble"
  elif log_s > -10:
    return "Poorly Soluble"
  else:
    return "Insoluble"


def check_structural_alerts(mol):
  """Ανίχνευση γνωστών τοξικών / δραστικών υποδομών μέσω SMARTS."""
  alerts = {
      "Epoxide (Alkylating)": "[OX1]1C2C1C2",
      "Michael Acceptor": "[C,c]=C-[C,c]=[O,N,S]",
      "Nitroaromatic (Mutagenic)": "c[N+](=O)[O-]",
      "Aromatic Amine": "c[NH2]",
      "Acyl Halide": "[CX3](=[OX1])[F,Cl,Br,I]",
      "Alkyl Halide": "[CX4][Cl,Br,I]",
      "Aldehyde (Reactive)": "[CX3H1](=O)[#6]",
      "Hydrazine": "[NX3][NX3]",
  }
  found_alerts = []
  for name, smarts in alerts.items():
    pat = Chem.MolFromSmarts(smarts)
    if pat and mol.HasSubstructMatch(pat):
      found_alerts.append(name)
  return found_alerts


def check_cyp450_inhibition(mol):
  """Εκτίμηση πιθανής αναστολής ηπατικών ενζύμων CYP450 μέσω υποδομών."""
  cyp_patterns = {
      "CYP3A4 Inhibitor Motif": "c1ccccc1[NX3](C)C",  # Ενδεικτικά αρωματικά αζώτα
      "CYP2D6 Inhibitor Motif": "[NX3](C)(C)CCC",
      "CYP2C9 Inhibitor Motif": "c1ccc(cc1)C(=O)O",
  }
  inhibited = []
  for name, smarts in cyp_patterns.items():
    pat = Chem.MolFromSmarts(smarts)
    if pat and mol.HasSubstructMatch(pat):
      inhibited.append(name.split()[0])  # Κρατάει CYP3A4, κλπ.
  return inhibited if inhibited else ["Low Risk"]


def mol_to_svg(mol, size=(120, 120)):
  """Μετατροπή μορίου RDKit σε καθαρό HTML SVG string για το dataframe."""
  if mol is None:
    return ""
  try:
    # Χρήση SVG για τέλεια ευκρίνεια στα DataFrames
    from rdkit.Chem.Draw import rdMolDraw2D
    import svgwrite  # fallback or direct drawing

    drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    # Κωδικοποίηση σε data URI για απευθείας εμφάνιση σε markdown/dataframe
    import base64

    b64 = base64.b64encode(svg.encode("utf-8")).decode("utf-8")
    return f'<img src="data:image/svg+xml;base64,{b64}" width="{size[0]}" height="{size[1]}"/>'
  except Exception:
    return ""


# ==========================================
# ΚΑΤΑΣΤΑΣΗ 1: ΕΙΣΑΓΩΓΗ ΔΕΔΟΜΕΝΩΝ (ΜΕ ΦΟΡΜΑ)
# ==========================================
if st.session_state.admet_results is None:
  with st.form("admet_input_form"):
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

    user_text_input = ""
    uploaded_file = None

    if input_option == "Δοκιμή με έτοιμα ονόματα":
      st.write(
          "Δοκιμαστική λίστα:"
          " Aspirin, Paracetamol, Caffeine, Glycidol, Emodin, Nitrobenzene"
      )

    elif input_option == "Πληκτρολόγηση ονομάτων":
      user_text_input = st.text_area(
          "Γράψτε τα ονόματα (ένα ανά γραμμή):",
          value="Aspirin\nParacetamol\nCaffeine\nGlycidol\nEmodin",
      )

    elif input_option == "Πληκτρολόγηση SMILES":
      user_text_input = st.text_area(
          "Γράψτε τα SMILES strings (ένα ανά γραμμή):",
          value=(
              "CC(=O)OC1=CC=CC=C1C(=O)O\n"
              "CC(=O)NC1=CC=C(C=C1)O\n"
              "C1CO1CO\n"
              "CC1=CC2=C(C(=C1)O)C(=O)C3=C(C2=O)C=C(C=C3O)O"
          ),
      )

    elif input_option == "Ανέβασμα αρχείου ονομάτων":
      uploaded_file = st.file_uploader(
          "Ανεβάστε αρχείο (CSV ή TXT):", type=["csv", "txt"]
      )

    else:  # Ανέβασμα αρχείου SMILES
      uploaded_file = st.file_uploader(
          "Ανεβάστε αρχείο SMILES (CSV ή TXT):", type=["csv", "txt"]
      )

    submitted = st.form_submit_button(
        "Εκτέλεση Enterprise ADMET Screening", type="primary"
    )

  # --- ΕΠΕΞΕΡΓΑΣΙΑ ΜΕΤΑ ΤΗΝ ΥΠΟΒΟΛΗ ---
  if submitted:
    input_data_list = []
    is_smiles_mode = False

    if input_option == "Δοκιμή με έτοιμα ονόματα":
      input_data_list = [
          "Aspirin",
          "Paracetamol",
          "Caffeine",
          "Glycidol",
          "Emodin",
          "Nitrobenzene",
      ]

    elif input_option == "Πληκτρολόγηση ονομάτων":
      input_data_list = [
          name.strip() for name in user_text_input.split("\n") if name.strip()
      ]

    elif input_option == "Πληκτρολόγηση SMILES":
      is_smiles_mode = True
      input_data_list = [
          smi.strip() for smi in user_text_input.split("\n") if smi.strip()
      ]

    elif input_option in [
        "Ανέβασμα αρχείου ονομάτων",
        "Ανέβασμα αρχείου SMILES",
    ]:
      if "SMILES" in input_option:
        is_smiles_mode = True

      if uploaded_file is not None:
        try:
          if uploaded_file.name.endswith(".csv"):
            df_upload = pd.read_csv(uploaded_file)
            col_name = (
                df_upload.columns[0]
                if len(df_upload.columns) > 0
                else "Name"
            )
            input_data_list = (
                df_upload[col_name].dropna().astype(str).tolist()
            )
          else:
            stringdata = uploaded_file.getvalue().decode("utf-8")
            input_data_list = [
                line.strip() for line in stringdata.splitlines() if line.strip()
            ]
        except Exception as e:
          st.error(f"Σφάλμα κατά την ανάγνωση του αρχείου: {e}")

    if not input_data_list:
      st.warning("Παρακαλώ εισάγετε δεδομένα πριν την εκτέλεση.")
    else:
      results = []
      not_found = []
      progress_bar = st.progress(0)

      params = FilterCatalogParams()
      params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
      catalog = FilterCatalog(params)

      for i, item in enumerate(input_data_list):
        if is_smiles_mode:
          smi = item
          display_name = item
        else:
          display_name = item
          smi = get_smiles_from_pubchem(item)

        if smi is None:
          not_found.append(display_name)
          progress_bar.progress((i + 1) / len(input_data_list))
          continue

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
          not_found.append(display_name)
          progress_bar.progress((i + 1) / len(input_data_list))
          continue

        # Υπολογισμοί
        esol = calculate_esol(mol)
        sol_class = classify_solubility(esol)
        pains_matches = catalog.GetMatches(mol)
        pains_descriptions = [m.GetDescription() for m in pains_matches]
        alerts = check_structural_alerts(mol)
        cyp_inhibitions = check_cyp450_inhibition(mol)

        # Υπολογισμός ADMET Score (0 - 100%)
        # Ξεκινάμε από 100% και αφαιρούμε ποσοστά ανάλογα με τα προβλήματα
        score = 100
        if pains_descriptions:
          score -= 40
        score -= len(alerts) * 25
        if "Insoluble" in sol_class or "Poorly" in sol_class:
          score -= 15
        score = max(0, score)

        # 2D Structure Image HTML
        svg_html = mol_to_svg(mol)

        failure_reasons = []
        if pains_descriptions:
          failure_reasons.append(f"PAINS: {', '.join(pains_descriptions)}")
        if alerts:
          failure_reasons.append(f"Alerts: {', '.join(alerts)}")

        if failure_reasons:
          toxicity_status = " | ".join(failure_reasons)
        else:
          toxicity_status = "Safe (No Alerts)"

        results.append({
            "Δομή (2D)": svg_html,
            "Όνομα / SMILES": display_name,
            "SMILES": smi,
            "ADMET Score (%)": score,
            "LogS (ESOL)": round(esol, 2),
            "Διαλυτότητα": sol_class,
            "CYP450 Profile": ", ".join(cyp_inhibitions),
            "PAINS Match": (
                ", ".join(pains_descriptions)
                if pains_descriptions
                else "None detected"
            ),
            "Reactive Alerts": (
                ", ".join(alerts) if alerts else "None detected"
            ),
            "Συνολική Κατάσταση": toxicity_status,
        })

        progress_bar.progress((i + 1) / len(input_data_list))

      if not_found:
        st.warning(
            f"Δεν βρέθηκαν ή δεν ήταν έγκυρα τα εξής: {', '.join(not_found)}"
        )

      if results:
        st.session_state.admet_results = results
        st.rerun()

# ==========================================
# ΚΑΤΑΣΤΑΣΗ 2: ΑΠΟΤΕΛΕΣΜΑΤΑ & ΛΗΨΗ
# ==========================================
else:
  if st.button("🔄 Νέα Αναζήτηση / Αλλαγή Δεδομένων"):
    st.session_state.admet_results = None
    st.rerun()

  df_results = pd.DataFrame(st.session_state.admet_results)

  # --- KPI METRICS ---
  total_mols = len(df_results)
  safe_mols = len(
      df_results[df_results["Συνολική Κατάσταση"].str.contains("Safe")]
  )
  alert_mols = total_mols - safe_mols
  avg_score = round(df_results["ADMET Score (%)"].mean(), 1)

  col1, col2, col3, col4 = st.columns(4)
  col1.metric("Συνολικά Μόρια", total_mols)
  col2.metric("Ασφαλή (No Alerts)", safe_mols)
  col3.metric("Alerts / PAINS", alert_mols)
  col4.metric("Μέσος ADMET Score", f"{avg_score}%")

  st.markdown("---")

  # --- VISUAL SCORECARD (PLOTLY BAR CHART) ---
  st.subheader("📈 ADMET Score & Safety Overview")
  fig = px.bar(
      df_results,
      x="Όνομα / SMILES",
      y="ADMET Score (%)",
      color="ADMET Score (%)",
      color_continuous_scale=["red", "orange", "green"],
      range_y=[0, 100],
      title="Συνολική Βαθμολογία Ασφάλειας & Φαρμακευτικότητας Μορίων",
  )
  fig.update_layout(
      plot_bgcolor="rgba(0,0,0,0)",
      paper_bgcolor="rgba(0,0,0,0)",
      font=dict(family="sans serif", size=12),
  )
  st.plotly_chart(fig, use_container_width=True)

  st.markdown("---")
  st.subheader("📊 Αναλυτικός Πίνακας ADMET")

  # Προετοιμασία DataFrame για εμφάνιση (κρύβουμε το RAW SMILES)
  df_display = df_results.drop(columns=["SMILES"])

  # Εμφάνιση πίνακα με υποστήριξη HTML (για τις εικόνες SVG)
  st.write(
      df_display.to_html(escape=False, index=False), unsafe_allow_html=True
  )

  # --- ΛΗΨΗ ΑΡΧΕΙΟΥ ---
  st.markdown("---")
  st.markdown("### 📥 Λήψη Αποτελεσμάτων ADMET")

  # Για το export αφαιρούμε το SVG tag ώστε να είναι καθαρό το αρχείο
  df_export = df_display.drop(columns=["Δομή (2D)"])

  file_format = st.radio(
      "Επιλέξτε μορφή αρχείου για λήψη:",
      options=["CSV (.csv)", "Text (.txt)", "Excel (.xlsx)"],
      horizontal=True,
      key="download_format_radio_admet",
  )

  if file_format == "CSV (.csv)":
    file_data = df_export.to_csv(index=False).encode("utf-8")
    file_name = "admet_results.csv"
    mime_type = "text/csv"

  elif file_format == "Text (.txt)":
    file_data = df_export.to_csv(sep="\t", index=False).encode("utf-8")
    file_name = "admet_results.txt"
    mime_type = "text/plain"

  elif file_format == "Excel (.xlsx)":
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
      df_export.to_excel(writer, index=False, sheet_name="ADMET_Results")
    file_data = buffer.getvalue()
    file_name = "admet_results.xlsx"
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
