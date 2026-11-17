import streamlit as st

st.set_page_config(
    page_title="biochemtools.tech | Enterprise Cheminformatics",
    page_icon="🧬",
    layout="wide",
)

from components import render_navbar

render_navbar()

# --- HERO SECTION ---
st.title("🧬 biochemtools.tech")
st.markdown(
    "### Η σύγχρονη, γρήγορη και επαγγελματική σουίτα εργαλείων Cheminformatics"
    " για In Silico Screening και Ανάλυση Μορίων."
)
st.markdown("---")


# --- GRID NAVIGATION ΜΕ ΤΑ 7 ΕΡΓΑΛΕΙΑ ---
st.subheader("🛠️ Διαθέσιμα Εργαλεία Σουίτας")

c1, c2, c3 = st.columns(3)

with c1:
  st.markdown("#### 1. Lipinski's Rule of 5")
  st.markdown("Έλεγχος κανόνων φαρμακευτικότητας και φυσικοχημικών ιδιοτήτων.")
  st.page_link("pages/1_lipinski.py", label="Άνοιγμα Lipinski", icon="💊")

  st.markdown("#### 4. Ghose Filter")
  st.markdown("Αξιολόγηση βάσει MW, LogP, βαρέων ατόμων και Molar Refractivity.")
  st.page_link("pages/4_ghose.py", label="Άνοιγμα Ghose", icon="🧪")

  st.markdown("#### 7. File Converter")
  st.markdown(
      "Universal μετατροπή μορφών αρχείων αναφορών (CSV ⇄ TXT ⇄ Excel)."
  )
  st.page_link("pages/7_file_converter.py", label="Άνοιγμα Converter", icon="📁")

with c2:
  st.markdown("#### 2. Veber's Rules")
  st.markdown(
      "Ανάλυση ευελιξίας δεσμών και πολικής επιφάνειας (TPSA) για"
      " βιοδιαθεσιμότητα."
  )
  st.page_link("pages/2_veber.py", label="Άνοιγμα Veber", icon="📐")

  st.markdown("#### 5. QED Evaluator")
  st.markdown(
      "Υπολογισμός ποσοτικού δείκτη ελκυστικότητας (0 έως 1) για φάρμακα."
  )
  st.page_link("pages/5_qed.py", label="Άνοιγμα QED", icon="📈")

with c3:
  st.markdown("#### 3. Enterprise ADMET")
  st.markdown(
      "Προηγμένος έλεγχος ADMET, ESOL, PAINS, CYP450 και τοξικότητας."
  )
  st.page_link("pages/3_admet.py", label="Άνοιγμα ADMET", icon="🛡️")

  st.markdown("#### 6. Name ⇄ SMILES")
  st.markdown(
      "Αμφίδρομη μαζική μετατροπή ονομάτων (PubChem) σε SMILES και αντίστροφα."
  )
  st.page_link(
      "pages/6_names_to_smiles.py", label="Άνοιγμα Name ⇄ SMILES", icon="🔤"
  )

st.markdown("---")
st.markdown(
    "💡 *Tip: Μπορείτε ανά πάσα στιγμή να πλοηγηθείτε σε οποιοδήποτε εργαλείο"
    " χρησιμοποιώντας το μενού (sidebar) στα αριστερά.*"
)
