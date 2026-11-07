import streamlit as st

# Ρύθμιση σελίδας με κλειστή πλαϊνή μπάρα για καθαρό dashboard look
st.set_page_config(
    page_title="biochemtools.tech | Cheminformatics Suite",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Custom CSS για μοντέρνα εμφάνιση και banners
st.markdown("""
<style>
    .main-title {
        font-size: 2.8rem;
        font-weight: 800;
        color: #1E3A8A;
        margin-bottom: 0px;
    }
    .subtitle {
        font-size: 1.2rem;
        color: #4B5563;
        margin-bottom: 30px;
    }
    .nav-bar {
        background-color: #F3F4F6;
        padding: 10px 20px;
        border-radius: 8px;
        margin-bottom: 30px;
        display: flex;
        justify-content: space-around;
        font-weight: 600;
    }
</style>
""", unsafe_allow_html=True)

# ----------------- HEADER & HORIZONTAL NAVIGATION -----------------
st.markdown('<p class="main-title">🧬 biochemtools.tech</p>', unsafe_allow_html=True)
st.markdown('<p class="subtitle">Ανοιχτή, ταχεία σουίτα υπολογιστικής χημείας και cheminformatics για την ανακάλυψη φαρμάκων.</p>', unsafe_allow_html=True)

# Οριζόντια μπάρα γρήγορων συνδέσμων
col_n1, col_n2, col_n3, col_n4 = st.columns(4)
with col_n1:
    st.page_link("Home.py", label="🏠 Αρχική", icon="🏠")
with col_n2:
    st.page_link("pages/1_lipinski.py", label="💊 Lipinski Rule of 5", icon="💊")
with col_n3:
    st.page_link("pages/3_admet.py", label="🧪 ADMET Profiler", icon="🧪")
with col_n4:
    st.page_link("pages/8_similarity_matrix.py", label="📈 Similarity Matrix", icon="📈")

st.markdown("---")

# ----------------- BANNERS / CARDS ΓΙΑ ΤΑ 8 ΕΡΓΑΛΕΙΑ -----------------
st.subheader("🚀 Εξερευνήστε τα Εργαλεία της Πλατφόρμας")

# Σειρά 1 (Εργαλεία 1-4)
col1, col2, col3, col4 = st.columns(4)

with col1:
    with st.container(border=True):
        st.markdown("### 💊 Lipinski Rule of 5")
        st.caption("Έλεγχος φαρμακευτικότητας (MW, LogP, H-Donors/Acceptors).")
        st.page_link("pages/1_lipinski.py", label="Άνοιγμα Εργαλείου", icon="👉")

with col2:
    with st.container(border=True):
        st.markdown("### 🔍 Tanimoto Search")
        st.caption("Αναζήτηση μοριακής ομοιότητας βάσει Morgan Fingerprints.")
        st.page_link("pages/2_tanimoto.py", label="Άνοιγμα Εργαλείου", icon="👉")

with col3:
    with st.container(border=True):
        st.markdown("### 💊 ADMET Profiler")
        st.caption("Υπολογισμός φυσικοχημικών ιδιοτήτων και προφίλ ADMET.")
        st.page_link("pages/3_admet.py", label="Άνοιγμα Εργαλείου", icon="👉")

with col4:
    with st.container(border=True):
        st.markdown("### 🧪 Functional Groups")
        st.caption("Αυτόματη ανίχνευση και καταμέτρηση βασικών ομάδων.")
        st.page_link("pages/4_functional_groups.py", label="Άνοιγμα Εργαλείου", icon="👉")

# Σειρά 2 (Εργαλεία 5-8)
col5, col6, col7, col8 = st.columns(4)

with col5:
    with st.container(border=True):
        st.markdown("### 🔬 Fragments (Ro3)")
        st.caption("Rule of 3 checker για Fragment-Based Drug Design.")
        st.page_link("pages/5_fragments.py", label="Άνοιγμα Εργαλείου", icon="👉")

with col6:
    with st.container(border=True):
        st.markdown("### 🔄 Format Converter")
        st.caption("Μετατροπή SMILES ⇄ InChI ⇄ SDF και λήψη αρχείων.")
        st.page_link("pages/6_converter.py", label="Άνοιγμα Εργαλείου", icon="👉")

with col7:
    with st.container(border=True):
        st.markdown("### 🧬 Scaffold Splitter")
        st.caption("Απομόνωση κεντρικού πυρήνα Bemis-Murcko.")
        st.page_link("pages/7_scaffold.py", label="Άνοιγμα Εργαλείου", icon="👉")

with col8:
    with st.container(border=True):
        st.markdown("### 📈 Similarity Matrix")
        st.caption("Υπολογισμός μήτρας ομοιότητας και Heatmap.")
        st.page_link("pages/8_similarity_matrix.py", label="Άνοιγμα Εργαλείου", icon="👉")

st.markdown("---")
st.markdown("<p style='text-align: center; color: #6B7280;'>biochemtools.tech © 2026 — Built with Streamlit & RDKit</p>", unsafe_allow_html=True)
