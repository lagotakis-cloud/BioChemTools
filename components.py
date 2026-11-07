import streamlit as st

def render_navbar():
    """Εμφανίζει την οριζόντια, σταθερή (fixed) μπάρα πλοήγησης στην κορυφή."""
    
    st.markdown("""
    <style>
        /* Στοχεύουμε αποκλειστικά τα πρώτα horizontal blocks της πλοήγησης 
           και τα κάνουμε fixed στην κορυφή του παραθύρου */
        .stMainBlockContainer {
            padding-top: 5rem !important; /* Προσθέτουμε κενό στην κορυφή για να μη κρύβεται το περιεχόμενο */
        }
        
        /* Κάνουμε τα δύο πρώτα blocks (τις 2 σειρές της μπάρας) σταθερά ψηλά */
        div.stHorizontalBlock:nth-of-type(1), 
        div.stHorizontalBlock:nth-of-type(2) {
            position: fixed !important;
            top: 0px !important;
            left: 0 !important;
            right: 0 !important;
            width: 100% !important;
            z-index: 999999 !important;
            background-color: var(--background-color, #090D16) !important;
            padding: 10px 20px !important;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.3);
        }

        /* Διορθώδεις αποστάσεις για να μην κολλάνε τα εικονίδια στις άκρες όταν γίνονται fixed */
        div.stHorizontalBlock:nth-of-type(1) {
            top: 0px !important;
        }
        div.stHorizontalBlock:nth-of-type(2) {
            top: 45px !important; /* Ακριβώς κάτω από την 1η σειρά */
        }
    </style>
    """, unsafe_allow_html=True)
    
    # --- 1η Σειρά Εργαλείων ---
    col1, col2, col3, col4, col5 = st.columns(5)
    with col1:
        st.page_link("app.py", label="🏠 Αρχική", icon="🏠")
    with col2:
        st.page_link("pages/1_lipinski.py", label="💊 Lipinski", icon="💊")
    with col3:
        st.page_link("pages/2_tanimoto.py", label="🔍 Tanimoto", icon="🔍")
    with col4:
        st.page_link("pages/3_admet.py", label="🧪 ADMET", icon="🧪")
    with col5:
        st.page_link("pages/4_functional_groups.py", label="⚛️ Func. Groups", icon="⚛️")

    # --- 2η Σειρά Εργαλείων ---
    col6, col7, col8, col9 = st.columns(4)
    with col6:
        st.page_link("pages/5_fragments.py", label="🔬 Fragments", icon="🔬")
    with col7:
        st.page_link("pages/6_converter.py", label="🔄 Converter", icon="🔄")
    with col8:
        st.page_link("pages/7_scaffold.py", label="🧬 Scaffold", icon="🧬")
    with col9:
        st.page_link("pages/8_similarity_matrix.py", label="📈 Similarity", icon="📈")
        
    # Κενό διάκενο από κάτω για να μη επικαλύπτονται τα επόμενα στοιχεία της σελίδας
    st.markdown("<div style='margin-bottom: 80px;'></div>", unsafe_allow_html=True)

