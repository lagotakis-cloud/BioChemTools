import streamlit as st

def render_navbar():
    """Εμφανίζει μια καθαρή, κομψή οριζόντια μπάρα πλοήγησης για όλα τα 8 εργαλεία."""
    
    st.markdown("""
    <style>
        /* Μικρό, καθαρό στυλ για τα links της μπάρας ώστε να δείχνουν ομοιόμορφα */
        .stPageLink {
            font-size: 0.9rem;
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
        
    st.markdown("---")
