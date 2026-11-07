import streamlit as st

def render_navbar():
    """Εμφανίζει την οριζόντια μπάρα πλοήγησης και τα γρήγορα links."""
    
    # Εδώ δηλώνονται κανονικά οι στήλες εντός της συνάρτησης
    col_n1, col_n2, col_n3, col_n4 = st.columns(4)
    
    with col_n1:
        st.page_link("app.py", label="🏠 Αρχική", icon="🏠")
    with col_n2:
        st.page_link("pages/1_lipinski.py", label="💊 Lipinski", icon="💊")
    with col_n3:
        st.page_link("pages/3_admet.py", label="🧪 ADMET", icon="🧪")
    with col_n4:
        st.page_link("pages/8_similarity_matrix.py", label="📈 Similarity", icon="📈")
        
    st.markdown("---")
