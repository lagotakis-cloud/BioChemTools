import streamlit as st

def render_navbar():
  col_n1, col_n2, col_n3, col_n4, col_n5, col_n6, col_n7, col_n8, col_n9 = st.columns(9)
with col_n1:
    st.page_link("app.py", label="Αρχική")
with col_n2:
    st.page_link("pages/1_lipinski.py", label="Lipinski Rule of 5")
with col_n3:
    st.page_link("pages/2_tanimoto.py", label="Tanimoto Similarity")
with col_n4:
    st.page_link("pages/3_admet.py", label="ADMET Profiler")
with col_n5:
    st.page_link("pages/4_functional_groups.py", label="Functional Groups")
with col_n6:
    st.page_link("pages/5_fragments.py", label="Molecular Fragments")
with col_n7:
    st.page_link("pages/6_converter.py", label="Format Converter")
with col_n8:
    st.page_link("pages/7_scaffold.py", label="Scaffold Splitter")
with col_n9:
    st.page_link("pages/8_similarity_matrix.py", label="Similarity Matrix")

st.markdown("---")

