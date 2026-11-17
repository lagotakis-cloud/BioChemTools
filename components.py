import streamlit as st

def render_navbar():
    """Εμφανίζει την οριζόντια μπάρα πλοήγησης, κολλημένη σταθερά ψηλά κατά το scroll."""
    
    st.markdown("""
    <style>
        /* Κάνουμε sticky ΜΟΝΟ τις δύο πρώτες σειρές blocks (δηλαδή το navbar) 
           στην κορυφή, χωρίς να επηρεάζουμε τα banners παρακάτω */
        div.stMainBlockContainer > div.stVerticalBlock > div.stHorizontalBlock:nth-of-type(1),
        div.stMainBlockContainer > div.stVerticalBlock > div.stHorizontalBlock:nth-of-type(2) {
            position: sticky !important;
            top: 0px !important;
            z-index: 99999 !important;
            background-color: var(--background-color, #090D16) !important;
            padding-top: 8px !important;
            padding-bottom: 8px !important;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.2);
        }
        
        /* Η 2η σειρά καθήμενή ακριβώς κάτω από την 1η */
        div.stMainBlockContainer > div.stVerticalBlock > div.stHorizontalBlock:nth-of-type(2) {
            top: 48px !important; 
        }
    </style>
    """, unsafe_allow_html=True)
    
    # --- 1η Σειρά Εργαλείων ---
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.page_link("app.py", label="Home Page")
    with col2:
        st.page_link("pages/1_lipinski.py", label="Lipinski Rules")
    with col3:
        st.page_link("pages/2_veber.py", label="Veber Rules")
    with col4:
        st.page_link("pages/3_admet.py", label="ADMET Profile")

    # --- 2η Σειρά Εργαλείων ---
    col5, col6, col7, col8 = st.columns(4)
    with col5:
        st.page_link("pages/4_ghose.py", label="Ghose Filter")
    with col6:
        st.page_link("pages/5_qed.py", label="QED Score")
    with col7:
        st.page_link("pages/6_names_to_smiles.py", label="Names to SMILES")
    with col8:
        st.page_link("pages/7_file_converter.py", label="File Converter")
        
    # Μικρό διακριτικό κενό μετά το navbar για να ξεκινήσει το περιεχόμενο της σελίδας
    st.markdown("<div style='margin-bottom: 25px;'></div>", unsafe_allow_html=True)
