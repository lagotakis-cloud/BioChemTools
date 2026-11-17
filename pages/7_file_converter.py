import io
import pandas as pd
import streamlit as st

st.set_page_config(
    page_title="Universal Report Exporter", page_icon="📊", layout="wide"
)

from components import render_navbar

render_navbar()

st.title("📊 Universal Report & Format Exporter")
st.markdown(
    "Μετατροπή μορφής αρχείων αναφορών και αποτελεσμάτων (CSV ⇄ TXT ⇄ Excel)"
    " με δυνατότητα προεπισκόπησης και άμεσου εξαγωγής."
)
st.markdown("---")

uploaded_report = st.file_uploader(
    "Ανεβάστε το αρχείο αποτελεσμάτων σας (.csv, .txt, .xlsx):",
    type=["csv", "txt", "xlsx"],
)

if uploaded_report is not None:
  try:
    # Ανίχνευση τύπου και ανάγνωση
    file_name = uploaded_report.name
    if file_name.endswith(".csv"):
      df = pd.read_csv(uploaded_report)
    elif file_name.endswith(".txt"):
      df = pd.read_csv(uploaded_report, sep=None, engine="python")
    else:
      df = pd.read_excel(uploaded_report)

    st.success(
        f"Επιτυχής φόρτωση αρχείου `{file_name}` ({len(df)} γραμμές,"
        f" {len(df.columns)} στήλες)."
    )

    st.markdown("### 👁️ Προεπισκόπηση Δεδομένων")
    st.dataframe(df, use_container_width=True)

    st.markdown("---")
    st.subheader("📥 Επιλογή Νέας Μορφής Εξαγωγής (Convert & Export)")

    target_format = st.selectbox(
        "Επιθυμητή μορφή αρχείου:",
        ("CSV (.csv)", "Tab-Separated Text (.txt)", "Excel Workbook (.xlsx)"),
    )

    new_file_name = st.text_input(
        "Όνομα νέου αρχείου (χωρίς επέκταση):", "converted_biochem_report"
    )

    if target_format == "CSV (.csv)":
      export_data = df.to_csv(index=False).encode("utf-8")
      export_mime = "text/csv"
      final_filename = f"{new_file_name}.csv"
    elif target_format == "Tab-Separated Text (.txt)":
      export_data = df.to_csv(sep="\t", index=False).encode("utf-8")
      export_mime = "text/plain"
      final_filename = f"{new_file_name}.txt"
    else:
      buffer = io.BytesIO()
      with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="Exported_Data")
      export_data = buffer.getvalue()
      export_mime = (
          "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
      )
      final_filename = f"{new_file_name}.xlsx"

    st.download_button(
        label=f"⬇️ Λήψη νέου αρχείου ως {final_filename}",
        data=export_data,
        file_name=final_filename,
        mime=export_mime,
        type="primary",
    )

  except Exception as e:
    st.error(f"Σφάλμα κατά την επεξεργασία του αρχείου: {e}")
else:
  st.info(
      "👆 Ανεβάστε κάποιο από τα αρχεία CSV/Excel που κατεβάσατε από τα"
      " προηγούμενα εργαλεία της σουίτας για να το μετατρέψετε."
  )
