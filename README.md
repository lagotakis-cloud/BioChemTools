# 🧬 BiochemTools (biochemtools.tech)

A modern, fast, and user-friendly Cheminformatics and Drug-Likeness Screening web suite, built in Python using Streamlit and RDKit.

---

## 🛠️ Included Tools

The platform is divided into independent pages, offering comprehensive molecular analysis:

1. **Lipinski Rule of 5 Checker (`1_lipinski.py`)**
   * Calculation of Molecular Weight (MW), LogP, Hydrogen Bond Donors (HBD), and Acceptors (HBA).
   * Rule violation check and automated status classification (Passed / Failed).
   
2. **Veber's Rules Checker (`2_veber.py`)**
   * Evaluation of Topological Polar Surface Area (TPSA) with a selectable threshold (≤ 140 Å² for cell membranes or ≤ 90 Å² for BBB).
   * Molecular flexibility assessment via Rotatable Bonds (≤ 10).

3. **ADMET & Toxicity Screener (`3_admet.py`)**
   * Estimation of aqueous solubility (ESOL / LogS).
   * PAINS filters check (Pan Assay Interference Compounds).
   * Detection of reactive structural alerts (e.g., epoxides, Michael acceptors, etc.).

---

## 🚀 Input & Export Capabilities
* **Multiple Input Methods:** Support for built-in sample molecules, bulk text input of Names or SMILES strings, and file uploads (`.csv` / `.txt`).
* **Batch Processing:** Rapid analysis of multiple compounds simultaneously.
* **Flexible Downloads:** Export clean results in **CSV**, **TXT**, or **Excel (.xlsx)** formats with built-in styling.

---

## 📂 Project Structure

```text
biochemtools/
│
├── pages/
│   ├── 1_lipinski.py       # Lipinski Rule of 5 tool
│   ├── 2_veber.py          # Veber's Rules tool
│   └── 3_admet.py          # ADMET & Toxicity Screener tool
│
├── assets/
│   └── logo.png            # Icons and graphic assets
│
├── requirements.txt        # Required Python libraries
└── README.md               # Project documentation
