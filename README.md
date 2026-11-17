# 🧬 biochemtools.tech | A Free Tool Suite for Drug Discovery

**biochemtools.tech** is a modern, fast, and free tool suite designed for researchers, medicinal chemists, and students. Built with **Streamlit**, **RDKit**, **PubChemPy**, and **Pandas**, it provides robust local tools for high-throughput in silico screening, molecular property evaluation, and chemical translation.

Designed with a **"First-Line Triage"** philosophy, the platform enables secure, local, and lightning-fast preliminary screening of chemical libraries without the security risks or latency of external third-party servers.

---

## The 6 Core Suite Tools

The platform consists of six specialized, standalone tools integrated via a seamless sticky navigation header:

1. **`1_lipinski.py` (Lipinski's Rule of 5 Checker)**
   * Evaluates core physicochemical properties (Molecular Weight, LogP, HBD, HBA) to assess oral bioavailability and drug-likeness.
2. **`2_veber.py` (Veber's Rules & TPSA Screener)**
   * Analyzes molecular flexibility (Rotatable Bonds) and polar surface area (TPSA) alongside Blood-Brain Barrier (BBB) potential.
3. **`3_admet.py` (Simple ADMET & Toxicity Screener)**
   * Computes aqueous solubility (ESOL/LogS), screens for Pan Assay Interference Compounds (PAINS), detects toxic structural alerts via SMARTS, and profiles hepatic CYP450 enzyme inhibition.
4. **`4_ghose.py` (Ghose Filter Checker)**
   * Validates molecules against comprehensive statistical ranges for MW, LogP, Molar Refractivity, and total atom count derived from clinical drug databases.
5. **`5_qed.py` (QED Drug-Likeness Evaluator)**
   * Calculates the continuous Quantitative Estimate of Drug-Likeness (0 to 1 score) based on desirability functions of eight fundamental properties.
6. **`6_names_to_smiles.py` (Name ⇄ SMILES Batch Converter)**
   * Provides high-throughput bidirectional batch translation between common/IUPAC chemical names and canonical SMILES strings via PubChem integration.

---

## Key Features

- **Local & Secure Processing:** Keeps your proprietary chemical structures and compound libraries strictly private on your local machine or secure deployment environment.
- **Flexible Data Input:** Supports single-compound name search (via PubChem), direct manual SMILES entry, and automated batch file uploads (`.csv` and `.txt`).
- **Interactive Visual Feedback:** Dynamic KPI metrics, automated error handling, and color-coded results tables for instant molecular prioritization.
- **Universal Export:** Instant export options (`.csv`, `.txt`, `.xlsx` via OpenPyXL) directly from generated results.

---

## Tech Stack

- **Frontend / UI:** [Streamlit](https://streamlit.io/) (with custom CSS components & sticky dual-row navigation)
- **Cheminformatics Engine:** [RDKit](https://www.rdkit.org/)
- **Database & Retrieval:** [PubChemPy](https://pypi.org/project/PubChemPy/)
- **Data Manipulation & Export:** Pandas, OpenPyXL

---

## Installation & Local Setup

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/your-username/biochemtools.git](https://github.com/your-username/biochemtools.git)
   cd biochemtools
