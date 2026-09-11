```python
# ============================================================
# PCOS TRANSCRIPTOMICS STREAMLIT APP
# Gene ID -> Gene Symbol -> Gene Description
# ============================================================

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from pathlib import Path

# ============================================================
# 1. APP CONFIGURATION
# ============================================================

st.set_page_config(
    page_title="PCOS Transcriptomics Portal",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("🧬 PCOS Transcriptomics Analysis Portal")
st.markdown(
    """
    **Integrated PCOS transcriptomic analysis**

    Differential Expression → GSEA → TF Regulation → Candidate Genes
    """
)

# ============================================================
# 2. FILE LOCATIONS
# ============================================================

BASE_DIR = Path(__file__).resolve().parent

DE_FILE = BASE_DIR / "pcos_deseq2_results.csv"
GSEA_FILE = BASE_DIR / "gsea_results_precalculated.csv"
TF_FILE = BASE_DIR / "tf_screening_results.csv"
ANNOTATION_FILE = BASE_DIR / "Human.GRCh38.p13.annot.tsv"


# ============================================================
# 3. LOAD FILES
# ============================================================

@st.cache_data
def load_csv(file_path):
    """Load CSV safely."""
    return pd.read_csv(file_path)


@st.cache_data
def load_annotation(file_path):
    """
    Load NCBI GRCh38 annotation file.

    Required columns:
    GeneID
    Symbol
    Description
    """
    annotation = pd.read_csv(
        file_path,
        sep="\t",
        dtype=str,
        low_memory=False
    )

    required_columns = [
        "GeneID",
        "Symbol",
        "Description"
    ]

    missing = [
        col for col in required_columns
        if col not in annotation.columns
    ]

    if missing:
        raise ValueError(
            f"Annotation file is missing columns: {missing}"
        )

    annotation = annotation[
        ["GeneID", "Symbol", "Description"]
    ].copy()

    annotation["GeneID"] = annotation["GeneID"].astype(str).str.strip()
    annotation["Symbol"] = annotation["Symbol"].fillna("")
    annotation["Description"] = annotation["Description"].fillna("")

    # Remove duplicate GeneIDs
    annotation = annotation.drop_duplicates(
        subset=["GeneID"],
        keep="first"
    )

    return annotation


# ============================================================
# 4. CHECK REQUIRED FILES
# ============================================================

missing_files = []

for file_path in [
    DE_FILE,
    GSEA_FILE,
    TF_FILE,
    ANNOTATION_FILE
]:
    if not file_path.exists():
        missing_files.append(file_path.name)

if missing_files:

    st.error("❌ The following required files are missing:")

    for file_name in missing_files:
        st.write(f"- `{file_name}`")

    st.info(
        """
        Put these files in the **same folder as appy.py**.
        """
    )

    st.stop()


# ============================================================
# 5. LOAD DATA
# ============================================================

try:

    de = load_csv(DE_FILE)
    gsea = load_csv(GSEA_FILE)
    tf = load_csv(TF_FILE)
    annotation = load_annotation(ANNOTATION_FILE)

except Exception as e:

    st.error("Error loading files.")
    st.exception(e)
    st.stop()


# ============================================================
# 6. STANDARDIZE GENE IDs
# ============================================================

if "GeneID" not in de.columns:

    st.error(
        "The DESeq2 file does not contain a `GeneID` column."
    )

    st.stop()


de["GeneID"] = (
    de["GeneID"]
    .astype(str)
    .str.replace(r"\.0$", "", regex=True)
    .str.strip()
)


# ============================================================
# 7. MERGE GENE ANNOTATION
# ============================================================

de = de.merge(
    annotation,
    on="GeneID",
    how="left"
)


# ============================================================
# 8. CLEAN ANNOTATION
# ============================================================

de["Symbol"] = de["Symbol"].fillna("")
de["Description"] = de["Description"].fillna("")


# If no Symbol is found, use GeneID
de["Display_Gene"] = np.where(
    de["Symbol"].str.strip() != "",
    de["Symbol"],
    "GeneID:" + de["GeneID"]
)


# ============================================================
# 9. CREATE USER-FRIENDLY GENE LABEL
# ============================================================

de["Gene_Label"] = np.where(
    de["Description"].str.strip() != "",
    de["Display_Gene"] + " — " + de["Description"],
    de["Display_Gene"]
)


# ============================================================
# 10. SIDEBAR
# ============================================================

st.sidebar.header("⚙️ Analysis Controls")

padj_cutoff = st.sidebar.number_input(
    "Adjusted P-value cutoff",
    min_value=0.001,
    max_value=0.20,
    value=0.05,
    step=0.001
)

log2fc_cutoff = st.sidebar.number_input(
    "Absolute log2 Fold Change cutoff",
    min_value=0.0,
    max_value=10.0,
    value=1.0,
    step=0.1
)


# ============================================================
# 11. PROCESS DESEQ2 RESULTS
# ============================================================

required_de_columns = [
    "padj",
    "log2FoldChange"
]

missing_de_columns = [
    col for col in required_de_columns
    if col not in de.columns
]

if missing_de_columns:

    st.error(
        f"DESeq2 file is missing: {missing_de_columns}"
    )

    st.stop()


de["padj"] = pd.to_numeric(
    de["padj"],
    errors="coerce"
)

de["log2FoldChange"] = pd.to_numeric(
    de["log2FoldChange"],
    errors="coerce"
)

de = de.replace(
    [np.inf, -np.inf],
    np.nan
)

de = de.dropna(
    subset=["padj", "log2FoldChange"]
)


# ============================================================
# 12. VOLCANO PLOT VARIABLES
# ============================================================

de["-log10(padj)"] = -np.log10(
    de["padj"].clip(lower=1e-300)
)

de["Significance"] = "Not significant"

de.loc[
    (
        (de["padj"] < padj_cutoff)
        &
        (de["log2FoldChange"] >= log2fc_cutoff)
    ),
    "Significance"
] = "Upregulated"

de.loc[
    (
        (de["padj"] < padj_cutoff)
        &
        (de["log2FoldChange"] <= -log2fc_cutoff)
    ),
    "Significance"
] = "Downregulated"


# ============================================================
# 13. SUMMARY COUNTS
# ============================================================

upregulated = int(
    (
        de["Significance"] == "Upregulated"
    ).sum()
)

downregulated = int(
    (
        de["Significance"] == "Downregulated"
    ).sum()
)

significant_genes = upregulated + downregulated


if "FDR q-val" in gsea.columns:

    gsea["FDR q-val"] = pd.to_numeric(
        gsea["FDR q-val"],
        errors="coerce"
    )

    significant_pathways = int(
        (
            gsea["FDR q-val"] < padj_cutoff
        ).sum()
    )

else:

    significant_pathways = 0


# ============================================================
# 14. MAIN TABS
# ============================================================

tab_overview, tab_deg, tab_gsea, tab_tf, tab_gene = st.tabs(
    [
        "📊 Overview",
        "🌋 Differential Expression",
        "🧬 GSEA",
        "🎯 TF Regulators",
        "🔎 Gene Explorer"
    ]
)


# ============================================================
# TAB 1 — OVERVIEW
# ============================================================

with tab_overview:

    st.header("Study Overview")

    col1, col2, col3, col4 = st.columns(4)

    col1.metric(
        "Genes analysed",
        f"{len(de):,}"
    )

    col2.metric(
        "Upregulated",
        f"{upregulated:,}"
    )

    col3.metric(
        "Downregulated",
        f"{downregulated:,}"
    )

    col4.metric(
        "Significant pathways",
        f"{significant_pathways:,}"
    )

    st.divider()

    st.subheader("🧬 Gene annotation")

    annotated_count = int(
        (
            de["Symbol"].str.strip() != ""
        ).sum()
    )

    annotation_percentage = (
        annotated_count / len(de) * 100
        if len(de) > 0
        else 0
    )

    st.write(
        f"""
        **Gene annotation coverage:**
        {annotated_count:,} / {len(de):,}
        genes ({annotation_percentage:.2f}%)
        """
    )

    st.info(
        """
        Gene identifiers are mapped using the supplied
        `Human.GRCh38.p13.annot.tsv` annotation file.
        """
    )

    st.subheader("Top differentially expressed genes")

    display_columns = [
        "GeneID",
        "Symbol",
        "Description",
        "log2FoldChange",
        "padj"
    ]

    display_columns = [
        col
        for col in display_columns
        if col in de.columns
    ]

    top_genes = (
        de
        .sort_values("padj")
        .head(25)
    )

    st.dataframe(
        top_genes[display_columns],
        use_container_width=True,
        hide_index=True
    )


# ============================================================
# TAB 2 — DIFFERENTIAL EXPRESSION
# ============================================================

with tab_deg:

    st.header("🌋 Differentially Expressed Genes")

    st.caption(
        f"""
        Thresholds:
        adjusted P-value < {padj_cutoff}
        and |log2FC| ≥ {log2fc_cutoff}
       
```
