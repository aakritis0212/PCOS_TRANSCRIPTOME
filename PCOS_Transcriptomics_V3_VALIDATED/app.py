import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA, RESULTS = ROOT/"data", ROOT/"results"
st.set_page_config(page_title="PCOS Transcriptomics", page_icon="🧬", layout="wide")

@st.cache_data
def read(path):
    return pd.read_csv(path) if path.exists() else pd.DataFrame()

de = read(DATA/"pcos_deseq2_results.csv")
gsea = read(DATA/"gsea_results_precalculated.csv")
tf = read(DATA/"tf_screening_results.csv")
cand = read(RESULTS/"current_leading_edge_candidates.csv")

st.title("🧬 PCOS Transcriptomics Research Portal")
st.caption("GSE168404 analysis with reproducible pathway/TF integration and independent-validation framework.")

if de.empty:
    st.error("DESeq2 result file is missing.")
    st.stop()

cut = st.sidebar.number_input("Adjusted P-value cutoff", .001, .20, .05, .001)
lfc = st.sidebar.number_input("|log2FC| cutoff", 0.0, 5.0, 1.0, .1)

de = de.dropna(subset=["padj","log2FoldChange"]).copy()
de["-log10(padj)"] = -np.log10(de.padj.clip(lower=1e-300))
de["Class"] = "Not significant"
de.loc[(de.padj < cut)&(de.log2FoldChange >= lfc),"Class"] = "Upregulated"
de.loc[(de.padj < cut)&(de.log2FoldChange <= -lfc),"Class"] = "Downregulated"

tabs = st.tabs(["Overview","DEGs","GSEA","TFs","Candidate module","Validation"])

with tabs[0]:
    a,b,c,d = st.columns(4)
    a.metric("Genes tested", f"{len(de):,}")
    b.metric("Upregulated", int((de.Class=="Upregulated").sum()))
    c.metric("Downregulated", int((de.Class=="Downregulated").sum()))
    d.metric("GSEA FDR<cut", int((gsea["FDR q-val"]<cut).sum()) if not gsea.empty else 0)
    st.info("This is a hypothesis-generation and validation framework; it does not establish causality.")

with tabs[1]:
    fig = px.scatter(de, x="log2FoldChange", y="-log10(padj)", color="Class",
                     hover_data=["GeneID","padj"], title="PCOS vs Control")
    fig.add_vline(x=lfc, line_dash="dash"); fig.add_vline(x=-lfc, line_dash="dash")
    fig.add_hline(y=-np.log10(cut), line_dash="dash")
    st.plotly_chart(fig, use_container_width=True)
    st.dataframe(de.sort_values("padj").head(100), use_container_width=True)

with tabs[2]:
    if gsea.empty:
        st.warning("No GSEA table found.")
    else:
        x = gsea.copy()
        x["FDR q-val"] = pd.to_numeric(x["FDR q-val"], errors="coerce")
        x["NES"] = pd.to_numeric(x["NES"], errors="coerce")
        x = x[x["FDR q-val"]<cut].copy()
        x["absNES"] = x.NES.abs()
        x = x.sort_values("absNES", ascending=False).head(25)
        st.plotly_chart(px.bar(x.sort_values("NES"), x="NES", y="Term", orientation="h",
                               color="FDR q-val"), use_container_width=True)
        st.dataframe(x[["Term","NES","FDR q-val","Lead_genes"]], use_container_width=True)

with tabs[3]:
    if tf.empty:
        st.warning("No TF results found.")
    else:
        x = tf[tf["FDR q-val"]<cut].sort_values("NES", key=np.abs, ascending=False)
        st.dataframe(x[["Term","NES","FDR q-val","Lead_genes"]], use_container_width=True)

with tabs[4]:
    if cand.empty:
        st.warning("No candidate module table.")
    else:
        st.dataframe(cand.head(100), use_container_width=True)
        st.download_button("Download candidate module",
                           cand.to_csv(index=False).encode(),
                           "current_leading_edge_candidates.csv", "text/csv")

with tabs[5]:
    st.subheader("Independent validation cohorts")
    st.markdown("""
    **GSE293353:** 9 PCOS + 9 controls, bulk granulosa-cell RNA-seq  
    **GSE193123:** 3 PCOS + 3 controls, bulk granulosa-cell RNA-seq  
    **GSE155489:** 4 PCOS + 4 controls, cumulus granulosa-cell RNA-seq  
    **GSE240688:** 3 PCOS + 3 controls, single-cell granulosa-cell RNA-seq
    """)
    st.warning("External count matrices must be downloaded and analyzed within each study before replication statistics are claimed.")
    st.code("python scripts/download_external_validation.py\n# then perform study-specific DE\n# then python scripts/validate_external.py")
