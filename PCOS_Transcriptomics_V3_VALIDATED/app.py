import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from pathlib import Path

st.set_page_config(page_title='PCOS Transcriptomics Portal', page_icon='🧬', layout='wide')

BASE = Path(__file__).resolve().parent
DE_FILE = BASE / 'pcos_deseq2_results.csv'
GSEA_FILE = BASE / 'gsea_results_precalculated.csv'
TF_FILE = BASE / 'tf_screening_results.csv'
ANNOT_FILE = BASE / 'Human.GRCh38.p13.annot.tsv'

@st.cache_data

def load_csv(path, sep=','):
    return pd.read_csv(path, sep=sep, dtype=str, low_memory=False)

for p in [DE_FILE, GSEA_FILE, TF_FILE, ANNOT_FILE]:
    if not p.exists():
        st.error(f'Missing required file: {p.name}')
        st.info('Put app.py and all four data files in the same Streamlit folder.')
        st.stop()

try:
    de = load_csv(DE_FILE)
    gsea = load_csv(GSEA_FILE)
    tf = load_csv(TF_FILE)
    annot = load_csv(ANNOT_FILE, sep='\t')
except Exception as exc:
    st.error('Could not read one or more input files.')
    st.exception(exc)
    st.stop()

# ---------- Annotation ----------
if 'GeneID' not in de.columns:
    st.error('pcos_deseq2_results.csv must contain a GeneID column.')
    st.stop()

if not {'GeneID', 'Symbol', 'Description'}.issubset(annot.columns):
    st.error('Human.GRCh38.p13.annot.tsv must contain GeneID, Symbol and Description columns.')
    st.stop()

def clean_id(series):
    return (series.astype(str).str.strip().str.replace(r'\.0$', '', regex=True))

de['GeneID'] = clean_id(de['GeneID'])
annot['GeneID'] = clean_id(annot['GeneID'])
annot = annot[['GeneID', 'Symbol', 'Description']].drop_duplicates('GeneID')

de = de.merge(annot, on='GeneID', how='left', suffixes=('', '_annotation'))
de['Symbol'] = de['Symbol'].fillna('').astype(str).str.strip()
de['Description'] = de['Description'].fillna('').astype(str).str.strip()
de['Display_Gene'] = np.where(de['Symbol'].ne(''), de['Symbol'], 'GeneID:' + de['GeneID'])

# ---------- Numeric fields ----------
for col in ['baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']:
    if col in de.columns:
        de[col] = pd.to_numeric(de[col], errors='coerce')

if 'padj' not in de.columns or 'log2FoldChange' not in de.columns:
    st.error('DESeq2 results must contain padj and log2FoldChange columns.')
    st.stop()

de = de.replace([np.inf, -np.inf], np.nan)
de = de.dropna(subset=['padj', 'log2FoldChange']).copy()

st.title('🧬 PCOS Transcriptomics Analysis Portal')
st.caption('Differential Expression → GSEA → TF Regulation → Gene Annotation → Validation')

# ---------- Sidebar ----------
st.sidebar.header('⚙️ Analysis Controls')
padj_cutoff = st.sidebar.number_input('Adjusted P-value cutoff', min_value=0.001, max_value=0.20, value=0.05, step=0.001)
lfc_cutoff = st.sidebar.number_input('|log2FC| cutoff', min_value=0.0, max_value=10.0, value=1.0, step=0.1)

# ---------- DEG classification ----------
de['-log10(padj)'] = -np.log10(de['padj'].clip(lower=1e-300))
de['Significance'] = 'Not significant'
de.loc[(de['padj'] < padj_cutoff) & (de['log2FoldChange'] >= lfc_cutoff), 'Significance'] = 'Upregulated'
de.loc[(de['padj'] < padj_cutoff) & (de['log2FoldChange'] <= -lfc_cutoff), 'Significance'] = 'Downregulated'

up = int((de['Significance'] == 'Upregulated').sum())
down = int((de['Significance'] == 'Downregulated').sum())

# ---------- Tabs ----------
tabs = st.tabs(['📊 Overview', '🌋 Differential Expression', '🧬 GSEA', '🎯 TF Regulators', '🔎 Gene Explorer', '🧪 Validation'])

# ---------- Overview ----------
with tabs[0]:
    c1, c2, c3, c4 = st.columns(4)
    c1.metric('Genes analysed', f'{len(de):,}')
    c2.metric('Upregulated', f'{up:,}')
    c3.metric('Downregulated', f'{down:,}')

    if 'FDR q-val' in gsea.columns:
        gsea_fdr = pd.to_numeric(gsea['FDR q-val'], errors='coerce')
        c4.metric('GSEA FDR < cutoff', f'{int((gsea_fdr < padj_cutoff).sum()):,}')
    else:
        c4.metric('GSEA pathways', 'N/A')

    st.divider()
    annotated = int(de['Symbol'].ne('').sum())
    coverage = annotated / len(de) * 100 if len(de) else 0
    st.subheader('Gene annotation')
    st.write(f'**Annotation coverage:** {annotated:,} / {len(de):,} genes ({coverage:.2f}%)')
    st.info('Gene IDs are mapped to Gene Symbol and Gene Description using Human.GRCh38.p13.annot.tsv.')

    cols = [c for c in ['GeneID', 'Symbol', 'Description', 'log2FoldChange', 'pvalue', 'padj'] if c in de.columns]
    st.subheader('Top genes by adjusted P-value')
    st.dataframe(de.sort_values('padj').head(25)[cols], use_container_width=True, hide_index=True)

# ---------- DEGs ----------
with tabs[1]:
    st.header('🌋 Differentially Expressed Genes')
    fig = px.scatter(
        de, x='log2FoldChange', y='-log10(padj)', color='Significance',
        hover_name='Display_Gene',
        hover_data={'GeneID': True, 'Symbol': True, 'Description': True, 'log2FoldChange': ':.3f', 'padj': ':.3e'},
        title='PCOS vs Control — Volcano Plot'
    )
    fig.add_vline(x=lfc_cutoff, line_dash='dash')
    fig.add_vline(x=-lfc_cutoff, line_dash='dash')
    fig.add_hline(y=-np.log10(padj_cutoff), line_dash='dash')
    fig.update_layout(height=650)
    st.plotly_chart(fig, use_container_width=True)

    sig = de[de['Significance'] != 'Not significant'].sort_values('padj')
    cols = [c for c in ['GeneID', 'Symbol', 'Description', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'] if c in sig.columns]
    st.subheader(f'Significant genes ({len(sig)})')
    st.dataframe(sig[cols], use_container_width=True, hide_index=True)
    st.download_button('⬇️ Download annotated DEG results', sig.to_csv(index=False).encode('utf-8'), 'PCOS_annotated_DEGs.csv', 'text/csv')

# ---------- GSEA ----------
with tabs[2]:
    st.header('🧬 Gene Set Enrichment Analysis')
    if gsea.empty or 'FDR q-val' not in gsea.columns:
        st.warning('GSEA results are not available.')
    else:
        x = gsea.copy()
        x['FDR q-val'] = pd.to_numeric(x['FDR q-val'], errors='coerce')
        if 'NES' in x.columns:
            x['NES'] = pd.to_numeric(x['NES'], errors='coerce')
        x = x[x['FDR q-val'] < padj_cutoff].copy()
        if 'NES' in x.columns and not x.empty:
            x['absNES'] = x['NES'].abs()
            top = x.sort_values('absNES', ascending=False).head(25)
            fig = px.bar(top.sort_values('NES'), x='NES', y='Term', orientation='h', color='FDR q-val', title='Top Significant Pathways')
            fig.update_layout(height=750)
            st.plotly_chart(fig, use_container_width=True)
        else:
            top = x.head(25)
        cols = [c for c in ['Term', 'ES', 'NES', 'FDR q-val', 'Lead_genes'] if c in top.columns]
        st.dataframe(top[cols], use_container_width=True, hide_index=True)

# ---------- TF ----------
with tabs[3]:
    st.header('🎯 Candidate Transcriptional Regulators')
    if tf.empty or 'FDR q-val' not in tf.columns:
        st.warning('TF screening results are not available.')
    else:
        x = tf.copy()
        x['FDR q-val'] = pd.to_numeric(x['FDR q-val'], errors='coerce')
        if 'NES' in x.columns:
            x['NES'] = pd.to_numeric(x['NES'], errors='coerce')
        x = x[x['FDR q-val'] < padj_cutoff].copy()
        if 'NES' in x.columns and not x.empty:
            x = x.sort_values('NES', key=lambda s: s.abs(), ascending=False)
            fig = px.bar(x, x='NES', y='Term', orientation='h', color='FDR q-val', title='Significant TF Signatures')
            fig.update_layout(height=550)
            st.plotly_chart(fig, use_container_width=True)
        cols = [c for c in ['Term', 'NES', 'FDR q-val', 'Lead_genes'] if c in x.columns]
        st.dataframe(x[cols], use_container_width=True, hide_index=True)

# ---------- Gene Explorer ----------
with tabs[4]:
    st.header('🔎 Gene Explorer')
    st.write('Search by **Gene ID**, **Gene Symbol**, or **Gene Description**.')
    query = st.text_input('Search gene', placeholder='Example: 205, FOXM1, kinase...')
    if query.strip():
        q = query.strip().lower()
        mask = (
            de['GeneID'].str.lower().str.contains(q, na=False) |
            de['Symbol'].str.lower().str.contains(q, na=False) |
            de['Description'].str.lower().str.contains(q, na=False)
        )
        result = de[mask].copy()
        if result.empty:
            st.warning('No matching genes found.')
        else:
            st.success(f'{len(result)} matching gene(s) found.')
            cols = [c for c in ['GeneID', 'Symbol', 'Description', 'log2FoldChange', 'pvalue', 'padj', 'Significance'] if c in result.columns]
            st.dataframe(result[cols], use_container_width=True, hide_index=True)
    st.divider()
    n = st.slider('Number of genes to browse', 10, 200, 50, 10)
    cols = [c for c in ['GeneID', 'Symbol', 'Description', 'log2FoldChange', 'padj'] if c in de.columns]
    st.dataframe(de.sort_values('padj').head(n)[cols], use_container_width=True, hide_index=True)

# ---------- Validation ----------
with tabs[5]:
    st.header('🧪 Independent Validation Framework')
    st.markdown('''
    **Selected external cohorts**

    - **GSE293353:** bulk granulosa-cell RNA-seq; 9 PCOS and 9 controls.
    - **GSE193123:** bulk granulosa-cell RNA-seq; small cohort, useful but lower statistical power.
    - **GSE155489:** cumulus granulosa-cell RNA-seq; useful as an orthogonal cell-context validation.
    - **GSE240688:** single-cell granulosa-cell RNA-seq; useful for cell-state/localization validation.
    ''')
    st.warning('External cohorts must be analysed within their own study before replication is claimed. Do not blindly merge raw count matrices across studies.')
    st.markdown('**Recommended validation:** compare predefined candidate genes for direction consistency, rank correlation, pathway replication, and TF-target enrichment in independent cohorts.')

st.divider()
st.caption('PCOS Transcriptomics Portal • Human GRCh38 annotation • Research / hypothesis-generation tool')
