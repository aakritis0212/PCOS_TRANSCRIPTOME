# PCOS Transcriptomics — Independent Validation Report

## 1. Current analysis

The supplied GSE168404-derived results contain:

- 21,597 genes in the DESeq2 table.
- 85 genes meeting padj < 0.05 and |log2FC| >= 1 in the supplied result table.
- 58 upregulated and 27 downregulated genes under that threshold.
- 65 GO terms with FDR q < 0.05 in the supplied GSEA result.
- 4 TF signatures with FDR q < 0.05: FOXM1 ENCODE, EZH2 CHEA, IRF8 CHEA and EZH2 ENCODE.

## 2. Independent cohorts selected

### GSE293353
Recent bulk granulosa-cell RNA-seq cohort: 9 PCOS and 9 controls. The associated 2025 study reported 199 consistently dysregulated genes and identified GPX3 as a convergent regulator of oxidative-stress, insulin-signaling, glucose-metabolism and mitochondrial programs.

### GSE193123
Bulk granulosa-cell RNA-seq: three control and three PCOS samples, with samples pooled from two individuals for each library. This is useful as a technically independent RNA-seq cohort, but the small sample size and pooling make it weaker than GSE293353.

### GSE155489
RNA-seq study of cumulus granulosa cells with four PCOS and four control samples. It is useful as an orthogonal granulosa-cell validation cohort, but the cell context differs from follicular granulosa cells.

### GSE240688
Single-cell granulosa-cell RNA-seq with three PCOS and three controls. It should be used for cell-state/localization validation rather than treated as a simple bulk replication cohort.

## 3. Novelty assessment

The following are NOT safe novelty claims:

- steroid/cholesterol biosynthesis dysregulation in PCOS granulosa cells
- FOXM1 involvement in PCOS granulosa-cell proliferation
- EZH2 involvement in ovarian/granulosa-cell epigenetic biology
- oxidative-stress/GPX3 involvement in PCOS granulosa cells

These areas have already been reported.

The defensible research opportunity is therefore an **integrated regulatory-module hypothesis**:

DESeq2 evidence + pathway leading-edge genes + TF target enrichment + PPI/network structure + independent cohort replication.

The candidate module should only be called novel if its *specific combination and reproducibility* survives external validation and a focused literature search.

## 4. Required next statistical test

For each independent cohort:

1. Reprocess within study; do not merge raw matrices across platforms.
2. Obtain gene-level DE statistics.
3. Map identifiers to a common gene namespace.
4. Compare candidate genes by direction and rank.
5. Test pathway-level replication using the same gene-set definitions.
6. Test TF-target enrichment independently.
7. Calculate module-level reproducibility.
8. Use an independent cohort only for validation, not for tuning the candidate list.

## 5. Decision rule

A candidate regulatory module becomes a strong research candidate when:

- direction is consistent in at least 2 independent cohorts,
- pathway-level enrichment replicates,
- regulator-target overlap remains significant after multiple-testing correction,
- the module is not explained solely by one previously published hub gene,
- and the result is robust to reasonable ranking/statistical choices.

## 6. Important limitation

The current execution environment can search and verify public GEO metadata, but the external compressed count matrices could not be downloaded into the analysis runtime. Therefore, this report does **not** falsely claim that GSE293353/GSE193123/GSE155489 have already been computationally reanalyzed here. The package includes the exact download and validation scripts needed to perform that analysis locally.
