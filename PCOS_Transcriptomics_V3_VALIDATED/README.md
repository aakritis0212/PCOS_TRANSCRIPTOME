# PCOS Transcriptomics — Validation Package

This package contains the real supplied GSE168404 analysis outputs plus an independent-validation framework.

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

For external validation:

```bash
python scripts/download_external_validation.py
```

Then analyze each cohort separately and provide DE tables to `scripts/validate_external.py`.

## Selected external datasets

- GSE293353 — bulk granulosa-cell RNA-seq, 9 PCOS/9 controls
- GSE193123 — bulk granulosa-cell RNA-seq, 3 PCOS/3 controls
- GSE155489 — cumulus granulosa-cell RNA-seq, 4 PCOS/4 controls
- GSE240688 — single-cell granulosa-cell RNA-seq, 3 PCOS/3 controls

## Scientific rule

Do not merge cohorts blindly. Validate within-study and compare gene direction/rank and pathway-level results across cohorts.
