"""
External validation framework.

This script is intentionally conservative. It does not merge raw cohorts blindly.
It expects each external dataset to be analyzed within its own study, then compares
the direction/rank of a predefined candidate module across studies.

Usage:
  1. Run download_external_validation.py
  2. Perform dataset-specific preprocessing/DE analysis.
  3. Save a table with columns: Gene, log2FC, padj for each cohort.
  4. Point this script to those tables.

This prevents cross-platform batch effects from being mistaken for biology.
"""
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CAND = pd.read_csv(ROOT/"results/current_leading_edge_candidates.csv")

tables = {
    # Fill these after dataset-specific DE analysis:
    # "GSE293353": ROOT/"external_validation/GSE293353_DE.csv",
    # "GSE193123": ROOT/"external_validation/GSE193123_DE.csv",
    # "GSE155489": ROOT/"external_validation/GSE155489_DE.csv",
}

available = []
for name, path in tables.items():
    if path.exists():
        d = pd.read_csv(path)
        d["Gene"] = d["Gene"].astype(str)
        available.append((name, d))

if not available:
    print("No external DE tables found yet.")
    print("Run the download script and perform study-specific DE analysis first.")
    raise SystemExit(0)

out = CAND[["Gene","Evidence_count"]].copy()
for name, d in available:
    keep = d[["Gene","log2FC","padj"]].copy()
    keep = keep.drop_duplicates("Gene")
    keep[name+"_log2FC"] = keep["log2FC"]
    keep[name+"_padj"] = keep["padj"]
    out = out.merge(keep[["Gene",name+"_log2FC",name+"_padj"]], on="Gene", how="left")

fc_cols = [c for c in out.columns if c.endswith("_log2FC")]
if fc_cols:
    out["direction_consistency"] = out[fc_cols].apply(
        lambda r: (np.sign(r.dropna()) > 0).sum() if len(r.dropna()) else np.nan, axis=1
    )

out.to_csv(ROOT/"results/external_validation_candidates.csv", index=False)
print("Saved external validation table.")
