import requests
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "external_validation"
OUT.mkdir(exist_ok=True)

FILES = {
    "GSE293353_counts.txt.gz":
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE293353&file=GSE293353_gene_count_matrix.txt.gz&format=file",
    "GSE293353_metadata.csv.gz":
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE293353&file=GSE293353_Sample_Information.csv.gz&format=file",
    "GSE193123_counts.txt.gz":
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE193123&file=GSE193123_gene_count.txt.gz&format=file",
    "GSE155489_gc_counts.csv.gz":
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE155489&file=GSE155489_gc_pcos_counts.csv.gz&format=file",
}

for name, url in FILES.items():
    target = OUT / name
    print("Downloading", name)
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    target.write_bytes(r.content)
    print("Saved", target, len(r.content), "bytes")

print("External validation datasets downloaded.")
