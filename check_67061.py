import os
import gzip
import numpy as np
import pandas as pd

dl = r"C:\Users\ethan\Downloads"

# Check HM450.hg38.manifest.gencode.v36.probeMap vs GPL13534-11288.txt
gpl = pd.read_csv(os.path.join(dl, "GPL13534-11288.txt"), sep="\t", skiprows=37, usecols=["ID", "UCSC_RefGene_Name", "UCSC_RefGene_Group"], low_memory=False).dropna(subset=["ID"])
def is_prom(s):
    if not isinstance(s, str): return False
    parts = set(s.split(";"))
    return ("TSS200" in parts) or ("TSS1500" in parts)

gpl["is_prom"] = gpl["UCSC_RefGene_Group"].apply(is_prom)
prom_gpl = set(gpl[gpl["is_prom"]]["ID"])
prom_gpl_cg = set(gpl[gpl["is_prom"] & gpl["ID"].str.startswith("cg")]["ID"])

X_raw = np.load("tcga_promoter_X_raw.npy")
probes_df = pd.read_csv("tcga_promoter_probes.csv")

# Wait: in extract_tcga_promoter.py, how many total_promoter_in_file were there?
# Let's check with cg-only and different SD / missing thresholds on X_raw:
is_cg = probes_df["probe"].str.startswith("cg").values
X_cg = X_raw[is_cg]
n_miss = np.isnan(X_cg).sum(axis=1)
n_obs = 456 - n_miss
sd0 = np.nanstd(X_cg, axis=1, ddof=0)
sd1 = np.nanstd(X_cg, axis=1, ddof=1)

for m_thresh in [91, 90, 92]:
    print(f"n_miss <= {m_thresh}, sd1 > 0.05 (cg only):", ((n_miss <= m_thresh) & (sd1 > 0.05)).sum())
    print(f"n_miss <= {m_thresh}, sd0 > 0.05 (cg only):", ((n_miss <= m_thresh) & (sd0 > 0.05)).sum())
    print(f"n_miss <= {m_thresh}, sd1 >= 0.05 (cg only):", ((n_miss <= m_thresh) & (sd1 >= 0.05)).sum())
    print(f"n_miss <= {m_thresh}, sd0 >= 0.05 (cg only):", ((n_miss <= m_thresh) & (sd0 >= 0.05)).sum())
