import os
import gzip
import numpy as np
import pandas as pd

dl = r"C:\Users\ethan\Downloads"

# 1. Check how many promoter probes in GPL13534 vs HM450.hg38.manifest
gpl = pd.read_csv(os.path.join(dl, "GPL13534-11288.txt"), sep="\t", skiprows=37, usecols=["ID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "CHR", "Probe_SNPs", "Probe_SNPs_10"], low_memory=False).dropna(subset=["ID"])

def is_prom(s):
    if not isinstance(s, str): return False
    parts = set(s.split(";"))
    return ("TSS200" in parts) or ("TSS1500" in parts)

gpl["is_prom"] = gpl["UCSC_RefGene_Group"].apply(is_prom)
print("GPL13534 total is_prom:", gpl["is_prom"].sum())
print("GPL13534 cg-only is_prom:", gpl[gpl["is_prom"] & gpl["ID"].str.startswith("cg")].shape[0])

# Check TCGA methylation file promoter probes
X_raw = np.load("tcga_promoter_X_raw.npy")
probes_df = pd.read_csv("tcga_promoter_probes.csv")
meta = pd.read_csv("tcga_456_meta.csv")

# Why did extract_tcga_promoter.py get 67,102 vs 67,061?
sd_ddof0 = np.nanstd(X_raw, axis=1, ddof=0)
sd_ddof1 = np.nanstd(X_raw, axis=1, ddof=1)
n_miss = np.isnan(X_raw).sum(axis=1)
print("Count with n_miss <= 91 and sd_ddof1 > 0.05:", ((n_miss <= 91) & (sd_ddof1 > 0.05)).sum())
print("Count with n_miss <= 91 and sd_ddof0 > 0.05:", ((n_miss <= 91) & (sd_ddof0 > 0.05)).sum())
print("Count with n_miss < 0.20*456 and sd_ddof0 > 0.05 and startswith('cg'):", ((n_miss <= 91) & (sd_ddof0 > 0.05) & probes_df["probe"].str.startswith("cg")).sum())

# 2. Check TCGA-SKCM.HiSeqV2.tsv.gz headers and sample overlap with meta
rnaseq_path = os.path.join(dl, "TCGA-SKCM.HiSeqV2.tsv.gz")
with gzip.open(rnaseq_path, "rt", encoding="utf-8") as f:
    rna_samples = f.readline().strip().split("\t")[1:]
print("Total RNA-seq samples:", len(rna_samples), "first 5:", rna_samples[:5])

# Note: HiSeqV2 uses 15-character sample IDs (e.g. TCGA-RP-A690-06) vs 16-character in methylation (TCGA-RP-A690-06A)!
meta["sample_15"] = meta["sample"].str.slice(0, 15)
rna_set = set(rna_samples)
print("Matched 456 survival cohort samples in RNA-seq (by 15-char barcode):", meta["sample_15"].isin(rna_set).sum())
