import os
import gzip
import time
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize
from sklearn.decomposition import PCA

dl = r"C:\Users\ethan\Downloads"

def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    out = np.empty_like(q)
    out[order] = q
    return out

# 1. Load GPL13534 promoter probes
print("Reading GPL13534...")
gpl_path = os.path.join(dl, "GPL13534-11288.txt")
gpl = pd.read_csv(gpl_path, sep="\t", skiprows=37, usecols=["ID", "UCSC_RefGene_Name", "UCSC_RefGene_Group"], low_memory=False)
gpl = gpl.dropna(subset=["ID"])

def is_promoter(group_str):
    if not isinstance(group_str, str):
        return False
    parts = set(group_str.split(";"))
    return ("TSS200" in parts) or ("TSS1500" in parts)

gpl["is_promoter"] = gpl["UCSC_RefGene_Group"].apply(is_promoter)
promoter_probes = set(gpl[gpl["is_promoter"]]["ID"])
probe_to_gene = dict(zip(gpl["ID"], gpl["UCSC_RefGene_Name"].fillna("")))
probe_to_group = dict(zip(gpl["ID"], gpl["UCSC_RefGene_Group"].fillna("")))

# 2. Prepare 456 patients
tcga_surv = pd.read_csv(os.path.join(dl, "TCGA-SKCM.survival.tsv.gz"), sep="\t")
tcga_clin = pd.read_csv(os.path.join(dl, "TCGA-SKCM.clinical.tsv.gz"), sep="\t", low_memory=False)
tcga_meth_path = os.path.join(dl, "TCGA-SKCM.methylation450.tsv.gz")
with gzip.open(tcga_meth_path, "rt", encoding="utf-8") as f:
    header_cols = f.readline().strip().split("\t")
    meth_samples = header_cols[1:]

meta = pd.DataFrame({"sample": meth_samples, "col_idx": range(1, len(header_cols))})
meta["sample_type_code"] = meta["sample"].apply(lambda s: s.split("-")[3] if len(s.split("-")) > 3 else "")
meta["patient"] = meta["sample"].apply(lambda s: "-".join(s.split("-")[:3]))
meta = meta[meta["sample_type_code"].isin(["01A", "06A", "06B", "07A"])].copy()
meta["is_primary"] = (meta["sample_type_code"] == "01A").astype(int)
meta = meta.merge(tcga_surv[["sample", "OS.time", "OS"]], on="sample", how="inner")
meta = meta[meta["OS.time"].notna() & (meta["OS.time"] > 0) & meta["OS"].notna()].copy()
clin_cols = ["sample", "age_at_index.demographic", "gender.demographic", "ajcc_pathologic_stage.diagnoses"]
meta = meta.merge(tcga_clin[clin_cols], on="sample", how="left")
meta = meta[meta["age_at_index.demographic"].notna()].copy()
meta = meta.sort_values(by=["patient", "is_primary", "sample"], ascending=[True, False, True])
meta = meta.drop_duplicates(subset=["patient"], keep="first").copy()
meta = meta.reset_index(drop=True)

print("Stage value counts:")
print(meta["ajcc_pathologic_stage.diagnoses"].value_counts(dropna=False))

def clean_stage(s):
    if not isinstance(s, str):
        return "Unknown"
    s = s.strip()
    if s in ["Stage IV"]:
        return "IV"
    if s.startswith("Stage III"):
        return "III"
    if s.startswith("Stage II"):
        return "II"
    if s.startswith("Stage I"):
        return "I"
    if s == "Stage 0":
        return "0"
    return "Unknown"

meta["stage_clean"] = meta["ajcc_pathologic_stage.diagnoses"].apply(clean_stage)
print("Cleaned stage counts:")
print(meta["stage_clean"].value_counts())

# Stream through TCGA methylation file and keep promoter probes for the 456 samples
col_indices = meta["col_idx"].values
n_samples = len(col_indices)
max_missing = int(0.20 * n_samples) # at most 20% missing (<= 91.2 -> 91)

t0 = time.time()
probe_ids = []
matrix_rows = []
total_promoter_in_file = 0

with gzip.open(tcga_meth_path, "rt", encoding="utf-8") as f:
    _ = f.readline() # skip header
    for line in f:
        tab_pos = line.find("\t")
        pid = line[:tab_pos]
        if pid not in promoter_probes:
            continue
        total_promoter_in_file += 1
        parts = line.rstrip("\r\n").split("\t")
        vals = [float(parts[i]) if parts[i] not in ("", "NA", "nan") else np.nan for i in col_indices]
        arr = np.array(vals, dtype=np.float64)
        n_miss = np.isnan(arr).sum()
        if n_miss <= 0.20 * n_samples:
            sd = np.nanstd(arr, ddof=1)
            if sd > 0.05:
                probe_ids.append(pid)
                matrix_rows.append(arr)

print(f"Total promoter CpGs in TCGA file: {total_promoter_in_file}")
print(f"Tested promoter CpGs (<=20% missing, SD(ddof=1)>0.05): {len(probe_ids)} in {time.time()-t0:.1f}s")

# Save intermediate filtered promoter matrix so we don't have to re-parse the 1.6GB file!
X_raw = np.vstack(matrix_rows) # shape: (n_probes, 456)
np.save("tcga_promoter_X_raw.npy", X_raw)
pd.DataFrame({"probe": probe_ids}).to_csv("tcga_promoter_probes.csv", index=False)
meta.to_csv("tcga_456_meta.csv", index=False)
print("Saved tcga_promoter_X_raw.npy and metadata!")
