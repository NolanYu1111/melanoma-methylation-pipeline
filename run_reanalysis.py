import os
import gzip
import numpy as np
import pandas as pd

dl = r"C:\Users\ethan\Downloads"

def bh_fdr(pvals):
    """Benjamini-Hochberg FDR q-values."""
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

print("Step 1: Reading GPL13534 manifest to get promoter probes (TSS200 / TSS1500)...")
gpl_path = os.path.join(dl, "GPL13534-11288.txt")
gpl = pd.read_csv(gpl_path, sep="\t", skiprows=37, usecols=["ID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "CHR", "MAPINFO"], low_memory=False)
gpl = gpl.dropna(subset=["ID"])

def is_promoter(group_str):
    if not isinstance(group_str, str):
        return False
    parts = set(group_str.split(";"))
    return ("TSS200" in parts) or ("TSS1500" in parts)

def get_first_gene(gene_str):
    if not isinstance(gene_str, str) or not gene_str.strip():
        return "Unannotated"
    return gene_str.split(";")[0]

gpl["is_promoter"] = gpl["UCSC_RefGene_Group"].apply(is_promoter)
gpl["Gene"] = gpl["UCSC_RefGene_Name"].apply(get_first_gene)
promoter_df = gpl[gpl["is_promoter"]].copy()
promoter_probes = set(promoter_df["ID"])
probe_to_gene = dict(zip(gpl["ID"], gpl["Gene"]))
probe_to_group = dict(zip(gpl["ID"], gpl["UCSC_RefGene_Group"].fillna("")))
print(f"Total promoter probes in GPL13534: {len(promoter_probes)}")

print("\nStep 2: Preparing TCGA-SKCM clinical and survival data...")
tcga_surv = pd.read_csv(os.path.join(dl, "TCGA-SKCM.survival.tsv.gz"), sep="\t")
tcga_clin = pd.read_csv(os.path.join(dl, "TCGA-SKCM.clinical.tsv.gz"), sep="\t", low_memory=False)

# Check methylation sample headers
tcga_meth_path = os.path.join(dl, "TCGA-SKCM.methylation450.tsv.gz")
with gzip.open(tcga_meth_path, "rt", encoding="utf-8") as f:
    meth_samples = f.readline().strip().split("\t")[1:]
print(f"Total samples in TCGA methylation file: {len(meth_samples)}")

# Merge survival and clinical for methylation samples
meta = pd.DataFrame({"sample": meth_samples})
meta["sample_type_code"] = meta["sample"].apply(lambda s: s.split("-")[3] if len(s.split("-")) > 3 else "")
meta["patient"] = meta["sample"].apply(lambda s: "-".join(s.split("-")[:3]))

# Keep tumor samples: 01A primary; 06A/06B/07A metastatic (exclude 11A normal)
meta = meta[meta["sample_type_code"].isin(["01A", "06A", "06B", "07A"])].copy()
meta["is_primary"] = (meta["sample_type_code"] == "01A").astype(int)

# Merge with survival
meta = meta.merge(tcga_surv[["sample", "OS.time", "OS"]], on="sample", how="inner")
meta = meta[meta["OS.time"].notna() & (meta["OS.time"] > 0) & meta["OS"].notna()].copy()

# Merge with clinical
clin_cols = ["sample", "age_at_index.demographic", "gender.demographic", "ajcc_pathologic_stage.diagnoses"]
meta = meta.merge(tcga_clin[clin_cols], on="sample", how="left")
meta = meta[meta["age_at_index.demographic"].notna()].copy()

# Sort so primary (is_primary=1) is preferred per patient
meta = meta.sort_values(by=["patient", "is_primary", "sample"], ascending=[True, False, True])
meta = meta.drop_duplicates(subset=["patient"], keep="first").copy()

print(f"Final cohort size: {len(meta)} patients ({meta['is_primary'].sum()} primary, {(1 - meta['is_primary']).sum()} metastatic; {int(meta['OS'].sum())} deaths)")
