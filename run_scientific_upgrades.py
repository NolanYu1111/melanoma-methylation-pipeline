import os
import gzip
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score, accuracy_score

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

# Load TCGA 456 cohort & promoter methylation matrix
X_raw = np.load("tcga_promoter_X_raw.npy")
probes_df = pd.read_csv("tcga_promoter_probes.csv")
meta = pd.read_csv("tcga_456_meta.csv")

# Filter to exact 67,061 cg-only promoter probes (exclude ch.* probes and keep top 67,061 by SD > 0.05)
gpl = pd.read_csv(os.path.join(dl, "GPL13534-11288.txt"), sep="\t", skiprows=37, usecols=["ID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "CHR", "Probe_SNPs", "Probe_SNPs_10"], low_memory=False).dropna(subset=["ID"])
probe_to_gene = dict(zip(gpl["ID"], gpl["UCSC_RefGene_Name"].fillna("").apply(lambda s: s.split(";")[0] if s else "Unannotated")))
probe_to_allgenes = dict(zip(gpl["ID"], gpl["UCSC_RefGene_Name"].fillna("")))
probe_has_snp = set(gpl[gpl["Probe_SNPs"].notna() | gpl["Probe_SNPs_10"].notna()]["ID"])

probes_df["Gene"] = probes_df["probe"].map(probe_to_gene)
probes_df["AllGenes"] = probes_df["probe"].map(probe_to_allgenes)
probes_df["has_snp"] = probes_df["probe"].isin(probe_has_snp)

# Impute missing values with row mean and standardize
row_means = np.nanmean(X_raw, axis=1, keepdims=True)
X_imp = np.where(np.isnan(X_raw), row_means, X_raw)
row_sds = np.std(X_imp, axis=1, ddof=1, keepdims=True)
X_std = (X_imp - row_means) / row_sds

# Compute top 5 PCs
pca = PCA(n_components=5, random_state=42)
pcs = pca.fit_transform(X_std.T) # (456, 5)

# Load TCGA-SKCM.HiSeqV2.tsv.gz (RNA-seq log2(norm_count+1))
print("1. Loading TCGA-SKCM.HiSeqV2.tsv.gz RNA-seq matrix...", flush=True)
rnaseq_df = pd.read_csv(os.path.join(dl, "TCGA-SKCM.HiSeqV2.tsv.gz"), sep="\t", index_col=0)
print("RNA-seq shape:", rnaseq_df.shape, flush=True)

# Match with the 456 patients via 15-character sample barcode (e.g., TCGA-RP-A690-06)
meta["sample_15"] = meta["sample"].str.slice(0, 15)
has_rna_mask = meta["sample_15"].isin(rnaseq_df.columns).values
rna_samples_matched = meta.loc[has_rna_mask, "sample_15"].values
print(f"Matched {has_rna_mask.sum()} of 456 patients with RNA-seq expression data.", flush=True)

rna_matched_df = rnaseq_df[rna_samples_matched] # (20530, 454)

# Compute Canonical Immune Infiltration Score from established T-cell / pan-leukocyte / IFN-gamma / cytolytic marker genes:
# PTPRC (CD45), CD8A, CD8B, CD3D, CD3E, GZMA, PRF1, IFNG, CXCL9, CXCL10, STAT1, PSMB9, TAP1, HLA-DRA
immune_markers = ["PTPRC", "CD8A", "CD8B", "CD3D", "CD3E", "GZMA", "PRF1", "IFNG", "CXCL9", "CXCL10", "STAT1", "PSMB9", "TAP1", "HLA-DRA"]
avail_markers = [g for g in immune_markers if g in rna_matched_df.index]
print("Available immune marker genes in RNA-seq:", avail_markers, flush=True)

immune_expr = rna_matched_df.loc[avail_markers].values # (14, 454)
# Standardize each marker gene across the 454 samples and average to form the composite RNA-seq Immune Score
immune_z = (immune_expr - immune_expr.mean(axis=1, keepdims=True)) / immune_expr.std(axis=1, ddof=1, keepdims=True)
immune_score_454 = immune_z.mean(axis=0)
immune_score_454_std = (immune_score_454 - np.mean(immune_score_454)) / np.std(immune_score_454, ddof=1)

# Also compute Cytolytic Score (mean of GZMA and PRF1, Rooney et al. Cell 2015) and Pan-Leukocyte PTPRC (CD45)
cytolytic_454 = rna_matched_df.loc[["GZMA", "PRF1"]].mean(axis=0).values

# Correlate Promoter Methylation PCs (PC1 - PC5) with RNA-seq Immune Score, Cytolytic Score, and PTPRC
print("\n--- Correlation of Methylation PCs with RNA-seq Immune Score (n = 454) ---", flush=True)
pcs_454 = pcs[has_rna_mask]
for k in range(5):
    r_imm, p_imm = stats.pearsonr(pcs_454[:, k], immune_score_454_std)
    r_cyt, p_cyt = stats.pearsonr(pcs_454[:, k], cytolytic_454)
    r_cd45, p_cd45 = stats.pearsonr(pcs_454[:, k], rna_matched_df.loc["PTPRC"].values)
    print(f"PC{k+1}: r(Immune Score) = {r_imm:+.4f} (p = {p_imm:.2e}), r(Cytolytic) = {r_cyt:+.4f}, r(PTPRC/CD45) = {r_cd45:+.4f}", flush=True)

# Multiple R^2 of PC1-3 and PC1-5 predicting the RNA-seq Immune Score
beta_ols3, _, _, _ = np.linalg.lstsq(np.column_stack([ np.ones(len(pcs_454)), pcs_454[:, :3] ]), immune_score_454_std, rcond=None)
pred3 = np.column_stack([ np.ones(len(pcs_454)), pcs_454[:, :3] ]) @ beta_ols3
r2_3 = np.corrcoef(pred3, immune_score_454_std)[0, 1] ** 2

beta_ols5, _, _, _ = np.linalg.lstsq(np.column_stack([ np.ones(len(pcs_454)), pcs_454[:, :5] ]), immune_score_454_std, rcond=None)
pred5 = np.column_stack([ np.ones(len(pcs_454)), pcs_454[:, :5] ]) @ beta_ols5
r2_5 = np.corrcoef(pred5, immune_score_454_std)[0, 1] ** 2
print(f"Multiple R of PC1-3 with RNA-seq Immune Score: R = {np.sqrt(r2_3):.4f} (R^2 = {r2_3:.4f})", flush=True)
print(f"Multiple R of PC1-5 with RNA-seq Immune Score: R = {np.sqrt(r2_5):.4f} (R^2 = {r2_5:.4f})", flush=True)

# Now test how many of the clinical-only FDR < 0.05 hits disappear when adjusting for RNA-seq Immune Score!
# Impute the 2 samples without RNA-seq with 0.0 (mean) or run on the 454 matched samples
immune_full = np.zeros(len(meta))
immune_full[has_rna_mask] = immune_score_454_std

# Also create a 3-component RNA-seq immune basis (T-cell CD8A/CD3D/GZMA/PRF1, IFN-gamma STAT1/CXCL9/CXCL10/PSMB9/TAP1, Leukocyte PTPRC/HLA-DRA)
tcell_score = immune_z[[avail_markers.index(g) for g in ["CD8A", "CD8B", "CD3D", "CD3E", "GZMA", "PRF1"]]].mean(axis=0)
ifn_score = immune_z[[avail_markers.index(g) for g in ["IFNG", "CXCL9", "CXCL10", "STAT1", "PSMB9", "TAP1"]]].mean(axis=0)
leuk_score = immune_z[[avail_markers.index(g) for g in ["PTPRC", "HLA-DRA"]]].mean(axis=0)
rna_imm3_full = np.zeros((len(meta), 3))
rna_imm3_full[has_rna_mask, 0] = (tcell_score - tcell_score.mean()) / tcell_score.std(ddof=1)
rna_imm3_full[has_rna_mask, 1] = (ifn_score - ifn_score.mean()) / ifn_score.std(ddof=1)
rna_imm3_full[has_rna_mask, 2] = (leuk_score - leuk_score.mean()) / leuk_score.std(ddof=1)

# Build clinical matrix Z_clin
meta["sex_male"] = (meta["gender.demographic"].str.lower() == "male").astype(float)
age = meta["age_at_index.demographic"].values.astype(float)
age_std = (age - np.mean(age)) / np.std(age, ddof=1)
is_prim = meta["is_primary"].values.astype(float)
stage_dummies = pd.get_dummies(meta["stage_clean"], drop_first=True, dtype=float)
Z_clin = np.column_stack([age_std, meta["sex_male"].values, is_prim, stage_dummies.values])

time_arr = meta["OS.time"].values.astype(float)
event_arr = meta["OS"].values.astype(int)
order = np.argsort(time_arr)
t_s = time_arr[order]
d_s = event_arr[order]
Z_clin_s = Z_clin[order]
pcs_s = pcs[order]
imm1_s = immune_full[order]
imm3_s = rna_imm3_full[order]
X_std_s = X_std[:, order]

from run_cox_and_ml import fit_breslow_null_and_score_test

# Check association of RNA-seq Immune Score itself with OS adjusting for clinical covariates:
z_imm, p_imm_os, hr_imm_os, _ = fit_breslow_null_and_score_test(Z_clin_s, imm1_s[None, :], t_s, d_s)
print(f"\nRNA-seq Immune Score association with OS (clinical-adjusted): HR = {hr_imm_os[0]:.4f}, z = {z_imm[0]:.2f}, p = {p_imm_os[0]:.4e}", flush=True)

# Run genome-wide Cox screen adjusting for Clinical + RNA-seq Immune Score (1 score and 3 sub-scores)
z_ci1, p_ci1, hr_ci1, _ = fit_breslow_null_and_score_test(np.column_stack([Z_clin_s, imm1_s]), X_std_s, t_s, d_s)
q_ci1 = bh_fdr(p_ci1)
z_ci3, p_ci3, hr_ci3, _ = fit_breslow_null_and_score_test(np.column_stack([Z_clin_s, imm3_s]), X_std_s, t_s, d_s)
q_ci3 = bh_fdr(p_ci3)
print(f"Genome-wide hits (q < 0.05) when adjusting for Clinical + Single RNA-seq Immune Score: {(q_ci1 < 0.05).sum()} (p < 0.05: {(p_ci1 < 0.05).sum()})", flush=True)
print(f"Genome-wide hits (q < 0.05) when adjusting for Clinical + 3 RNA-seq Immune Sub-scores: {(q_ci3 < 0.05).sum()} (p < 0.05: {(p_ci3 < 0.05).sum()})", flush=True)

# Item 2: Methylation-Expression Correlation for Top Candidate Loci (KRAS, ZCCHC11, APOL3, CRABP2, CPA6, PDE2A, PER2, IGFL4, GABRB3, CX3CR1)
print("\n--- Item 2: Promoter Methylation vs. Gene Expression Correlations (n = 454) ---", flush=True)
cand_probes = [
    ("cg06703343", "APOL3"),
    ("cg25053164", "LOC286094"),
    ("cg02363950", "CRABP2"),
    ("cg23152216", "KRAS"),
    ("cg21860285", "CPA6"),
    ("cg14790396", "PDE2A"),
    ("cg05664072", "PER2"),
    ("cg13095627", "IGFL4"),
    ("cg15245051", "APOL3"),
    ("cg14222314", "GABRB3"),
    ("cg24130739", "CX3CR1"),
    ("cg17720013", "IGFL4"),
    ("cg00754357", "ZCCHC11"),
]

meth_expr_rows = []
for pid, gsym in cand_probes:
    p_idx = probes_df[probes_df["probe"] == pid].index[0]
    m_vals = X_raw[p_idx, has_rna_mask]
    valid = ~np.isnan(m_vals)
    # Check if gsym is in rna_matched_df (note ZCCHC11 is also known as TUT4)
    rna_gene = gsym
    if rna_gene not in rna_matched_df.index and gsym == "ZCCHC11" and "TUT4" in rna_matched_df.index:
        rna_gene = "TUT4"
    if rna_gene in rna_matched_df.index:
        e_vals = rna_matched_df.loc[rna_gene].values[valid]
        r_p, p_p = stats.pearsonr(m_vals[valid], e_vals)
        r_s, p_s = stats.spearmanr(m_vals[valid], e_vals)
        meth_expr_rows.append({
            "probe": pid, "Gene": gsym, "RNA_Gene": rna_gene,
            "Pearson_r": r_p, "Pearson_p": p_p,
            "Spearman_rho": r_s, "Spearman_p": p_s,
            "Mean_Beta": np.mean(m_vals[valid]), "SD_Beta": np.std(m_vals[valid], ddof=1),
            "Mean_Log2Expr": np.mean(e_vals)
        })
    else:
        meth_expr_rows.append({
            "probe": pid, "Gene": gsym, "RNA_Gene": "Not in HiSeqV2",
            "Pearson_r": np.nan, "Pearson_p": np.nan,
            "Spearman_rho": np.nan, "Spearman_p": np.nan,
            "Mean_Beta": np.mean(m_vals[valid]), "SD_Beta": np.std(m_vals[valid], ddof=1),
            "Mean_Log2Expr": np.nan
        })

meth_expr_df = pd.DataFrame(meth_expr_rows)
print(meth_expr_df.to_string(index=False), flush=True)
meth_expr_df.to_csv("candidate_meth_expr_corr.tsv", sep="\t", index=False)
