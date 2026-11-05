import os
import gzip
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

dl = r"C:\Users\ethan\Downloads"

X_raw = np.load("tcga_promoter_X_raw.npy")
row_means = np.nanmean(X_raw, axis=1, keepdims=True)
X_imp = np.where(np.isnan(X_raw), row_means, X_raw)
row_sds = np.std(X_imp, axis=1, ddof=1, keepdims=True)
X_std = (X_imp - row_means) / row_sds

probes_df = pd.read_csv("tcga_promoter_probes.csv")
meta = pd.read_csv("tcga_456_meta.csv")
cox_pc5 = pd.read_csv("cox_pc5.tsv", sep="\t")

gpl = pd.read_csv(
    os.path.join(dl, "GPL13534-11288.txt"),
    sep="\t",
    skiprows=37,
    usecols=["ID", "CHR", "Probe_SNPs", "Probe_SNPs_10"],
    low_memory=False
).dropna(subset=["ID"])

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

cox_merged = cox_pc5.merge(gpl, left_on="probe", right_on="ID", how="left")
has_snp_main = cox_merged["Probe_SNPs"].notna() & (cox_merged["Probe_SNPs"].astype(str).str.strip() != "")
has_snp_10 = cox_merged["Probe_SNPs_10"].notna() & (cox_merged["Probe_SNPs_10"].astype(str).str.strip() != "")
clean_mask = (~has_snp_main & ~has_snp_10).values
cox_clean = cox_merged[clean_mask].copy()
cox_clean["q_pc5_clean"] = bh_fdr(cox_clean["p_pc5"].values)
clean_q_map = dict(zip(cox_clean["probe"], cox_clean["q_pc5_clean"]))
snp_map = dict(zip(cox_merged["probe"], cox_merged["Probe_SNPs"].fillna("None")))
snp10_map = dict(zip(cox_merged["probe"], cox_merged["Probe_SNPs_10"].fillna("None")))

pca = PCA(n_components=5, random_state=42)
pcs = pca.fit_transform(X_std.T)

rnaseq_path = os.path.join(dl, "TCGA-SKCM.HiSeqV2.tsv.gz")
immune_14 = ["PTPRC", "CD8A", "CD8B", "CD3D", "CD3E", "GZMA", "PRF1", "IFNG", "CXCL9", "CXCL10", "STAT1", "PSMB9", "TAP1", "HLA-DRA"]
target_genes = set(immune_14 + ["APOL3", "LOC286094", "CRABP2", "KRAS", "CPA6", "PDE2A", "PER2", "IGFL4", "GABRB3", "CX3CR1", "DYSFIP1", "ZCCHC11"])

expr_dict = {}
with gzip.open(rnaseq_path, "rt", encoding="utf-8") as f:
    rna_samples = f.readline().strip().split("\t")[1:]
    for line in f:
        gene = line[:line.find("\t")]
        if gene in target_genes:
            parts = line.rstrip("\r\n").split("\t")
            expr_dict[gene] = np.array([float(x) if x not in ("", "NA", "nan") else np.nan for x in parts[1:]])

rna_df = pd.DataFrame(expr_dict, index=rna_samples)
meta["sample_15"] = meta["sample"].str.slice(0, 15)
has_rna = meta["sample_15"].isin(rna_df.index).values
sub_samples = meta.loc[has_rna, "sample_15"].values

sub_imm = rna_df.loc[sub_samples, immune_14].values
z_imm = (sub_imm - np.mean(sub_imm, axis=0)) / np.std(sub_imm, axis=0, ddof=1)
score_14 = np.mean(z_imm, axis=1)
score_14 = (score_14 - np.mean(score_14)) / np.std(score_14, ddof=1)
pcs_sub = pcs[has_rna]

def partial_r(x, y, covs):
    C = np.column_stack([np.ones(len(x)), covs])
    rx = x - C @ np.linalg.lstsq(C, x, rcond=None)[0]
    ry = y - C @ np.linalg.lstsq(C, y, rcond=None)[0]
    return stats.pearsonr(rx, ry)

p_map = {p: i for i, p in enumerate(probes_df["probe"])}
t2_probes = ["cg06703343", "cg25053164", "cg02363950", "cg23152216", "cg21860285", "cg14790396", "cg05664072", "cg13095627", "cg15245051", "cg14222314", "cg24130739", "cg17720013", "cg03750525", "cg00754357"]

rows = []
for pid in t2_probes:
    r_cox = cox_pc5[cox_pc5["probe"] == pid].iloc[0]
    gname = r_cox["Gene"]
    idx = p_map[pid]
    m_vals = X_imp[idx, has_rna]
    if gname in rna_df.columns:
        e_vals = rna_df.loc[sub_samples, gname].values
        pr, pp = stats.pearsonr(m_vals, e_vals)
        sr, sp = stats.spearmanr(m_vals, e_vals)
        prim_r, prim_p = partial_r(m_vals, e_vals, score_14)
        ppc5_r, ppc5_p = partial_r(m_vals, e_vals, pcs_sub[:, :5])
    else:
        pr = pp = sr = sp = prim_r = prim_p = ppc5_r = ppc5_p = np.nan

    rows.append({
        "probe": pid,
        "Gene": gname,
        "HR_clin": r_cox["HR_clin"],
        "p_clin": r_cox["p_clin"],
        "q_clin": r_cox["q_clin"],
        "HR_pc5": r_cox["HR_pc5"],
        "p_pc5": r_cox["p_pc5"],
        "q_pc5": r_cox["q_pc5"],
        "q_pc5_noSNP": clean_q_map.get(pid, np.nan),
        "Probe_SNPs_10": snp10_map.get(pid, "None") if snp10_map.get(pid, "None") != "" else "None",
        "Pearson_r": pr,
        "Pearson_p": pp,
        "Spearman_rho": sr,
        "Spearman_p": sp,
        "Partial_r_Immune": prim_r,
        "Partial_p_Immune": prim_p,
        "Partial_r_PC5": ppc5_r,
        "Partial_p_PC5": ppc5_p
    })

out_df = pd.DataFrame(rows)
print(out_df.to_string(index=False))
out_df.to_csv("candidate_meth_expr_corr.tsv", sep="\t", index=False)
out_df.to_csv(os.path.join(dl, "candidate_meth_expr_corr.tsv"), sep="\t", index=False)
