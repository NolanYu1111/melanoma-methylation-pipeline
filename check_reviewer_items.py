import os
import gzip
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

dl = r"C:\Users\ethan\Downloads"

# 1. Load GPL13534 annotations including CHR, Probe_SNPs, Probe_SNPs_10
gpl = pd.read_csv(
    os.path.join(dl, "GPL13534-11288.txt"),
    sep="\t",
    skiprows=37,
    usecols=["ID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "CHR", "Probe_SNPs", "Probe_SNPs_10"],
    low_memory=False
).dropna(subset=["ID"])

X_raw = np.load("tcga_promoter_X_raw.npy")
probes_df = pd.read_csv("tcga_promoter_probes.csv")
meta = pd.read_csv("tcga_456_meta.csv")
cox_pc5 = pd.read_csv("cox_pc5.tsv", sep="\t")

# Merge GPL annotations onto cox_pc5 by probe ID
cox_merged = cox_pc5.merge(gpl, left_on="probe", right_on="ID", how="left")
has_snp_main = cox_merged["Probe_SNPs"].notna() & (cox_merged["Probe_SNPs"].astype(str).str.strip() != "")
has_snp_10 = cox_merged["Probe_SNPs_10"].notna() & (cox_merged["Probe_SNPs_10"].astype(str).str.strip() != "")
is_sex_chr = cox_merged["CHR"].astype(str).isin(["X", "Y"])

print("=== 1. PROBE FILTERING & GENOME-WIDE SNP CHECKS ===")
print(f"Total tested promoter CpGs (<=20% missing, observed sample SD > 0.05): {len(cox_merged)}")
row_means = np.nanmean(X_raw, axis=1, keepdims=True)
X_imp = np.where(np.isnan(X_raw), row_means, X_raw)
sd_after_imp_ddof1 = np.std(X_imp, axis=1, ddof=1)
sd_after_imp_ddof0 = np.std(X_imp, axis=1, ddof=0)
print(f"If SD > 0.05 is applied AFTER mean imputation: ddof=1 gives {(sd_after_imp_ddof1 > 0.05).sum()}, ddof=0 gives {(sd_after_imp_ddof0 > 0.05).sum()}")
print(f"Probes with primary Probe_SNPs flag: {has_snp_main.sum()}")
print(f"Probes with Probe_SNPs_10 flag: {has_snp_10.sum()}")
print(f"Probes with EITHER Probe_SNPs or Probe_SNPs_10: {(has_snp_main | has_snp_10).sum()}")
print(f"Probes on X/Y chromosomes: {is_sex_chr.sum()}")
print(f"Autosomal non-polymorphic probes (no Probe_SNPs, no Probe_SNPs_10, autosomal): {(~has_snp_main & ~has_snp_10 & ~is_sex_chr).sum()}")

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

clean_mask = (~has_snp_main & ~has_snp_10).values
cox_clean = cox_merged[clean_mask].copy()
cox_clean["q_clin_clean"] = bh_fdr(cox_clean["p_clin"].values)
cox_clean["q_pc3_clean"] = bh_fdr(cox_clean["p_pc3"].values)
cox_clean["q_pc5_clean"] = bh_fdr(cox_clean["p_pc5"].values)
print(f"After genome-wide removal of all Probe_SNPs & Probe_SNPs_10 ({clean_mask.sum()} probes):")
print(f"  Clinical-only q < 0.05: {(cox_clean['q_clin_clean'] < 0.05).sum()}")
print(f"  3-PC q < 0.05: {(cox_clean['q_pc3_clean'] < 0.05).sum()}")
print(f"  5-PC q < 0.05: {(cox_clean['q_pc5_clean'] < 0.05).sum()}")
print("  Top 5-PC hits after removing all Probe_SNPs and Probe_SNPs_10:")
print(cox_clean.sort_values("p_pc5")[["probe", "Gene", "HR_pc5", "p_pc5", "q_pc5_clean"]].head(12).to_string(index=False))

# 2. What does PC1 track?
print("\n=== 2. WHAT DOES PC1 TRACK? ===")
row_sds = np.std(X_imp, axis=1, ddof=1, keepdims=True)
X_std = (X_imp - row_means) / row_sds
pca = PCA(n_components=5, random_state=42)
pcs = pca.fit_transform(X_std.T)

global_mean_beta = np.mean(X_imp, axis=0)
global_median_beta = np.median(X_imp, axis=0)
for i in range(5):
    r_mean, p_mean = stats.pearsonr(pcs[:, i], global_mean_beta)
    r_prim, p_prim = stats.pointbiserialr(meta["is_primary"].values, pcs[:, i])
    r_age, p_age = stats.pearsonr(meta["age_at_index.demographic"].values, pcs[:, i])
    r_sex, p_sex = stats.pointbiserialr((meta["gender.demographic"].str.lower() == "male").astype(int), pcs[:, i])
    print(f"PC{i+1} ({pca.explained_variance_ratio_[i]*100:.2f}% var): r(global_mean_beta)={r_mean:+.4f} (p={p_mean:.2e}), r(primary)={r_prim:+.4f} (p={p_prim:.2e}), r(age)={r_age:+.4f} (p={p_age:.2e}), r(male)={r_sex:+.4f}")

# 3. Check RNA-seq immune score WITHOUT PSMB9 and TAP1 (and without any gene in Table 2)
print("\n=== 3. RNA-SEQ IMMUNE SCORE WITHOUT PSMB9 / TAP1 (NON-CIRCULARITY CHECK) ===")
rnaseq_path = os.path.join(dl, "TCGA-SKCM.HiSeqV2.tsv.gz")
immune_14 = ["PTPRC", "CD8A", "CD8B", "CD3D", "CD3E", "GZMA", "PRF1", "IFNG", "CXCL9", "CXCL10", "STAT1", "PSMB9", "TAP1", "HLA-DRA"]
immune_12_no_psmb9_tap1 = ["PTPRC", "CD8A", "CD8B", "CD3D", "CD3E", "GZMA", "PRF1", "IFNG", "CXCL9", "CXCL10", "STAT1", "HLA-DRA"]
# Also check pure T-cell / cytolytic 6-gene score (CD8A, CD8B,CD3D, CD3E, GZMA, PRF1)
immune_6_cytolytic = ["CD8A", "CD8B", "CD3D", "CD3E", "GZMA", "PRF1"]

expr_dict = {}
with gzip.open(rnaseq_path, "rt", encoding="utf-8") as f:
    rna_samples = f.readline().strip().split("\t")[1:]
    for line in f:
        gene = line[:line.find("\t")]
        if gene in set(immune_14 + ["CRABP2", "KRAS", "APOL3", "ZCCHC11"]):
            parts = line.rstrip("\r\n").split("\t")
            expr_dict[gene] = np.array([float(x) if x not in ("", "NA", "nan") else np.nan for x in parts[1:]])

rna_df = pd.DataFrame(expr_dict, index=rna_samples)
meta["sample_15"] = meta["sample"].str.slice(0, 15)
has_rna = meta["sample_15"].isin(rna_df.index).values
sub_samples = meta.loc[has_rna, "sample_15"].values

def make_score(genes):
    sub = rna_df.loc[sub_samples, genes].values
    z = (sub - np.mean(sub, axis=0)) / np.std(sub, axis=0, ddof=1)
    s = np.mean(z, axis=1)
    return (s - np.mean(s)) / np.std(s, ddof=1)

score_14 = make_score(immune_14)
score_12 = make_score(immune_12_no_psmb9_tap1)
score_6 = make_score(immune_6_cytolytic)

pcs_sub = pcs[has_rna]
for name, sc in [("14-gene", score_14), ("12-gene (no PSMB9/TAP1)", score_12), ("6-gene (CD3/CD8/GZMA/PRF1)", score_6)]:
    r_pc2, p_pc2 = stats.pearsonr(pcs_sub[:, 1], sc)
    r_pc1, p_pc1 = stats.pearsonr(pcs_sub[:, 0], sc)
    r_pc3, p_pc3 = stats.pearsonr(pcs_sub[:, 2], sc)
    X_ols = np.column_stack([np.ones(len(sc)), pcs_sub[:, :3]])
    beta_ols = np.linalg.lstsq(X_ols, sc, rcond=None)[0]
    y_hat = X_ols @ beta_ols
    r2 = 1.0 - np.sum((sc - y_hat)**2) / np.sum((sc - np.mean(sc))**2)
    print(f"{name}: PC2 r={r_pc2:+.4f} (p={p_pc2:.2e}), PC1 r={r_pc1:+.4f}, PC3 r={r_pc3:+.4f}, PC1-3 R={np.sqrt(r2):.4f} (R^2={r2:.4f})")

# 4. Partial correlation of CRABP2 (cg02363950) methylation and mRNA expression adjusting for immune score and PC1-5
print("\n=== 4. CRABP2 PARTIAL CORRELATION (ADJUSTED FOR IMMUNE SCORE AND PC1-5) ===")
for pid, gname in [("cg02363950", "CRABP2"), ("cg23152216", "KRAS"), ("cg06703343", "APOL3"), ("cg00754357", "ZCCHC11")]:
    idx = probes_df[probes_df["probe"] == pid].index[0]
    m_vals = X_imp[idx, has_rna]
    e_vals = rna_df.loc[sub_samples, gname].values
    r_raw, p_raw = stats.pearsonr(m_vals, e_vals)
    
    # Partial correlation adjusting for score_14
    def partial_r(x, y, covs):
        C = np.column_stack([np.ones(len(x)), covs])
        rx = x - C @ np.linalg.lstsq(C, x, rcond=None)[0]
        ry = y - C @ np.linalg.lstsq(C, y, rcond=None)[0]
        return stats.pearsonr(rx, ry)
    
    r_imm, p_imm = partial_r(m_vals, e_vals, score_14)
    r_pcs, p_pcs = partial_r(m_vals, e_vals, pcs_sub[:, :5])
    print(f"{pid} ({gname}): unadjusted r={r_raw:+.4f} (p={p_raw:.2e}) | adj for RNA immune score: r={r_imm:+.4f} (p={p_imm:.2e}) | adj for PC1-5: r={r_pcs:+.4f} (p={p_pcs:.2e})")

# 5. Check Table 4 Welch t-statistic, degrees of freedom, and p-value in GSE120878
print("\n=== 5. TABLE 4 WELCH T-TEST DEGREES OF FREEDOM CHECK ===")
X_120 = np.load("gse120878_prom_X.npy")
y_120 = np.load("gse120878_y.npy")
p_120 = pd.read_csv("gse120878_prom_probes.csv")["probe"].values
top10_df = pd.read_csv("top50_diagnostic_gse120878.tsv", sep="\t").head(10)
for _, row in top10_df.iterrows():
    pid = row["probe"]
    idx = np.where(p_120 == pid)[0][0]
    x1 = X_120[y_120 == 1, idx] # melanoma (n=89)
    x0 = X_120[y_120 == 0, idx] # nevus (n=73)
    s1, s0 = np.var(x1, ddof=1), np.var(x0, ddof=1)
    n1, n0 = len(x1), len(x0)
    se = np.sqrt(s1/n1 + s0/n0)
    t_stat = (np.mean(x1) - np.mean(x0)) / se
    # Welch-Satterthwaite df
    df_w = (s1/n1 + s0/n0)**2 / ((s1/n1)**2/(n1-1) + (s0/n0)**2/(n0-1))
    p_w = stats.t.sf(np.abs(t_stat), df=df_w) * 2
    p_stud = stats.ttest_ind(x1, x0, equal_var=True).pvalue
    print(f"{pid} ({row['gene']}): t={t_stat:+.2f}, Welch df={df_w:.1f}, Welch p={p_w:.2e}, Student p(df=160)={p_stud:.2e}, SD_mel={np.sqrt(s1):.3f}, SD_nev={np.sqrt(s0):.3f}")
