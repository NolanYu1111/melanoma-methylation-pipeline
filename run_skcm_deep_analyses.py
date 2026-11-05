import os
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedShuffleSplit
import statsmodels.api as sm
from statsmodels.duration.hazard_regression import PHReg
from run_cox_and_ml import fit_breslow_null_and_score_test

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

print("Loading TCGA-SKCM promoter methylation and RNA-seq data...", flush=True)
X_raw = np.load("tcga_promoter_X_raw.npy")
probes_df = pd.read_csv("tcga_promoter_probes.csv")
meta = pd.read_csv("tcga_456_meta.csv")

gpl = pd.read_csv(
    os.path.join(dl, "GPL13534-11288.txt"),
    sep="\t", skiprows=37,
    usecols=["ID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "CHR", "Probe_SNPs", "Probe_SNPs_10"],
    low_memory=False
).dropna(subset=["ID"])

probe_to_allgenes = dict(zip(gpl["ID"], gpl["UCSC_RefGene_Name"].fillna("")))
probe_to_gene = {
    k: (v.split(";")[0] if v else "Unannotated")
    for k, v in probe_to_allgenes.items()
}
probes_df["Gene"] = probes_df["probe"].map(probe_to_gene)
probes_df["AllGenes"] = probes_df["probe"].map(probe_to_allgenes)

# Mean-impute and standardize methylation
row_means = np.nanmean(X_raw, axis=1, keepdims=True)
X_imp = np.where(np.isnan(X_raw), row_means, X_raw)
row_sds = np.std(X_imp, axis=1, ddof=1, keepdims=True)
X_std = (X_imp - row_means) / row_sds

pca = PCA(n_components=5, random_state=42)
pcs = pca.fit_transform(X_std.T)

# Load RNA-seq
rnaseq_df = pd.read_csv(os.path.join(dl, "TCGA-SKCM.HiSeqV2.tsv.gz"), sep="\t", index_col=0)
meta["sample_15"] = meta["sample"].str.slice(0, 15)
has_rna_mask = meta["sample_15"].isin(rnaseq_df.columns).values
rna_samples = meta.loc[has_rna_mask, "sample_15"].values
rna_df = rnaseq_df[rna_samples] # (20530, 454)

# Restrict to the 454 patients with matched RNA-seq + methylation
meta_454 = meta.loc[has_rna_mask].copy().reset_index(drop=True)
X_std_454 = X_std[:, has_rna_mask]
pcs_454 = pcs[has_rna_mask]

# 14-gene immune score
immune_markers = ["PTPRC", "CD8A", "CD8B", "CD3D", "CD3E", "GZMA", "PRF1", "IFNG", "CXCL9", "CXCL10", "STAT1", "PSMB9", "TAP1", "HLA-DRA"]
immune_expr = rna_df.loc[immune_markers].values
immune_z = (immune_expr - immune_expr.mean(axis=1, keepdims=True)) / immune_expr.std(axis=1, ddof=1, keepdims=True)
imm_raw = immune_z.mean(axis=0)
imm_std = (imm_raw - imm_raw.mean()) / imm_raw.std(ddof=1)

# Clinical covariates for 454 cohort
sex_male = (meta_454["gender.demographic"].str.lower() == "male").astype(float).values
age = meta_454["age_at_index.demographic"].values.astype(float)
age_std = (age - age.mean()) / age.std(ddof=1)
is_prim = meta_454["is_primary"].values.astype(float)
stage_dummies = pd.get_dummies(meta_454["stage_clean"], drop_first=True, dtype=float).values
Z_clin_454 = np.column_stack([age_std, sex_male, is_prim, stage_dummies])

time_454 = meta_454["OS.time"].values.astype(float)
event_454 = meta_454["OS"].values.astype(int)

# ==============================================================================
# ITEM 4: CRABP2 (cg02363950) Deep Follow-Up (Expression Survival & Mediation)
# ==============================================================================
print("\n" + "="*80)
print("ITEM 4: CRABP2 (cg02363950) EXPRESSION SURVIVAL & MEDIATION ANALYSIS (n = 454)")
print("="*80, flush=True)

crabp2_idx = probes_df[probes_df["probe"] == "cg02363950"].index[0]
crabp2_meth = X_std_454[crabp2_idx] # standardized across 456; re-standardize across 454
crabp2_meth = (crabp2_meth - crabp2_meth.mean()) / crabp2_meth.std(ddof=1)

crabp2_expr_raw = rna_df.loc["CRABP2"].values.astype(float)
crabp2_expr = (crabp2_expr_raw - crabp2_expr_raw.mean()) / crabp2_expr_raw.std(ddof=1)

# Path a: linear regression of CRABP2 expression on cg02363950 methylation
# (1) Unadjusted, (2) Clinical-adjusted, (3) Clinical + Immune-adjusted, (4) Clinical + 5-PC-adjusted
for label, Z_cov in [
    ("Unadjusted", np.ones((454, 1))),
    ("Clinical-Adjusted", sm.add_constant(Z_clin_454)),
    ("Clinical + Immune-Adjusted", sm.add_constant(np.column_stack([Z_clin_454, imm_std]))),
    ("Clinical + 5-PC-Adjusted", sm.add_constant(np.column_stack([Z_clin_454, pcs_454]))),
]:
    ols_fit = sm.OLS(crabp2_expr, np.column_stack([crabp2_meth, Z_cov])).fit()
    a_coef = ols_fit.params[0]
    a_se = ols_fit.bse[0]
    a_p = ols_fit.pvalues[0]
    print(f"Path a (cg02363950 -> CRABP2 mRNA) [{label}]: beta_a = {a_coef:+.4f} (SE = {a_se:.4f}, p = {a_p:.4e})", flush=True)

# Cox models for (A) cg02363950 alone, (B) CRABP2 mRNA alone, (C) Joint cg02363950 + CRABP2 mRNA
cox_mediation_rows = []
for cov_name, Z_cov in [
    ("Clinical Only", Z_clin_454),
    ("Clinical + RNA Immune", np.column_stack([Z_clin_454, imm_std])),
    ("Clinical + 5 PCs", np.column_stack([Z_clin_454, pcs_454])),
]:
    # Fit with statsmodels PHReg (Breslow ties)
    res_m = PHReg(time_454, np.column_stack([crabp2_meth, Z_cov]), status=event_454, ties="breslow").fit()
    res_e = PHReg(time_454, np.column_stack([crabp2_expr, Z_cov]), status=event_454, ties="breslow").fit()
    res_j = PHReg(time_454, np.column_stack([crabp2_meth, crabp2_expr, Z_cov]), status=event_454, ties="breslow").fit()
    
    # Path a with same covariates
    ols_a = sm.OLS(crabp2_expr, np.column_stack([crabp2_meth, sm.add_constant(Z_cov)])).fit()
    a_val, a_se = ols_a.params[0], ols_a.bse[0]
    b_val, b_se = res_j.params[1], res_j.bse[1] # expression effect controlling for methylation
    c_val, c_se = res_m.params[0], res_m.bse[0] # total methylation effect
    cp_val, cp_se = res_j.params[0], res_j.bse[0] # direct methylation effect controlling for expression
    
    # Sobel test for indirect effect a * b
    indirect = a_val * b_val
    sobel_se = np.sqrt((b_val**2) * (a_se**2) + (a_val**2) * (b_se**2))
    sobel_z = indirect / sobel_se
    sobel_p = 2 * stats.norm.sf(abs(sobel_z))
    prop_mediated = (c_val - cp_val) / c_val if c_val != 0 else np.nan
    
    print(f"\n--- {cov_name} ---")
    print(f"  [Methylation alone] cg02363950: logHR = {c_val:+.4f}, HR = {np.exp(c_val):.4f} (95% CI: {np.exp(c_val-1.96*c_se):.4f}-{np.exp(c_val+1.96*c_se):.4f}), p = {res_m.pvalues[0]:.4e}")
    print(f"  [Expression alone]  CRABP2 mRNA: logHR = {res_e.params[0]:+.4f}, HR = {np.exp(res_e.params[0]):.4f} (95% CI: {np.exp(res_e.params[0]-1.96*res_e.bse[0]):.4f}-{np.exp(res_e.params[0]+1.96*res_e.bse[0]):.4f}), p = {res_e.pvalues[0]:.4e}")
    print(f"  [Joint Model]       cg02363950:  logHR = {cp_val:+.4f}, HR = {np.exp(cp_val):.4f}, p = {res_j.pvalues[0]:.4e}")
    print(f"                      CRABP2 mRNA: logHR = {b_val:+.4f}, HR = {np.exp(b_val):.4f}, p = {res_j.pvalues[1]:.4e}")
    print(f"  [Mediation]         Indirect (a*b) = {indirect:+.4f}, Sobel z = {sobel_z:+.3f}, p = {sobel_p:.4f}, Prop. change (c - c')/c = {prop_mediated*100:+.1f}%", flush=True)
    
    cox_mediation_rows.append({
        "Adjustment": cov_name,
        "Meth_Alone_HR": np.exp(c_val), "Meth_Alone_CI_L": np.exp(c_val-1.96*c_se), "Meth_Alone_CI_U": np.exp(c_val+1.96*c_se), "Meth_Alone_p": res_m.pvalues[0],
        "Expr_Alone_HR": np.exp(res_e.params[0]), "Expr_Alone_CI_L": np.exp(res_e.params[0]-1.96*res_e.bse[0]), "Expr_Alone_CI_U": np.exp(res_e.params[0]+1.96*res_e.bse[0]), "Expr_Alone_p": res_e.pvalues[0],
        "Joint_Meth_HR": np.exp(cp_val), "Joint_Meth_p": res_j.pvalues[0],
        "Joint_Expr_HR": np.exp(b_val), "Joint_Expr_p": res_j.pvalues[1],
        "Path_a_beta": a_val, "Path_a_p": ols_a.pvalues[0],
        "Indirect_ab": indirect, "Sobel_z": sobel_z, "Sobel_p": sobel_p,
        "Pct_Attenuated": prop_mediated * 100.0
    })

pd.DataFrame(cox_mediation_rows).to_csv("crabp2_mediation_results.tsv", sep="\t", index=False)

# ==============================================================================
# ITEM 3: Do survival-linked CpGs behave like real regulators?
# Genome-wide Promoter Methylation-Expression Correlation & Enrichment Test
# ==============================================================================
print("\n" + "="*80)
print("ITEM 3: GENOME-WIDE PROMOTER METHYLATION-EXPRESSION ENRICHMENT (n = 454)")
print("="*80, flush=True)

# Map each of the 67,102 promoter probes to its primary gene symbol in rna_df (if present)
rna_gene_set = set(rna_df.index)
# Also handle common aliases like ZCCHC11 -> TUT4
alias_map = {"ZCCHC11": "TUT4"}

matched_probe_indices = []
matched_rna_genes = []
for idx, all_g in enumerate(probes_df["AllGenes"].values):
    if not isinstance(all_g, str) or not all_g:
        continue
    found_g = None
    for g in all_g.split(";"):
        if g in rna_gene_set:
            found_g = g
            break
        elif g in alias_map and alias_map[g] in rna_gene_set:
            found_g = alias_map[g]
            break
    if found_g is not None:
        matched_probe_indices.append(idx)
        matched_rna_genes.append(found_g)

matched_probe_indices = np.array(matched_probe_indices, dtype=int)
print(f"Total tested promoter CpGs mapped to expressed genes in HiSeqV2: {len(matched_probe_indices)} / {len(probes_df)}", flush=True)

# Vectorized Pearson r and immune-adjusted partial r across all 56,000+ matched CpG-gene pairs!
M_mat = X_std_454[matched_probe_indices] # (N_matched, 454)
E_mat = rna_df.loc[matched_rna_genes].values.astype(float) # (N_matched, 454)

# Filter out any gene with 0 variance across the 454 samples
e_sd = E_mat.std(axis=1, ddof=1)
valid_e = e_sd > 1e-8
matched_probe_indices = matched_probe_indices[valid_e]
matched_rna_genes = np.array(matched_rna_genes)[valid_e]
M_mat = M_mat[valid_e]
E_mat = E_mat[valid_e]
print(f"Promoter CpGs with non-zero mRNA variance: {len(matched_probe_indices)}", flush=True)

# Center and normalize rows for fast Pearson correlation
M_c = M_mat - M_mat.mean(axis=1, keepdims=True)
M_n = M_c / np.linalg.norm(M_c, axis=1, keepdims=True)
E_c = E_mat - E_mat.mean(axis=1, keepdims=True)
E_n = E_c / np.linalg.norm(E_c, axis=1, keepdims=True)

r_unadj = np.sum(M_n * E_n, axis=1)
t_unadj = r_unadj * np.sqrt((454 - 2) / np.clip(1.0 - r_unadj**2, 1e-15, None))
p_unadj = 2 * stats.t.sf(np.abs(t_unadj), df=454 - 2)
q_unadj = bh_fdr(p_unadj)

# Immune-adjusted partial correlation: residualize both M_c and E_c on imm_std
imm_u = (imm_std - imm_std.mean())
imm_u = imm_u / np.linalg.norm(imm_u)
M_res = M_c - np.outer(M_c @ imm_u, imm_u)
E_res = E_c - np.outer(E_c @ imm_u, imm_u)
M_res_n = M_res / np.linalg.norm(M_res, axis=1, keepdims=True)
E_res_n = E_res / np.linalg.norm(E_res, axis=1, keepdims=True)

r_imm_adj = np.sum(M_res_n * E_res_n, axis=1)
t_imm_adj = r_imm_adj * np.sqrt((454 - 3) / np.clip(1.0 - r_imm_adj**2, 1e-15, None))
p_imm_adj = 2 * stats.t.sf(np.abs(t_imm_adj), df=454 - 3)
q_imm_adj = bh_fdr(p_imm_adj)

# Also 5-PC-adjusted partial correlation: residualize both M_c and E_c on pcs_454
Q_pc, _ = np.linalg.qr(pcs_454 - pcs_454.mean(axis=0, keepdims=True)) # (454, 5) orthonormal basis
M_res_pc = M_c - (M_c @ Q_pc) @ Q_pc.T
E_res_pc = E_c - (E_c @ Q_pc) @ Q_pc.T
r_pc5_adj = np.sum(
    (M_res_pc / np.linalg.norm(M_res_pc, axis=1, keepdims=True)) *
    (E_res_pc / np.linalg.norm(E_res_pc, axis=1, keepdims=True)),
    axis=1
)
t_pc5_adj = r_pc5_adj * np.sqrt((454 - 7) / np.clip(1.0 - r_pc5_adj**2, 1e-15, None))
p_pc5_adj = 2 * stats.t.sf(np.abs(t_pc5_adj), df=454 - 7)
q_pc5_adj = bh_fdr(p_pc5_adj)

# Load survival q-values for all 67,102 probes
cox_clin = pd.read_csv(os.path.join(dl, "cox_results.tsv"), sep="\t").set_index("probe").loc[probes_df["probe"]].reset_index()
cox_pc5 = pd.read_csv(os.path.join(dl, "cox_pc5.tsv"), sep="\t").set_index("probe").loc[probes_df["probe"]].reset_index()

# Recompute RNA-seq immune-adjusted Cox screen on the 454 matched samples
order_454 = np.argsort(time_454)
_, p_cox_imm, hr_cox_imm, _ = fit_breslow_null_and_score_test(
    np.column_stack([Z_clin_454[order_454], imm_std[order_454]]),
    X_std_454[:, order_454],
    time_454[order_454],
    event_454[order_454]
)
q_cox_imm = bh_fdr(p_cox_imm)

# Build DataFrame of all matched promoter probes
reg_df = pd.DataFrame({
    "probe_idx": matched_probe_indices,
    "probe": probes_df["probe"].values[matched_probe_indices],
    "Gene": matched_rna_genes,
    "r_unadj": r_unadj, "p_unadj": p_unadj, "q_unadj": q_unadj,
    "r_imm_adj": r_imm_adj, "p_imm_adj": p_imm_adj, "q_imm_adj": q_imm_adj,
    "r_pc5_adj": r_pc5_adj, "p_pc5_adj": p_pc5_adj, "q_pc5_adj": q_pc5_adj,
    "cox_clin_p": cox_clin["p_clin"].values[matched_probe_indices],
    "cox_clin_q": cox_clin["q_clin"].values[matched_probe_indices],
    "cox_pc5_p": cox_pc5["p_pc5"].values[matched_probe_indices],
    "cox_pc5_q": cox_pc5["q_pc5"].values[matched_probe_indices],
    "cox_imm_p": p_cox_imm[matched_probe_indices],
    "cox_imm_q": q_cox_imm[matched_probe_indices],
})
reg_df["beta_sd"] = np.nanstd(X_raw[matched_probe_indices][:, has_rna_mask], axis=1, ddof=1)
reg_df["beta_mean"] = np.nanmean(X_raw[matched_probe_indices][:, has_rna_mask], axis=1)

# Define canonical inverse correlation thresholds:
# (A) Any statistically significant inverse correlation after immune adjustment: r_imm_adj < 0 & q_imm_adj < 0.05
# (B) Moderate-to-strong inverse correlation after immune adjustment: r_imm_adj < -0.20 (and q_imm_adj < 0.05)
reg_df["is_inv_sig"] = (reg_df["r_imm_adj"] < 0) & (reg_df["q_imm_adj"] < 0.05)
reg_df["is_inv_mod"] = (reg_df["r_imm_adj"] < -0.20) & (reg_df["q_imm_adj"] < 0.05)
reg_df["is_pos_sig"] = (reg_df["r_imm_adj"] > 0) & (reg_df["q_imm_adj"] < 0.05)
reg_df.to_csv("promoter_meth_expr_survival_table.tsv", sep="\t", index=False)

# Matched null: bin all CpGs by SD decile x mean-beta decile; for each hit set draw 1,000 sets with identical bin composition
reg_df["sd_bin"] = pd.qcut(reg_df["beta_sd"], 10, labels=False, duplicates="drop")
reg_df["mean_bin"] = pd.qcut(reg_df["beta_mean"], 10, labels=False, duplicates="drop")
reg_df["bin"] = reg_df["sd_bin"].astype(int) * 10 + reg_df["mean_bin"].astype(int)
bin_members = {b: np.where(reg_df["bin"].values == b)[0] for b in reg_df["bin"].unique()}
inv_mod_arr = reg_df["is_inv_mod"].values
inv_sig_arr = reg_df["is_inv_sig"].values
rng_null = np.random.default_rng(2026)
print("\n--- SD- and mean-beta-matched null (1,000 draws per hit set) ---", flush=True)
matched_rows = []
for name, mask in [
    ("Clinical-Only Survival Hits (q < 0.05)", reg_df["cox_clin_q"].values < 0.05),
    ("5-PC-Adjusted Survival Hits (q < 0.05)", reg_df["cox_pc5_q"].values < 0.05),
    ("5-PC-Adjusted Nominal Hits (p < 0.001)", reg_df["cox_pc5_p"].values < 0.001),
    ("RNA-Immune-Adjusted Nominal Hits (p < 0.001)", reg_df["cox_imm_p"].values < 0.001),
]:
    hit_idx = np.where(mask)[0]
    bin_counts = pd.Series(reg_df["bin"].values[hit_idx]).value_counts()
    obs_mod = inv_mod_arr[hit_idx].mean()
    obs_sig = inv_sig_arr[hit_idx].mean()
    null_mod, null_sig = np.empty(1000), np.empty(1000)
    for it in range(1000):
        draw = np.concatenate([rng_null.choice(bin_members[b], size=c, replace=True) for b, c in bin_counts.items()])
        null_mod[it] = inv_mod_arr[draw].mean()
        null_sig[it] = inv_sig_arr[draw].mean()
    p_mod = (1 + np.sum(null_mod >= obs_mod)) / 1001
    p_sig = (1 + np.sum(null_sig >= obs_sig)) / 1001
    print(f"{name} (n = {len(hit_idx)}): inverse r<-0.20 observed {100*obs_mod:.1f}% vs matched null {100*null_mod.mean():.1f}% (95% {100*np.percentile(null_mod,2.5):.1f}-{100*np.percentile(null_mod,97.5):.1f}%), empirical p = {p_mod:.4f}; any inverse q<0.05 observed {100*obs_sig:.1f}% vs null {100*null_sig.mean():.1f}%, p = {p_sig:.4f}", flush=True)
    matched_rows.append({"Subset": name, "N": len(hit_idx),
                         "Obs_Pct_Inv_r20": 100*obs_mod, "Null_Mean_Pct_Inv_r20": 100*null_mod.mean(),
                         "Null_2.5": 100*np.percentile(null_mod, 2.5), "Null_97.5": 100*np.percentile(null_mod, 97.5), "Emp_p_Inv_r20": p_mod,
                         "Obs_Pct_Inv_q05": 100*obs_sig, "Null_Mean_Pct_Inv_q05": 100*null_sig.mean(), "Emp_p_Inv_q05": p_sig})
pd.DataFrame(matched_rows).to_csv("regulator_enrichment_matched_null.tsv", sep="\t", index=False)

# Also draw a fixed random sample of 10,000 promoter probes for comparison
rng = np.random.default_rng(42)
rand_10k_mask = np.zeros(len(reg_df), dtype=bool)
rand_10k_mask[rng.choice(len(reg_df), size=10000, replace=False)] = True

subsets = [
    ("All Mapped Promoter CpGs (Background)", np.ones(len(reg_df), dtype=bool)),
    ("Random 10,000 Promoter CpGs", rand_10k_mask),
    ("Clinical-Only Survival Hits (q < 0.05)", reg_df["cox_clin_q"].values < 0.05),
    ("5-PC-Adjusted Survival Hits (q < 0.05)", reg_df["cox_pc5_q"].values < 0.05),
    ("5-PC-Adjusted Nominal Hits (p < 0.001)", reg_df["cox_pc5_p"].values < 0.001),
    ("RNA-Immune-Adjusted Hits (q < 0.05)", reg_df["cox_imm_q"].values < 0.05),
    ("RNA-Immune-Adjusted Nominal Hits (p < 0.001)", reg_df["cox_imm_p"].values < 0.001),
]

bg_inv_sig = reg_df["is_inv_sig"].values
bg_inv_mod = reg_df["is_inv_mod"].values

enrichment_rows = []
for name, mask in subsets:
    sub = reg_df[mask]
    n_sub = len(sub)
    n_inv_sig = sub["is_inv_sig"].sum()
    n_inv_mod = sub["is_inv_mod"].sum()
    n_pos_sig = sub["is_pos_sig"].sum()
    mean_r_unadj = sub["r_unadj"].mean()
    mean_r_imm = sub["r_imm_adj"].mean()
    median_r_imm = sub["r_imm_adj"].median()
    
    # Fisher's exact test vs complement in background
    if mask.sum() < len(reg_df):
        comp = reg_df[~mask]
        tab_sig = [[n_inv_sig, n_sub - n_inv_sig], [comp["is_inv_sig"].sum(), len(comp) - comp["is_inv_sig"].sum()]]
        or_sig, p_fish_sig = stats.fisher_exact(tab_sig)
        tab_mod = [[n_inv_mod, n_sub - n_inv_mod], [comp["is_inv_mod"].sum(), len(comp) - comp["is_inv_mod"].sum()]]
        or_mod, p_fish_mod = stats.fisher_exact(tab_mod)
        _, p_mwu = stats.mannwhitneyu(sub["r_imm_adj"], comp["r_imm_adj"], alternative="two-sided")
    else:
        or_sig, p_fish_sig, or_mod, p_fish_mod, p_mwu = 1.0, 1.0, 1.0, 1.0, 1.0
        
    print(f"{name} (n = {n_sub}):")
    print(f"  Mean r_unadj = {mean_r_unadj:+.4f} | Mean r_imm_adj = {mean_r_imm:+.4f} (Median {median_r_imm:+.4f}, MWU p = {p_mwu:.3e})")
    print(f"  Inverse (r_imm < 0, q < 0.05):    {n_inv_sig}/{n_sub} ({100*n_inv_sig/n_sub:.1f}%) | OR = {or_sig:.2f}, Fisher p = {p_fish_sig:.3e}")
    print(f"  Inverse (r_imm < -0.20, q < 0.05): {n_inv_mod}/{n_sub} ({100*n_inv_mod/n_sub:.1f}%) | OR = {or_mod:.2f}, Fisher p = {p_fish_mod:.3e}")
    print(f"  Positive (r_imm > 0, q < 0.05):   {n_pos_sig}/{n_sub} ({100*n_pos_sig/n_sub:.1f}%)", flush=True)
    
    enrichment_rows.append({
        "Subset": name, "N_Mapped_Probes": n_sub,
        "Mean_r_Unadj": mean_r_unadj, "Mean_r_ImmAdj": mean_r_imm, "Median_r_ImmAdj": median_r_imm,
        "N_Inv_Sig_q05": n_inv_sig, "Pct_Inv_Sig_q05": 100.0 * n_inv_sig / n_sub, "OR_Inv_Sig": or_sig, "Fisher_p_Inv_Sig": p_fish_sig,
        "N_Inv_Mod_r20": n_inv_mod, "Pct_Inv_Mod_r20": 100.0 * n_inv_mod / n_sub, "OR_Inv_Mod": or_mod, "Fisher_p_Inv_Mod": p_fish_mod,
        "N_Pos_Sig_q05": n_pos_sig, "Pct_Pos_Sig_q05": 100.0 * n_pos_sig / n_sub,
        "MWU_p_vs_Background": p_mwu
    })

pd.DataFrame(enrichment_rows).to_csv("regulator_enrichment_results.tsv", sep="\t", index=False)

# List which of the 78 5-PC-adjusted hits (and 74 mapped to RNA-seq) have r_imm_adj < -0.20 & q_imm_adj < 0.05
pc5_mapped = reg_df[reg_df["cox_pc5_q"] < 0.05].sort_values("r_imm_adj")
print("\n5-PC-Adjusted Survival Hits (q < 0.05) with Strongest Inverse Immune-Adjusted mRNA Correlation (r_imm_adj < -0.20):")
print(pc5_mapped[pc5_mapped["is_inv_mod"]][["probe", "Gene", "cox_pc5_q", "r_unadj", "r_imm_adj", "p_imm_adj", "q_imm_adj", "r_pc5_adj"]].to_string(index=False), flush=True)

# ==============================================================================
# ITEM 2: Is methylation prognostic only in some tumors?
# Immune-Hot vs. Immune-Cold & Methylation x Immune Interaction with 50/50 Split
# ==============================================================================
print("\n" + "="*80)
print("ITEM 2: IMMUNE-HOT VS. IMMUNE-COLD & INTERACTION SCREEN (50/50 DISCOVERY/VALIDATION)")
print("="*80, flush=True)

# Pre-specified 50/50 stratified split by death indicator (OS)
sss = StratifiedShuffleSplit(n_splits=1, test_size=0.5, random_state=42)
disc_idx, val_idx = next(sss.split(np.zeros(454), event_454))
print(f"Discovery half: n = {len(disc_idx)} ({event_454[disc_idx].sum()} deaths)")
print(f"Validation half: n = {len(val_idx)} ({event_454[val_idx].sum()} deaths)", flush=True)

def run_subgroup_and_interaction(sub_indices, cohort_label):
    t_sub = time_454[sub_indices]
    d_sub = event_454[sub_indices]
    Z_sub = Z_clin_454[sub_indices]
    imm_sub = imm_std[sub_indices]
    X_sub = X_std_454[:, sub_indices]
    
    # 1. Immune-Cold (<= median) vs. Immune-Hot (> median) subgroup screens
    med_imm = np.median(imm_sub)
    cold_mask = imm_sub <= med_imm
    hot_mask = imm_sub > med_imm
    
    # Sort each subgroup by time
    ord_c = np.argsort(t_sub[cold_mask])
    z_cold, p_cold, hr_cold, _ = fit_breslow_null_and_score_test(
        np.column_stack([Z_sub[cold_mask][ord_c], imm_sub[cold_mask][ord_c]]),
        X_sub[:, cold_mask][:, ord_c],
        t_sub[cold_mask][ord_c],
        d_sub[cold_mask][ord_c]
    )
    q_cold = bh_fdr(p_cold)
    
    ord_h = np.argsort(t_sub[hot_mask])
    z_hot, p_hot, hr_hot, _ = fit_breslow_null_and_score_test(
        np.column_stack([Z_sub[hot_mask][ord_h], imm_sub[hot_mask][ord_h]]),
        X_sub[:, hot_mask][:, ord_h],
        t_sub[hot_mask][ord_h],
        d_sub[hot_mask][ord_h]
    )
    q_hot = bh_fdr(p_hot)
    
    # 2. Interaction term test: test difference in z-scores / slopes between Hot and Cold, AND score test for X * imm
    # First: fast score test for product term (X_sub * imm_sub) where each row X_j * imm is residualized on X_j!
    ord_all = np.argsort(t_sub)
    t_s = t_sub[ord_all]
    d_s = d_sub[ord_all]
    Z_s = np.column_stack([Z_sub[ord_all], imm_sub[ord_all]])
    X_s = X_sub[:, ord_all]
    imm_s = imm_sub[ord_all]
    
    # Construct interaction matrix I_s = X_s * imm_s, orthogonalized row-wise against X_s so main effect is removed
    I_raw = X_s * imm_s[None, :]
    x_norm2 = np.sum(X_s * X_s, axis=1, keepdims=True)
    proj = np.sum(I_raw * X_s, axis=1, keepdims=True) / np.clip(x_norm2, 1e-12, None)
    I_orth = I_raw - proj * X_s
    
    z_int, p_int, hr_int, _ = fit_breslow_null_and_score_test(Z_s, I_orth, t_s, d_s)
    q_int = bh_fdr(p_int)
    
    # For the top 20 interaction candidates by p_int, verify with exact full PHReg model (including main effect X_j and X_j * imm)
    top20_idx = np.argsort(p_int)[:20]
    exact_p_list = []
    for pidx in top20_idx:
        xj = X_s[pidx]
        ij = xj * imm_s
        try:
            res_full = PHReg(t_s, np.column_stack([xj, ij, Z_s]), status=d_s, ties="breslow").fit()
            exact_p_list.append(res_full.pvalues[1])
        except Exception:
            exact_p_list.append(p_int[pidx])
    min_exact_p = np.min(exact_p_list)
    
    print(f"\n[{cohort_label}] (n = {len(sub_indices)}; Cold n = {cold_mask.sum()}, Hot n = {hot_mask.sum()}):")
    print(f"  Immune-Cold subgroup (q < 0.05): {(q_cold < 0.05).sum()} hits (p < 0.001: {(p_cold < 0.001).sum()}, min p = {p_cold.min():.3e})")
    print(f"  Immune-Hot subgroup  (q < 0.05): {(q_hot < 0.05).sum()} hits (p < 0.001: {(p_hot < 0.001).sum()}, min p = {p_hot.min():.3e})")
    print(f"  CpG x Immune Interaction (q < 0.05): {(q_int < 0.05).sum()} hits (p < 0.001: {(p_int < 0.001).sum()}, min score p = {p_int.min():.3e}, min exact Wald p = {min_exact_p:.3e})", flush=True)
    
    return {
        "z_cold": z_cold, "p_cold": p_cold, "q_cold": q_cold, "hr_cold": hr_cold,
        "z_hot": z_hot, "p_hot": p_hot, "q_hot": q_hot, "hr_hot": hr_hot,
        "z_int": z_int, "p_int": p_int, "q_int": q_int, "min_exact_p": min_exact_p
    }

res_disc = run_subgroup_and_interaction(disc_idx, "Discovery Half")
res_val = run_subgroup_and_interaction(val_idx, "Validation Half")
res_full = run_subgroup_and_interaction(np.arange(454), "Full Cohort (n = 454)")

# Check Discovery -> Validation replication for nominal p < 0.001 candidates in Cold, Hot, and Interaction!
for mode, k_p, k_z in [
    ("Immune-Cold Subgroup", "p_cold", "z_cold"),
    ("Immune-Hot Subgroup", "p_hot", "z_hot"),
    ("CpG x Immune Interaction", "p_int", "z_int"),
]:
    disc_cand = np.where(res_disc[k_p] < 0.001)[0]
    if len(disc_cand) > 0:
        # Check how many replicate in Validation at p < 0.05 (or Bonferroni p < 0.05 / len(disc_cand)) with concordant sign
        val_p = res_val[k_p][disc_cand]
        concordant = np.sign(res_disc[k_z][disc_cand]) == np.sign(res_val[k_z][disc_cand])
        n_rep_nom = ((val_p < 0.05) & concordant).sum()
        n_rep_bonf = ((val_p < (0.05 / len(disc_cand))) & concordant).sum()
        n_rep_fdr = ((bh_fdr(val_p) < 0.05) & concordant).sum()
        print(f"Replication [{mode}]: {len(disc_cand)} Discovery candidates at p < 1e-3 -> {n_rep_nom} nominal p < 0.05 concordant in Validation ({100*n_rep_nom/len(disc_cand):.1f}%; ~{0.025*len(disc_cand):.1f} expected by chance), {n_rep_fdr} at Validation FDR < 0.05, {n_rep_bonf} at Bonferroni.", flush=True)

# Also check CRABP2 (cg02363950) specifically in Immune-Cold vs Immune-Hot and interaction in Full Cohort
cidx = crabp2_idx
print(f"\nCRABP2 (cg02363950) in Immune-Cold vs Immune-Hot (Full Cohort n = 454):")
print(f"  Immune-Cold (n = 227): HR = {res_full['hr_cold'][cidx]:.3f}, p = {res_full['p_cold'][cidx]:.4e}")
print(f"  Immune-Hot  (n = 227): HR = {res_full['hr_hot'][cidx]:.3f}, p = {res_full['p_hot'][cidx]:.4e}")
print(f"  Interaction p = {res_full['p_int'][cidx]:.4f}", flush=True)
