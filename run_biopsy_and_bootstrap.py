import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize
import statsmodels.api as sm
from statsmodels.duration.hazard_regression import PHReg
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

def fit_breslow_null_and_score_test(Z_cov, X_genes, t_s, d_s):
    N, p = Z_cov.shape
    utimes = np.unique(t_s[d_s == 1])
    n_events_at_t = np.array([np.sum((t_s == ut) & (d_s == 1)) for ut in utimes], dtype=float)
    Z_events_sum = np.array([np.sum(Z_cov[(t_s == ut) & (d_s == 1)], axis=0) for ut in utimes])
    risk_masks = np.array([t_s >= ut for ut in utimes])

    def neg_log_lik_and_grad(beta):
        eta = np.clip(Z_cov @ beta, -50, 50)
        w = np.exp(eta)
        S0 = risk_masks @ w
        S1 = (risk_masks * w[None, :]) @ Z_cov
        ll = np.sum(Z_events_sum * beta[None, :]) - np.sum(n_events_at_t * np.log(S0))
        grad = np.sum(Z_events_sum, axis=0) - np.sum(n_events_at_t[:, None] * (S1 / S0[:, None]), axis=0)
        return -ll, -grad

    res = minimize(neg_log_lik_and_grad, np.zeros(p), jac=True, method="BFGS")
    beta_hat = res.x
    eta = np.clip(Z_cov @ beta_hat, -50, 50)
    w = np.exp(eta)
    S0 = risk_masks @ w
    E_Z = ((risk_masks * w[None, :]) @ Z_cov) / S0[:, None]

    I_ZZ = np.zeros((p, p))
    for k in range(len(utimes)):
        mask_k = risk_masks[k]
        wk = w[mask_k] / S0[k]
        Zk_centered = Z_cov[mask_k] - E_Z[k]
        I_ZZ += n_events_at_t[k] * (Zk_centered.T * wk) @ Zk_centered
    I_ZZ_inv = np.linalg.pinv(I_ZZ)

    cum_hazard_factor = (risk_masks.T @ (n_events_at_t / S0))
    martingale_resid = d_s - w * cum_hazard_factor
    U_x = X_genes @ martingale_resid

    weighted_risk = risk_masks.T * w[:, None]
    E_X = (X_genes @ weighted_risk) / S0[None, :]
    X_sq = X_genes ** 2
    E_X2 = (X_sq @ weighted_risk) / S0[None, :]
    I_XX = np.sum(n_events_at_t[None, :] * (E_X2 - E_X ** 2), axis=1)

    h_w = w * cum_hazard_factor
    term1_XZ = X_genes @ (Z_cov * h_w[:, None])
    term2_XZ = (E_X * n_events_at_t[None, :]) @ E_Z
    I_XZ = term1_XZ - term2_XZ

    V_x = I_XX - np.sum((I_XZ @ I_ZZ_inv) * I_XZ, axis=1)
    z_stat = U_x / np.sqrt(np.maximum(V_x, 1e-12))
    p_vals = 2.0 * stats.norm.sf(np.abs(z_stat))
    beta_1step = np.clip(U_x / np.maximum(V_x, 1e-12), -10, 10)
    hr_est = np.exp(beta_1step)
    return z_stat, p_vals, hr_est, beta_hat

X_raw = np.load("tcga_promoter_X_raw.npy")
probes_df = pd.read_csv("tcga_promoter_probes.csv")
meta = pd.read_csv("tcga_456_meta.csv")
clin = pd.read_csv(os.path.join(dl, "TCGA-SKCM.clinical.tsv.gz"), sep="\t", low_memory=False)
meta = meta.merge(clin[["sample", "site_of_resection_or_biopsy.diagnoses"]], on="sample", how="left")

def classify_site(s):
    s_low = str(s).lower()
    if "lymph node" in s_low:
        return "Lymph_Node"
    elif "skin" in s_low:
        return "Skin"
    elif "connective" in s_low or "subcutaneous" in s_low or "soft tissue" in s_low:
        return "Soft_Tissue_Subcut"
    else:
        return "Visceral_Other"

meta["biopsy_group"] = meta["site_of_resection_or_biopsy.diagnoses"].apply(classify_site)

row_means = np.nanmean(X_raw, axis=1, keepdims=True)
X_imp = np.where(np.isnan(X_raw), row_means, X_raw)
row_sds = np.std(X_imp, axis=1, ddof=1, keepdims=True)
X_std = (X_imp - row_means) / row_sds
pcs = PCA(n_components=5, random_state=42).fit_transform(X_std.T)

rnaseq_df = pd.read_csv(os.path.join(dl, "TCGA-SKCM.HiSeqV2.tsv.gz"), sep="\t", index_col=0)
meta["sample_15"] = meta["sample"].str.slice(0, 15)
has_rna = meta["sample_15"].isin(rnaseq_df.columns).values
rna_df = rnaseq_df[meta.loc[has_rna, "sample_15"].values]
immune_markers = ["PTPRC", "CD8A", "CD8B", "CD3D", "CD3E", "GZMA", "PRF1", "IFNG", "CXCL9", "CXCL10", "STAT1", "PSMB9", "TAP1", "HLA-DRA"]
imm_expr = rna_df.loc[immune_markers].values
imm_z = (imm_expr - imm_expr.mean(axis=1, keepdims=True)) / imm_expr.std(axis=1, ddof=1, keepdims=True)
imm_raw = imm_z.mean(axis=0)
imm_std = (imm_raw - imm_raw.mean()) / imm_raw.std(ddof=1)

meta_454 = meta.loc[has_rna].copy().reset_index(drop=True)

sex = (meta["gender.demographic"].str.lower() == "male").astype(float).values
age = meta["age_at_index.demographic"].values.astype(float)
age_std = (age - age.mean()) / age.std(ddof=1)
is_prim = meta["is_primary"].values.astype(float)
stage_d = pd.get_dummies(meta["stage_clean"], drop_first=True, dtype=float).values
site_d = pd.get_dummies(meta["biopsy_group"], drop_first=True, dtype=float).values

Z_clin = np.column_stack([age_std, sex, is_prim, stage_d])
Z_clin_site = np.column_stack([age_std, sex, is_prim, stage_d, site_d])

time_v = meta["OS.time"].values.astype(float)
event_v = meta["OS"].values.astype(int)

# 1. Clinical + Biopsy Site screen (n = 456)
_, p_site, hr_site, _ = fit_breslow_null_and_score_test(Z_clin_site, X_std, time_v, event_v)
q_site = bh_fdr(p_site)
print(f"Clinical + Biopsy Site screen (n=456): q < 0.05 = {(q_site < 0.05).sum()}, p < 0.05 = {(p_site < 0.05).sum()}")

# 2. Clinical + Biopsy Site + 5 PCs (n = 456)
_, p_site_pc5, _, _ = fit_breslow_null_and_score_test(np.column_stack([Z_clin_site, pcs]), X_std, time_v, event_v)
q_site_pc5 = bh_fdr(p_site_pc5)
print(f"Clinical + Biopsy Site + 5 PCs screen (n=456): q < 0.05 = {(q_site_pc5 < 0.05).sum()}")

# 3. Clinical + Biopsy Site + RNA Immune (n = 454)
Z_clin_454 = Z_clin[has_rna]
Z_clin_site_454 = Z_clin_site[has_rna]
res_imm_site = PHReg(time_v[has_rna], np.column_stack([imm_std, Z_clin_site_454]), status=event_v[has_rna], ties="breslow").fit()
print(f"Immune score -> OS adjusted for Clinical + Biopsy Site (n=454): HR = {np.exp(res_imm_site.params[0]):.4f}, p = {res_imm_site.pvalues[0]:.4e}")

_, p_site_imm, _, _ = fit_breslow_null_and_score_test(np.column_stack([Z_clin_site_454, imm_std]), X_std[:, has_rna], time_v[has_rna], event_v[has_rna])
q_site_imm = bh_fdr(p_site_imm)
print(f"Clinical + Biopsy Site + RNA Immune screen (n=454): q < 0.05 = {(q_site_imm < 0.05).sum()}")

# 4. Within non-lymph-node tumors (n=246 RNA-matched) and within lymph-node tumors (n=208 RNA-matched)
non_ln = (meta_454["biopsy_group"] != "Lymph_Node").values
ln_only = (meta_454["biopsy_group"] == "Lymph_Node").values
r_nonln, p_nonln = stats.pearsonr(pcs[has_rna][non_ln, 1], imm_std[non_ln])
r_ln, p_ln = stats.pearsonr(pcs[has_rna][ln_only, 1], imm_std[ln_only])
print(f"PC2 vs RNA immune within non-lymph-node tumors (n={non_ln.sum()}): r = {r_nonln:+.4f}, p = {p_nonln:.4e}")
print(f"PC2 vs RNA immune within lymph-node tumors (n={ln_only.sum()}): r = {r_ln:+.4f}, p = {p_ln:.4e}")

# 5. Nonparametric bootstrap 95% CI (1,000 resamples) for CRABP2 mediation (a * b)
crabp2_idx = probes_df[probes_df["probe"] == "cg02363950"].index[0]
crabp2_meth = X_std[crabp2_idx, has_rna]
crabp2_meth = (crabp2_meth - crabp2_meth.mean()) / crabp2_meth.std(ddof=1)
crabp2_expr_raw = rna_df.loc["CRABP2"].values.astype(float)
crabp2_expr = (crabp2_expr_raw - crabp2_expr_raw.mean()) / crabp2_expr_raw.std(ddof=1)

time_454 = time_v[has_rna]
event_454 = event_v[has_rna]
pcs_454 = pcs[has_rna]

rng = np.random.default_rng(42)
n_boot = 1000
med_df = pd.read_csv("crabp2_mediation_results.tsv", sep="\t")
boot_ci_l = []
boot_ci_u = []

for cov_name, Z_cov in [
    ("Clinical Only", Z_clin_454),
    ("Clinical + RNA Immune", np.column_stack([Z_clin_454, imm_std])),
    ("Clinical + 5 PCs", np.column_stack([Z_clin_454, pcs_454])),
]:
    ab_boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, 454, size=454)
        m_b = crabp2_meth[idx]
        e_b = crabp2_expr[idx]
        Z_b = Z_cov[idx]
        t_b = time_454[idx]
        d_b = event_454[idx]
        try:
            X_ols = np.column_stack([m_b, np.ones(454), Z_b])
            a_b = np.linalg.lstsq(X_ols, e_b, rcond=None)[0][0]
            res_j = PHReg(t_b, np.column_stack([m_b, e_b, Z_b]), status=d_b, ties="breslow").fit()
            b_b = res_j.params[1]
            if np.isfinite(a_b) and np.isfinite(b_b):
                ab_boots.append(a_b * b_b)
        except Exception:
            continue
    ci_l, ci_u = np.percentile(ab_boots, [2.5, 97.5])
    boot_ci_l.append(ci_l)
    boot_ci_u.append(ci_u)
    print(f"Bootstrap 95% CI for CRABP2 indirect effect a*b [{cov_name}]: [{ci_l:+.4f}, {ci_u:+.4f}] (n_valid={len(ab_boots)})")

med_df["Boot_Indirect_CI_L"] = boot_ci_l
med_df["Boot_Indirect_CI_U"] = boot_ci_u
med_df.to_csv("crabp2_mediation_results.tsv", sep="\t", index=False)
