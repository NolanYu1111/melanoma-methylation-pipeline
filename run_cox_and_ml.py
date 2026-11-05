import os
import gzip
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, accuracy_score, roc_curve

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

X_raw = np.load("tcga_promoter_X_raw.npy")
probes_df = pd.read_csv("tcga_promoter_probes.csv")
meta = pd.read_csv("tcga_456_meta.csv")
print("Loaded X_raw shape:", X_raw.shape)
print("Loaded meta shape:", meta.shape)

# Map probes to genes using GPL13534
gpl = pd.read_csv(os.path.join(dl, "GPL13534-11288.txt"), sep="\t", skiprows=37, usecols=["ID", "UCSC_RefGene_Name", "UCSC_RefGene_Group"], low_memory=False).dropna(subset=["ID"])
probe_to_gene = dict(zip(gpl["ID"], gpl["UCSC_RefGene_Name"].fillna("").apply(lambda s: s.split(";")[0] if s else "")))
probe_to_allgenes = dict(zip(gpl["ID"], gpl["UCSC_RefGene_Name"].fillna("")))
probes_df["Gene"] = probes_df["probe"].map(probe_to_gene)
probes_df["AllGenes"] = probes_df["probe"].map(probe_to_allgenes)

# Check KRAS cg23152216 and ZCCHC11 cg00754357
for pid in ["cg23152216", "cg00754357"]:
    idx = probes_df[probes_df["probe"] == pid].index
    print(f"Probe {pid} in probes_df:", len(idx), probes_df.loc[idx].to_dict(orient="records") if len(idx) else "None")

# Impute missing values in X_raw with row mean and standardize each CpG (mean 0, SD 1)
row_means = np.nanmean(X_raw, axis=1, keepdims=True)
X_imp = np.where(np.isnan(X_raw), row_means, X_raw)
row_sds = np.std(X_imp, axis=1, ddof=1, keepdims=True)
X_std = (X_imp - row_means) / row_sds # (67061, 456)

# Compute Top 5 PCs of standardized promoter methylation matrix (samples x probes)
pca = PCA(n_components=5, random_state=42)
pcs = pca.fit_transform(X_std.T) # (456, 5)
print("Explained variance ratio of top 5 PCs:", pca.explained_variance_ratio_)

# Check correlation of PC1 with immune genes (PSMB9, TAP1, PTPRC)
for gname in ["PSMB9", "TAP1", "PTPRC"]:
    g_idxs = probes_df[probes_df["AllGenes"].str.contains(gname, na=False)].index
    if len(g_idxs) > 0:
        corrs = [np.corrcoef(pcs[:, 0], X_std[i, :])[0, 1] for i in g_idxs]
        print(f"Max abs corr of PC1 with {gname} promoter CpGs:", np.max(np.abs(corrs)))

# Build clinical covariates matrix Z: age, sex, AJCC stage (0, I, II, III, IV, unknown), primary vs metastatic
meta["sex_male"] = (meta["gender.demographic"].str.lower() == "male").astype(float)
age = meta["age_at_index.demographic"].values.astype(float)
age_std = (age - np.mean(age)) / np.std(age, ddof=1)
is_prim = meta["is_primary"].values.astype(float)

# Stage dummies (reference: Stage II or Stage I, doesn't affect partial likelihood score test!)
stage_dummies = pd.get_dummies(meta["stage_clean"], drop_first=True, dtype=float)
Z_clin = np.column_stack([age_std, meta["sex_male"].values, is_prim, stage_dummies.values])
print("Z_clin shape:", Z_clin.shape)

# Implement exact Breslow Cox fit for null model and vectorized score test for all 67,061 CpGs
time_arr = meta["OS.time"].values.astype(float)
event_arr = meta["OS"].values.astype(int)

# Sort samples by time ascending
order = np.argsort(time_arr)
t_sorted = time_arr[order]
d_sorted = event_arr[order]
Z_clin_s = Z_clin[order]
pcs_s = pcs[order]
X_std_s = X_std[:, order] # (67061, 456)

# Unique event times for Breslow ties
uniq_times, inv_idx = np.unique(t_sorted[d_sorted == 1], return_inverse=True)

def fit_breslow_null_and_score_test(Z_cov, X_genes, t_s, d_s):
    """
    Fits Cox proportional hazards model with Breslow ties for covariates Z_cov (N, p),
    then performs vectorized score test for each row of X_genes (M, N) adjusting for Z_cov.
    Note: Can handle missing values in X_genes or use standardized X_genes.
    """
    N, p = Z_cov.shape
    # Group by unique failure times
    utimes = np.unique(t_s[d_s == 1])
    n_events_at_t = np.array([np.sum((t_s == ut) & (d_s == 1)) for ut in utimes], dtype=float)
    Z_events_sum = np.array([np.sum(Z_cov[(t_s == ut) & (d_s == 1)], axis=0) for ut in utimes]) # (K, p)
    risk_masks = np.array([t_s >= ut for ut in utimes]) # (K, N)
    
    def neg_log_lik_and_grad(beta):
        eta = np.clip(Z_cov @ beta, -50, 50)
        w = np.exp(eta) # (N,)
        # For each unique event time k, sum(w) over risk set
        S0 = risk_masks @ w # (K,)
        S1 = (risk_masks * w[None, :]) @ Z_cov # (K, p)
        ll = np.sum(Z_events_sum * beta[None, :]) - np.sum(n_events_at_t * np.log(S0))
        grad = np.sum(Z_events_sum, axis=0) - np.sum(n_events_at_t[:, None] * (S1 / S0[:, None]), axis=0)
        return -ll, -grad

    res = minimize(neg_log_lik_and_grad, np.zeros(p), jac=True, method="BFGS")
    beta_hat = res.x
    eta = np.clip(Z_cov @ beta_hat, -50, 50)
    w = np.exp(eta) # (N,)
    
    S0 = risk_masks @ w # (K,)
    E_Z = ((risk_masks * w[None, :]) @ Z_cov) / S0[:, None] # (K, p)
    
    # Information matrix I_ZZ (p, p)
    I_ZZ = np.zeros((p, p))
    for k in range(len(utimes)):
        mask_k = risk_masks[k]
        wk = w[mask_k] / S0[k]
        Zk_centered = Z_cov[mask_k] - E_Z[k]
        I_ZZ += n_events_at_t[k] * (Zk_centered.T * wk) @ Zk_centered
    I_ZZ_inv = np.linalg.inv(I_ZZ)
    
    # Vectorized score statistic for X_genes (M, N)
    # Score U_x = sum_{i: d_i=1} x_i - sum_k d_k * S1_x(k) / S0(k)
    # Using martingale residuals: r_i = d_i - w_i * sum_{k: t_k <= t_i} (d_k / S0_k)
    cum_hazard_factor = (risk_masks.T @ (n_events_at_t / S0)) # (N,)
    martingale_resid = d_s - w * cum_hazard_factor # (N,)
    
    U_x = X_genes @ martingale_resid # (M,)
    
    # Expected X at each failure time: E_X = (X_genes @ (risk_masks.T * w[:, None])) / S0[None, :] -> (M, K)
    weighted_risk = risk_masks.T * w[:, None] # (N, K)
    E_X = (X_genes @ weighted_risk) / S0[None, :] # (M, K)
    
    # I_XX = sum_k d_k * ( sum_{i in R_k} w_i x_i^2 / S0_k - E_X_k^2 )
    X_sq = X_genes ** 2
    E_X2 = (X_sq @ weighted_risk) / S0[None, :] # (M, K)
    I_XX = np.sum(n_events_at_t[None, :] * (E_X2 - E_X ** 2), axis=1) # (M,)
    
    # Cross information I_XZ (M, p):
    # sum_k d_k * ( (sum_{i in R_k} w_i x_i Z_i / S0_k) - E_X_k E_Z_k )
    # Note: sum_k d_k * (sum_{i in R_k} w_i x_i Z_i / S0_k) = X_genes @ (w[:, None] * sum_{k: t_k <= t_i} (d_k / S0_k) * Z_i - wait: Z_i is independent of k!)
    # Yes! Because Z_i and x_i do not depend on k!
    # So sum_k d_k * sum_{i in R_k} (w_i / S0_k) x_i Z_i = sum_i x_i (w_i * cum_hazard_factor_i) Z_i !
    h_w = w * cum_hazard_factor # (N,)
    term1_XZ = X_genes @ (Z_cov * h_w[:, None]) # (M, p)
    term2_XZ = (E_X * n_events_at_t[None, :]) @ E_Z # (M, p)
    I_XZ = term1_XZ - term2_XZ # (M, p)
    
    # Variance of score U_x adjusted for Z: V_x = I_XX - diag(I_XZ @ I_ZZ_inv @ I_XZ.T)
    V_x = I_XX - np.sum((I_XZ @ I_ZZ_inv) * I_XZ, axis=1)
    z_stat = U_x / np.sqrt(np.maximum(V_x, 1e-12))
    p_vals = 2.0 * stats.norm.sf(np.abs(z_stat))
    # One-step HR approximation: exp(U_x / V_x)
    beta_1step = U_x / np.maximum(V_x, 1e-12)
    hr_est = np.exp(beta_1step)
    return z_stat, p_vals, hr_est, beta_hat

if __name__ == "__main__":
    print("\nRunning Cox score tests...")
    z0, p0, hr0, _ = fit_breslow_null_and_score_test(Z_clin_s, X_std_s, t_sorted, d_sorted)
    q0 = bh_fdr(p0)
    print(f"Unadjusted for PCs (Clinical only): p<0.05 = {(p0 < 0.05).sum()}, q<0.05 = {(q0 < 0.05).sum()} (protective: {((q0 < 0.05) & (hr0 < 1)).sum()}, hazardous: {((q0 < 0.05) & (hr0 > 1)).sum()})")

    z3, p3, hr3, _ = fit_breslow_null_and_score_test(np.column_stack([Z_clin_s, pcs_s[:, :3]]), X_std_s, t_sorted, d_sorted)
    q3 = bh_fdr(p3)
    print(f"Adjusted for 3 PCs: p<0.05 = {(p3 < 0.05).sum()}, q<0.05 = {(q3 < 0.05).sum()}, q<0.01 = {(q3 < 0.01).sum()}")

    z5, p5, hr5, _ = fit_breslow_null_and_score_test(np.column_stack([Z_clin_s, pcs_s[:, :5]]), X_std_s, t_sorted, d_sorted)
    q5 = bh_fdr(p5)
    print(f"Adjusted for 5 PCs: p<0.05 = {(p5 < 0.05).sum()}, q<0.05 = {(q5 < 0.05).sum()}, q<0.01 = {(q5 < 0.01).sum()}")
