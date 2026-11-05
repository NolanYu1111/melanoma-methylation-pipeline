"""
Permutation FDR calibration check across all 5 TCGA cohorts (SKCM, HNSC, LUAD, BLCA, LIHC).
Tests whether the covariate-adjusted Cox score test + Benjamini-Hochberg FDR procedure
produces false-positive inflation when conditioning on 3 PCs, 5 PCs, or the RNA-seq immune score.
Specifically, we permute survival outcomes (time, event) across patients (50 permutations per model
and cohort) while keeping covariates (clinical + PCs / immune score) and the promoter methylation
matrix fixed (preserving all CpG-PC and CpG-immune collinearity), and record:
  1. Mean percentage of probes with p < 0.05 (nominal 5.0% expected under calibration)
  2. Mean, median, and proportion of permutations yielding 0 hits at Benjamini-Hochberg q < 0.05.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize
from sklearn.decomposition import PCA

WORK_DIR = r"C:\Users\ethan\.gemini\antigravity\scratch\melanoma_methylation_paper"
DL_DIR = r"C:\Users\ethan\Downloads"

IMMUNE_MARKERS = [
    "PTPRC", "CD8A", "CD8B", "CD3D", "CD3E", "GZMA", "PRF1",
    "IFNG", "CXCL9", "CXCL10", "STAT1", "PSMB9", "TAP1", "HLA-DRA"
]

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
    return p_vals

def run_calibration():
    os.chdir(WORK_DIR)
    rng = np.random.RandomState(2026)
    N_PERMS = 50
    rows = []

    for code in ["SKCM", "HNSC", "LUAD", "BLCA", "LIHC"]:
        if code == "SKCM":
            X_raw = np.load("tcga_promoter_X_raw.npy")
            meta = pd.read_csv("tcga_456_meta.csv")
        else:
            X_raw = np.load(f"tcga_{code.lower()}_promoter_X_raw.npy")
            meta = pd.read_csv(f"tcga_{code.lower()}_meta.csv")

        row_means = np.nanmean(X_raw, axis=1, keepdims=True)
        X_imp = np.where(np.isnan(X_raw), row_means, X_raw)
        row_sds = np.std(X_imp, axis=1, ddof=1, keepdims=True)
        X_std = (X_imp - row_means) / row_sds

        pca = PCA(n_components=5, random_state=42)
        pcs = pca.fit_transform(X_std.T)

        rna_path = os.path.join(DL_DIR, f"TCGA-{code}.HiSeqV2.tsv.gz")
        rnaseq_df = pd.read_csv(rna_path, sep="\t", index_col=0)
        meta["sample_15"] = meta["sample"].str.slice(0, 15)
        has_rna = meta["sample_15"].isin(rnaseq_df.columns).values
        rna_samples = meta.loc[has_rna, "sample_15"].values
        rna_sub = rnaseq_df[rna_samples]

        avail_markers = [g for g in IMMUNE_MARKERS if g in rna_sub.index]
        imm_expr = rna_sub.loc[avail_markers].values.astype(float)
        imm_z = (imm_expr - imm_expr.mean(axis=1, keepdims=True)) / imm_expr.std(axis=1, ddof=1, keepdims=True)
        imm_score = imm_z.mean(axis=0)
        imm_score_std = (imm_score - imm_score.mean()) / imm_score.std(ddof=1)

        sex_male = (meta["gender.demographic"].astype(str).str.lower() == "male").astype(float).values
        age = meta["age_at_index.demographic"].values.astype(float)
        age_std = (age - age.mean()) / age.std(ddof=1)
        stage_dummies = pd.get_dummies(meta["stage_clean"], drop_first=True, dtype=float).values
        cov_list = [age_std, sex_male]
        if meta["is_primary"].nunique() > 1:
            cov_list.append(meta["is_primary"].values.astype(float))
        if stage_dummies.shape[1] > 0:
            cov_list.append(stage_dummies)
        Z_clin = np.column_stack(cov_list)

        time_arr = meta["OS.time"].values.astype(float)
        event_arr = meta["OS"].values.astype(int)

        models = [
            ("Clinical Only", Z_clin, X_std, time_arr, event_arr),
            ("Clinical + 3 PCs", np.column_stack([Z_clin, pcs[:, :3]]), X_std, time_arr, event_arr),
            ("Clinical + 5 PCs", np.column_stack([Z_clin, pcs[:, :5]]), X_std, time_arr, event_arr),
            ("Clinical + RNA Immune", np.column_stack([Z_clin[has_rna], imm_score_std]), X_std[:, has_rna], time_arr[has_rna], event_arr[has_rna]),
        ]

        for model_name, Z_mat, X_mat, t_vec, d_vec in models:
            n_s = len(t_vec)
            n_p = X_mat.shape[0]
            p05_pcts = []
            q05_hits = []
            for _ in range(N_PERMS):
                perm = rng.permutation(n_s)
                t_p = t_vec[perm]
                d_p = d_vec[perm]
                ord_p = np.argsort(t_p)
                pvals = fit_breslow_null_and_score_test(Z_mat[ord_p], X_mat[:, ord_p], t_p[ord_p], d_p[ord_p])
                qvals = bh_fdr(pvals)
                p05_pcts.append(100.0 * np.mean(pvals < 0.05))
                q05_hits.append(int(np.sum(qvals < 0.05)))
            mean_p05 = float(np.mean(p05_pcts))
            mean_q05 = float(np.mean(q05_hits))
            med_q05 = float(np.median(q05_hits))
            zero_frac = int(np.sum(np.array(q05_hits) == 0))
            ge1_frac = int(np.sum(np.array(q05_hits) >= 1))
            pct_ge1 = 100.0 * ge1_frac / N_PERMS
            max_q05 = int(np.max(q05_hits))
            print(f"TCGA-{code} | {model_name:21s} | mean(p<0.05)={mean_p05:.2f}% | P(>=1 hit at q<0.05)={ge1_frac}/{N_PERMS} ({pct_ge1:.0f}%) | median(q<0.05)={med_q05:.0f} | mean(q<0.05)={mean_q05:.2f}")
            rows.append({
                "cohort": f"TCGA-{code}",
                "model": model_name,
                "n_patients": n_s,
                "n_probes": n_p,
                "n_perms": N_PERMS,
                "mean_pct_p_lt_005": round(mean_p05, 3),
                "ge1_hit_perms": f"{ge1_frac}/{N_PERMS}",
                "pct_perms_ge1_q05": round(pct_ge1, 1),
                "median_q05_hits": int(med_q05),
                "mean_q05_hits": round(mean_q05, 2),
                "zero_hit_perms": f"{zero_frac}/{N_PERMS}",
                "max_q05_hits": max_q05,
            })

    out_df = pd.DataFrame(rows)
    out_df.to_csv("permutation_fdr_calibration.tsv", sep="\t", index=False)
    print("Saved permutation_fdr_calibration.tsv")

if __name__ == "__main__":
    run_calibration()
