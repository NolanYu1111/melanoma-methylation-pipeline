import os
import gzip
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

X_raw = np.load("tcga_promoter_X_raw.npy")
probes_df = pd.read_csv("tcga_promoter_probes.csv")
meta = pd.read_csv("tcga_456_meta.csv")

# Test standardizing observed values BEFORE imputing 0 (equivalent to mean imputation before standardizing, except SD is over non-missing vs all 456!)
obs_mean = np.nanmean(X_raw, axis=1, keepdims=True)
obs_sd = np.nanstd(X_raw, axis=1, ddof=1, keepdims=True)
X_std_obs = np.where(np.isnan(X_raw), 0.0, (X_raw - obs_mean) / obs_sd)

# Also test ddof=0
obs_sd0 = np.nanstd(X_raw, axis=1, ddof=0, keepdims=True)
X_std_obs0 = np.where(np.isnan(X_raw), 0.0, (X_raw - obs_mean) / obs_sd0)

pca = PCA(n_components=5, random_state=42)
pcs = pca.fit_transform(X_std_obs.T)

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
X_s = X_std_obs[:, order]

from run_cox_and_ml import fit_breslow_null_and_score_test
z0, p0, hr0, _ = fit_breslow_null_and_score_test(Z_clin_s, X_s, t_s, d_s)
q0 = bh_fdr(p0)
print(f"X_std_obs: p<0.05 = {(p0 < 0.05).sum()}, q<0.05 = {(q0 < 0.05).sum()} (prot: {((q0<0.05)&(hr0<1)).sum()}, haz: {((q0<0.05)&(hr0>1)).sum()})")
for pid in ["cg06703343", "cg02363950", "cg15245051", "cg23152216", "cg00754357"]:
    idx = probes_df[probes_df["probe"] == pid].index[0]
    print(f"  {pid}: HR={hr0[idx]:.4f}, p={p0[idx]:.4e}, q={q0[idx]:.4e}")
