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

# Filter to cg-only and check how 67,061 is obtained
is_cg = probes_df["probe"].str.startswith("cg").values
X_cg = X_raw[is_cg]
probes_cg = probes_df[is_cg].reset_index(drop=True)

# What if we impute missing with median or mean FIRST, then take std?
row_means = np.nanmean(X_cg, axis=1, keepdims=True)
X_imp_mean = np.where(np.isnan(X_cg), row_means, X_cg)
sd_imp_ddof0 = np.std(X_imp_mean, axis=1, ddof=0)
sd_imp_ddof1 = np.std(X_imp_mean, axis=1, ddof=1)
print("Count sd_imp_ddof0 > 0.05:", (sd_imp_ddof0 > 0.05).sum())
print("Count sd_imp_ddof1 > 0.05:", (sd_imp_ddof1 > 0.05).sum())

# Check why APOL3 cg06703343 had HR 0.72 vs HR 0.79 in clinical-only run!
# Let's check what covariates were in Z_clin!
# Wait: in our previous script:
# stage_dummies = pd.get_dummies(meta["stage_clean"], drop_first=True, dtype=float)
# Wait: how did the original script clean ajcc_pathologic_stage.diagnoses or age?
print("\nUnique values of ajcc_pathologic_stage.diagnoses:")
print(meta["ajcc_pathologic_stage.diagnoses"].value_counts(dropna=False))

# What if PCs were computed BEFORE or AFTER standardizing, or with signs flipped?
# Wait: PCA signs don't change linear regression / Cox adjustment!
# Then why did APOL3 cg06703343 have HR 0.72 in our clinical-only run vs HR 0.79 in the user's clinical-only run?
# Wait: look at Table 2 in manuscript.md that I wrote earlier!
# Let's check what `run_cox_and_ml.py` ACTUALLY computed for `cg06703343` in `cox_pc5.tsv`!
cox_pc5 = pd.read_csv("cox_pc5.tsv", sep="\t")
print("\nActual values in cox_pc5.tsv for the top 12 rows + ZCCHC11:")
top_probes = ["cg06703343", "cg25053164", "cg02363950", "cg23152216", "cg21860285", "cg14790396", "cg05664072", "cg13095627", "cg15245051", "cg14222314", "cg24130739", "cg17720013", "cg00754357"]
print(cox_pc5[cox_pc5["probe"].isin(top_probes)][["probe", "Gene", "HR_clin", "p_clin", "q_clin", "HR_pc5", "p_pc5", "q_pc5"]].to_string(index=False))
