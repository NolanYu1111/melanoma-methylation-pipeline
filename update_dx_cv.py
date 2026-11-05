import os
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, accuracy_score
from run_diag_ml import load_geo_fast

dl = r"C:\Users\ethan\Downloads"

X_raw = np.load("tcga_promoter_X_raw.npy")
row_means = np.nanmean(X_raw, axis=1, keepdims=True)
X_imp = np.where(np.isnan(X_raw), row_means, X_raw)
probes_df = pd.read_csv("tcga_promoter_probes.csv")
meta = pd.read_csv("tcga_456_meta.csv")
prim_mask = (meta["is_primary"].values == 1)
met_mask = (meta["is_primary"].values == 0)
p_tcga_map = {p: i for i, p in enumerate(probes_df["probe"])}

X122, probes122, sids122, chars122 = load_geo_fast(os.path.join(dl, "GSE122909_series_matrix.txt.gz"))
is_cell_line = np.array(["cell line" in x.lower() for x in chars122[0]])
is_melanoma = np.array(["melanoma" in x.lower() for x in chars122[1]])
keep_15 = ~is_cell_line
y122_15 = is_melanoma[keep_15].astype(int)
p122_map = {p: i for i, p in enumerate(probes122)}

X120_prom = np.load("gse120878_prom_X.npy")
y120 = np.load("gse120878_y.npy")
prom_probes120 = pd.read_csv("gse120878_prom_probes.csv")["probe"].values

common_all3_mask = np.array([(p in p122_map) and (p in p_tcga_map) for p in prom_probes120])
probes_all3 = prom_probes120[common_all3_mask]
X120_all3 = X120_prom[common_all3_mask].T
X_tcga_all3 = X_imp[[p_tcga_map[p] for p in probes_all3], :].T
X122_all3 = X122[[p122_map[p] for p in probes_all3]][:, keep_15].T
tr_means_all3 = X120_all3.mean(axis=0, keepdims=True)
X122_all3 = np.where(np.isnan(X122_all3), tr_means_all3, X122_all3)

def bootstrap_auc_ci(y_true, y_score, n_boot=2000, seed=42):
    rng = np.random.default_rng(seed)
    pos = np.where(y_true == 1)[0]
    neg = np.where(y_true == 0)[0]
    aucs = []
    for _ in range(n_boot):
        b_pos = rng.choice(pos, size=len(pos), replace=True)
        b_neg = rng.choice(neg, size=len(neg), replace=True)
        b_idx = np.concatenate([b_pos, b_neg])
        aucs.append(roc_auc_score(y_true[b_idx], y_score[b_idx]))
    return np.percentile(aucs, [2.5, 97.5])

rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=2, random_state=42)
fold_top200_prom = []
for train_idx, test_idx in rskf.split(X120_prom.T, y120):
    X_tr, y_tr = X120_prom.T[train_idx], y120[train_idx]
    t_fold, _ = stats.ttest_ind(X_tr[y_tr == 1], X_tr[y_tr == 0], axis=0, equal_var=False)
    top200 = np.argsort(np.abs(t_fold))[::-1][:200]
    fold_top200_prom.append((train_idx, test_idx, top200))

t_full_all3, _ = stats.ttest_ind(X120_all3[y120 == 1], X120_all3[y120 == 0], axis=0, equal_var=False)

dx_rows = []
configs = [
    (10, "L1 Logistic Regression", lambda: LogisticRegression(l1_ratio=1.0, solver="liblinear", C=1.0, random_state=42)),
    (10, "Random Forest", lambda: RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1)),
    (40, "ElasticNet (Conway-style)", lambda: LogisticRegression(l1_ratio=0.5, solver="saga", C=1.0, max_iter=2000, random_state=42)),
    (50, "L1 Logistic Regression", lambda: LogisticRegression(l1_ratio=1.0, solver="liblinear", C=1.0, random_state=42)),
    (50, "Random Forest", lambda: RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1)),
    (200, "L1 Logistic Regression", lambda: LogisticRegression(l1_ratio=1.0, solver="liblinear", C=1.0, random_state=42)),
    (200, "Random Forest", lambda: RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1)),
]

for k_feat, model_name, clf_fn in configs:
    aucs, accs, sens_cv, spec_cv = [], [], [], []
    for train_idx, test_idx, top200 in fold_top200_prom:
        top_idx = top200[:k_feat]
        X_tr, y_tr = X120_prom.T[train_idx][:, top_idx], y120[train_idx]
        X_te, y_te = X120_prom.T[test_idx][:, top_idx], y120[test_idx]
        clf = clf_fn()
        clf.fit(X_tr, y_tr)
        probs = clf.predict_proba(X_te)[:, 1]
        preds = (probs >= 0.5).astype(int)
        aucs.append(roc_auc_score(y_te, probs))
        accs.append(accuracy_score(y_te, preds) * 100.0)
        sens_cv.append(np.mean(preds[y_te == 1]) * 100.0)
        spec_cv.append(np.mean(preds[y_te == 0] == 0) * 100.0)

    top_ext_idx = np.argsort(np.abs(t_full_all3))[::-1][:k_feat]
    clf_ext = clf_fn()
    clf_ext.fit(X120_all3[:, top_ext_idx], y120)

    preds_tcga = (clf_ext.predict_proba(X_tcga_all3[:, top_ext_idx])[:, 1] >= 0.5).astype(int)
    sens_all = np.mean(preds_tcga) * 100.0
    sens_prim = np.mean(preds_tcga[prim_mask]) * 100.0
    sens_met = np.mean(preds_tcga[met_mask]) * 100.0

    probs_122 = clf_ext.predict_proba(X122_all3[:, top_ext_idx])[:, 1]
    auc_122 = roc_auc_score(y122_15, probs_122)
    ci_lo, ci_hi = bootstrap_auc_ci(y122_15, probs_122)
    mcyte_fp = np.mean(probs_122[y122_15 == 0] >= 0.5) * 100.0

    dx_rows.append({
        "Probe_Set": "Promoter-Only",
        "k_Features": k_feat,
        "Model": model_name,
        "CV_AUC_Mean": np.mean(aucs),
        "CV_AUC_SD": np.std(aucs),
        "CV_Acc_Mean": np.mean(accs),
        "CV_Acc_SD": np.std(accs),
        "CV_Sens_Mean": np.mean(sens_cv),
        "CV_Spec_Mean": np.mean(spec_cv),
        "TCGA_Sens_All": sens_all,
        "TCGA_Sens_Primary": sens_prim,
        "TCGA_Sens_Metastatic": sens_met,
        "GSE122909_AUC": auc_122,
        "GSE122909_AUC_CI_Low": ci_lo,
        "GSE122909_AUC_CI_High": ci_hi,
        "GSE122909_Melanocyte_FP_Pct": mcyte_fp
    })

cv_prev = pd.read_csv("gse120878_cv_summary.tsv", sep="\t")
for _, r in cv_prev[cv_prev["Probe_Set"] == "All-Probes (60k)"].iterrows():
    dx_rows.append({
        "Probe_Set": r["Probe_Set"],
        "k_Features": int(r["k_Features"]),
        "Model": r["Model"],
        "CV_AUC_Mean": r["Mean_AUC"],
        "CV_AUC_SD": r["SD_AUC"],
        "CV_Acc_Mean": r["Mean_Acc"],
        "CV_Acc_SD": r["SD_Acc"],
        "CV_Sens_Mean": np.nan,
        "CV_Spec_Mean": np.nan,
        "TCGA_Sens_All": np.nan,
        "TCGA_Sens_Primary": np.nan,
        "TCGA_Sens_Metastatic": np.nan,
        "GSE122909_AUC": np.nan,
        "GSE122909_AUC_CI_Low": np.nan,
        "GSE122909_AUC_CI_High": np.nan,
        "GSE122909_Melanocyte_FP_Pct": np.nan
    })

dx_df = pd.DataFrame(dx_rows)
print(dx_df.to_string(index=False))
dx_df.to_csv("dx_cv.tsv", sep="\t", index=False)
dx_df.to_csv(os.path.join(dl, "dx_cv.tsv"), sep="\t", index=False)
