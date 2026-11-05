import os
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, accuracy_score
from run_diag_ml import load_geo_fast

dl = r"C:\Users\ethan\Downloads"

X122, probes122, sids122, chars122 = load_geo_fast(os.path.join(dl, "GSE122909_series_matrix.txt.gz"))
print("=== GSE122909 SAMPLE CHARACTERISTICS ===")
for i in range(len(sids122)):
    c_str = " | ".join([c[i] for c in chars122])
    print(f"  {i+1:2d}. {sids122[i]}: {c_str}")

is_cell_line = np.array(["cell line" in x.lower() for x in chars122[0]])
is_melanoma = np.array(["melanoma" in x.lower() for x in chars122[1]])
mel_tissue_mask = (~is_cell_line) & is_melanoma
melanocyte_mask = (~is_cell_line) & (~is_melanoma)

p122_map = {p: i for i, p in enumerate(probes122)}
X120_prom = np.load("gse120878_prom_X.npy")
y120 = np.load("gse120878_y.npy")
prom_probes120 = pd.read_csv("gse120878_prom_probes.csv")["probe"].values

common_mask = np.array([p in p122_map for p in prom_probes120])
common_p = prom_probes120[common_mask]
X120_c = X120_prom[common_mask]
X122_c = X122[[p122_map[p] for p in common_p]]

t120, p120 = stats.ttest_ind(X120_c[:, y120==1].T, X120_c[:, y120==0].T, axis=0, equal_var=False)
d120 = np.nanmean(X120_c[:, y120==1], axis=1) - np.nanmean(X120_c[:, y120==0], axis=1)

t122, p122 = stats.ttest_ind(X122_c[:, mel_tissue_mask].T, X122_c[:, melanocyte_mask].T, axis=0, equal_var=False)
d122 = np.nanmean(X122_c[:, mel_tissue_mask], axis=1) - np.nanmean(X122_c[:, melanocyte_mask], axis=1)

pass120_1e10 = (p120 < 1e-10) & (np.abs(d120) > 0.20)
pass120_001 = (p120 < 0.01) & (np.abs(d120) > 0.20)
pass122_001 = (p122 < 0.01) & (np.abs(d122) > 0.20)

print("\n=== FAM174B EXPLORATORY OVERLAP & EXPECTED BY CHANCE ===")
print(f"Total common promoter probes tested (SD > 0.15 in GSE120878): {len(common_p)}")
print(f"GSE120878 hits with |dBeta|>0.20 and p<0.01: {pass120_001.sum()} ({pass120_001.mean()*100:.2f}%)")
print(f"GSE120878 hits with |dBeta|>0.20 and p<1e-10: {pass120_1e10.sum()} ({pass120_1e10.mean()*100:.2f}%)")
print(f"GSE122909 hits with |dBeta|>0.20 and p<0.01: {pass122_001.sum()} ({pass122_001.mean()*100:.2f}%)")

both_001_conc = pass120_001 & pass122_001 & (np.sign(d120) == np.sign(d122))
both_001_disc = pass120_001 & pass122_001 & (np.sign(d120) != np.sign(d122))
exp_001_total = pass120_001.sum() * pass122_001.mean()

both_1e10_conc = pass120_1e10 & pass122_001 & (np.sign(d120) == np.sign(d122))
both_1e10_disc = pass120_1e10 & pass122_001 & (np.sign(d120) != np.sign(d122))
exp_1e10_total = pass120_1e10.sum() * pass122_001.mean()

print(f"At p<0.01 in both: observed {both_001_conc.sum()} concordant, {both_001_disc.sum()} discordant (expected by chance under independence: {exp_001_total:.2f} total, ~{0.5*exp_001_total:.2f} concordant)")
print(f"At p<1e-10 in GSE120878 & p<0.01 in GSE122909: observed {both_1e10_conc.sum()} concordant, {both_1e10_disc.sum()} discordant (expected by chance: {exp_1e10_total:.2f} total, ~{0.5*exp_1e10_total:.2f} concordant)")
for idx in np.where(both_001_conc)[0]:
    print(f"  Concordant probe: {common_p[idx]}: d120={d120[idx]:+.3f} (p={p120[idx]:.2e}), d122={d122[idx]:+.3f} (p={p122[idx]:.2e})")

# Bootstrap 95% CI for GSE122909 ROC-AUC (n=15: 7 melanomas, 8 melanocytes)
print("\n=== GSE122909 95% CONFIDENCE INTERVALS & CONWAY-STYLE ELASTICNET 40-CpG COMPARISON ===")
keep_15 = ~is_cell_line
y122_15 = is_melanoma[keep_15].astype(int)
X122_15 = X122_c[:, keep_15].T # (15, 13243)
X120_T = X120_c.T # (162, 13243)
X122_15 = np.where(np.isnan(X122_15), X120_T.mean(axis=0, keepdims=True), X122_15)

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

t_full, _ = stats.ttest_ind(X120_T[y120 == 1], X120_T[y120 == 0], axis=0, equal_var=False)
for k_feat in [10, 40, 50, 200]:
    top_idx = np.argsort(np.abs(t_full))[::-1][:k_feat]
    clf_l1 = LogisticRegression(penalty="l1", solver="liblinear", C=1.0, random_state=42)
    clf_l1.fit(X120_T[:, top_idx], y120)
    probs = clf_l1.predict_proba(X122_15[:, top_idx])[:, 1]
    auc_val = roc_auc_score(y122_15, probs)
    ci_lo, ci_hi = bootstrap_auc_ci(y122_15, probs)
    print(f"Promoter L1 Logistic (k={k_feat}): GSE122909 AUC = {auc_val:.3f} (95% CI: {ci_lo:.3f}-{ci_hi:.3f})")

# Also test ElasticNet (l1_ratio=0.5, Conway et al. method) with k=40 features on the same 5-fold CV (repeated twice) and GSE122909
rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=2, random_state=42)
enet_aucs, enet_accs, enet_sens, enet_spec = [], [], [], []
for tr_idx, te_idx in rskf.split(X120_T, y120):
    X_tr, y_tr = X120_T[tr_idx], y120[tr_idx]
    X_te, y_te = X120_T[te_idx], y120[te_idx]
    tf, _ = stats.ttest_ind(X_tr[y_tr == 1], X_tr[y_tr == 0], axis=0, equal_var=False)
    top40 = np.argsort(np.abs(tf))[::-1][:40]
    clf_en = LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5, C=1.0, max_iter=2000, random_state=42)
    clf_en.fit(X_tr[:, top40], y_tr)
    pr = clf_en.predict_proba(X_te[:, top40])[:, 1]
    pd_lbl = (pr >= 0.5).astype(int)
    enet_aucs.append(roc_auc_score(y_te, pr))
    enet_accs.append(accuracy_score(y_te, pd_lbl) * 100.0)
    enet_sens.append(np.mean(pd_lbl[y_te == 1]) * 100.0)
    enet_spec.append(np.mean(pd_lbl[y_te == 0] == 0) * 100.0)

top40_full = np.argsort(np.abs(t_full))[::-1][:40]
clf_en_full = LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5, C=1.0, max_iter=2000, random_state=42)
clf_en_full.fit(X120_T[:, top40_full], y120)
pr_122_en = clf_en_full.predict_proba(X122_15[:, top40_full])[:, 1]
auc_122_en = roc_auc_score(y122_15, pr_122_en)
ci_en = bootstrap_auc_ci(y122_15, pr_122_en)
print(f"ElasticNet 40-CpG (Conway-style): GSE120878 CV AUC={np.mean(enet_aucs):.3f}+/-{np.std(enet_aucs):.3f}, Acc={np.mean(enet_accs):.1f}%+/-{np.std(enet_accs):.1f}%, Sens={np.mean(enet_sens):.1f}%, Spec={np.mean(enet_spec):.1f}% | GSE122909 AUC={auc_122_en:.3f} (95% CI: {ci_en[0]:.3f}-{ci_en[1]:.3f})")
