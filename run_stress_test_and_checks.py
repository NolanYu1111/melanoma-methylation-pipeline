import os
import gzip
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, accuracy_score

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

def concordance_index(event_times, predicted_risk, event_observed):
    """Exact Harrell's C-index (higher predicted_risk = shorter survival)."""
    n = len(event_times)
    num = 0.0
    den = 0.0
    for i in range(n):
        if event_observed[i] == 1:
            for j in range(n):
                if event_times[j] > event_times[i]:
                    den += 1.0
                    if predicted_risk[i] > predicted_risk[j]:
                        num += 1.0
                    elif predicted_risk[i] == predicted_risk[j]:
                        num += 0.5
                elif event_times[j] == event_times[i] and event_observed[j] == 0:
                    den += 1.0
                    if predicted_risk[i] > predicted_risk[j]:
                        num += 1.0
                    elif predicted_risk[i] == predicted_risk[j]:
                        num += 0.5
    return num / den if den > 0 else 0.5

# 1. Load Cox results and check the 2 CpGs that remain at q < 0.05 after Clinical + RNA-seq Immune Score
X_raw = np.load("tcga_promoter_X_raw.npy")
probes_df = pd.read_csv("tcga_promoter_probes.csv")
meta = pd.read_csv("tcga_456_meta.csv")
cox_pc5 = pd.read_csv("cox_pc5.tsv", sep="\t")

# Load RNA-seq immune score from run_scientific_upgrades logic
rnaseq_df = pd.read_csv(os.path.join(dl, "TCGA-SKCM.HiSeqV2.tsv.gz"), sep="\t", index_col=0)
meta["sample_15"] = meta["sample"].str.slice(0, 15)
has_rna_mask = meta["sample_15"].isin(rnaseq_df.columns).values
rna_samples_matched = meta.loc[has_rna_mask, "sample_15"].values
rna_matched_df = rnaseq_df[rna_samples_matched]

immune_markers = ["PTPRC", "CD8A", "CD8B", "CD3D", "CD3E", "GZMA", "PRF1", "IFNG", "CXCL9", "CXCL10", "STAT1", "PSMB9", "TAP1", "HLA-DRA"]
immune_expr = rna_matched_df.loc[immune_markers].values
immune_z = (immune_expr - immune_expr.mean(axis=1, keepdims=True)) / immune_expr.std(axis=1, ddof=1, keepdims=True)
immune_score_454 = immune_z.mean(axis=0)
immune_score_454_std = (immune_score_454 - np.mean(immune_score_454)) / np.std(immune_score_454, ddof=1)
immune_full = np.zeros(len(meta))
immune_full[has_rna_mask] = immune_score_454_std

row_means = np.nanmean(X_raw, axis=1, keepdims=True)
X_imp = np.where(np.isnan(X_raw), row_means, X_raw)
row_sds = np.std(X_imp, axis=1, ddof=1, keepdims=True)
X_std = (X_imp - row_means) / row_sds

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
imm1_s = immune_full[order]
X_std_s = X_std[:, order]

from run_cox_and_ml import fit_breslow_null_and_score_test
z_ci1, p_ci1, hr_ci1, _ = fit_breslow_null_and_score_test(np.column_stack([Z_clin_s, imm1_s]), X_std_s, t_s, d_s)
q_ci1 = bh_fdr(p_ci1)

gpl = pd.read_csv(os.path.join(dl, "GPL13534-11288.txt"), sep="\t", skiprows=37, usecols=["ID", "UCSC_RefGene_Name", "Probe_SNPs", "Probe_SNPs_10"], low_memory=False).dropna(subset=["ID"])
probe_to_gene = dict(zip(gpl["ID"], gpl["UCSC_RefGene_Name"].fillna("").apply(lambda s: s.split(";")[0] if s else "Unannotated")))
probe_to_snp = dict(zip(gpl["ID"], gpl["Probe_SNPs"].fillna("")))
probe_to_snp10 = dict(zip(gpl["ID"], gpl["Probe_SNPs_10"].fillna("")))

probes_df["Gene"] = probes_df["probe"].map(probe_to_gene)
probes_df["HR_imm"] = hr_ci1
probes_df["p_imm"] = p_ci1
probes_df["q_imm"] = q_ci1
print("Top 10 CpGs after adjusting for Clinical + RNA-seq Immune Score:")
print(probes_df.sort_values("p_imm")[["probe", "Gene", "HR_imm", "p_imm", "q_imm"]].head(10).to_string(index=False), flush=True)

# 2. Check Stratified Cox (stratified by primary vs metastatic) and PH interaction with log(time) for Table 2 candidates
print("\n--- Stratified Cox (primary vs metastatic strata) & Proportional Hazards Check ---", flush=True)
Z_clin_nostrat = np.column_stack([age_std, meta["sex_male"].values, stage_dummies.values])
prim_mask = (is_prim == 1.0)
met_mask = (is_prim == 0.0)

# Score test separately in primary and metastatic and combine via inverse-variance / score addition
def stratified_cox_hits(Z_cov, X_mat, t_vec, d_vec, strata_bool):
    # Stratum 1
    m1 = strata_bool
    o1 = np.argsort(t_vec[m1])
    z1, p1, hr1, _ = fit_breslow_null_and_score_test(Z_cov[m1][o1], X_mat[:, m1][:, o1], t_vec[m1][o1], d_vec[m1][o1])
    # Stratum 2
    m2 = ~strata_bool
    o2 = np.argsort(t_vec[m2])
    z2, p2, hr2, _ = fit_breslow_null_and_score_test(Z_cov[m2][o2], X_mat[:, m2][:, o2], t_vec[m2][o2], d_vec[m2][o2])
    # Meta-analytic / stratified score combination weighted by sqrt(events)
    w1 = np.sqrt(d_vec[m1].sum())
    w2 = np.sqrt(d_vec[m2].sum())
    z_strat = (w1 * z1 + w2 * z2) / np.sqrt(w1**2 + w2**2)
    p_strat = 2.0 * stats.norm.sf(np.abs(z_strat))
    q_strat = bh_fdr(p_strat)
    return z_strat, p_strat, q_strat

z_strat_clin, p_strat_clin, q_strat_clin = stratified_cox_hits(Z_clin_nostrat, X_std, time_arr, event_arr, prim_mask)
print(f"Stratified Clinical-only Cox screen: {(p_strat_clin < 0.05).sum()} at p<0.05, {(q_strat_clin < 0.05).sum()} at q<0.05", flush=True)

# Check how many of the 67,061 tested probes have a probe SNP in GPL13534
has_snp = probes_df["probe"].map(lambda x: bool(probe_to_snp.get(x, "")) or bool(probe_to_snp10.get(x, ""))).values
print(f"Probes with annotated SNPs (Probe_SNPs or Probe_SNPs_10): {has_snp.sum()} / {len(probes_df)} ({has_snp.mean()*100:.1f}%)", flush=True)

# Check Table 2 candidates for SNP status and PH violation (correlation of CpG * log(time) interaction)
t2_probes = ["cg06703343", "cg25053164", "cg02363950", "cg23152216", "cg21860285", "cg14790396", "cg05664072", "cg13095627", "cg15245051", "cg14222314", "cg24130739", "cg17720013", "cg00754357"]
p_tcga_map = {p: i for i, p in enumerate(probes_df["probe"])}
log_t_s = np.log(t_s + 1.0)
log_t_s = (log_t_s - np.mean(log_t_s)) / np.std(log_t_s, ddof=1)

print("\nTable 2 Candidates SNP & Proportional-Hazards (CpG x log(time)) Check:")
for pid in t2_probes:
    idx = p_tcga_map[pid]
    x_vec = X_std_s[idx]
    inter_vec = (x_vec * log_t_s).reshape(1, -1)
    Z_base = np.column_stack([Z_clin_s, x_vec])
    _, p_ph, _, _ = fit_breslow_null_and_score_test(Z_base, inter_vec, t_s, d_s)
    print(f"  {pid} ({probe_to_gene.get(pid, '')}): SNP={probe_to_snp.get(pid, 'None') or 'None'}, SNP10={probe_to_snp10.get(pid, 'None') or 'None'}, PH_p={p_ph[0]:.4f}", flush=True)

# 3. 5-Fold Cross-Validation C-Index & 95% CI on Differences
from sklearn.decomposition import PCA
pca5 = PCA(n_components=5, random_state=42)
pcs5 = pca5.fit_transform(X_std.T)
pcs5_std = (pcs5 - pcs5.mean(axis=0, keepdims=True)) / pcs5.std(axis=0, ddof=1, keepdims=True)

skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
c_clin_folds, c_imm_folds, c_top50_folds, c_imm_top50_folds = [], [], [], []

for tr_idx, te_idx in skf.split(Z_clin, event_arr):
    ord_tr = np.argsort(time_arr[tr_idx])
    t_tr, d_tr = time_arr[tr_idx][ord_tr], event_arr[tr_idx][ord_tr]
    Z_tr = Z_clin[tr_idx][ord_tr]
    _, _, _, beta_clin = fit_breslow_null_and_score_test(Z_tr, np.zeros((1, len(tr_idx))), t_tr, d_tr)
    risk_clin = Z_clin[te_idx] @ beta_clin
    c_clin_folds.append(concordance_index(time_arr[te_idx], risk_clin, event_arr[te_idx]))

    Z_tr_imm = np.column_stack([Z_clin[tr_idx], immune_full[tr_idx]])[ord_tr]
    _, _, _, beta_imm = fit_breslow_null_and_score_test(Z_tr_imm, np.zeros((1, len(tr_idx))), t_tr, d_tr)
    risk_imm = np.column_stack([Z_clin[te_idx], immune_full[te_idx]]) @ beta_imm
    c_imm_folds.append(concordance_index(time_arr[te_idx], risk_imm, event_arr[te_idx]))

    # In-fold top 50 CpGs adjusted for clinical + 5 PCs
    Z_tr_pc5 = np.column_stack([Z_clin[tr_idx], pcs5_std[tr_idx]])[ord_tr]
    X_tr_s = X_std[:, tr_idx][:, ord_tr]
    z_f, p_f, _, _ = fit_breslow_null_and_score_test(Z_tr_pc5, X_tr_s, t_tr, d_tr)
    top50_idx = np.argsort(p_f)[:50]
    # Penalized Cox (L2 / elastic-net style shrinkage on CpGs + clinical)
    # Using score-weighted composite risk of top 50 CpGs + clinical
    comp_score_tr = (X_std[top50_idx][:, tr_idx].T @ (z_f[top50_idx] / np.sqrt(50)))
    comp_score_te = (X_std[top50_idx][:, te_idx].T @ (z_f[top50_idx] / np.sqrt(50)))
    c_mean, c_sd = np.mean(comp_score_tr), np.std(comp_score_tr, ddof=1)
    comp_score_tr = (comp_score_tr - c_mean) / c_sd
    comp_score_te = (comp_score_te - c_mean) / c_sd

    Z_tr_c50 = np.column_stack([Z_clin[tr_idx], comp_score_tr])[ord_tr]
    _, _, _, beta_c50 = fit_breslow_null_and_score_test(Z_tr_c50, np.zeros((1, len(tr_idx))), t_tr, d_tr)
    risk_c50 = np.column_stack([Z_clin[te_idx], comp_score_te]) @ beta_c50
    c_top50_folds.append(concordance_index(time_arr[te_idx], risk_c50, event_arr[te_idx]))

    Z_tr_imm_c50 = np.column_stack([Z_clin[tr_idx], immune_full[tr_idx], comp_score_tr])[ord_tr]
    _, _, _, beta_imm_c50 = fit_breslow_null_and_score_test(Z_tr_imm_c50, np.zeros((1, len(tr_idx))), t_tr, d_tr)
    risk_imm_c50 = np.column_stack([Z_clin[te_idx], immune_full[te_idx], comp_score_te]) @ beta_imm_c50
    c_imm_top50_folds.append(concordance_index(time_arr[te_idx], risk_imm_c50, event_arr[te_idx]))

print(f"\n5-Fold CV C-Index Clinical Only: {np.mean(c_clin_folds):.3f} +/- {np.std(c_clin_folds):.3f}", flush=True)
print(f"5-Fold CV C-Index Clinical + RNA-seq Immune Score: {np.mean(c_imm_folds):.3f} +/- {np.std(c_imm_folds):.3f}", flush=True)
print(f"5-Fold CV C-Index Clinical + Top 50 5-PC CpGs: {np.mean(c_top50_folds):.3f} +/- {np.std(c_top50_folds):.3f}", flush=True)
print(f"5-Fold CV C-Index Clinical + Immune Score + Top 50 5-PC CpGs: {np.mean(c_imm_top50_folds):.3f} +/- {np.std(c_imm_top50_folds):.3f}", flush=True)

diff_imm = np.array(c_imm_folds) - np.array(c_clin_folds)
diff_c50 = np.array(c_top50_folds) - np.array(c_clin_folds)
diff_imm_c50 = np.array(c_imm_top50_folds) - np.array(c_imm_folds)
print(f"Delta C-Index (Immune Score vs Clinical): {np.mean(diff_imm):+.3f} (95% CI: {np.mean(diff_imm)-1.96*np.std(diff_imm, ddof=1)/np.sqrt(5):+.3f} to {np.mean(diff_imm)+1.96*np.std(diff_imm, ddof=1)/np.sqrt(5):+.3f})", flush=True)
print(f"Delta C-Index (Top 50 CpGs vs Clinical): {np.mean(diff_c50):+.3f} (95% CI: {np.mean(diff_c50)-1.96*np.std(diff_c50, ddof=1)/np.sqrt(5):+.3f} to {np.mean(diff_c50)+1.96*np.std(diff_c50, ddof=1)/np.sqrt(5):+.3f})", flush=True)
print(f"Delta C-Index (Top 50 CpGs + Immune vs Immune alone): {np.mean(diff_imm_c50):+.3f} (95% CI: {np.mean(diff_imm_c50)-1.96*np.std(diff_imm_c50, ddof=1)/np.sqrt(5):+.3f} to {np.mean(diff_imm_c50)+1.96*np.std(diff_imm_c50, ddof=1)/np.sqrt(5):+.3f})", flush=True)

# 4. Item 3 & Reviewer Point 7: GSE120878 vs GSE122909 Beta-Value Comparison Table for Top 10 Diagnostic CpGs
print("\n--- Top 10 Diagnostic CpGs Beta-Values Across GSE120878, GSE122909, and TCGA-SKCM ---", flush=True)
top_diag = pd.read_csv("top50_diagnostic_gse120878.tsv", sep="\t").head(10)
from run_diag_ml import load_geo_fast
X122, probes122, sids122, chars122 = load_geo_fast(os.path.join(dl, "GSE122909_series_matrix.txt.gz"))
is_cell_line = np.array(["cell line" in x.lower() for x in chars122[0]])
is_melanoma = np.array(["melanoma" in x.lower() for x in chars122[1]])
mel_tissue_mask = (~is_cell_line) & is_melanoma
melanocyte_mask = (~is_cell_line) & (~is_melanoma)
cell_line_mask = is_cell_line

p122_map = {p: i for i, p in enumerate(probes122)}

comp_rows = []
for _, r in top_diag.iterrows():
    pid = r["probe"]
    i122 = p122_map.get(pid, None)
    itcga = p_tcga_map.get(pid, None)
    b_mel_122 = np.nanmean(X122[i122, mel_tissue_mask]) if i122 is not None else np.nan
    b_mcyte_122 = np.nanmean(X122[i122, melanocyte_mask]) if i122 is not None else np.nan
    b_cline_122 = np.nanmean(X122[i122, cell_line_mask]) if i122 is not None else np.nan
    b_tcga = np.nanmean(X_raw[itcga, :]) if itcga is not None else np.nan
    comp_rows.append({
        "probe": pid, "Gene": r["Gene"],
        "GSE120878_Nevus": r["mean_nevus"],
        "GSE120878_Melanoma": r["mean_melanoma"],
        "TCGA_SKCM_Tumor": b_tcga,
        "GSE122909_TissueMel": b_mel_122,
        "GSE122909_Melanocyte": b_mcyte_122,
        "GSE122909_CellLine": b_cline_122
    })
comp_df = pd.DataFrame(comp_rows)
print(comp_df.to_string(index=False), flush=True)
comp_df.to_csv("top10_cross_cohort_betas.tsv", sep="\t", index=False)

hyper9 = comp_df[comp_df["probe"] != "cg27030541"]
print("\nMean across the 9 hypermethylated probes:")
print(hyper9[["GSE120878_Nevus", "GSE120878_Melanoma", "TCGA_SKCM_Tumor", "GSE122909_TissueMel", "GSE122909_Melanocyte", "GSE122909_CellLine"]].mean().to_string(), flush=True)

# 5. Check which promoter CpGs separate melanoma from BOTH nevi (in GSE120878) AND melanocytes (in GSE122909)!
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

concordant_mask = (p120 < 1e-10) & (p122 < 0.01) & (np.sign(d120) == np.sign(d122)) & (np.abs(d120) > 0.20) & (np.abs(d122) > 0.20)
conc_df = pd.DataFrame({
    "probe": common_p[concordant_mask],
    "Gene": [probe_to_gene.get(p, "") for p in common_p[concordant_mask]],
    "delta_120878": d120[concordant_mask], "p_120878": p120[concordant_mask],
    "delta_122909": d122[concordant_mask], "p_122909": p122[concordant_mask],
    "min_abs_t": np.minimum(np.abs(t120[concordant_mask]), np.abs(t122[concordant_mask]))
}).sort_values("min_abs_t", ascending=False)
print(f"\nPromoter CpGs separating melanoma from BOTH benign nevi (p<1e-10, |d|>0.2) AND cultured melanocytes (p<0.01, |d|>0.2): {len(conc_df)}", flush=True)
print(conc_df.head(15).to_string(index=False), flush=True)
conc_df.to_csv("concordant_melanoma_vs_nevus_and_melanocyte.tsv", sep="\t", index=False)

# 6. Complete Table 3 Re-generation (dx_cv.tsv) with IN-FOLD variance filtering AND external testing on TCGA-SKCM (overall, primary, metastatic) and GSE122909
print("\n--- Generating Complete Reproducible dx_cv.tsv for Table 3 ---", flush=True)
# Common promoter probes present across GSE120878, TCGA-SKCM (tested), and GSE122909
common_all3_mask = np.array([(p in p122_map) and (p in p_tcga_map) for p in prom_probes120])
probes_all3 = prom_probes120[common_all3_mask]
X120_all3 = X120_prom[common_all3_mask].T # (162, n_common)
X_tcga_all3 = X_imp[[p_tcga_map[p] for p in probes_all3], :].T # (456, n_common)
keep_15 = ~is_cell_line
y122_15 = is_melanoma[keep_15].astype(int)
X122_all3 = X122[[p122_map[p] for p in probes_all3]][:, keep_15].T # (15, n_common)
tr_means_all3 = X120_all3.mean(axis=0, keepdims=True)
X122_all3 = np.where(np.isnan(X122_all3), tr_means_all3, X122_all3)

rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=2, random_state=42)
fold_top200_prom = []
for train_idx, test_idx in rskf.split(X120_prom.T, y120):
    X_tr, y_tr = X120_prom.T[train_idx], y120[train_idx]
    t_fold, _ = stats.ttest_ind(X_tr[y_tr == 1], X_tr[y_tr == 0], axis=0, equal_var=False)
    top200 = np.argsort(np.abs(t_fold))[::-1][:200]
    fold_top200_prom.append((train_idx, test_idx, top200))

t_full_all3, _ = stats.ttest_ind(X120_all3[y120 == 1], X120_all3[y120 == 0], axis=0, equal_var=False)

dx_rows = []
for k_feat in [10, 50, 200]:
    for model_name, clf_fn in [
        ("L1 Logistic Regression", lambda: LogisticRegression(penalty="l1", solver="liblinear", C=1.0, random_state=42)),
        ("Random Forest", lambda: RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1))
    ]:
        aucs, accs = [], []
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

        # Fit on all GSE120878 for external testing on TCGA-SKCM and GSE122909
        top_ext_idx = np.argsort(np.abs(t_full_all3))[::-1][:k_feat]
        clf_ext = clf_fn()
        clf_ext.fit(X120_all3[:, top_ext_idx], y120)

        # TCGA-SKCM sensitivity
        preds_tcga = (clf_ext.predict_proba(X_tcga_all3[:, top_ext_idx])[:, 1] >= 0.5).astype(int)
        sens_all = np.mean(preds_tcga) * 100.0
        sens_prim = np.mean(preds_tcga[prim_mask]) * 100.0
        sens_met = np.mean(preds_tcga[met_mask]) * 100.0

        # GSE122909 test AUC
        probs_122 = clf_ext.predict_proba(X122_all3[:, top_ext_idx])[:, 1]
        auc_122 = roc_auc_score(y122_15, probs_122)
        mcyte_fp = np.mean(probs_122[y122_15 == 0] >= 0.5) * 100.0

        dx_rows.append({
            "Probe_Set": "Promoter-Only",
            "k_Features": k_feat,
            "Model": model_name,
            "CV_AUC_Mean": np.mean(aucs),
            "CV_AUC_SD": np.std(aucs),
            "CV_Acc_Mean": np.mean(accs),
            "CV_Acc_SD": np.std(accs),
            "TCGA_Sens_All": sens_all,
            "TCGA_Sens_Primary": sens_prim,
            "TCGA_Sens_Metastatic": sens_met,
            "GSE122909_AUC": auc_122,
            "GSE122909_Melanocyte_FP_Pct": mcyte_fp
        })

# Append the All-Probes (60k) rows from gse120878_cv_summary.tsv
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
        "TCGA_Sens_All": np.nan,
        "TCGA_Sens_Primary": np.nan,
        "TCGA_Sens_Metastatic": np.nan,
        "GSE122909_AUC": np.nan,
        "GSE122909_Melanocyte_FP_Pct": np.nan
    })

dx_df = pd.DataFrame(dx_rows)
print(dx_df.to_string(index=False), flush=True)
dx_df.to_csv("dx_cv.tsv", sep="\t", index=False)
dx_df.to_csv(os.path.join(dl, "dx_cv.tsv"), sep="\t", index=False)
