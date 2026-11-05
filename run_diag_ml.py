import os
import gzip
import time
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score, accuracy_score

dl = r"C:\Users\ethan\Downloads"

def load_geo_fast(gz_path):
    t0 = time.time()
    sample_ids = []
    characteristics = []
    probes = []
    rows = []
    with gzip.open(gz_path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("!Sample_geo_accession"):
                sample_ids = [x.strip('"') for x in line.strip().split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                characteristics.append([x.strip('"') for x in line.strip().split("\t")[1:]])
            elif line.startswith("!series_matrix_table_begin"):
                break
        header = f.readline() # skip header
        n_cols = len(sample_ids)
        for line in f:
            if line.startswith("!series_matrix_table_end") or not line.strip():
                break
            tab_pos = line.find("\t")
            pid = line[:tab_pos].strip('"')
            if not (pid.startswith("cg") or pid.startswith("ch")):
                continue
            parts = line[tab_pos+1:].rstrip("\r\n").split("\t")
            if len(parts) == n_cols:
                arr = np.array([float(x) if x not in ("", "NA", "null", "NULL", "NaN", "nan", '""') else np.nan for x in parts], dtype=np.float32)
                probes.append(pid)
                rows.append(arr)
    X = np.vstack(rows)
    print(f"Loaded {gz_path} ({X.shape}) in {time.time()-t0:.1f}s", flush=True)
    return X, np.array(probes), sample_ids, characteristics

if __name__ == "__main__":
    print("1. Loading GPL13534 promoter map...", flush=True)
    gpl = pd.read_csv(os.path.join(dl, "GPL13534-11288.txt"), sep="\t", skiprows=37, usecols=["ID", "UCSC_RefGene_Name", "UCSC_RefGene_Group"], low_memory=False).dropna(subset=["ID"])
    def is_promoter(group_str):
        if not isinstance(group_str, str):
            return False
        parts = set(group_str.split(";"))
        return ("TSS200" in parts) or ("TSS1500" in parts)
    gpl["is_promoter"] = gpl["UCSC_RefGene_Group"].apply(is_promoter)
    promoter_probes = set(gpl[gpl["is_promoter"]]["ID"])
    probe_to_gene = dict(zip(gpl["ID"], gpl["UCSC_RefGene_Name"].fillna("").apply(lambda s: s.split(";")[0] if s else "Unannotated")))

    print("2. Loading GSE120878...", flush=True)
    X120, probes120, sids120, chars120 = load_geo_fast(os.path.join(dl, "GSE120878_series_matrix.txt.gz"))
    disease_states = [c for c in chars120 if any("disease state" in x.lower() or "tissue" in x.lower() for x in c)][0]
    y120 = np.array([1 if "melanoma" in x.lower() else 0 for x in disease_states])
    print(f"GSE120878 class counts: {np.sum(y120==1)} melanomas, {np.sum(y120==0)} nevi", flush=True)

    # Remove probes with >= 5% missing values and keep top 60,000 highest variance probes
    miss_frac = np.isnan(X120).mean(axis=1)
    keep_mask = miss_frac < 0.05
    X120_clean = X120[keep_mask].copy()
    probes120_clean = probes120[keep_mask]

    row_means = np.nanmean(X120_clean, axis=1, keepdims=True)
    X120_clean = np.where(np.isnan(X120_clean), row_means, X120_clean)

    variances = np.var(X120_clean, axis=1, ddof=1)
    top60k_idx = np.argsort(variances)[::-1][:60000]
    X120_60k = X120_clean[top60k_idx]
    probes120_60k = probes120_clean[top60k_idx]

    prom_mask = np.array([p in promoter_probes for p in probes120_60k])
    X120_prom = X120_60k[prom_mask]
    probes120_prom = probes120_60k[prom_mask]
    print(f"Top 60k probes: {len(probes120_60k)}, promoter subset: {len(probes120_prom)}", flush=True)

    # Save GSE120878 promoter subset for fast figure plotting
    np.save("gse120878_prom_X.npy", X120_prom)
    np.save("gse120878_y.npy", y120)
    pd.DataFrame({"probe": probes120_prom}).to_csv("gse120878_prom_probes.csv", index=False)

    # Welch t-statistic across all GSE120878 promoter probes
    X_prom_T = X120_prom.T # (162, n_prom)
    mel_mask = (y120 == 1)
    nev_mask = (y120 == 0)
    t_stat, p_val = stats.ttest_ind(X_prom_T[mel_mask], X_prom_T[nev_mask], axis=0, equal_var=False)
    delta_beta = X_prom_T[mel_mask].mean(axis=0) - X_prom_T[nev_mask].mean(axis=0)

    top_diag_df = pd.DataFrame({
        "probe": probes120_prom,
        "Gene": [probe_to_gene.get(p, "") for p in probes120_prom],
        "t_stat": t_stat,
        "abs_t": np.abs(t_stat),
        "p_val": p_val,
        "mean_melanoma": X_prom_T[mel_mask].mean(axis=0),
        "mean_nevus": X_prom_T[nev_mask].mean(axis=0),
        "delta_beta": delta_beta
    }).sort_values("abs_t", ascending=False)

    print("\nTop 15 Diagnostic Promoter CpGs in GSE120878:", flush=True)
    print(top_diag_df.head(15).to_string(index=False), flush=True)
    top_diag_df.head(50).to_csv("top50_diagnostic_gse120878.tsv", sep="\t", index=False)

    # 5-fold CV repeated twice
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=2, random_state=42)
    cv_results = []

    for probe_set_name, X_sub_mat in [("Promoter-Only", X120_prom.T), ("All-Probes (60k)", X120_60k.T)]:
        # Precompute fold t-stats once per fold!
        fold_top200 = []
        for train_idx, test_idx in rskf.split(X_sub_mat, y120):
            X_tr, y_tr = X_sub_mat[train_idx], y120[train_idx]
            t_fold, _ = stats.ttest_ind(X_tr[y_tr == 1], X_tr[y_tr == 0], axis=0, equal_var=False)
            top200 = np.argsort(np.abs(t_fold))[::-1][:200]
            fold_top200.append((train_idx, test_idx, top200))

        for k_feat in [10, 50, 200]:
            for model_name, clf_fn in [
                ("L1 Logistic Regression", lambda: LogisticRegression(penalty="l1", solver="liblinear", C=1.0, random_state=42)),
                ("Random Forest", lambda: RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1))
            ]:
                aucs, accs = [], []
                for train_idx, test_idx, top200 in fold_top200:
                    top_idx = top200[:k_feat]
                    X_tr, y_tr = X_sub_mat[train_idx][:, top_idx], y120[train_idx]
                    X_te, y_te = X_sub_mat[test_idx][:, top_idx], y120[test_idx]
                    clf = clf_fn()
                    clf.fit(X_tr, y_tr)
                    probs = clf.predict_proba(X_te)[:, 1]
                    preds = (probs >= 0.5).astype(int)
                    aucs.append(roc_auc_score(y_te, probs))
                    accs.append(accuracy_score(y_te, preds) * 100)
                cv_results.append({
                    "Probe_Set": probe_set_name,
                    "k_Features": k_feat,
                    "Model": model_name,
                    "Mean_AUC": np.mean(aucs),
                    "SD_AUC": np.std(aucs),
                    "Mean_Acc": np.mean(accs),
                    "SD_Acc": np.std(accs)
                })

    cv_df = pd.DataFrame(cv_results)
    print("\nGSE120878 Cross-Validation Summary:", flush=True)
    print(cv_df.to_string(index=False), flush=True)
    cv_df.to_csv("gse120878_cv_summary.tsv", sep="\t", index=False)

    # 3. External Testing on GSE122909
    print("\n3. Loading GSE122909...", flush=True)
    X122, probes122, sids122, chars122 = load_geo_fast(os.path.join(dl, "GSE122909_series_matrix.txt.gz"))
    is_cell_line = np.array(["cell line" in x.lower() for x in chars122[0]])
    is_melanoma = np.array(["melanoma" in x.lower() for x in chars122[1]])
    keep_15 = ~is_cell_line
    y122_15 = is_melanoma[keep_15].astype(int)
    print(f"GSE122909 testing subset (n={keep_15.sum()}): {y122_15.sum()} tissue melanomas, {(1-y122_15).sum()} melanocytes", flush=True)

    p122_map = {p: i for i, p in enumerate(probes122)}
    common_mask = np.array([p in p122_map for p in probes120_prom])
    X_tr_all = X120_prom[common_mask].T
    idx_in_122 = [p122_map[p] for p in probes120_prom[common_mask]]
    X_te_122 = X122[idx_in_122][:, keep_15].T
    # Impute any NaNs in X_te_122 with training mean
    tr_means = X_tr_all.mean(axis=0, keepdims=True)
    X_te_122 = np.where(np.isnan(X_te_122), tr_means, X_te_122)

    t_all, _ = stats.ttest_ind(X_tr_all[y120 == 1], X_tr_all[y120 == 0], axis=0, equal_var=False)
    for k_feat in [10, 50, 200]:
        top_idx = np.argsort(np.abs(t_all))[::-1][:k_feat]
        for model_name, clf in [
            ("L1 Logistic Regression", LogisticRegression(penalty="l1", solver="liblinear", C=1.0, random_state=42)),
            ("Random Forest", RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1))
        ]:
            clf.fit(X_tr_all[:, top_idx], y120)
            p122 = clf.predict_proba(X_te_122[:, top_idx])[:, 1]
            auc122 = roc_auc_score(y122_15, p122)
            melanocyte_called_mel = np.mean(p122[y122_15 == 0] >= 0.5) * 100
            print(f"GSE122909 (k={k_feat}, {model_name}): AUC = {auc122:.3f}, melanocytes called melanoma = {melanocyte_called_mel:.1f}%", flush=True)
