import os
import gzip
import time
import urllib.request
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

def clean_stage(s):
    if not isinstance(s, str):
        return "Unknown"
    s = s.strip()
    if "Stage IV" in s:
        return "IV"
    if "Stage III" in s:
        return "III"
    if "Stage II" in s:
        return "II"
    if "Stage I" in s:
        return "I"
    if "Stage 0" in s:
        return "0"
    return "Unknown"

# 1. Load GPL13534 promoter probes
print("Loading GPL13534 promoter manifest...", flush=True)
gpl_path = os.path.join(dl, "GPL13534-11288.txt")
gpl = pd.read_csv(gpl_path, sep="\t", skiprows=37, usecols=["ID", "UCSC_RefGene_Name", "UCSC_RefGene_Group"], low_memory=False).dropna(subset=["ID"])
def is_promoter(group_str):
    if not isinstance(group_str, str):
        return False
    parts = set(group_str.split(";"))
    return ("TSS200" in parts) or ("TSS1500" in parts)

gpl["is_promoter"] = gpl["UCSC_RefGene_Group"].apply(is_promoter)
promoter_probes = set(gpl[gpl["is_promoter"]]["ID"])
probe_to_gene = dict(zip(gpl["ID"], gpl["UCSC_RefGene_Name"].fillna("").apply(lambda s: s.split(";")[0] if s else "Unannotated")))

immune_markers = ["PTPRC", "CD8A", "CD8B", "CD3D", "CD3E", "GZMA", "PRF1", "IFNG", "CXCL9", "CXCL10", "STAT1", "PSMB9", "TAP1", "HLA-DRA"]

def download_if_needed(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        return
    print(f"Downloading {url} -> {dest} ...", flush=True)
    urllib.request.urlretrieve(url, dest)

def run_cohort_pipeline(cancer_code, delete_raw_meth=True):
    print("\n" + "="*80)
    print(f"RUNNING PIPELINE FOR TCGA-{cancer_code}")
    print("="*80, flush=True)
    
    surv_path = os.path.join(dl, f"TCGA-{cancer_code}.survival.tsv.gz")
    clin_path = os.path.join(dl, f"TCGA-{cancer_code}.clinical.tsv.gz")
    rna_path  = os.path.join(dl, f"TCGA-{cancer_code}.HiSeqV2.tsv.gz")
    meth_path = os.path.join(dl, f"TCGA-{cancer_code}.methylation450.tsv.gz")
    npy_path  = f"tcga_{cancer_code.lower()}_promoter_X_raw.npy"
    prb_path  = f"tcga_{cancer_code.lower()}_promoter_probes.csv"
    meta_path = f"tcga_{cancer_code.lower()}_meta.csv"

    if cancer_code == "SKCM":
        X_raw = np.load("tcga_promoter_X_raw.npy")
        probes_df = pd.read_csv("tcga_promoter_probes.csv")
        meta = pd.read_csv("tcga_456_meta.csv")
    elif os.path.exists(npy_path) and os.path.exists(prb_path) and os.path.exists(meta_path):
        X_raw = np.load(npy_path)
        probes_df = pd.read_csv(prb_path)
        meta = pd.read_csv(meta_path)
    else:
        download_if_needed(f"https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-{cancer_code}.survival.tsv.gz", surv_path)
        download_if_needed(f"https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-{cancer_code}.clinical.tsv.gz", clin_path)
        download_if_needed(f"https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.{cancer_code}.sampleMap%2FHiSeqV2.gz", rna_path)
        download_if_needed(f"https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-{cancer_code}.methylation450.tsv.gz", meth_path)

        tcga_surv = pd.read_csv(surv_path, sep="\t")
        tcga_clin = pd.read_csv(clin_path, sep="\t", low_memory=False)
        with gzip.open(meth_path, "rt", encoding="utf-8") as f:
            header_cols = f.readline().strip().split("\t")
            meth_samples = header_cols[1:]

        meta = pd.DataFrame({"sample": meth_samples, "col_idx": range(1, len(header_cols))})
        meta["sample_type_code"] = meta["sample"].apply(lambda s: s.split("-")[3] if len(s.split("-")) > 3 else "")
        meta["patient"] = meta["sample"].apply(lambda s: "-".join(s.split("-")[:3]))
        # Retain primary and metastatic solid tumor samples (01A/01B/06A/06B/07A), exclude 11A normal
        meta = meta[meta["sample_type_code"].str.startswith(("01", "06", "07"))].copy()
        meta["is_primary"] = meta["sample_type_code"].str.startswith("01").astype(int)
        meta = meta.merge(tcga_surv[["sample", "OS.time", "OS"]], on="sample", how="inner")
        meta = meta[meta["OS.time"].notna() & (meta["OS.time"] > 0) & meta["OS"].notna()].copy()
        
        # Age, sex, stage columns in GDC clinical file
        age_col = "age_at_index.demographic" if "age_at_index.demographic" in tcga_clin.columns else [c for c in tcga_clin.columns if "age_at_index" in c][0]
        sex_col = "gender.demographic" if "gender.demographic" in tcga_clin.columns else [c for c in tcga_clin.columns if "gender" in c][0]
        stage_cols = [c for c in tcga_clin.columns if "ajcc_pathologic_stage" in c]
        stage_col = stage_cols[0] if stage_cols else None

        keep_c = ["sample", age_col, sex_col] + ([stage_col] if stage_col else [])
        meta = meta.merge(tcga_clin[keep_c], on="sample", how="left")
        meta = meta.rename(columns={age_col: "age_at_index.demographic", sex_col: "gender.demographic"})
        if stage_col:
            meta = meta.rename(columns={stage_col: "ajcc_pathologic_stage.diagnoses"})
        else:
            meta["ajcc_pathologic_stage.diagnoses"] = "Unknown"

        meta = meta[meta["age_at_index.demographic"].notna()].copy()
        meta = meta.sort_values(by=["patient", "is_primary", "sample"], ascending=[True, False, True])
        meta = meta.drop_duplicates(subset=["patient"], keep="first").copy().reset_index(drop=True)
        meta["stage_clean"] = meta["ajcc_pathologic_stage.diagnoses"].apply(clean_stage)

        col_indices = meta["col_idx"].values
        n_samples = len(col_indices)
        t0 = time.time()
        probe_ids = []
        matrix_rows = []
        with gzip.open(meth_path, "rt", encoding="utf-8") as f:
            _ = f.readline()
            for line in f:
                tab_pos = line.find("\t")
                pid = line[:tab_pos]
                if pid not in promoter_probes:
                    continue
                parts = line.rstrip("\r\n").split("\t")
                vals = [float(parts[i]) if (i < len(parts) and parts[i] not in ("", "NA", "nan")) else np.nan for i in col_indices]
                arr = np.array(vals, dtype=np.float64)
                if np.isnan(arr).sum() <= 0.20 * n_samples:
                    if np.nanstd(arr, ddof=1) > 0.05:
                        probe_ids.append(pid)
                        matrix_rows.append(arr)

        print(f"TCGA-{cancer_code}: Extracted {len(probe_ids)} promoter CpGs across {n_samples} patients in {time.time()-t0:.1f}s", flush=True)
        X_raw = np.vstack(matrix_rows)
        probes_df = pd.DataFrame({"probe": probe_ids})
        np.save(npy_path, X_raw)
        probes_df.to_csv(prb_path, index=False)
        meta.to_csv(meta_path, index=False)

        if delete_raw_meth and cancer_code != "SKCM" and os.path.exists(meth_path):
            os.remove(meth_path)
            print(f"Deleted raw {meth_path} to conserve disk space.", flush=True)

    # Impute & standardize
    row_means = np.nanmean(X_raw, axis=1, keepdims=True)
    X_imp = np.where(np.isnan(X_raw), row_means, X_raw)
    row_sds = np.std(X_imp, axis=1, ddof=1, keepdims=True)
    X_std = (X_imp - row_means) / row_sds

    # Global mean promoter beta per sample
    global_beta = np.mean(X_imp, axis=0)

    # Top 5 PCs
    pca = PCA(n_components=5, random_state=42)
    pcs = pca.fit_transform(X_std.T)
    evr = pca.explained_variance_ratio_

    # Load RNA-seq & match by 15-character sample barcode
    rnaseq_df = pd.read_csv(rna_path, sep="\t", index_col=0)
    meta["sample_15"] = meta["sample"].str.slice(0, 15)
    has_rna = meta["sample_15"].isin(rnaseq_df.columns).values
    rna_samples = meta.loc[has_rna, "sample_15"].values
    rna_sub = rnaseq_df[rna_samples]

    avail_markers = [g for g in immune_markers if g in rna_sub.index]
    imm_expr = rna_sub.loc[avail_markers].values.astype(float)
    imm_z = (imm_expr - imm_expr.mean(axis=1, keepdims=True)) / imm_expr.std(axis=1, ddof=1, keepdims=True)
    imm_score = imm_z.mean(axis=0)
    imm_score_std = (imm_score - imm_score.mean()) / imm_score.std(ddof=1)

    # Correlate PCs with global_beta and imm_score_std
    r_pc1_glob, _ = stats.pearsonr(pcs[:, 0], global_beta)
    r_pc2_glob, _ = stats.pearsonr(pcs[:, 1], global_beta)
    pcs_rna = pcs[has_rna]
    pc_imm_corrs = [stats.pearsonr(pcs_rna[:, k], imm_score_std)[0] for k in range(5)]
    pc_imm_pvals = [stats.pearsonr(pcs_rna[:, k], imm_score_std)[1] for k in range(5)]
    best_pc_idx = int(np.argmax(np.abs(pc_imm_corrs)))

    # Multiple R of PC1-3 and PC1-5 predicting RNA-seq immune score
    X_ols3 = np.column_stack([np.ones(len(pcs_rna)), pcs_rna[:, :3]])
    pred3 = X_ols3 @ np.linalg.lstsq(X_ols3, imm_score_std, rcond=None)[0]
    mult_r3 = abs(stats.pearsonr(pred3, imm_score_std)[0])

    X_ols5 = np.column_stack([np.ones(len(pcs_rna)), pcs_rna[:, :5]])
    pred5 = X_ols5 @ np.linalg.lstsq(X_ols5, imm_score_std, rcond=None)[0]
    mult_r5 = abs(stats.pearsonr(pred5, imm_score_std)[0])

    # Build clinical covariate matrix Z_clin
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
    order = np.argsort(time_arr)
    t_s = time_arr[order]
    d_s = event_arr[order]
    Z_clin_s = Z_clin[order]
    pcs_s = pcs[order]
    X_std_s = X_std[:, order]

    # 1. Clinical-only Cox screen
    z_clin, p_clin, hr_clin, _ = fit_breslow_null_and_score_test(Z_clin_s, X_std_s, t_s, d_s)
    q_clin = bh_fdr(p_clin)

    # 2. Clinical + 3 PCs Cox screen
    z_pc3, p_pc3, hr_pc3, _ = fit_breslow_null_and_score_test(np.column_stack([Z_clin_s, pcs_s[:, :3]]), X_std_s, t_s, d_s)
    q_pc3 = bh_fdr(p_pc3)

    # 3. Clinical + 5 PCs Cox screen
    z_pc5, p_pc5, hr_pc5, _ = fit_breslow_null_and_score_test(np.column_stack([Z_clin_s, pcs_s[:, :5]]), X_std_s, t_s, d_s)
    q_pc5 = bh_fdr(p_pc5)

    # 4. Clinical + RNA-seq Immune Score Cox screen (on matched RNA-seq patients)
    meta_rna = meta.loc[has_rna].reset_index(drop=True)
    t_r = meta_rna["OS.time"].values.astype(float)
    d_r = meta_rna["OS"].values.astype(int)
    ord_r = np.argsort(t_r)
    Z_clin_r = Z_clin[has_rna][ord_r]
    imm_r = imm_score_std[ord_r]
    X_std_r = X_std[:, has_rna][:, ord_r]

    # Test immune score itself on OS
    z_imm_os, p_imm_os, hr_imm_os, _ = fit_breslow_null_and_score_test(Z_clin_r, imm_r[None, :], t_r[ord_r], d_r[ord_r])
    # Clinical-only on exact RNA-matched subset
    _, p_clin_r, _, _ = fit_breslow_null_and_score_test(Z_clin_r, X_std_r, t_r[ord_r], d_r[ord_r])
    q_clin_r = bh_fdr(p_clin_r)
    # Clinical + RNA immune on exact RNA-matched subset
    z_imm, p_imm, hr_imm, _ = fit_breslow_null_and_score_test(np.column_stack([Z_clin_r, imm_r]), X_std_r, t_r[ord_r], d_r[ord_r])
    q_imm = bh_fdr(p_imm)

    # Also check CRABP2 (cg02363950) methylation-expression correlation in this cancer
    crabp2_r, crabp2_p = np.nan, np.nan
    if "cg02363950" in probes_df["probe"].values and "CRABP2" in rna_sub.index:
        cidx = probes_df[probes_df["probe"] == "cg02363950"].index[0]
        m_c = X_std[cidx, has_rna]
        e_c = rna_sub.loc["CRABP2"].values.astype(float)
        crabp2_r, crabp2_p = stats.pearsonr(m_c, e_c)

    n_clin_q05 = int((q_clin < 0.05).sum())
    n_clin_p05 = int((p_clin < 0.05).sum())
    n_clin_r_q05 = int((q_clin_r < 0.05).sum())
    n_pc3_q05  = int((q_pc3 < 0.05).sum())
    n_pc5_q05  = int((q_pc5 < 0.05).sum())
    n_imm_q05  = int((q_imm < 0.05).sum())
    n_imm_p05  = int((p_imm < 0.05).sum())

    ret_pc3 = int(((q_clin < 0.05) & (q_pc3 < 0.05)).sum())
    ret_pc5 = int(((q_clin < 0.05) & (q_pc5 < 0.05)).sum())
    ret_imm = int(((q_clin_r < 0.05) & (q_imm < 0.05)).sum())

    pct_att_pc5 = 100.0 * (1.0 - ret_pc5 / n_clin_q05) if n_clin_q05 > 0 else 0.0
    pct_att_imm = 100.0 * (1.0 - ret_imm / n_clin_r_q05) if n_clin_r_q05 > 0 else 0.0

    n_prim = int((meta["is_primary"] == 1).sum())
    n_met  = int((meta["is_primary"] == 0).sum())

    print(f"TCGA-{cancer_code} Summary:")
    print(f"  Patients: n = {len(meta)} ({n_prim} primary, {n_met} metastatic; {int(event_arr.sum())} deaths; { int(has_rna.sum()) } with matched RNA-seq)")
    print(f"  Tested Promoter CpGs: {len(probes_df)}")
    print(f"  PC1 ({100*evr[0]:.1f}% var): r(Global Beta) = {r_pc1_glob:+.3f}, r(Immune) = {pc_imm_corrs[0]:+.3f}")
    print(f"  PC2 ({100*evr[1]:.1f}% var): r(Global Beta) = {r_pc2_glob:+.3f}, r(Immune) = {pc_imm_corrs[1]:+.3f}")
    print(f"  Strongest Immune PC: PC{best_pc_idx+1} ({100*evr[best_pc_idx]:.1f}% var), r = {pc_imm_corrs[best_pc_idx]:+.3f} (p = {pc_imm_pvals[best_pc_idx]:.2e}) | PC1-3 Multiple R = {mult_r3:.3f}, PC1-5 R = {mult_r5:.3f}")
    print(f"  Immune Score -> OS (clinical-adj): HR = {hr_imm_os[0]:.3f}, p = {p_imm_os[0]:.3e}")
    print(f"  Clinical-Only Hits: p < 0.05 = {n_clin_p05}, q < 0.05 = {n_clin_q05} (RNA subset q < 0.05 = {n_clin_r_q05})")
    print(f"  + 3 PCs Hits (q < 0.05):       {n_pc3_q05} (retained {ret_pc3}/{n_clin_q05})")
    print(f"  + 5 PCs Hits (q < 0.05):       {n_pc5_q05} (retained {ret_pc5}/{n_clin_q05}; {pct_att_pc5:.1f}% clin hits attenuated)")
    print(f"  + RNA Immune Hits (q < 0.05):  {n_imm_q05} (retained {ret_imm}/{n_clin_r_q05}; {pct_att_imm:.1f}% clin hits attenuated)")
    print(f"  CRABP2 (cg02363950) vs CRABP2 mRNA: r = {crabp2_r:+.3f} (p = {crabp2_p:.2e})", flush=True)

    return {
        "Cohort": f"TCGA-{cancer_code}",
        "N_Patients": len(meta),
        "N_Primary": n_prim,
        "N_Metastatic": n_met,
        "N_Deaths": int(event_arr.sum()),
        "N_RNA_Matched": int(has_rna.sum()),
        "N_Promoter_CpGs": len(probes_df),
        "PC1_Var_Pct": 100 * evr[0],
        "PC1_r_GlobalBeta": r_pc1_glob,
        "PC1_r_Immune": pc_imm_corrs[0],
        "PC2_Var_Pct": 100 * evr[1],
        "PC2_r_GlobalBeta": r_pc2_glob,
        "PC2_r_Immune": pc_imm_corrs[1],
        "Best_Immune_PC": f"PC{best_pc_idx+1}",
        "Best_Immune_PC_r": pc_imm_corrs[best_pc_idx],
        "Best_Immune_PC_p": pc_imm_pvals[best_pc_idx],
        "PC1_3_Multiple_R_Immune": mult_r3,
        "PC1_5_Multiple_R_Immune": mult_r5,
        "Immune_OS_HR": hr_imm_os[0],
        "Immune_OS_p": p_imm_os[0],
        "Clin_p05": n_clin_p05,
        "Clin_q05": n_clin_q05,
        "Clin_RNA_q05": n_clin_r_q05,
        "PC3_p05": int((p_pc3 < 0.05).sum()),
        "PC3_q05": n_pc3_q05,
        "Retained_PC3_q05": ret_pc3,
        "PC5_p05": int((p_pc5 < 0.05).sum()),
        "PC5_q05": n_pc5_q05,
        "Retained_PC5_q05": ret_pc5,
        "Pct_Clin_Att_PC5": pct_att_pc5,
        "RNA_Imm_p05": n_imm_p05,
        "RNA_Imm_q05": n_imm_q05,
        "Retained_Imm_q05": ret_imm,
        "Pct_Clin_Att_Imm": pct_att_imm,
        "CRABP2_meth_expr_r": crabp2_r,
        "CRABP2_meth_expr_p": crabp2_p,
    }

rows = []
for code in ["SKCM", "HNSC", "LUAD", "BLCA", "LIHC"]:
    rows.append(run_cohort_pipeline(code, delete_raw_meth=True))

cross_df = pd.DataFrame(rows)
print("\n" + "="*80)
print("CROSS-CANCER SUMMARY TABLE (5 COHORTS)")
print("="*80)
print(cross_df.to_string(index=False), flush=True)
cross_df.to_csv("cross_cancer_summary.tsv", sep="\t", index=False)

