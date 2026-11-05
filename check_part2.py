import os
import re
import gzip
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, accuracy_score

dl = r"C:\Users\ethan\Downloads"

# 1. Table 4 Welch t-test degrees of freedom check
print("=== 1. TABLE 4 WELCH T-TEST DEGREES OF FREEDOM CHECK ===")
X_120 = np.load("gse120878_prom_X.npy") # shape (13243, 162)
y_120 = np.load("gse120878_y.npy")
p_120 = pd.read_csv("gse120878_prom_probes.csv")["probe"].values
top10_df = pd.read_csv("top50_diagnostic_gse120878.tsv", sep="\t").head(10)
for _, row in top10_df.iterrows():
    pid = row["probe"]
    idx = np.where(p_120 == pid)[0][0]
    x1 = X_120[idx, y_120 == 1] # melanoma (n=89)
    x0 = X_120[idx, y_120 == 0] # nevus (n=73)
    s1, s0 = np.var(x1, ddof=1), np.var(x0, ddof=1)
    n1, n0 = len(x1), len(x0)
    se = np.sqrt(s1/n1 + s0/n0)
    t_stat = (np.mean(x1) - np.mean(x0)) / se
    df_w = (s1/n1 + s0/n0)**2 / ((s1/n1)**2/(n1-1) + (s0/n0)**2/(n0-1))
    p_w = stats.t.sf(np.abs(t_stat), df=df_w) * 2
    p_stud = stats.ttest_ind(x1, x0, equal_var=True).pvalue
    print(f"{pid} ({row['Gene']:8s}): t={t_stat:+6.2f}, Welch df={df_w:5.1f}, Welch p={p_w:.2e}, Student p(df=160)={p_stud:.2e}, SD_mel={np.sqrt(s1):.3f}, SD_nev={np.sqrt(s0):.3f}")

# 2. Extract Conway et al. (2019) 40 CpG probes from PMC6531328 BioC JSON
print("\n=== 2. CONWAY ET AL. (2019) 40-CpG PANEL EXTRACTION & HEAD-TO-HEAD EVALUATION ===")
bioc_path = r"C:\Users\ethan\.gemini\antigravity\brain\28ec116b-b931-4f66-87c4-029fabff2fd8\.system_generated\steps\505\content.md"
with open(bioc_path, "r", encoding="utf-8") as f:
    bioc_text = f.read()

cgs_in_paper = []
for m in re.finditer(r"cg\d{8}", bioc_text):
    cg = m.group(0)
    if cg not in cgs_in_paper:
        cgs_in_paper.append(cg)
print(f"Found {len(cgs_in_paper)} unique cg IDs in Conway et al. (2019) PMC6531328: {cgs_in_paper}")

# Also search for sensitivity, specificity, AUC sentences in Conway et al.
for m in re.finditer(r"[^.]*(?:sensitivity|specificity|AUC|area under|classifier|40 CpG)[^.]*\.", bioc_text, re.IGNORECASE):
    s = m.group(0).strip()
    if any(k in s.lower() for k in ["96", "100", "0.99", "40", "validation", "elastic"]):
        if len(s) < 350:
            print("  [Conway text]:", s)

# 3. Check GSE122909 sample metadata & 95% CI on ROC-AUC + FAM174B expected by chance
print("\n=== 3. GSE122909 METADATA, BOOTSTRAP 95% CI ON AUC, AND FAM174B CHANCE EXPECTATION ===")
gse122909_path = os.path.join(dl, "GSE122909_series_matrix.txt.gz")
with gzip.open(gse122909_path, "rt", encoding="utf-8") as f:
    for line in f:
        if line.startswith("!Sample_title") or line.startswith("!Sample_source_name_ch1") or line.startswith("!Sample_characteristics_ch1"):
            print(line.strip()[:250])
        if line.startswith("!series_matrix_table_begin"):
            break

# Check how many common promoter probes exist between GSE120878 and GSE122909 and expected number passing |dBeta|>0.20, p<0.01 in both
gpl = pd.read_csv(os.path.join(dl, "GPL13534-11288.txt"), sep="\t", skiprows=37, usecols=["ID", "UCSC_RefGene_Group"], low_memory=False).dropna(subset=["ID"])
def is_prom(s):
    if not isinstance(s, str): return False
    parts = set(s.split(";"))
    return ("TSS200" in parts) or ("TSS1500" in parts)
prom_set = set(gpl[gpl["UCSC_RefGene_Group"].apply(is_prom)]["ID"])

# Parse GSE122909 promoter probes
p122_list = []
X122_rows = []
with gzip.open(gse122909_path, "rt", encoding="utf-8") as f:
    in_table = False
    for line in f:
        if line.startswith("!Sample_title"):
            titles = [x.strip('"') for x in line.strip().split("\t")[1:]]
        if line.startswith("!series_matrix_table_begin"):
            in_table = True
            header = f.readline()
            continue
        if in_table:
            if line.startswith("!series_matrix_table_end"):
                break
            parts = line.rstrip("\r\n").split("\t")
            pid = parts[0].strip('"')
            if pid in prom_set:
                vals = [float(x) if x not in ("", "NA", "null") else np.nan for x in parts[1:]]
                p122_list.append(pid)
                X122_rows.append(vals)

X122 = np.array(X122_rows)
p122 = np.array(p122_list)
print("GSE122909 sample titles:", titles)
is_mel_tissue = np.array([t.startswith("MM") for t in titles])
is_melanocyte = np.array([t.startswith("MC") or t.startswith("NHEM") for t in titles])

# Match p_120 (13,243 variance-filtered promoter probes) or all promoter probes in GSE120878
# Wait: how many total promoter probes in GSE120878 and GSE122909?
# In run_stress_test_and_checks.py, how many common promoter probes were tested for concordant |dBeta|>0.20, p<0.01?
common_probes, idx120, idx122 = np.intersect1d(p_120, p122, return_indices=True)
print(f"Common promoter probes between GSE120878 filtered set and GSE122909: {len(common_probes)}")

X120_c = X_120[idx120]
X122_c = X122[idx122]

db120 = np.nanmean(X120_c[:, y_120 == 1], axis=1) - np.nanmean(X120_c[:, y_120 == 0], axis=1)
p_val120 = np.array([stats.ttest_ind(X120_c[i, y_120 == 1], X120_c[i, y_120 == 0], equal_var=False).pvalue for i in range(len(common_probes))])

db122 = np.nanmean(X122_c[:, is_mel_tissue], axis=1) - np.nanmean(X122_c[:, is_melanocyte], axis=1)
p_val122 = np.array([stats.ttest_ind(X122_c[i, is_mel_tissue], X122_c[i, is_melanocyte], equal_var=False, nan_policy="omit").pvalue for i in range(len(common_probes))])

pass120 = (np.abs(db120) > 0.20) & (p_val120 < 0.01)
pass122 = (np.abs(db122) > 0.20) & (p_val122 < 0.01)
concordant = pass120 & pass122 & (np.sign(db120) == np.sign(db122))
discordant = pass120 & pass122 & (np.sign(db120) != np.sign(db122))
print(f"Pass in GSE120878 (|dBeta|>0.20, p<0.01): {pass120.sum()} / {len(common_probes)} ({pass120.mean()*100:.2f}%)")
print(f"Pass in GSE122909 (|dBeta|>0.20, p<0.01): {pass122.sum()} / {len(common_probes)} ({pass122.mean()*100:.2f}%)")
print(f"Expected overlap by chance (if independent): {pass120.sum() * pass122.mean():.2f} total ({0.5 * pass120.sum() * pass122.mean():.2f} concordant)")
print(f"Observed concordant: {concordant.sum()}, discordant: {discordant.sum()}")
