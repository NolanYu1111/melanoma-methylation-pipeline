import os
import gzip
import pandas as pd
import numpy as np

dl = r"C:\Users\ethan\Downloads"

print("1. Checking survival_results.xlsx...")
surv_xlsx = pd.read_excel(os.path.join(dl, "survival_results.xlsx"))
print("survival_results.xlsx shape:", surv_xlsx.shape)
print("survival_results.xlsx columns:", surv_xlsx.columns.tolist())
print(surv_xlsx.head(25))

print("\n2. Checking TCGA-SKCM.survival.tsv.gz...")
tcga_surv = pd.read_csv(os.path.join(dl, "TCGA-SKCM.survival.tsv.gz"), sep="\t")
print("tcga_surv shape:", tcga_surv.shape)
print("tcga_surv columns:", tcga_surv.columns.tolist())
print(tcga_surv.head(5))

print("\n3. Checking TCGA-SKCM.clinical.tsv.gz...")
tcga_clin = pd.read_csv(os.path.join(dl, "TCGA-SKCM.clinical.tsv.gz"), sep="\t", low_memory=False)
print("tcga_clin shape:", tcga_clin.shape)
stage_cols = [c for c in tcga_clin.columns if "stage" in c.lower() or "age" in c.lower() or "sex" in c.lower() or "gender" in c.lower() or "sample" in c.lower()]
print("Relevant clinical cols:", stage_cols)
