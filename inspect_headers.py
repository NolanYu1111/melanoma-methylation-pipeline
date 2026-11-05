import os
import gzip
import pandas as pd
import numpy as np

dl = r"C:\Users\ethan\Downloads"

# 1. Inspect GPL13534-11288.txt header
gpl_path = os.path.join(dl, "GPL13534-11288.txt")
with open(gpl_path, "r", encoding="utf-8", errors="ignore") as f:
    for i in range(50):
        line = f.readline()
        if line.startswith("ID\t") or line.startswith("IlmnID"):
            print(f"GPL header at line {i}:", line.strip().split("\t")[:25])
            break

# 2. Inspect GSE120878_series_matrix.txt.gz metadata
gse120_path = os.path.join(dl, "GSE120878_series_matrix.txt.gz")
with gzip.open(gse120_path, "rt", encoding="utf-8", errors="ignore") as f:
    for i in range(100):
        line = f.readline()
        if line.startswith("!Sample_title") or line.startswith("!Sample_Geo_accession") or line.startswith("!Sample_characteristics_ch1") or line.startswith("!Sample_source_name_ch1"):
            print("GSE120878:", line[:120].strip())
        if line.startswith("!series_matrix_table_begin"):
            header = f.readline().strip().split("\t")
            first_row = f.readline().strip().split("\t")
            print("GSE120878 table cols:", len(header), "first row:", first_row[:5])
            break

# 3. Inspect GSE122909_series_matrix.txt.gz metadata
gse122_path = os.path.join(dl, "GSE122909_series_matrix.txt.gz")
with gzip.open(gse122_path, "rt", encoding="utf-8", errors="ignore") as f:
    for i in range(100):
        line = f.readline()
        if line.startswith("!Sample_title") or line.startswith("!Sample_characteristics_ch1") or line.startswith("!Sample_source_name_ch1"):
            print("GSE122909:", line.strip())
        if line.startswith("!series_matrix_table_begin"):
            header = f.readline().strip().split("\t")
            first_row = f.readline().strip().split("\t")
            print("GSE122909 table cols:", len(header), "first row:", first_row[:5])
            break

# 4. Inspect TCGA-SKCM.methylation450.tsv.gz first 2 lines
tcga_meth_path = os.path.join(dl, "TCGA-SKCM.methylation450.tsv.gz")
with gzip.open(tcga_meth_path, "rt", encoding="utf-8", errors="ignore") as f:
    h = f.readline().strip().split("\t")
    r1 = f.readline().strip().split("\t")
    print("TCGA meth cols:", len(h), "first 5:", h[:5], "row1:", r1[:5])
