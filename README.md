# Transcriptomic Immune Infiltration Confounds Bulk Promoter Methylation–Survival Screens in Cutaneous Melanoma but Not in Four Other Solid Cancers, with a Promoter Classifier Separating Melanoma from Benign Nevi

**Authors:** Nolan Yu, Valentina L. Kouznetsova, Santosh Kesari, and Igor F. Tsigelny

This repository contains the complete Python analysis scripts, figure generation code, and output tables (`.tsv`) for the study:

> **Transcriptomic Immune Infiltration Confounds Bulk Promoter Methylation–Survival Screens in Cutaneous Melanoma but Not in Four Other Solid Cancers, with a Promoter Classifier Separating Melanoma from Benign Nevi**

## Overview of Scripts

1. `extract_tcga_promoter.py` — Extracts and deduplicates TCGA-SKCM Illumina 450K promoter methylation data ($n = 456$ patients, 221 deaths; 67,102 promoter CpGs with $\le 20\%$ missingness and $\text{SD}(\beta) > 0.05$).
2. `run_cox_and_ml.py` — Genome-wide clinical covariate-adjusted and 5-PC composition-adjusted Cox proportional hazards screening (`cox_results.tsv`, `cox_pc5.tsv`, `top15_clin.tsv`) and five-fold cross-validated concordance (Harrell's C-index) evaluation in TCGA-SKCM.
3. `run_diag_ml.py` — Five-fold cross-validation repeated twice for L1 logistic regression and random forest melanoma-vs.-nevus classifiers in GSE120878 ($n = 162$: 89 primary melanomas, 73 benign nevi), plus external sensitivity in TCGA-SKCM and external stress testing in GSE122909 ($n = 15$: 7 tissue melanomas, 8 normal cultured melanocytes; `dx_cv.tsv`, `top50_diagnostic_gse120878.tsv`).
4. `run_stress_test_and_checks.py` — Principal component biological characterization, 3-PC and 5-PC attenuation checks, proportional-hazards checks, Conway-style 40-CpG ElasticNet evaluation (`gse120878_cv_summary.tsv`), and cross-cohort $\beta$-value profiling (`top10_cross_cohort_betas.tsv`, `concordant_melanoma_vs_nevus_and_melanocyte.tsv`).
5. `run_scientific_upgrades.py` — Integration of matched TCGA-SKCM HiSeqV2 RNA-seq data ($n = 454$), computation of the 14-gene transcriptomic immune infiltration score (plus 12-gene non-overlapping and 6-gene T-cell/cytolytic sensitivity scores), RNA-seq immune-adjusted Cox screening, and candidate promoter methylation vs. mRNA expression correlations.
6. `update_candidate_table.py` — Genome-wide manifest SNP filtering (`Probe_SNPs` and `Probe_SNPs_10`, $n = 57,423$ SNP-free promoter CpGs) and immune-adjusted partial mRNA correlation table generation (`candidate_meth_expr_corr.tsv`, including *DYSFIP1* `cg03750525`).
7. `run_skcm_deep_analyses.py` — Pre-specified 50/50 discovery/validation ($n = 227 / 227$) immune-cold vs. immune-hot subgroup Cox screens and $\text{CpG} \times \text{immune-score}$ interaction testing, genome-wide regulator enrichment across all 65,227 mRNA-mapped promoter CpGs (including 1,000-draw SD- and mean-$\beta$-matched null distributions), and *CRABP2* (`cg02363950`) mRNA survival and Sobel mediation models (`statsmodels` `PHReg`).
8. `run_biopsy_and_bootstrap.py` — Anatomic biopsy/resection site adjustment (`site_of_resection_or_biopsy.diagnoses`: lymph node, skin, soft tissue, visceral) and within-compartment PC2–immune correlations in TCGA-SKCM ($n = 456$), plus 1,000-resample nonparametric bootstrap 95% confidence intervals for *CRABP2* (`cg02363950`) mRNA mediation (`crabp2_mediation_results.tsv`).
9. `run_cross_cancer_generalization.py` — Cross-cancer generalization of promoter methylation PCA, 14-gene RNA-seq immune infiltration scoring, and clinical/PC/immune-adjusted Cox survival screening across five TCGA cohorts ($n = 2,198$ total): TCGA-SKCM ($n = 456$), TCGA-HNSC ($n = 527$, 224 deaths; 64,044 promoter CpGs), TCGA-LUAD ($n = 437$, 158 deaths; 62,060 promoter CpGs), TCGA-BLCA ($n = 407$, 180 deaths; 63,201 promoter CpGs), and TCGA-LIHC ($n = 371$, 133 deaths; 64,673 promoter CpGs; `cross_cancer_summary.tsv`).
10. `run_permutation_calibration.py` — Empirical FDR calibration across all five TCGA cohorts (`SKCM`, `HNSC`, `LUAD`, `BLCA`, `LIHC`) and four covariate specifications (`Clinical Only`, `Clinical + 3 PCs`, `Clinical + 5 PCs`, `Clinical + RNA Immune`) using 50 survival-time permutations per model with covariates and methylation fixed (`permutation_fdr_calibration.tsv`, reporting nominal $p < 0.05$ fraction, fraction of permutations with $\ge 1$ discovery at $q < 0.05$, and median/mean $q < 0.05$ false positives).
11. `generate_new_figures.py` — Generates all publication figures (`pipeline_flowchart.png`, `km_curve.png`, `ml_classifier_accuracy.png`, and `methylation_boxplots_combined.png`).
12. `reproduce_all_results.py` — Master script that executes the complete pipeline sequentially.

## Output Tables Included

- `cox_results.tsv` — Clinical covariate-adjusted Cox proportional hazards results across 67,102 promoter CpGs in TCGA-SKCM ($n = 456$).
- `cox_pc5.tsv` — Clinical + 5-PC composition-adjusted Cox results across 67,102 promoter CpGs in TCGA-SKCM ($n = 456$).
- `top15_clin.tsv` — Top 15 clinical-only survival hits with Breslow vs. Efron tie-handling comparison.
- `candidate_meth_expr_corr.tsv` — Summary table of leading candidate promoter CpGs (Table 2), including clinical, 5-PC, SNP-free, and RNA-seq immune-adjusted Cox statistics alongside unadjusted and immune-adjusted matched mRNA expression correlations.
- `crabp2_mediation_results.tsv` — *CRABP2* (`cg02363950`) and *CRABP2* mRNA survival, Sobel mediation, and 1,000-resample bootstrap 95% confidence interval results in TCGA-SKCM ($n = 454$).
- `regulator_enrichment_results.tsv` and `regulator_enrichment_matched_null.tsv` — Genome-wide enrichment of survival-associated promoter CpGs for inverse immune-adjusted mRNA correlations across 65,227 mapped probes and 1,000 SD/mean-$\beta$-matched null draws.
- `promoter_meth_expr_survival_table.tsv` — Full table of all 65,227 mRNA-mapped promoter CpGs in TCGA-SKCM with unadjusted, immune-adjusted, and 5-PC-adjusted mRNA correlations, plus discovery/validation subgroup and $\text{CpG} \times \text{immune}$ interaction statistics.
- `cross_cancer_summary.tsv` — Cross-cancer comparison across five TCGA cohorts: TCGA-SKCM ($n = 456$), TCGA-HNSC ($n = 527$), TCGA-LUAD ($n = 437$), TCGA-BLCA ($n = 407$), and TCGA-LIHC ($n = 371$) (Table 3).
- `permutation_fdr_calibration.tsv` — Empirical false-discovery calibration across 50 survival-time permutations per model and cohort across all five TCGA cohorts (`SKCM`, `HNSC`, `LUAD`, `BLCA`, `LIHC`), including `ge1_hit_perms` and `pct_perms_ge1_q05` (Table 3).
- `dx_cv.tsv` and `gse120878_cv_summary.tsv` — Five-fold cross-validation (repeated twice) ROC-AUC and accuracy across promoter-only and all-probe diagnostic classifiers in GSE120878, plus external TCGA-SKCM sensitivity and GSE122909 stress-test ROC-AUC (Table 4).
- `top50_diagnostic_gse120878.tsv` — Welch's $t$-test statistics, Satterthwaite degrees of freedom, and mean $\beta$-values for the top 50 differentially methylated promoter CpGs in GSE120878 (Figure 3B).
- `top10_cross_cohort_betas.tsv` and `concordant_melanoma_vs_nevus_and_melanocyte.tsv` — Cross-cohort mean $\beta$-values for the top 10 diagnostic promoter CpGs across GSE120878, TCGA-SKCM, and GSE122909 (Table 4 note) and concordant probe check.

## Software Requirements

- Python $\ge 3.10$
- `numpy`, `scipy`, `pandas`, `scikit-learn`, `statsmodels`, `matplotlib`, `seaborn`
