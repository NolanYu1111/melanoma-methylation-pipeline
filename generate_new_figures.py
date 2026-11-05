import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns

dl = r"C:\Users\ethan\Downloads"

def plot_km_with_ci(ax, times, events, color, label):
    order = np.argsort(times)
    t_s = times[order]
    e_s = events[order]
    utimes = np.unique(t_s[e_s == 1])
    S = [1.0]
    var_sum = [0.0]
    t_pts = [0.0]
    cur_s = 1.0
    cur_v = 0.0
    for ut in utimes:
        n_i = np.sum(t_s >= ut)
        d_i = np.sum((t_s == ut) & (e_s == 1))
        if n_i > 0:
            cur_s *= (1.0 - d_i / n_i)
            if n_i > d_i:
                cur_v += d_i / (n_i * (n_i - d_i))
        t_pts.append(ut)
        S.append(cur_s)
        var_sum.append(cur_v)
    if t_s.max() > t_pts[-1]:
        t_pts.append(t_s.max())
        S.append(cur_s)
        var_sum.append(cur_v)
    t_pts = np.array(t_pts)
    S = np.array(S)
    se = S * np.sqrt(np.array(var_sum))
    ci_low = np.clip(S - 1.96 * se, 0.0, 1.0)
    ci_high = np.clip(S + 1.96 * se, 0.0, 1.0)
    ax.step(t_pts, S, where='post', color=color, linewidth=2.2, label=label)
    ax.fill_between(t_pts, ci_low, ci_high, step='post', color=color, alpha=0.18)

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8

def generate_figure1_flowchart():
    print("Generating Figure 1: Dual-Branch Study Design & Pipeline Flowchart...")
    fig, ax = plt.subplots(figsize=(16, 9.2), dpi=300)
    ax.axis('off')

    # Row 1: Prognostic Survival Branch (5 TCGA Cohorts, led by TCGA-SKCM n = 456)
    boxes_r1 = [
        {"x": 0.02, "y": 0.56, "w": 0.19, "h": 0.33, "title": "5 TCGA Cohorts", "sub": "SKCM (n=456), HNSC (527)\nLUAD (437), BLCA (407), LIHC (371)\n62k-67k Promoter CpGs"},
        {"x": 0.27, "y": 0.56, "w": 0.19, "h": 0.33, "title": "Clinical Cox Screen", "sub": "Age, Sex, Stage (+Type)\nSKCM: 4,956 Hits (5,279 RNA)\nOther 4: 15-85 Hits (11-91 RNA)"},
        {"x": 0.52, "y": 0.56, "w": 0.19, "h": 0.33, "title": "Immune vs. PC Adjust", "sub": "RNA Imm: -99.96% SKCM (0% Others)\n5 PCs: -28% to -100% Across All 5\nSKCM Biopsy Site: 4,281 Hits"},
        {"x": 0.77, "y": 0.56, "w": 0.19, "h": 0.33, "title": "CV & Mediation Check", "sub": "No C-Index Gain (0.647 vs 0.653)\nCRABP2 Partial r = -0.32\nNo mRNA Mediation (p = 0.91)"}
    ]

    # Row 2: Diagnostic Classifier Branch (GSE120878 -> GSE122909 & TCGA-SKCM)
    boxes_r2 = [
        {"x": 0.02, "y": 0.07, "w": 0.19, "h": 0.33, "title": "GSE120878 Training", "sub": "89 Primary Melanomas\nvs. 73 Benign Nevi\n(FFPE Biopsies)"},
        {"x": 0.27, "y": 0.07, "w": 0.19, "h": 0.33, "title": "In-Fold Welch t-Test", "sub": "Top 10, 40, 50, 200 CpGs\nL1, ElasticNet & RF\n(5-Fold CV, 2 Repeats)"},
        {"x": 0.52, "y": 0.07, "w": 0.19, "h": 0.33, "title": "Held-Out CV Results", "sub": "Promoter AUC 0.995-0.998\nAccuracy 95.7-96.9%\n(HOXD13, PAX3, PAX6)"},
        {"x": 0.77, "y": 0.07, "w": 0.19, "h": 0.33, "title": "External Stress Test", "sub": "TCGA: 88.8-99.6% Sensitivity\nGSE122909: AUC 0.45-0.50\n(Melanoma vs. Nevus Scope)"}
    ]

    for b in boxes_r1:
        rect = patches.FancyBboxPatch(
            (b["x"], b["y"]), b["w"], b["h"],
            boxstyle="round,pad=0.008", fc="#e8f4fc", ec="#0055b3", lw=3.0
        )
        ax.add_patch(rect)
        ax.text(b["x"] + b["w"]/2, b["y"] + b["h"]*0.72, b["title"],
                ha='center', va='center', fontsize=16.5, fontweight='bold', color='#00264d')
        ax.text(b["x"] + b["w"]/2, b["y"] + b["h"]*0.32, b["sub"],
                ha='center', va='center', fontsize=12.0, fontstyle='italic', fontweight='bold', color='#004080')

    for b in boxes_r2:
        rect = patches.FancyBboxPatch(
            (b["x"], b["y"]), b["w"], b["h"],
            boxstyle="round,pad=0.008", fc="#eefaf0", ec="#1e7b34", lw=3.0
        )
        ax.add_patch(rect)
        ax.text(b["x"] + b["w"]/2, b["y"] + b["h"]*0.72, b["title"],
                ha='center', va='center', fontsize=16.5, fontweight='bold', color='#0d471a')
        ax.text(b["x"] + b["w"]/2, b["y"] + b["h"]*0.32, b["sub"],
                ha='center', va='center', fontsize=12.0, fontstyle='italic', fontweight='bold', color='#19662b')

    arrow_blue = dict(arrowstyle="-|>", lw=4.0, color='#0055b3', mutation_scale=24)
    arrow_green = dict(arrowstyle="-|>", lw=4.0, color='#1e7b34', mutation_scale=24)

    for x_start, x_end in [(0.222, 0.258), (0.472, 0.508), (0.722, 0.758)]:
        ax.annotate("", xy=(x_end, 0.725), xytext=(x_start, 0.725), arrowprops=arrow_blue)
        ax.annotate("", xy=(x_end, 0.235), xytext=(x_start, 0.235), arrowprops=arrow_green)

    ax.text(0.49, 0.945, "(A) Cross-Cancer Survival Pipeline, Immune Score Attenuation & Mediation Analysis (5 TCGA Cohorts, n = 2,198)",
            ha='center', va='center', fontsize=17.5, fontweight='bold', color='#00264d')
    ax.text(0.49, 0.455, "(B) Diagnostic Melanoma vs. Nevus Classification & External Stress Testing (GSE120878 / GSE122909)",
            ha='center', va='center', fontsize=17.5, fontweight='bold', color='#0d471a')

    plt.tight_layout()
    plt.savefig('pipeline_flowchart.png', bbox_inches='tight', dpi=300)
    plt.close()

def generate_figure2_survival_and_composition():
    print("Generating Figure 2: Tumor Composition Attenuation, RNA-Seq Immune Correlation, Cross-Cancer Generalization & Candidate Survival Curves...")
    X_raw = np.load("tcga_promoter_X_raw.npy")
    probes_df = pd.read_csv("tcga_promoter_probes.csv")
    meta = pd.read_csv("tcga_456_meta.csv")

    from sklearn.decomposition import PCA
    row_means = np.nanmean(X_raw, axis=1, keepdims=True)
    X_imp = np.where(np.isnan(X_raw), row_means, X_raw)
    row_sds = np.std(X_imp, axis=1, ddof=1, keepdims=True)
    X_std = (X_imp - row_means) / row_sds
    pca = PCA(n_components=5, random_state=42)
    pcs = pca.fit_transform(X_std.T)

    rnaseq_df = pd.read_csv(os.path.join(dl, "TCGA-SKCM.HiSeqV2.tsv.gz"), sep="\t", index_col=0)
    meta["sample_15"] = meta["sample"].str.slice(0, 15)
    has_rna_mask = meta["sample_15"].isin(rnaseq_df.columns).values
    rna_samples_matched = meta.loc[has_rna_mask, "sample_15"].values
    rna_matched_df = rnaseq_df[rna_samples_matched]
    immune_markers = ["PTPRC", "CD8A", "CD8B", "CD3D", "CD3E", "GZMA", "PRF1", "IFNG", "CXCL9", "CXCL10", "STAT1", "PSMB9", "TAP1", "HLA-DRA"]
    immune_expr = rna_matched_df.loc[immune_markers].values
    immune_z = (immune_expr - immune_expr.mean(axis=1, keepdims=True)) / immune_expr.std(axis=1, ddof=1, keepdims=True)
    immune_score_454 = immune_z.mean(axis=0)
    pc2_454 = pcs[has_rna_mask, 1]

    fig = plt.figure(figsize=(13.5, 9.0), dpi=300)
    gs = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.28)

    # Panel A: Bar chart of FDR-significant CpGs across adjustment levels
    ax0 = fig.add_subplot(gs[0, 0])
    stages = ["Clinical\nOnly", "+ Biopsy\nSite", "+ Top 3\nPCs", "+ Top 5\nPCs", "+ RNA\nImmune"]
    counts_q05 = [4956, 4281, 49, 78, 2]
    bars = ax0.bar(stages, counts_q05, color=['#3b82f6', '#64748b', '#f59e0b', '#10b981', '#8b5cf6'], edgecolor='black', width=0.58)
    ax0.set_yscale('log')
    ax0.set_ylim(1, 25000)
    ax0.set_ylabel('FDR-Significant CpGs (q < 0.05)', fontsize=11.0, fontweight='bold')
    ax0.set_title('(A) SKCM Survival Hits After Adjustment', fontsize=12.0, fontweight='bold', pad=8)
    ax0.tick_params(axis='both', labelsize=9.8)
    ax0.grid(axis='y', linestyle='--', alpha=0.4)
    for bar, c in zip(bars, counts_q05):
        ax0.text(bar.get_x() + bar.get_width()/2, c * 1.30, f"{c:,}", ha='center', va='bottom', fontsize=10.0, fontweight='bold')

    # Panel B: Promoter Methylation PC2 vs RNA-seq Immune Infiltration Score (n = 454)
    ax1 = fig.add_subplot(gs[0, 1])
    ax1.scatter(pc2_454, immune_score_454, alpha=0.45, s=20, color='#2563eb', edgecolors='none')
    m_fit, b_fit = np.polyfit(pc2_454, immune_score_454, 1)
    x_line = np.linspace(pc2_454.min(), pc2_454.max(), 100)
    ax1.plot(x_line, m_fit * x_line + b_fit, color='#dc2626', lw=2.4, label='PC2 r = +0.715 (p = 2.6e-72)\nNon-LN r = 0.692 | LN r = 0.722')
    ax1.set_xlabel('Promoter Methylation PC2 Score', fontsize=11.0, fontweight='bold')
    ax1.set_ylabel('RNA-Seq Immune Score (Z)', fontsize=11.0, fontweight='bold')
    ax1.set_title('(B) PC2 Tracks Immune Infiltration', fontsize=12.0, fontweight='bold', pad=8)
    ax1.tick_params(axis='both', labelsize=10.0)
    ax1.legend(loc='upper left', fontsize=9.2, framealpha=0.92)
    ax1.grid(True, linestyle='--', alpha=0.4)
    ax1.text(0.97, 0.05, "CV C-Index: Clin 0.653\n+Imm: 0.675 | +50 CpGs: 0.647",
             transform=ax1.transAxes, ha='right', va='bottom', fontsize=9.2, fontweight='bold',
             bbox=dict(boxstyle="round,pad=0.3", fc="#f8fafc", ec="#cbd5e1", lw=1.0))

    # Panel C: Cross-Cancer Generalization across 5 TCGA Cohorts (Attenuation % vs Immune Score HR)
    ax2 = fig.add_subplot(gs[0, 2])
    cc_df = pd.read_csv("cross_cancer_summary.tsv", sep="\t")
    c_names = [c.replace("TCGA-", "") for c in cc_df["Cohort"]]
    imm_hrs = cc_df["Immune_OS_HR"].values
    pct_att_5pc = cc_df["Pct_Clin_Att_PC5"].values
    pct_att_imm = cc_df["Pct_Clin_Att_Imm"].values

    ax2.scatter(imm_hrs, pct_att_imm, s=95, color='#dc2626', edgecolor='black', zorder=4, label='+ RNA Immune Score')
    ax2.scatter(imm_hrs, pct_att_5pc, s=90, marker='s', color='#10b981', edgecolor='black', zorder=3, label='+ Top 5 Promoter PCs')
    ax2.axvline(1.0, color='#64748b', linestyle=':', lw=1.2)

    offsets_imm = {
        "SKCM": (0.012, -10.0),
        "HNSC": (-0.082, 8.0),
        "LIHC": (-0.032, 33.0),
        "BLCA": (-0.025, 13.0),
        "LUAD": (0.012, -4.0)
    }
    for name, hr, p_imm, p_5pc, hits_full, hits_rna in zip(
        c_names, imm_hrs, pct_att_imm, pct_att_5pc, cc_df["Clin_q05"], cc_df["Clin_RNA_q05"]
    ):
        dx, dy = offsets_imm.get(name, (0.01, 5.0))
        if name == "SKCM":
            ax2.annotate(f"SKCM ({hits_rna:,}/{hits_full:,})\nImm HR=0.68 (p=7e-8)", xy=(hr, p_imm),
                         xytext=(hr + 0.015, p_imm - 22), fontsize=8.4, fontweight='bold', color='#991b1b',
                         arrowprops=dict(arrowstyle="->", color='#991b1b', lw=1.0))
        else:
            ax2.text(hr + dx, p_imm + dy, f"{name}\n({hits_rna}/{hits_full})", fontsize=8.0, fontweight='bold', color='#1e293b')

    ax2.set_xlim(0.64, 1.01)
    ax2.set_ylim(-14, 112)
    ax2.set_xlabel('14-Gene Immune Score Survival HR', fontsize=11.0, fontweight='bold')
    ax2.set_ylabel('% of Clinical Hits Attenuated', fontsize=11.0, fontweight='bold')
    ax2.set_title('(C) RNA Immune vs. 5-PC Attenuation', fontsize=12.0, fontweight='bold', pad=8)
    ax2.tick_params(axis='both', labelsize=10.0)
    ax2.legend(loc='center left', fontsize=8.6, framealpha=0.92)
    ax2.grid(True, linestyle='--', alpha=0.4)

    # Panels D, E, F: Real Kaplan-Meier plots from TCGA-SKCM (n = 456) for KRAS, APOL3, CRABP2
    km_loci = [
        {"probe": "cg23152216", "gene": "KRAS", "panel": "D", "sub": "5-PC HR=0.66, q=0.0054\nmRNA r=-0.004 (10bp SNP)", "grid_pos": (1, 0)},
        {"probe": "cg06703343", "gene": "APOL3", "panel": "E", "sub": "5-PC HR=0.60, q=0.0019\nClin q=0.031 | mRNA r=+0.10", "grid_pos": (1, 1)},
        {"probe": "cg02363950", "gene": "CRABP2", "panel": "F", "sub": "5-PC HR=0.64, q=0.0054\nmRNA r=-0.31 (No med. p=0.91)", "grid_pos": (1, 2)},
    ]

    t_years = meta["OS.time"].values / 365.25
    events = meta["OS"].values.astype(int)

    for loc in km_loci:
        ax = fig.add_subplot(gs[loc["grid_pos"][0], loc["grid_pos"][1]])
        p_idx = probes_df[probes_df["probe"] == loc["probe"]].index[0]
        vals = X_raw[p_idx, :]
        valid = ~np.isnan(vals)
        med = np.median(vals[valid])
        high_mask = valid & (vals >= med)
        low_mask = valid & (vals < med)

        plot_km_with_ci(ax, t_years[high_mask], events[high_mask], '#1f77b4', f"High Meth. (n={high_mask.sum()})")
        plot_km_with_ci(ax, t_years[low_mask], events[low_mask], '#d62728', f"Low Meth. (n={low_mask.sum()})")

        ax.set_title(f'({loc["panel"]}) {loc["gene"]} ({loc["probe"]})', fontsize=12.0, fontweight='bold', pad=6)
        ax.set_xlabel('Overall Survival Time (Years)', fontsize=11.0, fontweight='bold')
        ax.set_ylabel('Survival Probability', fontsize=11.0, fontweight='bold')
        ax.tick_params(axis='both', labelsize=10.0)
        ax.set_ylim(-0.02, 1.03)
        ax.set_xlim(0, 20)
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.legend(loc='upper right', fontsize=9.4, framealpha=0.92)
        ax.text(0.03, 0.05, loc["sub"], transform=ax.transAxes, fontsize=9.2, fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.28", fc="#ffffff", ec="#cccccc", lw=1.0))

    plt.suptitle('Cross-Cancer Immune Attenuation (5 TCGA Cohorts) and Candidate Survival Loci in TCGA-SKCM', fontsize=14.0, fontweight='bold', y=0.98)
    plt.savefig('km_curve.png', bbox_inches='tight', dpi=300)
    plt.close()

def generate_figure3_boxplots():
    print("Generating Supplementary Figure S1: Real Promoter Methylation Boxplots in GSE120878...")
    top_diag = pd.read_csv("top50_diagnostic_gse120878.tsv", sep="\t")
    X120_prom = np.load("gse120878_prom_X.npy")
    y120 = np.load("gse120878_y.npy")
    prom_probes = pd.read_csv("gse120878_prom_probes.csv")["probe"].values
    probe_to_idx = {p: i for i, p in enumerate(prom_probes)}
    groups = ["Melanoma (n=89)" if y == 1 else "Benign Nevus (n=73)" for y in y120]

    exact_probes = [
        "cg10418524", # HOXD13 (#1)
        "cg19352038", # PAX3 (#2)
        "cg16529477", # PAX3 (#3)
        "cg09656389", # PAX6 (#5)
        "cg05845376", # SLC25A2 (#7)
        "cg27030541", # REXO1L2P (#8)
        "cg26149275", # EVX2 (#9)
        "cg05787556", # TLX3 (#10)
    ]
    selected_rows = [top_diag[top_diag["probe"] == p].iloc[0] for p in exact_probes]

    fig, axes = plt.subplots(2, 4, figsize=(15, 7.8), dpi=300)
    axes = axes.flatten()

    for idx, row in enumerate(selected_rows):
        ax = axes[idx]
        probe = row["probe"]
        gene = row["Gene"]
        vals = X120_prom[probe_to_idx[probe]].astype(float)
        plot_df = pd.DataFrame({"Group": groups, "Beta Value": vals})

        sns.boxplot(
            x="Group", y="Beta Value", data=plot_df, ax=ax, hue="Group",
            palette={"Benign Nevus (n=73)": "#93c5fd", "Melanoma (n=89)": "#fca5a5"},
            order=["Benign Nevus (n=73)", "Melanoma (n=89)"],
            width=0.48, showfliers=False, linewidth=1.4, legend=False
        )
        sns.stripplot(
            x="Group", y="Beta Value", data=plot_df, ax=ax,
            order=["Benign Nevus (n=73)", "Melanoma (n=89)"],
            color="#1e293b", alpha=0.45, jitter=0.18, size=3.5
        )
        ax.set_title(f"{gene} ({probe})\nt = {row['t_stat']:.1f}, p = {row['p_val']:.1e}", fontsize=10.5, fontweight='bold', pad=6)
        ax.set_xlabel("")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Nevus\n(n=73)", "Melanoma\n(n=89)"], fontsize=9.5, fontweight='bold')
        if idx % 4 == 0:
            ax.set_ylabel("Promoter Methylation (Beta)", fontsize=10, fontweight='bold')
        else:
            ax.set_ylabel("")
        ax.set_ylim(-0.04, 1.04)
        ax.grid(True, linestyle='--', alpha=0.35)

    plt.suptitle("Top Differentially Methylated Promoter CpGs in GSE120878 (89 Primary Melanomas vs. 73 Benign Nevi)", fontsize=13.5, fontweight='bold', y=0.99)
    plt.tight_layout()
    plt.savefig("methylation_boxplots_combined.png", bbox_inches='tight', dpi=300)
    plt.close()

def generate_figure4_ml():
    print("Generating Figure 3 (Diagnostic CV Performance and Top Features)...")
    dx_df = pd.read_csv("dx_cv.tsv", sep="\t")
    top_diag = pd.read_csv("top50_diagnostic_gse120878.tsv", sep="\t")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.0, 6.2), dpi=300)

    prom_df = dx_df[dx_df["Probe_Set"] == "Promoter-Only"].copy()
    prom_df["Label"] = prom_df["Model"] + "\n(" + prom_df["k_Features"].astype(str) + " CpGs)"

    y_pos = np.arange(len(prom_df))
    ax1.barh(y_pos + 0.17, prom_df["CV_AUC_Mean"], xerr=prom_df["CV_AUC_SD"], height=0.32,
             color='#2563eb', edgecolor='black', capsize=3, label='GSE120878 CV AUC (Melanoma vs. Nevus)')
    ax1.barh(y_pos - 0.17, prom_df["GSE122909_AUC"], height=0.32,
             color='#f59e0b', edgecolor='black', label='GSE122909 Test AUC (Melanoma vs. Melanocyte)')
    ax1.axvline(0.50, color='#dc2626', linestyle='--', lw=1.4, label='Chance (AUC = 0.50)')
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(prom_df["Label"], fontsize=9.2, fontweight='bold')
    ax1.set_xlim(0.0, 1.22)
    ax1.set_ylim(-1.35, len(prom_df) - 0.35)
    ax1.set_xlabel("ROC-AUC (5-Fold CV Repeated Twice & External Stress Test)", fontsize=10, fontweight='bold')
    ax1.set_title("(A) Classifier ROC-AUC: Nevi vs. Cultured Melanocytes", fontsize=11, fontweight='bold', pad=10)
    ax1.legend(loc='lower right', fontsize=8.3, framealpha=0.95)
    ax1.grid(axis='x', linestyle='--', alpha=0.5)
    for i, (auc_cv, acc_cv, auc_ext) in enumerate(zip(prom_df["CV_AUC_Mean"], prom_df["CV_Acc_Mean"], prom_df["GSE122909_AUC"])):
        ax1.text(auc_cv + 0.012, i + 0.17, f"{auc_cv:.3f} ({acc_cv:.1f}%)", va='center', fontsize=8.2, fontweight='bold', color='#1e3a8a')
        ax1.text(auc_ext + 0.012, i - 0.17, f"{auc_ext:.3f}", va='center', fontsize=8.2, fontweight='bold', color='#92400e')

    # Panel B: Top 10 Promoter CpGs by |Welch t-statistic| in GSE120878
    top10 = top_diag.head(10).iloc[::-1]
    labels_10 = [f"{r['probe']} ({r['Gene']})" for _, r in top10.iterrows()]
    colors_10 = ['#dc2626' if d > 0 else '#2563eb' for d in top10["delta_beta"]]
    ax2.barh(np.arange(10), top10["abs_t"], color=colors_10, edgecolor='black', height=0.6)
    ax2.set_yticks(np.arange(10))
    ax2.set_yticklabels(labels_10, fontsize=9.5, fontweight='bold')
    ax2.set_xlabel("|Welch t-Statistic| (Melanoma vs. Benign Nevus)", fontsize=10, fontweight='bold')
    ax2.set_title("(B) Top 10 Diagnostic Promoter CpGs (GSE120878)", fontsize=11, fontweight='bold', pad=10)
    ax2.grid(axis='x', linestyle='--', alpha=0.5)
    for i, (_, r) in enumerate(top10.iterrows()):
        direction = "+Δβ" if r["delta_beta"] > 0 else "-Δβ"
        ax2.text(r["abs_t"] + 0.25, i, f"t={r['t_stat']:.1f} ({direction}={r['delta_beta']:+.2f})", va='center', fontsize=8.6, fontweight='bold')
    ax2.set_xlim(0, top10["abs_t"].max() * 1.30)

    plt.tight_layout()
    plt.savefig("ml_classifier_accuracy.png", bbox_inches='tight', dpi=300)
    plt.savefig("ml_roc_and_feature_importance.png", bbox_inches='tight', dpi=300)
    plt.close()

if __name__ == "__main__":
    generate_figure1_flowchart()
    generate_figure2_survival_and_composition()
    generate_figure3_boxplots()
    generate_figure4_ml()
    print("All 4 figures generated successfully!")
