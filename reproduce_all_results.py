"""
reproduce_all_results.py
Self-contained verification script for:
'Transcriptomic Immune Infiltration Confounds Bulk Promoter Methylation–Survival Screens in Cutaneous Melanoma
 but Not in Four Other Solid Cancers, with a Promoter Classifier Separating Melanoma from Benign Nevi'

Run from the directory containing the pre-extracted matrices or raw data in C:/Users/ethan/Downloads:
    python reproduce_all_results.py
"""
import os
import subprocess
import sys

SCRIPTS = [
    "run_cox_and_ml.py",
    "run_scientific_upgrades.py",
    "run_stress_test_and_checks.py",
    "check_reviewer_items.py",
    "check_gse122909_and_conway.py",
    "update_candidate_table.py",
    "update_dx_cv.py",
    "run_skcm_deep_analyses.py",
    "run_biopsy_and_bootstrap.py",
    "run_cross_cancer_generalization.py",
    "run_permutation_calibration.py",
    "generate_new_figures.py",
]

if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    for script in SCRIPTS:
        path = os.path.join(here, script)
        if os.path.exists(path):
            print(f"\n==================== Running {script} ====================")
            subprocess.run([sys.executable, path], check=True)
