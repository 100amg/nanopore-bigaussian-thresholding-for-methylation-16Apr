import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import matthews_corrcoef

# ============================================================
# CONFIG
# ============================================================

TSV_ROOT = Path("/Volumes/Amishi_SSD/bio_data/16Apr/modkit_tsv_outputs")

# ============================================================
# IMPORTANT
# ============================================================

DATASET_TYPE = "epic"

OUTPUT_DIR = Path("final_results_epic")
OUTPUT_DIR.mkdir(exist_ok=True)

# ============================================================
# EpiC positions (32 positions only)
# ============================================================

EPIC_POSITIONS = list(range(60, 805, 24))

# ============================================================
# EpiC Ground Truth
# ============================================================

GROUND_TRUTH_EPIC = [
    0,1,0,0,0,1,0,1,
    0,1,1,1,0,0,0,0,
    0,1,1,0,1,0,0,1,
    0,1,0,0,0,0,1,1
]

# ============================================================
# LOAD THRESHOLDS
# ============================================================

THRESHOLD_FILE = (
    "/Volumes/Amishi_SSD/bio_data/16Apr/"
    "gaussian_outputs_EpiC/"
    "gaussian_threshold_EpiC.csv"
)

th_df = pd.read_csv(THRESHOLD_FILE)

threshold_map = dict(
    zip(th_df["group"], th_df["threshold"])
)

# ============================================================
# LOAD ONLY EPIC FILES
# ============================================================

def load_data(folder):

    dfs = []

    for tsv in folder.glob("*.tsv"):

        name = tsv.name.lower()

        # ====================================================
        # ONLY LOAD EPIC FILES
        # ====================================================

        if "epic" not in name:
            continue

        df = pd.read_csv(tsv, sep="\t")

        # ====================================================
        # KEEP ONLY METHYLATION CALLS
        # ====================================================

        df = df[df["mod_code"] == "m"].copy()

        df = df.rename(columns={
            "ref_position": "position",
            "mod_qual": "prob"
        })

        df = df[["read_id", "position", "prob"]]

        # ====================================================
        # POSITION OFFSET
        # ====================================================

        df["position"] += 1

        dfs.append(df)

    if not dfs:
        return None

    return pd.concat(dfs, ignore_index=True)

# ============================================================
# COMPUTE MCC
# ============================================================

def compute_mcc(df, threshold):

    pred_bits = []

    for pos in EPIC_POSITIONS:

        site_df = df[df["position"] == pos]

        # ====================================================
        # MINIMUM COVERAGE FILTER
        # ====================================================

        if len(site_df) < 5:
            pred_bits.append(0)
            continue

        # ====================================================
        # APPLY GAUSSIAN THRESHOLD
        # ====================================================

        calls = (
            site_df["prob"] > threshold
        ).astype(int)

        # ====================================================
        # FRACTION OF METHYLATED READS
        # ====================================================

        frac_high = calls.mean()

        # ====================================================
        # SITE DECISION
        # ====================================================

        site_call = int(frac_high > 0.05)

        pred_bits.append(site_call)

    y_true = np.array(GROUND_TRUTH_EPIC)
    y_pred = np.array(pred_bits)

    return matthews_corrcoef(y_true, y_pred)

# ============================================================
# MAIN LOOP
# ============================================================

results = []

for folder in TSV_ROOT.glob("*/*"):

    if not folder.is_dir():
        continue

    group_name = "/".join(folder.parts[-2:])

    print(f"\nProcessing: {group_name} (epic)")

    # ========================================================
    # THRESHOLD CHECK
    # ========================================================

    if group_name not in threshold_map:
        print("  No threshold found, skipping")
        continue

    threshold = threshold_map[group_name]

    df = load_data(folder)

    if df is None:
        print("  No EpiC TSV files")
        continue

    # ========================================================
    # KEEP ONLY VALID EPIC POSITIONS
    # ========================================================

    df = df[df["position"].isin(EPIC_POSITIONS)]

    print(f"  Data points: {len(df)}")
    print(f"  Threshold: {threshold:.4f}")
    print(f"  Positions: {len(EPIC_POSITIONS)}")

    mcc = compute_mcc(df, threshold)

    print(f"  MCC: {mcc:.4f}")

    results.append({
        "group": group_name,
        "threshold": round(threshold, 4),
        "MCC": round(mcc, 4),
        "n_points": len(df)
    })

# ============================================================
# SAVE RESULTS
# ============================================================

df_out = pd.DataFrame(results)

out_file = (
    OUTPUT_DIR /
    "final_results_5_reads_filter_EpiC.csv"
)

df_out.to_csv(out_file, index=False)

print(f"\nSaved → {out_file}")
