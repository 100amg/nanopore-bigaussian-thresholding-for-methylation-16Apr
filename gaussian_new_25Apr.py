import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
import matplotlib.pyplot as plt

# ============================================================
# OUTPUT SETUP
# ============================================================

DATASET_TYPE = "EpiC"   # change to "AHEAD" when needed

OUTPUT_DIR = Path(f"gaussian_outputs_{DATASET_TYPE}")
OUTPUT_DIR.mkdir(exist_ok=True)

PLOT_DIR = OUTPUT_DIR / "plots"
PLOT_DIR.mkdir(exist_ok=True)

DATA_DIR = OUTPUT_DIR / "distributions"
DATA_DIR.mkdir(exist_ok=True)

# ============================================================
# CONFIG
# ============================================================

TSV_ROOT = Path("modkit_tsv_outputs")
T960_POSITIONS = list(range(60, 901, 24))

# ============================================================
# FILTER FUNCTION
# ============================================================

def is_correct_sample(tsv_path):
    path_str = str(tsv_path).lower()
    return DATASET_TYPE.lower() in path_str

# ============================================================
# LOAD DATA
# ============================================================

def load_group_data(folder):

    dfs = []

    for tsv in folder.glob("*.tsv"):

        if not is_correct_sample(tsv):
            continue

        df = pd.read_csv(tsv, sep="\t")

        df = df[df["mod_code"] == "m"].copy()

        df = df.rename(columns={
            "ref_position": "position",
            "mod_qual": "prob"
        })

        df = df[["read_id", "position", "prob"]]
        df["position"] += 1

        df = df[df["position"].isin(T960_POSITIONS)]

        dfs.append(df)

    if not dfs:
        return None

    return pd.concat(dfs, ignore_index=True)

# ============================================================
# GAUSSIAN FIT
# ============================================================

def fit_gaussian(df):

    site_params = {}

    for pos in sorted(df["position"].unique()):

        vals = df[df["position"] == pos]["prob"].values.reshape(-1,1)

        if len(vals) < 20:
            continue

        gmm = GaussianMixture(n_components=2, random_state=0)
        gmm.fit(vals)

        means = gmm.means_.flatten()
        stds = np.sqrt(gmm.covariances_.flatten())
        weights = gmm.weights_.flatten()

        order = np.argsort(means)

        mu0, mu1 = means[order]
        sigma0, sigma1 = stds[order]
        alpha0, alpha1 = weights[order]

        site_params[pos] = (mu0, sigma0, alpha0, mu1, sigma1, alpha1)

    return site_params

# ============================================================
# COMPUTE THRESHOLD + OBJECTIVE CURVE
# ============================================================

def compute_threshold(site_params):

    def obj(th):
        total = 0
        for (mu0, s0, a0, mu1, s1, a1) in site_params.values():
            fp = a0 * (1 - norm.cdf(th, mu0, s0))
            fn = a1 * norm.cdf(th, mu1, s1)
            total += (fp + fn)
        return -total

    ths = np.linspace(0,1,500)
    vals = [obj(t) for t in ths]

    best_idx = np.argmax(vals)
    threshold = ths[best_idx]

    return threshold, ths, vals

# ============================================================
# MAIN LOOP
# ============================================================

results = []

for folder in TSV_ROOT.glob("*/*"):

    if not folder.is_dir():   # skip .DS_Store
        continue

    group_name = folder.as_posix().replace("modkit_tsv_outputs/", "")
    safe_name = group_name.replace("/", "_")

    print(f"\nProcessing: {group_name} ({DATASET_TYPE})")

    df = load_group_data(folder)

    if df is None:
        print("  No matching samples")
        continue

    print(f"  Total reads: {len(df)}")

    site_params = fit_gaussian(df)

    if not site_params:
        print("  No valid sites")
        continue

    threshold, ths, vals = compute_threshold(site_params)

    print(f"  Threshold: {threshold:.4f}")

    # ========================================================
    # SAVE RESULTS
    # ========================================================

    results.append({
        "group": group_name,
        "dataset": DATASET_TYPE,
        "threshold": round(threshold,4),
        "n_reads": len(df),
        "n_sites": len(site_params)
    })

    # ========================================================
    # SAVE RAW DISTRIBUTION
    # ========================================================

    df.to_csv(DATA_DIR / f"{safe_name}_{DATASET_TYPE}.csv", index=False)

    # ========================================================
    # SAVE HISTOGRAM PLOT
    # ========================================================

    plt.figure(figsize=(6,4))
    plt.hist(df["prob"], bins=50, alpha=0.6)

    plt.axvline(threshold, color='black', linestyle='--',
                label=f"Threshold = {threshold:.3f}")

    plt.title(f"{group_name} ({DATASET_TYPE})")
    plt.xlabel("Methylation probability")
    plt.ylabel("Counts")
    plt.legend()

    plt.tight_layout()
    plt.savefig(PLOT_DIR / f"{safe_name}_{DATASET_TYPE}.png")
    plt.close()

    # ========================================================
    # SAVE OBJECTIVE CURVE
    # ========================================================

    plt.figure(figsize=(6,4))
    plt.plot(ths, vals)

    plt.axvline(threshold, color='red', linestyle='--',
                label=f"chosen = {threshold:.3f}")

    plt.title(f"{group_name} ({DATASET_TYPE})")
    plt.xlabel("Threshold")
    plt.ylabel("Objective")
    plt.legend()

    plt.tight_layout()
    plt.savefig(PLOT_DIR / f"{safe_name}_objective_{DATASET_TYPE}.png")
    plt.close()

# ============================================================
# SAVE FINAL CSV
# ============================================================

df_out = pd.DataFrame(results)

df_out.to_csv(
    OUTPUT_DIR / f"gaussian_threshold_{DATASET_TYPE}.csv",
    index=False
)

print(f"\nSaved → {OUTPUT_DIR}")