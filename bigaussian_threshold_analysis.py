#!/usr/bin/env python3
"""
Bi-Gaussian threshold analysis for T960 methylation calling.

Builds on the AUC analysis by finding the optimal decoding threshold per
(basecaller, model, qscore, sample_type) combination using a 2-component
Gaussian Mixture Model fit on pooled fraction_modified values — the same
approach Zhang et al. used to derive their 0.345 global threshold on T960.

For each condition:
  1. Pool fraction_modified from BED cache across all replicates of one
     sample type (e.g. 3 AHEAD samples x 36 sites = 108 values).
  2. Fit 2-component GMM (n_init=10) to the pooled distribution.
  3. Validate the fit — reject clearly unimodal distributions rather than
     reporting a meaningless threshold.
  4. If valid, compute the analytical crossover between the two weighted
     Gaussian PDFs — this is the bi-Gaussian threshold.
  5. Apply that single threshold to every replicate of that sample type
     and compute decoding accuracy against ground truth.
  6. Also compute a pooled-substrate threshold (AHEAD+EpiC together) as a
     deployment realism check.

Usage:
    python bigaussian_threshold_analysis.py <auc_results_dir>

    where <auc_results_dir> is the directory containing bed_cache/ from
    the AUC analysis script.

Outputs (written to <auc_results_dir>/bigaussian/):
    bigaussian_thresholds.csv           one row per condition x sample_type
    bigaussian_per_sample_decoding.csv  one row per individual sample
    bigaussian_ranking.csv              conditions ranked by mean accuracy
    diagnostic_plots/                   histogram + GMM fit per condition
"""

import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.mixture import GaussianMixture

# ============================================================
# GROUND TRUTH (matches auc_model_qscore_analysis.py)
# ============================================================

T960_CPG_POSITIONS = list(range(60, 901, 24))  # 36 sites

AHEAD_BINARY = "100000110010001000101100000110001000"
AHEAD_GROUND_TRUTH = {pos: int(bit) for pos, bit in zip(T960_CPG_POSITIONS, AHEAD_BINARY)}

EPIC_BINARY = "010001010111000001101001010000110000"
EPIC_GROUND_TRUTH = {pos: int(bit) for pos, bit in zip(T960_CPG_POSITIONS, EPIC_BINARY)}

GROUND_TRUTH = {"AHEAD": AHEAD_GROUND_TRUTH, "EpiC": EPIC_GROUND_TRUTH}

# ============================================================
# FIT VALIDITY THRESHOLDS
# ============================================================
# A 2-component GMM will always "fit" any data, but if the two components
# are essentially the same distribution (or one has almost no mass), the
# crossover is meaningless. These checks catch that.

MIN_COMPONENT_WEIGHT  = 0.05   # each Gaussian must hold >=5% of the mass
                               # (low methylation efficiency often produces
                               #  weight splits like 92%/8%)
MIN_MEAN_SEPARATION   = 0.25   # means must differ by >=0.25
MIN_HIGHER_MEAN       = 0.30   # methylated component must sit at >=0.30
                               # (Zhang et al.'s 0.345 threshold came from a
                               #  methylated peak around ~0.5; some substrates
                               #  show a methylated peak even lower)
MAX_LOWER_MEAN        = 0.225  # unmethylated component must sit at <=0.225
                               # (Zhang et al. Fig 2c shows epi-bit 0 sites
                               #  capping below 0.20; 0.225 gives a small
                               #  margin for real experimental variation)

# Threshold (crossover) sanity range: a biologically sensible global
# threshold on fraction_modified should sit in this band. Outside it, the
# GMM crossover is a math artifact of a degenerate fit (e.g. crossover at
# 0.05 means "call everything methylated", crossover at 0.99 means "call
# nothing methylated" — neither is a real decoding threshold).
MIN_THRESHOLD = 0.15
MAX_THRESHOLD = 0.75

# ============================================================
# PARSING THE BED CACHE
# ============================================================

# BED filenames look like AHEAD_10_11_25_hac_methyl_q0.bed,
# AHEAD_10_11_25_guppy_hac_methyl_q0.bed, etc.
# We parse out: sample_type, basecaller, model, qscore.

BED_RE = re.compile(
    r"^(AHEAD|EpiC)_.+?_"                  # sample type + date
    r"(dorado_|guppy_|)"                   # optional explicit basecaller prefix in name
    r"(fast|hac|sup)_"                     # model
    r"methyl_q(\d+)"                       # qscore
    r"\.bed$"
)

def parse_bed_filename(fname: str):
    """
    Parse sample_type, basecaller, model, qscore from a cached BED filename.

    The AUC script cached BEDs using the raw sample_name (from bam.stem),
    so we infer basecaller from whether 'guppy' appears in the stem.
    Returns (sample_type, basecaller, model, qscore, sample_stem) or None.
    """
    stem = fname[:-4] if fname.endswith(".bed") else fname

    # Sample type from prefix
    if stem.startswith("AHEAD"):
        sample_type = "AHEAD"
    elif stem.startswith("EpiC"):
        sample_type = "EpiC"
    else:
        return None

    # Basecaller from presence of "guppy" in stem
    basecaller = "guppy" if "guppy" in stem else "dorado"

    # Model
    m_model = re.search(r"_(fast|hac|sup)_methyl_q(\d+)$", stem)
    if not m_model:
        return None
    model  = m_model.group(1)
    qscore = int(m_model.group(2))

    return sample_type, basecaller, model, qscore, stem


def extract_fractions_from_bed(bed_path: Path, target_positions: set) -> dict:
    """Read {pos: fraction_modified} for modkit-style BED, 5mC channel only.

    modkit pileup outputs column 11 as a PERCENTAGE (0-100), not a fraction.
    We divide by 100 to normalize to the [0, 1] range expected by the rest
    of the pipeline (AUC, GMM fitting, threshold validity checks).
    """
    fractions = {}
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 11 or cols[3] != "m":
                continue
            try:
                pos  = int(cols[2])
                frac = float(cols[10]) / 100.0   # percent -> fraction
            except ValueError:
                continue
            if pos in target_positions:
                fractions[pos] = frac
    return fractions


# ============================================================
# BI-GAUSSIAN FITTING
# ============================================================

def fit_bigaussian(values: np.ndarray, random_state: int = 0):
    """
    Fit a 2-component GaussianMixture and return component parameters
    sorted so that component 0 has the lower mean.

    Returns dict with: means, stds, weights (all as np arrays of length 2),
    plus 'converged' and the fitted GMM object.
    """
    X = values.reshape(-1, 1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gmm = GaussianMixture(
            n_components=2,
            n_init=10,
            random_state=random_state,
            covariance_type="full",
        )
        gmm.fit(X)

    means   = gmm.means_.flatten()
    stds    = np.sqrt(gmm.covariances_.flatten())
    weights = gmm.weights_.flatten()

    # Sort so component 0 is the "unmethylated" (lower mean) one
    order   = np.argsort(means)
    means   = means[order]
    stds    = stds[order]
    weights = weights[order]

    return {
        "means":    means,
        "stds":     stds,
        "weights":  weights,
        "converged": gmm.converged_,
        "gmm":      gmm,
    }


def check_fit_validity(fit: dict) -> tuple[bool, str]:
    """
    Returns (is_valid, reason_if_invalid). Applies the three sanity
    checks defined at the top of the module.
    """
    means, weights = fit["means"], fit["weights"]

    if not fit["converged"]:
        return False, "GMM did not converge"
    if weights.min() < MIN_COMPONENT_WEIGHT:
        return False, f"component weight {weights.min():.3f} < {MIN_COMPONENT_WEIGHT}"
    if (means[1] - means[0]) < MIN_MEAN_SEPARATION:
        return False, f"mean separation {means[1] - means[0]:.3f} < {MIN_MEAN_SEPARATION}"
    if means[1] < MIN_HIGHER_MEAN:
        return False, f"higher mean {means[1]:.3f} < {MIN_HIGHER_MEAN} (no clear methylated peak)"
    if means[0] > MAX_LOWER_MEAN:
        return False, f"lower mean {means[0]:.3f} > {MAX_LOWER_MEAN} (no clear unmethylated peak)"
    return True, ""


def find_crossover(fit: dict) -> float | None:
    """
    Analytical crossover of two weighted Gaussians. Solves
      w0 * N(x; mu0, s0) = w1 * N(x; mu1, s1)
    which reduces to a quadratic in x. Returns the crossover that falls
    BETWEEN the two means AND within a sensible methylation range.
    Returns None if no such crossover exists.
    """
    mu0, mu1   = fit["means"]
    s0, s1     = fit["stds"]
    w0, w1     = fit["weights"]

    # Coefficients of a*x^2 + b*x + c = 0 from
    # log(w0/s0) - (x-mu0)^2/(2*s0^2)  =  log(w1/s1) - (x-mu1)^2/(2*s1^2)
    a = 1.0/(2*s0**2) - 1.0/(2*s1**2)
    b = mu1/(s1**2) - mu0/(s0**2)
    c = (mu0**2)/(2*s0**2) - (mu1**2)/(2*s1**2) - np.log((w0*s1)/(w1*s0))

    if abs(a) < 1e-12:
        # Equal variances -> linear equation
        if abs(b) < 1e-12:
            return None
        x = -c / b
        if mu0 <= x <= mu1 and MIN_THRESHOLD <= x <= MAX_THRESHOLD:
            return float(x)
        return None

    disc = b*b - 4*a*c
    if disc < 0:
        return None

    sqrt_disc = np.sqrt(disc)
    roots = [(-b + sqrt_disc) / (2*a), (-b - sqrt_disc) / (2*a)]

    # Accept only crossovers that both (a) lie between the two means and
    # (b) fall within the biologically sensible threshold range. (b) rejects
    # degenerate fits where the crossover math technically converges but
    # lands at a value that can't function as a decoding threshold (e.g.
    # 0.05 or 0.99).
    between = [
        r for r in roots
        if mu0 <= r <= mu1 and MIN_THRESHOLD <= r <= MAX_THRESHOLD
    ]

    if not between:
        return None
    return float(between[0])


# ============================================================
# DECODING WITH A FIXED THRESHOLD
# ============================================================

def decode_with_threshold(fractions: dict, ground_truth: dict, threshold: float) -> dict:
    """Apply a single threshold and score against ground truth."""
    positions = sorted(ground_truth.keys())
    n_matches = n_1to0 = n_0to1 = 0
    n = 0
    per_site = []
    for pos in positions:
        if pos not in fractions:
            continue
        n += 1
        score = fractions[pos]
        called = 1 if score > threshold else 0
        truth  = ground_truth[pos]
        correct = called == truth
        if correct:
            n_matches += 1
        elif truth == 1:
            n_1to0 += 1
        else:
            n_0to1 += 1
        per_site.append({
            "position": pos, "score": score, "truth": truth,
            "called": called, "correct": correct,
        })

    return {
        "n_sites":      n,
        "n_matches":    n_matches,
        "n_1to0_flips": n_1to0,
        "n_0to1_flips": n_0to1,
        "accuracy":     n_matches / n if n else 0.0,
        "per_site":     per_site,
    }


# ============================================================
# DIAGNOSTIC PLOT
# ============================================================

def find_oracle_threshold(scores, labels):
    """
    Given fraction_modified values and their ground-truth labels,
    find the threshold that maximizes decoding accuracy. This is the
    best threshold achievable *with hindsight* — a ceiling to compare
    the unsupervised bi-Gaussian crossover against.

    Sweeps candidate thresholds at every midpoint between sorted scores
    (plus the endpoints) and returns the one with the highest accuracy.

    Returns (threshold, accuracy). If multiple thresholds tie, returns
    the one closest to 0.5 for stability.
    """
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels, dtype=int)

    if len(scores) == 0:
        return None, 0.0

    # Candidate thresholds: midpoints between consecutive sorted unique
    # scores, plus a bit below the min and a bit above the max.
    sorted_scores = np.sort(np.unique(scores))
    if len(sorted_scores) < 2:
        # Degenerate: all scores identical
        return float(sorted_scores[0]), (labels == 0).mean()

    midpoints = (sorted_scores[:-1] + sorted_scores[1:]) / 2.0
    candidates = np.concatenate([[sorted_scores[0] - 0.001],
                                 midpoints,
                                 [sorted_scores[-1] + 0.001]])

    best_acc = -1.0
    best_thresh = None
    for t in candidates:
        called = (scores > t).astype(int)
        acc = (called == labels).mean()
        # Prefer the threshold closest to 0.5 on ties for stability
        if acc > best_acc or (acc == best_acc and
                              abs(t - 0.5) < abs(best_thresh - 0.5)):
            best_acc = acc
            best_thresh = t

    return float(best_thresh), float(best_acc)


def plot_diagnostic(values, fit, threshold, valid, reason, title, out_path):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(values, bins=40, density=True, alpha=0.5, color="steelblue",
            edgecolor="white")

    x = np.linspace(0, 1, 500)
    colors = ["tab:blue", "tab:orange"]
    for i in range(2):
        pdf = fit["weights"][i] * norm.pdf(x, fit["means"][i], fit["stds"][i])
        ax.plot(x, pdf, color=colors[i], linewidth=2,
                label=f"comp{i}: μ={fit['means'][i]:.3f}, "
                      f"σ={fit['stds'][i]:.3f}, w={fit['weights'][i]:.2f}")

    if threshold is not None:
        ax.axvline(threshold, linestyle="--", color="crimson",
                   linewidth=2, label=f"crossover = {threshold:.3f}")

    status = "VALID fit" if valid else f"INVALID fit ({reason})"
    ax.set_title(f"{title}\n{status}", fontsize=11)
    ax.set_xlabel("fraction_modified")
    ax.set_ylabel("density")
    ax.set_xlim(0, 1)
    ax.legend(fontsize=8, loc="upper right")
    plt.tight_layout()
    plt.savefig(out_path, dpi=120)
    plt.close(fig)


def plot_truth_colored(scores, labels, bigaussian_threshold, title, out_path):
    """
    Paper Fig 4g style: stacked histogram colored by ground truth, with
    the oracle threshold and the bi-Gaussian crossover both overlaid.

    - Blue bars = sites whose ground truth is 0 (unmethylated)
    - Red bars  = sites whose ground truth is 1 (methylated)
    - Dashed black line = oracle threshold (max-accuracy given truth)
    - Dashed red line = bi-Gaussian crossover (if valid)
    """
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels, dtype=int)

    scores_0 = scores[labels == 0]
    scores_1 = scores[labels == 1]

    oracle_thresh, oracle_acc = find_oracle_threshold(scores, labels)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    bins = np.linspace(0, 1, 41)
    ax.hist(scores_0, bins=bins, alpha=0.6, color="steelblue",
            edgecolor="white", label=f"Epi-bit 0  (n={len(scores_0)})")
    ax.hist(scores_1, bins=bins, alpha=0.6, color="crimson",
            edgecolor="white", label=f"Epi-bit 1  (n={len(scores_1)})")

    if oracle_thresh is not None:
        ax.axvline(oracle_thresh, linestyle="--", color="black",
                   linewidth=2,
                   label=f"oracle = {oracle_thresh:.3f}  "
                         f"(acc={oracle_acc:.3f})")

    if bigaussian_threshold is not None:
        ax.axvline(bigaussian_threshold, linestyle="--", color="darkred",
                   linewidth=2,
                   label=f"bi-Gaussian = {bigaussian_threshold:.3f}")

    ax.set_title(f"{title}  (ground-truth coloured)", fontsize=11)
    ax.set_xlabel("fraction_modified")
    ax.set_ylabel("count")
    ax.set_xlim(0, 1)
    ax.legend(fontsize=8, loc="upper right")
    plt.tight_layout()
    plt.savefig(out_path, dpi=120)
    plt.close(fig)

    return oracle_thresh, oracle_acc


# ============================================================
# MAIN
# ============================================================

def main(auc_results_dir: str):
    auc_dir   = Path(auc_results_dir)
    bed_cache = auc_dir / "bed_cache"
    if not bed_cache.is_dir():
        print(f"ERROR: bed_cache/ not found under {auc_dir}"); sys.exit(1)

    out_dir   = auc_dir / "bigaussian"
    plot_dir  = out_dir / "diagnostic_plots"
    truth_dir = out_dir / "truth_colored_plots"
    out_dir.mkdir(exist_ok=True)
    plot_dir.mkdir(exist_ok=True)
    truth_dir.mkdir(exist_ok=True)

    # ---- Step 1: load all cached BEDs into a long-form DataFrame ----
    # One row per (sample_stem, position): fraction_modified + context
    rows = []
    for bed_path in sorted(bed_cache.glob("*.bed")):
        parsed = parse_bed_filename(bed_path.name)
        if parsed is None:
            print(f"[SKIP] Unparseable BED filename: {bed_path.name}")
            continue
        sample_type, basecaller, model, qscore, sample_stem = parsed

        gt = GROUND_TRUTH[sample_type]
        fractions = extract_fractions_from_bed(bed_path, set(gt.keys()))
        if not fractions:
            print(f"[WARN] {bed_path.name}: no CpG fractions extracted")
            continue

        for pos, frac in fractions.items():
            rows.append({
                "basecaller":  basecaller,
                "model":       model,
                "qscore":      qscore,
                "sample_type": sample_type,
                "sample_stem": sample_stem,
                "position":    pos,
                "fraction":    frac,
                "truth":       gt[pos],
            })

    if not rows:
        print("No data loaded from BED cache — nothing to do.")
        return

    df = pd.DataFrame(rows)
    print(f"Loaded {len(df)} rows across {df['sample_stem'].nunique()} samples "
          f"and {df.groupby(['basecaller','model','qscore']).ngroups} conditions.")

    # ---- Step 2: per (basecaller, model, qscore, sample_type) bi-Gaussian fit ----
    threshold_records = []
    per_sample_records = []

    group_keys = ["basecaller", "model", "qscore", "sample_type"]
    for (basecaller, model, qscore, sample_type), gdf in df.groupby(group_keys):
        values = gdf["fraction"].values.astype(float)

        # Need enough data to even attempt a fit
        if len(values) < 20 or len(np.unique(values)) < 10:
            threshold_records.append({
                "basecaller": basecaller, "model": model, "qscore": qscore,
                "sample_type": sample_type,
                "n_values": len(values),
                "fit_valid": False, "fit_reason": "too few unique values",
                "threshold": None,
                "mean0": None, "std0": None, "weight0": None,
                "mean1": None, "std1": None, "weight1": None,
                "n_samples": gdf["sample_stem"].nunique(),
                "mean_accuracy": None, "std_accuracy": None,
                "mean_1to0_flips": None, "mean_0to1_flips": None,
                "oracle_threshold": None, "oracle_accuracy": None,
                "threshold_gap": None,
            })
            continue

        fit = fit_bigaussian(values)
        valid, reason = check_fit_validity(fit)
        threshold = find_crossover(fit) if valid else None
        if valid and threshold is None:
            valid, reason = False, "no crossover between means"

        # Apply threshold to each replicate (only if valid)
        per_sample_accuracies = []
        per_sample_1to0 = []
        per_sample_0to1 = []
        if valid:
            for sample_stem, sdf in gdf.groupby("sample_stem"):
                fractions = dict(zip(sdf["position"], sdf["fraction"]))
                result = decode_with_threshold(
                    fractions, GROUND_TRUTH[sample_type], threshold
                )
                per_sample_accuracies.append(result["accuracy"])
                per_sample_1to0.append(result["n_1to0_flips"])
                per_sample_0to1.append(result["n_0to1_flips"])
                per_sample_records.append({
                    "basecaller": basecaller, "model": model, "qscore": qscore,
                    "sample_type": sample_type, "sample_stem": sample_stem,
                    "threshold": round(threshold, 4),
                    "n_sites":      result["n_sites"],
                    "n_matches":    result["n_matches"],
                    "n_1to0_flips": result["n_1to0_flips"],
                    "n_0to1_flips": result["n_0to1_flips"],
                    "accuracy":     round(result["accuracy"], 4),
                })

        threshold_records.append({
            "basecaller": basecaller, "model": model, "qscore": qscore,
            "sample_type": sample_type,
            "n_values": len(values),
            "n_samples": gdf["sample_stem"].nunique(),
            "fit_valid":  valid,
            "fit_reason": reason if not valid else "",
            "threshold":  round(threshold, 4) if threshold is not None else None,
            "mean0":   round(float(fit["means"][0]), 4),
            "std0":    round(float(fit["stds"][0]),  4),
            "weight0": round(float(fit["weights"][0]), 4),
            "mean1":   round(float(fit["means"][1]), 4),
            "std1":    round(float(fit["stds"][1]),  4),
            "weight1": round(float(fit["weights"][1]), 4),
            "mean_accuracy": round(float(np.mean(per_sample_accuracies)), 4) if per_sample_accuracies else None,
            "std_accuracy":  round(float(np.std(per_sample_accuracies,  ddof=1)), 4) if len(per_sample_accuracies) > 1 else None,
            "mean_1to0_flips": round(float(np.mean(per_sample_1to0)), 2) if per_sample_1to0 else None,
            "mean_0to1_flips": round(float(np.mean(per_sample_0to1)), 2) if per_sample_0to1 else None,
        })

        # Diagnostic plot (unlabeled)
        safe_name = f"{basecaller}_{model}_q{qscore}_{sample_type}.png"
        plot_diagnostic(
            values, fit, threshold, valid, reason,
            title=f"{basecaller} {model} Q{qscore} — {sample_type}  (n={len(values)})",
            out_path=plot_dir / safe_name,
        )

        # Truth-coloured plot (Fig 4g style) — uses ground truth labels
        truth_labels = gdf["truth"].values.astype(int)
        oracle_thresh, oracle_acc = plot_truth_colored(
            values, truth_labels, threshold,
            title=f"{basecaller} {model} Q{qscore} — {sample_type}",
            out_path=truth_dir / safe_name,
        )

        # Record oracle results on the last-appended threshold_record
        if threshold_records:
            threshold_records[-1]["oracle_threshold"] = (
                round(float(oracle_thresh), 4) if oracle_thresh is not None else None
            )
            threshold_records[-1]["oracle_accuracy"] = (
                round(float(oracle_acc), 4) if oracle_thresh is not None else None
            )
            # Gap between the unsupervised crossover and the supervised ceiling.
            # Large gap -> bi-Gaussian is leaving accuracy on the table.
            if threshold is not None and oracle_thresh is not None:
                threshold_records[-1]["threshold_gap"] = round(
                    float(threshold) - float(oracle_thresh), 4
                )
            else:
                threshold_records[-1]["threshold_gap"] = None

    # ---- Step 3: pooled-substrate threshold (realism check) ----
    pooled_records = []
    for (basecaller, model, qscore), gdf in df.groupby(
        ["basecaller", "model", "qscore"]
    ):
        values = gdf["fraction"].values.astype(float)
        if len(values) < 20 or len(np.unique(values)) < 10:
            continue

        fit = fit_bigaussian(values, random_state=1)
        valid, reason = check_fit_validity(fit)
        threshold = find_crossover(fit) if valid else None
        if valid and threshold is None:
            valid, reason = False, "no crossover between means"

        pooled_accs = []
        if valid:
            for sample_stem, sdf in gdf.groupby("sample_stem"):
                sample_type = sdf["sample_type"].iloc[0]
                fractions = dict(zip(sdf["position"], sdf["fraction"]))
                result = decode_with_threshold(
                    fractions, GROUND_TRUTH[sample_type], threshold
                )
                pooled_accs.append(result["accuracy"])

        pooled_records.append({
            "basecaller": basecaller, "model": model, "qscore": qscore,
            "pooled_fit_valid": valid,
            "pooled_fit_reason": reason if not valid else "",
            "pooled_threshold": round(threshold, 4) if threshold is not None else None,
            "pooled_mean_accuracy": round(float(np.mean(pooled_accs)), 4) if pooled_accs else None,
        })

    # ---- Step 4: save outputs ----
    thresh_df = pd.DataFrame(threshold_records).sort_values(
        ["basecaller", "model", "qscore", "sample_type"]
    )
    thresh_path = out_dir / "bigaussian_thresholds.csv"
    thresh_df.to_csv(thresh_path, index=False)
    print(f"\nPer-condition thresholds saved -> {thresh_path}")

    if per_sample_records:
        per_sample_df = pd.DataFrame(per_sample_records).sort_values(
            ["basecaller", "model", "qscore", "sample_type", "sample_stem"]
        )
        per_sample_path = out_dir / "bigaussian_per_sample_decoding.csv"
        per_sample_df.to_csv(per_sample_path, index=False)
        print(f"Per-sample decoding saved    -> {per_sample_path}")

    pooled_df = pd.DataFrame(pooled_records).sort_values(
        ["basecaller", "model", "qscore"]
    )
    pooled_path = out_dir / "bigaussian_pooled_threshold.csv"
    pooled_df.to_csv(pooled_path, index=False)
    print(f"Pooled-substrate thresholds  -> {pooled_path}")

    # ---- Step 5: ranking ----
    # Rank by mean_accuracy averaged across AHEAD and EpiC
    ranking = (
        thresh_df[thresh_df["fit_valid"] == True]
        .groupby(["basecaller", "model", "qscore"])
        .agg(
            mean_accuracy_across_substrates=("mean_accuracy", "mean"),
            ahead_accuracy=("mean_accuracy",
                            lambda s: s[thresh_df.loc[s.index, "sample_type"] == "AHEAD"].mean()),
            epic_accuracy=("mean_accuracy",
                            lambda s: s[thresh_df.loc[s.index, "sample_type"] == "EpiC"].mean()),
            ahead_threshold=("threshold",
                            lambda s: s[thresh_df.loc[s.index, "sample_type"] == "AHEAD"].mean()),
            epic_threshold=("threshold",
                            lambda s: s[thresh_df.loc[s.index, "sample_type"] == "EpiC"].mean()),
            n_valid_substrates=("fit_valid", "sum"),
        )
        .reset_index()
        .sort_values("mean_accuracy_across_substrates", ascending=False)
        .reset_index(drop=True)
    )
    ranking["rank"] = ranking.index + 1

    ranking = ranking.merge(pooled_df, on=["basecaller", "model", "qscore"], how="left")

    ranking_path = out_dir / "bigaussian_ranking.csv"
    ranking.to_csv(ranking_path, index=False)
    print(f"Ranking saved                -> {ranking_path}")

    # ---- Print summary ----
    print("\n" + "="*80)
    print("  TOP 10 BASECALLER x MODEL x Q-SCORE BY BI-GAUSSIAN DECODING ACCURACY")
    print("="*80)
    cols = ["rank", "basecaller", "model", "qscore",
            "mean_accuracy_across_substrates",
            "ahead_accuracy", "ahead_threshold",
            "epic_accuracy",  "epic_threshold",
            "pooled_mean_accuracy", "pooled_threshold",
            "n_valid_substrates"]
    cols = [c for c in cols if c in ranking.columns]
    print(ranking[cols].head(10).to_string(index=False))

    print("\n" + "="*80)
    print(f"  Diagnostic plots written to: {plot_dir}")
    print("="*80)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python bigaussian_threshold_analysis.py <auc_results_dir>")
        print("       where <auc_results_dir> contains bed_cache/ from the AUC analysis.")
        sys.exit(1)
    main(sys.argv[1])
