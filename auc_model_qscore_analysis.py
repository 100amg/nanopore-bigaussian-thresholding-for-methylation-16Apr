#!/usr/bin/env python3
"""
AUC-based model x Q-score comparison for T960 methylation calling.

Supports both Dorado and Guppy aligned output folders.

Expected directory structure (pass the comparison-level parent as bam_root_dir):
    <bam_root_dir>/
        dorado_fast_comparison/
            dorado_fast_methyl_q0/
                aligned/
                    AHEAD_10_11_25_fast_methyl_q0.aligned.sorted.bam
                    AHEAD_10_11_25_fast_methyl_q0.aligned.sorted.bam.bai
                unaligned/                    (ignored)
            dorado_fast_basecall_q0/          (ignored — no MM/ML tags)
                ...
        dorado_hac_comparison/
        dorado_sup_comparison/
        guppy_comparison/
            guppy_hac_methyl_q0/
                aligned/
                    AHEAD_10_11_25_guppy_hac_methyl_q0.aligned.sorted.bam
                    ...

Also supports the legacy flat structure (*_sorted.bam directly in the
q-score folder) for backward compatibility.

Usage:
    python auc_model_qscore_analysis.py <bam_root_dir> <reference_fasta>

Outputs (written to ./auc_results/):
    per_sample_auc.csv        AUC for every individual sample
    model_qscore_summary.csv  Mean AUC per basecaller x model x Q-score
    top5_combinations.csv     Top 5 ranked combinations
"""

import re
import sys
import subprocess
from pathlib import Path

import pandas as pd
from sklearn.metrics import roc_auc_score

# ============================================================
# GROUND TRUTH CONFIGURATION
# ============================================================

T960_CPG_POSITIONS = list(range(60, 901, 24))  # [60, 84, ..., 900], 36 sites

# AHEAD: 5 letters x 7 bits = 35 bits + 1 trailing padding zero = 36 bits total.
# All 36 CpG sites are present on T960, so we include all 36 in the AUC.
AHEAD_BINARY = "100000110010001000101100000110001000"  # 36 chars (35 encoded + 1 padding 0)
AHEAD_GROUND_TRUTH = {
    pos: int(bit)
    for pos, bit in zip(T960_CPG_POSITIONS[:36], AHEAD_BINARY)
}

# EpiC: 4 letters x 8 bits = 32 bits + 4 trailing padding zeros = 36 bits.
EPIC_BINARY = "010001010111000001101001010000110000"
EPIC_GROUND_TRUTH = {
    pos: int(bit)
    for pos, bit in zip(T960_CPG_POSITIONS[:36], EPIC_BINARY)
}

GROUND_TRUTH = {
    "AHEAD": AHEAD_GROUND_TRUTH,
    "EpiC":  EPIC_GROUND_TRUTH,
}

assert len(AHEAD_GROUND_TRUTH) == 36
assert len(EPIC_GROUND_TRUTH)  == 36

# ============================================================
# FOLDER PARSING
# ============================================================

# Matches methyl q-score folders for both Dorado and Guppy:
#   dorado_fast_methyl_q10  ->  basecaller=dorado, model=fast,     qscore=10
#   dorado_hac_methyl_q5   ->  basecaller=dorado, model=hac,      qscore=5
#   guppy_hac_methyl_q0    ->  basecaller=guppy,  model=hac,      qscore=0
QFOLDER_RE = re.compile(r"^(dorado|guppy)_(\w+?)_methyl_q(\d+)$")

def parse_qfolder(folder_name: str):
    """
    Returns (basecaller, model, qscore) or (None, None, None).
    e.g. 'guppy_hac_methyl_q10'      -> ('guppy',  'hac',  10)
         'dorado_fast_methyl_q5'      -> ('dorado', 'fast',  5)
         'dorado_sup_methyl_q15'      -> ('dorado', 'sup',  15)
    """
    m = QFOLDER_RE.match(folder_name)
    if m:
        return m.group(1), m.group(2), int(m.group(3))
    return None, None, None


def get_sample_type(filename: str):
    if filename.startswith("AHEAD"):
        return "AHEAD"
    if filename.startswith("EpiC"):
        return "EpiC"
    return None

# ============================================================
# MODKIT PILEUP
# ============================================================

def run_modkit_pileup(bam_path: Path, ref_fasta: Path, out_bed: Path) -> bool:
    cmd = [
        "modkit", "pileup",
        "--cpg",
        "--mod-thresholds", "C:0.0",
        "--ref", str(ref_fasta),
        str(bam_path),
        str(out_bed),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"    [modkit ERROR] {result.stderr[:400]}")
        return False
    return True

# ============================================================
# BED PARSING
# ============================================================

def extract_fractions(bed_path: Path, target_positions: set) -> dict:
    """
    Returns {chromEnd_position: fraction_modified} for T960 CpG sites.

    modkit pileup BED columns (0-indexed):
        2   chromEnd          <- position key
        3   mod_code          'm' = 5mC
        10  fraction_modified (0.0 to 1.0)
    """
    fractions = {}
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) < 11:
                continue
            if cols[3] != "m":
                continue
            try:
                pos  = int(cols[2])
                frac = float(cols[10])
            except ValueError:
                continue
            if pos in target_positions:
                fractions[pos] = frac
    return fractions

# ============================================================
# AUC COMPUTATION
# ============================================================

def compute_metrics(fractions: dict, ground_truth: dict):
    """
    Computes AUC and accuracy metrics for one sample.

    Threshold for binary calling = mean of all fraction_modified values
    at the T960 CpG sites. This matches the existing pipeline (get_stats.py)
    which uses the per-sample mean as its threshold, correctly handling
    background DNMT1 methylation that pushes all sites above 0.5.

        fraction_modified > mean  ->  called methylated (1)
        fraction_modified <= mean ->  called unmethylated (0)

    Returns a dict with keys:
        auc, n_sites, threshold_used,
        n_matches, n_1to0_flips, n_0to1_flips, accuracy
    """
    positions = sorted(ground_truth.keys())
    scores, labels = [], []
    for pos in positions:
        if pos in fractions:
            scores.append(fractions[pos])
            labels.append(ground_truth[pos])

    n = len(scores)
    result = dict(auc=None, n_sites=n, threshold_used=None,
                  n_matches=0, n_1to0_flips=0, n_0to1_flips=0, accuracy=0.0)

    if n < 2:
        return result
    if len(set(labels)) < 2:
        print(f"    [WARN] Only one class present — AUC undefined")
        return result

    result["auc"] = roc_auc_score(labels, scores)

    # Use mean of all site fractions as threshold — matches get_stats.py behaviour
    threshold = sum(scores) / len(scores)
    result["threshold_used"] = round(threshold, 4)

    n_matches = n_1to0 = n_0to1 = 0
    for score, label in zip(scores, labels):
        called = 1 if score > threshold else 0
        if called == label:
            n_matches += 1
        elif label == 1 and called == 0:
            n_1to0 += 1
        else:
            n_0to1 += 1

    result["n_matches"]    = n_matches
    result["n_1to0_flips"] = n_1to0
    result["n_0to1_flips"] = n_0to1
    result["accuracy"]     = round(n_matches / n, 4)
    return result

# ============================================================
# MAIN
# ============================================================

def main(bam_root: str, ref_fasta: str):
    bam_root  = Path(bam_root)
    ref_fasta = Path(ref_fasta)

    if not bam_root.is_dir():
        print(f"ERROR: BAM root not found: {bam_root}"); sys.exit(1)
    if not ref_fasta.is_file():
        print(f"ERROR: Reference FASTA not found: {ref_fasta}"); sys.exit(1)

    results_dir   = Path("auc_results")
    bed_cache_dir = results_dir / "bed_cache"
    bed_cache_dir.mkdir(parents=True, exist_ok=True)

    records = []

    # ---- Level 1: *_comparison folders (dorado_* and guppy_comparison) ----
    comp_folders = sorted(
        f for f in bam_root.iterdir()
        if f.is_dir() and f.name.endswith("_comparison")
    )

    if not comp_folders:
        print(f"No '*_comparison' folders found under {bam_root}")
        sys.exit(1)

    for comp_folder in comp_folders:
        print(f"\n{'#'*60}")
        print(f"  Comparison folder: {comp_folder.name}")
        print(f"{'#'*60}")

        # ---- Level 2: *_methyl_q* folders ----
        qscore_folders = sorted(
            f for f in comp_folder.iterdir()
            if f.is_dir() and parse_qfolder(f.name)[0] is not None
        )

        if not qscore_folders:
            print(f"  No methyl q-score folders found — skipping")
            continue

        for qfolder in qscore_folders:
            basecaller, model, qscore = parse_qfolder(qfolder.name)

            print(f"\n  {'='*52}")
            print(f"  Basecaller: {basecaller}  |  Model: {model}  |  Q: {qscore}")
            print(f"  {'='*52}")

            # Locate BAMs. Support both new and legacy layouts:
            #   new:    <qfolder>/aligned/*.aligned.sorted.bam
            #   legacy: <qfolder>/*_sorted.bam
            aligned_subdir = qfolder / "aligned"
            if aligned_subdir.is_dir():
                bam_files = sorted(
                    f for f in aligned_subdir.glob("*.aligned.sorted.bam")
                    if not f.name.endswith(".bai")
                )
            else:
                bam_files = sorted(
                    f for f in qfolder.glob("*_sorted.bam")
                    if not f.name.endswith(".bai")
                )

            if not bam_files:
                print(f"  No aligned BAM files found — skipping")
                print(f"  (Expected either {qfolder.name}/aligned/*.aligned.sorted.bam"
                      f" or {qfolder.name}/*_sorted.bam)")
                continue

            for bam_file in bam_files:
                # Strip suffix for the sample name / cache key:
                #   AHEAD_10_11_25_fast_methyl_q0.aligned.sorted.bam -> AHEAD_10_11_25_fast_methyl_q0
                #   AHEAD_10_11_25_fast_methyl_q0_sorted.bam         -> AHEAD_10_11_25_fast_methyl_q0
                sample_name = bam_file.name
                for suffix in (".aligned.sorted.bam", "_sorted.bam", ".bam"):
                    if sample_name.endswith(suffix):
                        sample_name = sample_name[: -len(suffix)]
                        break
                sample_type = get_sample_type(bam_file.name)

                if sample_type not in GROUND_TRUTH:
                    print(f"  [SKIP] {bam_file.name} — unrecognised sample type")
                    continue

                print(f"\n  Sample : {sample_name}  ({sample_type})")

                if not Path(str(bam_file) + ".bai").exists():
                    print(f"  [WARN] No .bai index — skipping")
                    continue

                # BED cache — keyed by unique sample name
                bed_path = bed_cache_dir / f"{sample_name}.bed"
                if not bed_path.exists():
                    print(f"  Running modkit pileup...")
                    if not run_modkit_pileup(bam_file, ref_fasta, bed_path):
                        print(f"  [ERROR] modkit failed — skipping")
                        continue
                else:
                    print(f"  Using cached BED")

                gt        = GROUND_TRUTH[sample_type]
                fractions = extract_fractions(bed_path, set(gt.keys()))
                coverage  = len(fractions)
                expected  = len(gt)

                print(f"  CpG sites found : {coverage} / {expected}")
                if coverage < expected * 0.8:
                    print(f"  [WARN] <80% of expected CpG sites — results may be unreliable")

                m = compute_metrics(fractions, gt)
                if m["auc"] is None:
                    print(f"  AUC : could not compute")
                    continue

                print(f"  AUC      : {m['auc']:.4f}  (n={m['n_sites']} sites)")
                print(f"  Threshold: {m['threshold_used']:.4f} (mean of site fractions)")
                print(f"  Accuracy : {m['accuracy']:.4f}  "
                      f"| Matches: {m['n_matches']}  "
                      f"| 1->0 flips: {m['n_1to0_flips']}  "
                      f"| 0->1 flips: {m['n_0to1_flips']}")

                records.append({
                    "basecaller":  basecaller,
                    "model":       model,
                    "qscore":      qscore,
                    "sample_type": sample_type,
                    "sample_name": sample_name,
                    "n_sites":       m["n_sites"],
                    "auc":           round(m["auc"], 6),
                    "threshold_used": m["threshold_used"],
                    "accuracy":      m["accuracy"],
                    "n_matches":   m["n_matches"],
                    "n_1to0_flips": m["n_1to0_flips"],
                    "n_0to1_flips": m["n_0to1_flips"],
                })

    # ---- Save and display results ----
    if not records:
        print("\nNo AUC results computed. Check that align_methyl_bams.sh ran successfully.")
        return

    df = pd.DataFrame(records)
    per_sample_path = results_dir / "per_sample_auc.csv"
    df.to_csv(per_sample_path, index=False)
    print(f"\nPer-sample results saved -> {per_sample_path}")

    # Summary: mean AUC and accuracy per basecaller x model x qscore
    summary = (
        df.groupby(["basecaller", "model", "qscore"])
        .agg(mean_auc=("auc", "mean"),
             std_auc=("auc", "std"),
             mean_accuracy=("accuracy", "mean"),
             mean_1to0_flips=("n_1to0_flips", "mean"),
             mean_0to1_flips=("n_0to1_flips", "mean"),
             n_samples=("auc", "count"))
        .reset_index()
    )

    # Separate AHEAD and EpiC columns
    for stype in ["AHEAD", "EpiC"]:
        sub = (
            df[df["sample_type"] == stype]
            .groupby(["basecaller", "model", "qscore"])
            .agg(**{f"{stype.lower()}_mean_auc": ("auc", "mean")})
            .reset_index()
        )
        summary = summary.merge(sub, on=["basecaller", "model", "qscore"], how="left")

    summary = summary.sort_values("mean_auc", ascending=False).reset_index(drop=True)
    summary["rank"] = summary.index + 1

    summary_path = results_dir / "model_qscore_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"Summary saved -> {summary_path}")

    top5 = summary.head(5)
    top5_path = results_dir / "top5_combinations.csv"
    top5.to_csv(top5_path, index=False)

    cols = ["rank", "basecaller", "model", "qscore", "mean_auc", "std_auc",
            "mean_accuracy", "mean_1to0_flips", "mean_0to1_flips",
            "ahead_mean_auc", "epic_mean_auc", "n_samples"]
    cols = [c for c in cols if c in summary.columns]

    print("\n" + "="*65)
    print("  TOP 5 BASECALLER x MODEL x Q-SCORE COMBINATIONS")
    print("="*65)
    print(top5[cols].to_string(index=False))

    print("\n" + "="*65)
    print("  FULL RANKING")
    print("="*65)
    print(summary[cols].to_string(index=False))

    print(f"\nTop-5 saved -> {top5_path}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python auc_model_qscore_analysis.py <bam_root_dir> <reference_fasta>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
