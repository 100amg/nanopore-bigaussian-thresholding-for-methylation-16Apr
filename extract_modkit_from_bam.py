'''
#!/usr/bin/env python3

import re
import subprocess
from pathlib import Path

# ============================================================
# CONFIG
# ============================================================

QFOLDER_RE = re.compile(r"^(dorado|guppy)_(\w+?)_methyl_q(\d+)$")

def parse_qfolder(folder_name):
    m = QFOLDER_RE.match(folder_name)
    if m:
        return m.group(1), m.group(2), int(m.group(3))
    return None, None, None

# ============================================================
# RUN MODKIT EXTRACT
# ============================================================

def run_modkit_extract(bam_path: Path, out_tsv: Path):
    cmd = [
        "modkit", "extract",
        "--cpg",
        "--mapped-only",
        str(bam_path),
        str(out_tsv)
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"    [ERROR] modkit failed:\n{result.stderr[:300]}")
        return False

    return True

# ============================================================
# MAIN
# ============================================================

def main(bam_root_dir):

    bam_root = Path(bam_root_dir)
    output_root = Path("modkit_extract_results")
    output_root.mkdir(exist_ok=True)

    comp_folders = sorted(
        f for f in bam_root.iterdir()
        if f.is_dir() and f.name.endswith("_comparison")
    )

    for comp_folder in comp_folders:
        print(f"\n{'#'*60}")
        print(f"Processing: {comp_folder.name}")
        print(f"{'#'*60}")

        qfolders = sorted(
            f for f in comp_folder.iterdir()
            if f.is_dir() and parse_qfolder(f.name)[0] is not None
        )

        for qfolder in qfolders:
            basecaller, model, qscore = parse_qfolder(qfolder.name)

            print(f"\n  {basecaller} | {model} | Q{qscore}")

            aligned_dir = qfolder / "aligned"

            if aligned_dir.is_dir():
                bam_files = sorted(aligned_dir.glob("*.aligned.sorted.bam"))
            else:
                bam_files = sorted(qfolder.glob("*_sorted.bam"))

            if not bam_files:
                print("  No BAM files found")
                continue

            # Create mirrored output folder
            out_dir = output_root / comp_folder.name / qfolder.name
            out_dir.mkdir(parents=True, exist_ok=True)

            for bam in bam_files:

                sample_name = bam.name.replace(".aligned.sorted.bam", "") \
                                      .replace("_sorted.bam", "") \
                                      .replace(".bam", "")

                out_tsv = out_dir / f"{sample_name}.tsv"

                print(f"\n    Sample: {sample_name}")

                if out_tsv.exists():
                    print("    [SKIP] Already extracted")
                    continue

                if not Path(str(bam) + ".bai").exists():
                    print("    [WARN] Missing index (.bai), skipping")
                    continue

                print("    Running modkit extract...")

                success = run_modkit_extract(bam, out_tsv)

                if success:
                    print(f"    Saved → {out_tsv}")
                else:
                    print("    Failed")

    print("\nDONE.")

# ============================================================

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python modkit_extract_all.py <bam_root_dir>")
        exit(1)

    main(sys.argv[1])
    '''

    #!/usr/bin/env python3

import re
import subprocess
from pathlib import Path

# ============================================================
# CONFIG
# ============================================================

QFOLDER_RE = re.compile(r"^(dorado|guppy)_(\w+?)_methyl_q(\d+)$")

def parse_qfolder(name):
    m = QFOLDER_RE.match(name)
    if m:
        return m.group(1), m.group(2), int(m.group(3))
    return None, None, None

# ============================================================
# MODKIT RUNNER
# ============================================================

def run_modkit_extract(bam_path, out_tsv):
    cmd = [
        "modkit", "extract", "full",
        str(bam_path),
        str(out_tsv)
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print("      ❌ ERROR running modkit")
        print(result.stderr[:300])
        return False

    return True

# ============================================================
# MAIN
# ============================================================

def main(bam_root_dir):

    bam_root = Path(bam_root_dir)

    if not bam_root.exists():
        print(f"ERROR: {bam_root} not found")
        return

    output_root = Path("modkit_tsv_outputs")
    output_root.mkdir(exist_ok=True)

    total_bams = 0
    total_success = 0

    # ========================================================
    # LOOP OVER COMPARISON FOLDERS
    # ========================================================

    comp_folders = sorted(
        f for f in bam_root.iterdir()
        if f.is_dir() and f.name.endswith("_comparison")
    )

    print(f"\nFound {len(comp_folders)} comparison folders")

    for comp in comp_folders:
        print(f"\n{'='*60}")
        print(f"Processing: {comp.name}")
        print(f"{'='*60}")

        qfolders = sorted(
            f for f in comp.iterdir()
            if f.is_dir() and parse_qfolder(f.name)[0] is not None
        )

        print(f"  Found {len(qfolders)} q-score folders")

        for qf in qfolders:
            basecaller, model, q = parse_qfolder(qf.name)

            print(f"\n  → {basecaller} | {model} | Q{q}")

            aligned_dir = qf / "aligned"

            if aligned_dir.exists():
                bam_files = sorted(aligned_dir.glob("*.bam"))
            else:
                bam_files = sorted(qf.glob("*.bam"))

            print(f"    Found {len(bam_files)} BAM files")

            if not bam_files:
                continue

            # output folder mirrors structure
            out_dir = output_root / comp.name / qf.name
            out_dir.mkdir(parents=True, exist_ok=True)

            for bam in bam_files:

                if bam.name.endswith(".bai"):
                    continue

                total_bams += 1

                sample_name = bam.name.replace(".aligned.sorted.bam", "") \
                                      .replace("_sorted.bam", "") \
                                      .replace(".bam", "")

                out_tsv = out_dir / f"{sample_name}.tsv"

                print(f"\n    Sample: {sample_name}")

                if out_tsv.exists():
                    print("      ⚠ Already exists → skipping")
                    continue

                if not Path(str(bam) + ".bai").exists():
                    print("      ⚠ Missing .bai → skipping")
                    continue

                print("      Running modkit extract...")

                success = run_modkit_extract(bam, out_tsv)

                if success:
                    print(f"      ✅ Saved: {out_tsv}")
                    total_success += 1
                else:
                    print("      ❌ Failed")

    print("\n" + "#"*60)
    print(f"TOTAL BAM FILES SEEN: {total_bams}")
    print(f"TOTAL TSV CREATED:   {total_success}")
    print("#"*60)


# ============================================================

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python modkit_extract_clean.py <bam_root_dir>")
        exit(1)

    main(sys.argv[1])