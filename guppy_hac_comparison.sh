#!/bin/bash
set -euo pipefail

# ─────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────

FAST5_DIR="/DATA4/amishi/5Apr/fast5_data"
BASE_OUTPUT_DIR="/DATA4/amishi/16Apr/guppy_hac_comparison"

# Reference genome (FASTA, FASTQ, or .mmi)
REFERENCE="/DATA4/amishi/reference.fasta"

GUPPY_BIN="/DATA4/amishi/ont-guppy/bin/guppy_basecaller"

# Models
BASE_MODEL="/DATA4/amishi/ont-guppy/data/dna_r10.4.1_e8.2_400bps_5khz_hac_prom.cfg"
METHYL_MODEL="/DATA4/amishi/ont-guppy/data/dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac_prom.cfg"
# ^^ verify this filename matches your actual methyl config; edit if different.

# Q-score thresholds
Q_SCORES=(0 5 10 15 20)

# Threads for sorting / alignment
THREADS=8

# ─────────────────────────────────────────────
# SETUP
# ─────────────────────────────────────────────

mkdir -p "$BASE_OUTPUT_DIR"
LOG="$BASE_OUTPUT_DIR/run_log.txt"
echo "Run started: $(date)" | tee "$LOG"

# Sanity checks
if [ ! -f "$REFERENCE" ]; then
    echo "ERROR: Reference file not found: $REFERENCE" | tee -a "$LOG"
    exit 1
fi
echo "Reference: $REFERENCE" | tee -a "$LOG"

if ! command -v minimap2 &> /dev/null; then
    echo "ERROR: minimap2 not found on PATH" | tee -a "$LOG"
    exit 1
fi

# ─────────────────────────────────────────────
# STEP 1: BASECALLING ONLY
#         Outputs: unaligned BAM (with mv tags) + aligned/sorted BAM
# ─────────────────────────────────────────────

echo "" | tee -a "$LOG"
echo "==============================" | tee -a "$LOG"
echo "  STEP 1: BASECALLING ONLY" | tee -a "$LOG"
echo "==============================" | tee -a "$LOG"

for Q in "${Q_SCORES[@]}"; do

    FOLDER_NAME="guppy_hac_basecall_q${Q}"
    OUTPUT_DIR="$BASE_OUTPUT_DIR/$FOLDER_NAME"
    UNALIGNED_DIR="$OUTPUT_DIR/unaligned"
    ALIGNED_DIR="$OUTPUT_DIR/aligned"
    mkdir -p "$UNALIGNED_DIR" "$ALIGNED_DIR"

    echo "" | tee -a "$LOG"
    echo "-> Basecalling | Q${Q}" | tee -a "$LOG"

    for sample_dir in "$FAST5_DIR"/*/; do
        [ -d "$sample_dir" ] || continue
        sample=$(basename "$sample_dir")

        UNALIGNED_BAM="$UNALIGNED_DIR/${sample}_guppy_hac_basecall_q${Q}.unaligned.bam"
        ALIGNED_BAM="$ALIGNED_DIR/${sample}_guppy_hac_basecall_q${Q}.aligned.sorted.bam"

        # ─── Unaligned BAM (with move tables) ───
        if [ -f "$UNALIGNED_BAM" ]; then
            echo "  [SKIP unaligned] $sample already processed" | tee -a "$LOG"
        else
            echo "  Basecalling: $sample (Q${Q})" | tee -a "$LOG"

            GUPPY_OUT="$OUTPUT_DIR/${sample}_guppy_raw"
            mkdir -p "$GUPPY_OUT"

            "$GUPPY_BIN" \
                -i "$sample_dir" \
                -s "$GUPPY_OUT" \
                -c "$BASE_MODEL" \
                --bam_out \
                --moves_out \
                --min_qscore "$Q" \
                --device cuda:0 \
                2>> "$LOG"

            # Merge pass BAMs into a single unaligned, unsorted BAM
            if ls "$GUPPY_OUT"/pass/*.bam 1>/dev/null 2>&1; then
                samtools merge -f "$UNALIGNED_BAM" "$GUPPY_OUT"/pass/*.bam

                READ_COUNT=$(samtools view -c "$UNALIGNED_BAM")
                echo "  ✓ $READ_COUNT reads → $(basename $UNALIGNED_BAM)" | tee -a "$LOG"

                # Verify mv (move table) tags are present
                MV_CHECK=$(samtools view "$UNALIGNED_BAM" | head -5 | grep -c "mv:B:" || true)
                if [ "$MV_CHECK" -gt 0 ]; then
                    echo "  ✓ Move table (mv) tags confirmed present" | tee -a "$LOG"
                else
                    echo "  ✗ WARNING: mv tags not found in $sample" | tee -a "$LOG"
                fi
            else
                echo "  ✗ WARNING: No pass BAMs for $sample at Q${Q}" | tee -a "$LOG"
                rm -rf "$GUPPY_OUT"
                continue
            fi

            rm -rf "$GUPPY_OUT"
        fi

        # ─── Aligned + coordinate-sorted BAM ───
        if [ -f "$ALIGNED_BAM" ]; then
            echo "  [SKIP aligned] $sample already aligned" | tee -a "$LOG"
            continue
        fi

        echo "  Aligning to reference: $sample" | tee -a "$LOG"

        # minimap2 -y preserves all tags (mv, MM, ML) from the unaligned BAM
        samtools fastq -T "*" "$UNALIGNED_BAM" 2>> "$LOG" \
        | minimap2 -ax lr:hq -y -t "$THREADS" "$REFERENCE" - 2>> "$LOG" \
        | samtools sort -@ "$THREADS" -o "$ALIGNED_BAM"

        samtools index "$ALIGNED_BAM"

        MAPPED=$(samtools view -c -F 4 "$ALIGNED_BAM")
        TOTAL=$(samtools view -c "$ALIGNED_BAM")
        echo "  ✓ $MAPPED / $TOTAL mapped → $(basename $ALIGNED_BAM)" | tee -a "$LOG"

    done
done

# ─────────────────────────────────────────────
# STEP 2: METHYLATION CALLING (5mCG + 5hmCG)
#         Outputs: unaligned BAM (with mv + MM/ML tags) + aligned/sorted BAM
# ─────────────────────────────────────────────

echo "" | tee -a "$LOG"
echo "==============================" | tee -a "$LOG"
echo "  STEP 2: METHYLATION CALLING" | tee -a "$LOG"
echo "==============================" | tee -a "$LOG"

for Q in "${Q_SCORES[@]}"; do

    FOLDER_NAME="guppy_hac_methyl_q${Q}"
    OUTPUT_DIR="$BASE_OUTPUT_DIR/$FOLDER_NAME"
    UNALIGNED_DIR="$OUTPUT_DIR/unaligned"
    ALIGNED_DIR="$OUTPUT_DIR/aligned"
    mkdir -p "$UNALIGNED_DIR" "$ALIGNED_DIR"

    echo "" | tee -a "$LOG"
    echo "-> Methylation calling | Q${Q}" | tee -a "$LOG"

    for sample_dir in "$FAST5_DIR"/*/; do
        [ -d "$sample_dir" ] || continue
        sample=$(basename "$sample_dir")

        UNALIGNED_BAM="$UNALIGNED_DIR/${sample}_guppy_hac_methyl_q${Q}.unaligned.bam"
        ALIGNED_BAM="$ALIGNED_DIR/${sample}_guppy_hac_methyl_q${Q}.aligned.sorted.bam"

        # ─── Unaligned BAM (with move tables + MM/ML) ───
        if [ -f "$UNALIGNED_BAM" ]; then
            echo "  [SKIP unaligned] $sample already processed" | tee -a "$LOG"
        else
            echo "  Basecalling + mod calling: $sample (Q${Q})" | tee -a "$LOG"

            GUPPY_OUT="$OUTPUT_DIR/${sample}_guppy_raw"
            mkdir -p "$GUPPY_OUT"

            "$GUPPY_BIN" \
                -i "$sample_dir" \
                -s "$GUPPY_OUT" \
                -c "$METHYL_MODEL" \
                --bam_out \
                --moves_out \
                --min_qscore "$Q" \
                --device cuda:0 \
                2>> "$LOG"

            if ls "$GUPPY_OUT"/pass/*.bam 1>/dev/null 2>&1; then
                samtools merge -f "$UNALIGNED_BAM" "$GUPPY_OUT"/pass/*.bam

                READ_COUNT=$(samtools view -c "$UNALIGNED_BAM")
                echo "  ✓ $READ_COUNT reads → $(basename $UNALIGNED_BAM)" | tee -a "$LOG"

                # Verify MM/ML tags
                MM_CHECK=$(samtools view "$UNALIGNED_BAM" | head -5 | grep -c "MM:Z:" || true)
                if [ "$MM_CHECK" -gt 0 ]; then
                    echo "  ✓ MM/ML methylation tags confirmed present" | tee -a "$LOG"
                else
                    echo "  ✗ WARNING: MM/ML tags not found in $sample" | tee -a "$LOG"
                fi

                # Verify mv (move table) tags
                MV_CHECK=$(samtools view "$UNALIGNED_BAM" | head -5 | grep -c "mv:B:" || true)
                if [ "$MV_CHECK" -gt 0 ]; then
                    echo "  ✓ Move table (mv) tags confirmed present" | tee -a "$LOG"
                else
                    echo "  ✗ WARNING: mv tags not found in $sample" | tee -a "$LOG"
                fi
            else
                echo "  ✗ WARNING: No pass BAMs for $sample at Q${Q}" | tee -a "$LOG"
                rm -rf "$GUPPY_OUT"
                continue
            fi

            rm -rf "$GUPPY_OUT"
        fi

        # ─── Aligned + coordinate-sorted BAM ───
        if [ -f "$ALIGNED_BAM" ]; then
            echo "  [SKIP aligned] $sample already aligned" | tee -a "$LOG"
            continue
        fi

        echo "  Aligning to reference: $sample" | tee -a "$LOG"

        samtools fastq -T "*" "$UNALIGNED_BAM" 2>> "$LOG" \
        | minimap2 -ax lr:hq -y -t "$THREADS" "$REFERENCE" - 2>> "$LOG" \
        | samtools sort -@ "$THREADS" -o "$ALIGNED_BAM"

        samtools index "$ALIGNED_BAM"

        MAPPED=$(samtools view -c -F 4 "$ALIGNED_BAM")
        TOTAL=$(samtools view -c "$ALIGNED_BAM")
        echo "  ✓ $MAPPED / $TOTAL mapped → $(basename $ALIGNED_BAM)" | tee -a "$LOG"

    done
done

# ─────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────

echo "" | tee -a "$LOG"
echo "==============================" | tee -a "$LOG"
echo "  FINAL SUMMARY" | tee -a "$LOG"
echo "==============================" | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "--- Basecall outputs (unaligned) ---" | tee -a "$LOG"
for bam in "$BASE_OUTPUT_DIR"/guppy_hac_basecall_q*/unaligned/*.bam; do
    [ -f "$bam" ] || continue
    count=$(samtools view -c "$bam")
    printf "  %-10s reads → %s\n" "$count" "$(basename $bam)" | tee -a "$LOG"
done

echo "" | tee -a "$LOG"
echo "--- Basecall outputs (aligned) ---" | tee -a "$LOG"
for bam in "$BASE_OUTPUT_DIR"/guppy_hac_basecall_q*/aligned/*.bam; do
    [ -f "$bam" ] || continue
    mapped=$(samtools view -c -F 4 "$bam")
    total=$(samtools view -c "$bam")
    printf "  %-10s mapped / %-10s total → %s\n" "$mapped" "$total" "$(basename $bam)" | tee -a "$LOG"
done

echo "" | tee -a "$LOG"
echo "--- Methylation outputs (unaligned) ---" | tee -a "$LOG"
for bam in "$BASE_OUTPUT_DIR"/guppy_hac_methyl_q*/unaligned/*.bam; do
    [ -f "$bam" ] || continue
    count=$(samtools view -c "$bam")
    printf "  %-10s reads → %s\n" "$count" "$(basename $bam)" | tee -a "$LOG"
done

echo "" | tee -a "$LOG"
echo "--- Methylation outputs (aligned) ---" | tee -a "$LOG"
for bam in "$BASE_OUTPUT_DIR"/guppy_hac_methyl_q*/aligned/*.bam; do
    [ -f "$bam" ] || continue
    mapped=$(samtools view -c -F 4 "$bam")
    total=$(samtools view -c "$bam")
    printf "  %-10s mapped / %-10s total → %s\n" "$mapped" "$total" "$(basename $bam)" | tee -a "$LOG"
done

echo "" | tee -a "$LOG"
echo "All done: $(date)" | tee -a "$LOG"
