#!/bin/bash
set -euo pipefail

# ─────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────

POD5_DIR="/DATA4/amishi/16Apr/pod5_data"
BASE_OUTPUT_DIR="/DATA4/amishi/16Apr/dorado_sup_comparison"

# Reference genome (FASTA). Can be plain .fa/.fasta or bgzipped .fa.gz
# Change this to point to your reference file.
REFERENCE="/DATA4/amishi/reference.fasta"

DORADO_BIN="/usr/bin/dorado"

# Models
BASE_MODEL="/DATA4/amishi/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
METHYL_MODEL="/DATA4/amishi/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0_5mCG_5hmCG@v3"

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

# Sanity check: reference exists
if [ ! -f "$REFERENCE" ]; then
    echo "ERROR: Reference file not found: $REFERENCE" | tee -a "$LOG"
    exit 1
fi
echo "Reference: $REFERENCE" | tee -a "$LOG"

# ─────────────────────────────────────────────
# STEP 1: BASECALLING ONLY (no methylation)
#         Outputs: unaligned BAM (with mv tags) + aligned/sorted BAM
# ─────────────────────────────────────────────

echo "" | tee -a "$LOG"
echo "==============================" | tee -a "$LOG"
echo "  STEP 1: BASECALLING ONLY" | tee -a "$LOG"
echo "==============================" | tee -a "$LOG"

for Q in "${Q_SCORES[@]}"; do

    FOLDER_NAME="dorado_sup_basecall_q${Q}"
    OUTPUT_DIR="$BASE_OUTPUT_DIR/$FOLDER_NAME"
    UNALIGNED_DIR="$OUTPUT_DIR/unaligned"
    ALIGNED_DIR="$OUTPUT_DIR/aligned"
    mkdir -p "$UNALIGNED_DIR" "$ALIGNED_DIR"

    echo "" | tee -a "$LOG"
    echo "→ Basecalling | Q${Q}" | tee -a "$LOG"

    for pod5 in "$POD5_DIR"/*.pod5; do
        [ -e "$pod5" ] || continue
        name=$(basename "$pod5" .pod5)

        UNALIGNED_BAM="$UNALIGNED_DIR/${name}_sup_basecall_q${Q}.unaligned.bam"
        ALIGNED_BAM="$ALIGNED_DIR/${name}_sup_basecall_q${Q}.aligned.sorted.bam"

        # ─── Unaligned BAM (with move tables) ───
        if [ -f "$UNALIGNED_BAM" ]; then
            echo "  [SKIP unaligned] $name already processed" | tee -a "$LOG"
        else
            echo "  Basecalling: $name" | tee -a "$LOG"

            "$DORADO_BIN" basecaller \
                "$BASE_MODEL" \
                "$pod5" \
                --min-qscore "$Q" \
                --emit-moves \
                > "$UNALIGNED_BAM" \
                2>> "$LOG"

            READ_COUNT=$(samtools view -c "$UNALIGNED_BAM")
            echo "  ✓ $READ_COUNT reads → $(basename $UNALIGNED_BAM)" | tee -a "$LOG"

            # Verify mv (move table) tags are present
            MV_CHECK=$(samtools view "$UNALIGNED_BAM" | head -5 | grep -c "mv:B:" || true)
            if [ "$MV_CHECK" -gt 0 ]; then
                echo "  ✓ Move table (mv) tags confirmed present" | tee -a "$LOG"
            else
                echo "  ✗ WARNING: mv tags not found in $name" | tee -a "$LOG"
            fi
        fi

        # ─── Aligned + coordinate-sorted BAM ───
        if [ -f "$ALIGNED_BAM" ]; then
            echo "  [SKIP aligned] $name already aligned" | tee -a "$LOG"
            continue
        fi

        echo "  Aligning to reference: $name" | tee -a "$LOG"

        "$DORADO_BIN" aligner \
            "$REFERENCE" \
            "$UNALIGNED_BAM" \
            2>> "$LOG" \
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

    FOLDER_NAME="dorado_sup_methyl_q${Q}"
    OUTPUT_DIR="$BASE_OUTPUT_DIR/$FOLDER_NAME"
    UNALIGNED_DIR="$OUTPUT_DIR/unaligned"
    ALIGNED_DIR="$OUTPUT_DIR/aligned"
    mkdir -p "$UNALIGNED_DIR" "$ALIGNED_DIR"

    echo "" | tee -a "$LOG"
    echo "→ Methylation calling | Q${Q}" | tee -a "$LOG"

    for pod5 in "$POD5_DIR"/*.pod5; do
        [ -e "$pod5" ] || continue
        name=$(basename "$pod5" .pod5)

        UNALIGNED_BAM="$UNALIGNED_DIR/${name}_sup_methyl_q${Q}.unaligned.bam"
        ALIGNED_BAM="$ALIGNED_DIR/${name}_sup_methyl_q${Q}.aligned.sorted.bam"

        # ─── Unaligned BAM (with move tables + MM/ML) ───
        if [ -f "$UNALIGNED_BAM" ]; then
            echo "  [SKIP unaligned] $name already processed" | tee -a "$LOG"
        else
            echo "  Basecalling + mod calling: $name" | tee -a "$LOG"

            "$DORADO_BIN" basecaller \
                "$BASE_MODEL" \
                "$pod5" \
                --modified-bases-models "$METHYL_MODEL" \
                --min-qscore "$Q" \
                --emit-moves \
                > "$UNALIGNED_BAM" \
                2>> "$LOG"

            READ_COUNT=$(samtools view -c "$UNALIGNED_BAM")
            echo "  ✓ $READ_COUNT reads → $(basename $UNALIGNED_BAM)" | tee -a "$LOG"

            # Verify MM/ML tags are present
            MM_CHECK=$(samtools view "$UNALIGNED_BAM" | head -5 | grep -c "MM:Z:" || true)
            if [ "$MM_CHECK" -gt 0 ]; then
                echo "  ✓ MM/ML methylation tags confirmed present" | tee -a "$LOG"
            else
                echo "  ✗ WARNING: MM/ML tags not found in $name" | tee -a "$LOG"
            fi

            # Verify mv (move table) tags are present
            MV_CHECK=$(samtools view "$UNALIGNED_BAM" | head -5 | grep -c "mv:B:" || true)
            if [ "$MV_CHECK" -gt 0 ]; then
                echo "  ✓ Move table (mv) tags confirmed present" | tee -a "$LOG"
            else
                echo "  ✗ WARNING: mv tags not found in $name" | tee -a "$LOG"
            fi
        fi

        # ─── Aligned + coordinate-sorted BAM ───
        if [ -f "$ALIGNED_BAM" ]; then
            echo "  [SKIP aligned] $name already aligned" | tee -a "$LOG"
            continue
        fi

        echo "  Aligning to reference: $name" | tee -a "$LOG"

        "$DORADO_BIN" aligner \
            "$REFERENCE" \
            "$UNALIGNED_BAM" \
            2>> "$LOG" \
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
for bam in "$BASE_OUTPUT_DIR"/dorado_sup_basecall_q*/unaligned/*.bam; do
    [ -f "$bam" ] || continue
    count=$(samtools view -c "$bam")
    printf "  %-10s reads → %s\n" "$count" "$(basename $bam)" | tee -a "$LOG"
done

echo "" | tee -a "$LOG"
echo "--- Basecall outputs (aligned) ---" | tee -a "$LOG"
for bam in "$BASE_OUTPUT_DIR"/dorado_sup_basecall_q*/aligned/*.bam; do
    [ -f "$bam" ] || continue
    mapped=$(samtools view -c -F 4 "$bam")
    total=$(samtools view -c "$bam")
    printf "  %-10s mapped / %-10s total → %s\n" "$mapped" "$total" "$(basename $bam)" | tee -a "$LOG"
done

echo "" | tee -a "$LOG"
echo "--- Methylation outputs (unaligned) ---" | tee -a "$LOG"
for bam in "$BASE_OUTPUT_DIR"/dorado_sup_methyl_q*/unaligned/*.bam; do
    [ -f "$bam" ] || continue
    count=$(samtools view -c "$bam")
    printf "  %-10s reads → %s\n" "$count" "$(basename $bam)" | tee -a "$LOG"
done

echo "" | tee -a "$LOG"
echo "--- Methylation outputs (aligned) ---" | tee -a "$LOG"
for bam in "$BASE_OUTPUT_DIR"/dorado_sup_methyl_q*/aligned/*.bam; do
    [ -f "$bam" ] || continue
    mapped=$(samtools view -c -F 4 "$bam")
    total=$(samtools view -c "$bam")
    printf "  %-10s mapped / %-10s total → %s\n" "$mapped" "$total" "$(basename $bam)" | tee -a "$LOG"
done

echo "" | tee -a "$LOG"
echo "All done: $(date)" | tee -a "$LOG"
