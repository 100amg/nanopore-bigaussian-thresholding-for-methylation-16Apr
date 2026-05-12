# Nanopore Bigaussian Methylation Thresholding Analysis Pipeline

Pipeline for Oxford Nanopore methylation analysis, threshold optimisation, and model/q-score benchmarking using AUC, GMM (bi-Gaussian), and MCC analyses. 

## Repository Structure

```bash
repo/
├── dorado_fast_comparison.sh
├── dorado_hac_comparison.sh
├── dorado_sup_comparison.sh
├── guppy_hac_comparison.sh
├── extract_modkit_from_bam.py
├── auc_model_qscore_analysis.py
├── bigaussian_threshold_analysis.py
├── gaussian_new_EpiC_5_reads_filter.py
├── gaussian_new_25Apr.py
├── plot_mcc_vs_threshold.py
└── README.md
```

## Scripts

| Script                                | Purpose                                     |
| ------------------------------------- | ------------------------------------------- |
| `dorado_fast_comparison.sh`           | Dorado FAST methylation basecalling         |
| `dorado_hac_comparison.sh`            | Dorado HAC methylation basecalling          |
| `dorado_sup_comparison.sh`            | Dorado SUP methylation basecalling          |
| `guppy_hac_comparison.sh`             | Guppy HAC methylation basecalling           |
| `extract_modkit_from_bam.py`          | Extract per-read methylation probabilities  |
| `auc_model_qscore_analysis.py`        | BED aggregation, AUC analysis, thresholding |
| `bigaussian_threshold_analysis.py`    | Bi-Gaussian / GMM threshold analysis        |
| `gaussian_new_EpiC_5_reads_filter.py` | MCC analysis with ≥5 read filter            |
| `gaussian_new_25Apr.py`               | MCC threshold evaluation                    |
| `plot_mcc_vs_threshold.py`            | MCC comparison plots                        |

## Workflow

### 1. Basecalling

```bash
bash dorado_fast_comparison.sh
bash dorado_hac_comparison.sh
bash dorado_sup_comparison.sh
bash guppy_hac_comparison.sh
```

### 2. Extract Modkit Data

```bash
python extract_modkit_from_bam.py
```

### 3. Run AUC + BED Analysis

```bash
python auc_model_qscore_analysis.py
```

### 4. Run Bi-Gaussian Threshold Analysis

```bash
python bigaussian_threshold_analysis.py
```

### 5. Run MCC Threshold Analysis

```bash
python gaussian_new_25Apr.py
python gaussian_new_EpiC_5_reads_filter.py
```

### 6. Generate MCC Plots

```bash
python plot_mcc_vs_threshold.py
```

## Main Outputs

* BAM files
* Modkit TSV outputs
* BED methylation summaries
* AUC rankings
* Gaussian thresholds
* MCC evaluation tables
* Diagnostic plots

## Full Documentation

Detailed pipeline documentation:

[Google Docs Documentation](https://docs.google.com/document/d/1cf8tUt95NHWAcoEajb62hgJrUrc3dvQ7wSXJtKJ2hVc/edit?tab=t.y62x5ccfrnsn)
