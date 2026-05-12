import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ============================================================
# LOAD DATA
# ============================================================

df = pd.read_csv("/Volumes/Amishi_SSD/bio_data/16Apr/final_results_epic/final_results_epic.csv")

# ============================================================
# PARSE MODEL + QSCORE
# ============================================================

def parse_group(group):
    parts = group.split("/")
    model_part = parts[1]  # dorado_fast_methyl_q10

    tokens = model_part.split("_")

    basecaller = tokens[0]   # dorado / guppy
    model = tokens[1]        # fast / hac / sup
    qscore = int(tokens[-1].replace("q", ""))

    return f"{basecaller}_{model}", qscore

df[["model", "qscore"]] = df["group"].apply(
    lambda x: pd.Series(parse_group(x))
)

# ============================================================
# PREP DATA
# ============================================================

# convert MCC to percentage
df["MCC_percent"] = df["MCC"] * 100

models = sorted(df["model"].unique())
qscores = sorted(df["qscore"].unique())

# ============================================================
# PLOT
# ============================================================

x = np.arange(len(models))
width = 0.15

plt.figure(figsize=(10,6))

for i, q in enumerate(qscores):
    sub = df[df["qscore"] == q]

    # align bars by model
    y = []
    for m in models:
        val = sub[sub["model"] == m]["MCC_percent"]
        y.append(val.values[0] if len(val) > 0 else 0)

    plt.bar(
        x + i*width,
        y,
        width=width,
        label=f"Q{q}"
    )

    # add value labels
    for j, v in enumerate(y):
        plt.text(
            x[j] + i*width,
            v + 1,
            f"{v:.1f}",
            ha='center',
            fontsize=8
        )

# ============================================================
# FINAL TOUCHES
# ============================================================

plt.xticks(x + width*len(qscores)/2, models, rotation=20)
plt.ylabel("MCC (%)")
plt.title("MCC across Models and Q-scores for EpiC with at least 5 reads at each site")
plt.legend(title="Q-score")
plt.grid(axis='y', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.show()
