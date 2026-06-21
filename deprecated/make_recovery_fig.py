# -*- coding: utf-8 -*-
"""Bar chart: severity-classification balanced accuracy, all genes vs candidate
biomarkers, for the best-performing scale-tolerant metric. Reads
knn_paper_summary.csv (written by distance_knn.py --paper-tables)."""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

df = pd.read_csv("knn_paper_summary.csv")

# scale-tolerant metrics we report; pick the one with the strongest recovery
KEEP = ["cosine", "spearman"]
df = df[df["metric"].isin(KEEP)]
best_metric = (df[df.gene_set == "biomarker"]
               .groupby("metric")["balanced_acc"].max().idxmax())

def best(gene_set):
    sub = df[(df.gene_set == gene_set) & (df.metric == best_metric)]
    return float(sub["balanced_acc"].max())

vals = [best("all"), best("biomarker")]
labels = ["All genes\n(512)", "Candidate\nbiomarkers (22)"]

fig, ax = plt.subplots(figsize=(3.6, 3.2))
bars = ax.bar(labels, vals, width=0.6, color=["#9aa0a6", "#3b6fb6"])
ax.axhline(0.50, ls="--", c="k", lw=1)
ax.axhline(0.65, ls=":", c="0.45", lw=1)
ax.text(1.46, 0.505, "chance (0.50)", fontsize=7.5, ha="right", va="bottom")
ax.text(1.46, 0.655, "majority baseline (0.65)", fontsize=7.5, ha="right", va="bottom")
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width() / 2, v + 0.02, f"{v:.2f}",
            ha="center", fontsize=10)
ax.set_ylim(0, 1.0)
ax.set_ylabel("Balanced accuracy")
ax.set_title(f"Severity recovery ({best_metric})")
plt.tight_layout()
plt.savefig("Manuscript/fig_recovery.png", dpi=200)
print(f"best_metric={best_metric}  all={vals[0]:.3f}  biomarker={vals[1]:.3f}")
print("wrote Manuscript/fig_recovery.png")
