"""
build_paper_figures.py - generate the manuscript result figures from the paper_results/
CSVs (the SNV-derived population-aware analysis; NOT the rejected GEO-classification plots).

  fig_divergence.png  (IV.E)  per-gene severity lean (high-rate - low-rate) in AA vs EA
                              as diverging horizontal bars -> shows PTEN consistent-low
                              and TP53/PIK3CA/SMAD4/GNAS population-divergent at a glance.
  fig_burden.png      (IV.B)  recurrent SNV counts per cohort (IN/EA/AA) across the three
                              annotation buckets, log scale -> the targeted-panel (IN)
                              vs broad-sequencing (AA/EA) orders-of-magnitude difference.

Run (matplotlib lives in the WSL py3 env):
  wsl -d Ubuntu -- bash -lc "cd /mnt/c/Users/s1banerj/Develop/vcf && \
      ~/miniconda3/bin/conda run -n py3 python build_paper_figures.py"
Outputs PNGs into manuscript/.
"""

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = "."   # write figures into the current folder (manuscript is maintained separately)
RESULTS = "paper_results"

# Severity cohort sizes per population at PT510_low (Table II): AA = 2 high / 3 low,
# EA = 2 high / 3 low. Rates are taken over these denominators.
H = {"AA": 2, "EA": 2}
L = {"AA": 3, "EA": 3}


def parse_pair(s):
    a, b = s.split("/")
    return int(a), int(b)


def fig_divergence():
    rows = {}
    with open(f"{RESULTS}/tableIVE_cross_pop.csv", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows[r["gene"]] = r
    # consistent signal, then the four divergent genes, then a flat reference
    genes = [g for g in ["PTEN", "TP53", "PIK3CA", "SMAD4", "GNAS", "APC"] if g in rows]
    aa_gap, ea_gap = [], []
    for g in genes:
        nha, nla = parse_pair(rows[g]["AA_hi/lo"])
        nhe, nle = parse_pair(rows[g]["EA_hi/lo"])
        aa_gap.append(nha / H["AA"] - nla / L["AA"])
        ea_gap.append(nhe / H["EA"] - nle / L["EA"])

    y = np.arange(len(genes))[::-1]          # first listed gene at the top
    bh = 0.38
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.barh(y + bh / 2, aa_gap, height=bh, color="#4C72B0", label="African American")
    ax.barh(y - bh / 2, ea_gap, height=bh, color="#C44E52", label="European American")
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(genes, fontstyle="italic")
    ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel("Severity lean  (carrier-rate gap: high $-$ low)")
    ax.text(-1.0, len(genes) - 0.4, "$\\leftarrow$ low-severity", ha="left",
            va="center", fontsize=9, color="gray")
    ax.text(1.0, len(genes) - 0.4, "high-severity $\\rightarrow$", ha="right",
            va="center", fontsize=9, color="gray")
    ax.set_title("Population-specific severity lean of candidate genes")
    ax.legend(loc="lower right", frameon=True, fontsize=9)
    ax.grid(axis="x", ls=":", alpha=0.5)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_divergence.png", dpi=200)
    plt.close(fig)
    print("wrote fig_divergence.png")


def fig_burden():
    counts = {}
    with open(f"{RESULTS}/tableIVB_bucket_counts.csv", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            counts.setdefault(r["bucket"], {})[r["group"]] = int(r["n_snv"])
    buckets = ["pathogenic/known", "benign/likely-benign", "unknown/unannotated"]
    blabels = ["Pathogenic/\nknown", "Benign/\nlikely-benign", "Unknown/\nunannotated"]
    groups = ["IN", "EA", "AA"]
    colors = {"IN": "#55A868", "EA": "#C44E52", "AA": "#4C72B0"}

    x = np.arange(len(buckets))
    w = 0.26
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for i, grp in enumerate(groups):
        vals = [counts[b].get(grp, 0) for b in buckets]
        bars = ax.bar(x + (i - 1) * w, vals, width=w, label=grp, color=colors[grp])
        ax.bar_label(bars, labels=[f"{v:,}" for v in vals], padding=2, fontsize=7,
                     rotation=90)
    ax.set_yscale("log")
    ax.set_ylim(1, 2e5)
    ax.set_xticks(x)
    ax.set_xticklabels(blabels)
    ax.set_ylabel("Recurrent SNVs (log scale)")
    ax.set_title("Variant burden by cohort and annotation bucket")
    ax.legend(title="Cohort", fontsize=9)
    ax.grid(axis="y", ls=":", alpha=0.5, which="both")
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_burden.png", dpi=200)
    plt.close(fig)
    print("wrote fig_burden.png")


def fig_candidates():
    """Integrated candidate-gene landscape: each gene's computational signal (left dot)
    and experimental evidence (right dot), unifying Tables VIII (experiment-anchored)
    and IX (computation-anchored) into one view."""
    from matplotlib.lines import Line2D

    genes = ["APC", "ATM", "CDH1", "TP53", "PTEN", "PIK3CA", "SMAD4", "GNAS"]
    comp = {
        "APC": "#7f7f7f", "ATM": "#55A868", "CDH1": "#55A868", "TP53": "#DD8452",
        "PTEN": "#4C72B0", "PIK3CA": "#DD8452", "SMAD4": "#DD8452", "GNAS": "#DD8452",
    }
    expr = {
        "APC": "#2f2f2f", "ATM": "#2f2f2f", "CDH1": "#9e9e9e", "TP53": "#9e9e9e",
        "PTEN": "white", "PIK3CA": "white", "SMAD4": "white", "GNAS": "white",
    }
    y = np.arange(len(genes))[::-1]
    fig, ax = plt.subplots(figsize=(7, 4.3))
    for gene, yi in zip(genes, y):
        ax.scatter(0, yi, s=340, color=comp[gene], edgecolor="black", zorder=3)
        ax.scatter(1, yi, s=340, color=expr[gene], edgecolor="black", zorder=3)
    ax.axhline(3.5, color="gray", ls="--", lw=0.8)   # both-arm (top) vs computation-only
    ax.set_yticks(y)
    ax.set_yticklabels(genes, fontstyle="italic")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Computational\nsignal", "Experimental\nevidence"])
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(-0.6, len(genes) - 0.4)

    comp_leg = [
        Line2D([], [], marker="o", ls="", markersize=11, markeredgecolor="black",
               markerfacecolor=c, label=l)
        for l, c in [("Population-divergent", "#DD8452"),
                     ("Consistent low-severity", "#4C72B0"),
                     ("Recurrent in IN panel", "#55A868"),
                     ("Recurrent (flat)", "#7f7f7f")]]
    expr_leg = [
        Line2D([], [], marker="o", ls="", markersize=11, markeredgecolor="black",
               markerfacecolor=c, label=l)
        for l, c in [("Recurrent detection", "#2f2f2f"),
                     ("Context-dependent", "#9e9e9e"),
                     ("Not tested", "white")]]
    leg1 = ax.legend(handles=comp_leg, title="Computational signal", loc="upper left",
                     bbox_to_anchor=(1.04, 1.0), fontsize=8, title_fontsize=9)
    ax.add_artist(leg1)
    ax.legend(handles=expr_leg, title="Experimental evidence", loc="lower left",
              bbox_to_anchor=(1.04, 0.0), fontsize=8, title_fontsize=9)
    ax.set_title("Integrated candidate-gene landscape")
    ax.grid(axis="y", ls=":", alpha=0.3)
    fig.savefig(f"{OUT}/fig_candidates.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote fig_candidates.png")


def fig_recovery():
    """Figure 7 = a direct visual of Table I vs Table X. Every bar is a number that
    appears in one of those two tables: the genome-wide cross-modality comparison
    (Table I: cosine 0.50, spearman 0.60) versus the comparison restricted to the
    clinically-prioritized candidate panel
    (Table X: cosine 0.81, spearman 0.78), per similarity metric, with permutation p
    above each bar. Both metrics recover under the panel restriction. No permutation
    interval is drawn: the p-values carry the significance, and a null band that the
    observed value escapes reads misleadingly as an error bar."""
    metrics = ["Cosine", "Spearman"]
    genome = [0.50, 0.60]          # Table I  (tab:distance_knn)
    biomk = [0.81, 0.78]           # Table X  (tab:distance_knn_bm)
    p_genome = ["p = 1.00", "p = 0.27"]
    p_biomk = ["p = 0.014", "p = 0.016"]

    x = np.arange(len(metrics))
    w = 0.36
    fig, ax = plt.subplots(figsize=(5.0, 3.7))
    b1 = ax.bar(x - w / 2, genome, width=w, color="#8C8C8C", edgecolor="black",
                linewidth=0.5, label="Genome-wide (Table I)")
    b2 = ax.bar(x + w / 2, biomk, width=w, color="#4C72B0", edgecolor="black",
                linewidth=0.5, label="Panel-restricted (Table X)")
    for bars, vals, ps in ((b1, genome, p_genome), (b2, biomk, p_biomk)):
        for rect, v, p in zip(bars, vals, ps):
            ax.text(rect.get_x() + rect.get_width() / 2, v + 0.015,
                    f"{v:.2f}\n{p}", ha="center", va="bottom", fontsize=8)
    ax.axhline(0.50, color="black", ls="--", lw=0.8, label="chance (0.50)")
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.set_ylabel("Balanced accuracy")
    ax.set_ylim(0, 1.16)
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_recovery.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("wrote fig_recovery.png")


def main():
    os.makedirs(OUT, exist_ok=True)
    fig_divergence()
    fig_burden()
    fig_candidates()
    fig_recovery()
    print("done")


if __name__ == "__main__":
    main()
