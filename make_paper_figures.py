# -*- coding: utf-8 -*-
"""make_paper_figures.py -- manuscript Fig. 2 and Fig. 3a, drawn from the result CSVs.

Run build_paper_results.py first; this script reads its outputs in paper_results/:
  fig_burden.png  Fig. 2   variants found in >=2 patients, by annotation category and
                           cohort (log scale)            <- tableIVB_bucket_counts.csv
  fig_lean.png    Fig. 3a  per gene and population: % of high-grade patients with a
                           variant minus % of low-grade patients with one
                                                          <- tableIVE_cross_pop.csv
Black-and-white safe: AA plain, EA diagonal stripes, IN crosshatch.

Usage:  python make_paper_figures.py [outdir]      (default: current folder)
"""
import csv
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

RESULTS = "paper_results"
out = sys.argv[1] if len(sys.argv) > 1 else "."
CM = 1 / 2.54
AA, EA, IN = "#2A78D6", "#EB6834", "#1BAF7A"
STYLE = {"AA": dict(color=AA, hatch=""), "EA": dict(color=EA, hatch="////"), "IN": dict(color=IN, hatch="xxxx")}
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["TeX Gyre Termes", "Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 7,
    "hatch.linewidth": 0.7, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
})


def legend_handles(keys, labels):
    return [Patch(facecolor=STYLE[k]["color"], hatch=STYLE[k]["hatch"], edgecolor="black", linewidth=0.5, label=l)
            for k, l in zip(keys, labels)]


def read(name):
    with open(f"{RESULTS}/{name}", newline="") as fh:
        return list(csv.DictReader(fh))


# ---- Fig. 2: variants found in >=2 patients, by cohort and annotation category (log scale)
BUCKETS = [("pathogenic/known", "Catalogued,\nnot benign"),
           ("benign/likely-benign", "Catalogued\nbenign/likely benign"),
           ("unknown/unannotated", "Not catalogued")]
counts = {(r["bucket"], r["group"]): int(r["n_snv"]) for r in read("tableIVB_bucket_counts.csv")}
data = {k: [counts[(b, k)] for b, _ in BUCKETS] for k in ("AA", "EA", "IN")}
fig, ax = plt.subplots(figsize=(8.3 * CM, 4.4 * CM))
w = 0.22
for i, k in enumerate(["AA", "EA", "IN"]):
    xs = [b + (i - 1) * w for b in range(3)]
    ax.bar(xs, [v - 1 for v in data[k]], w * 0.92, bottom=1, color=STYLE[k]["color"], hatch=STYLE[k]["hatch"],
           edgecolor="black", linewidth=0.5)
    for x, v in zip(xs, data[k]):
        ax.text(x, v * 1.6, f"{v:,}", rotation=90, ha="center", va="bottom", fontsize=5.5)
ax.set_yscale("log"); ax.set_ylim(1, 2e6); ax.minorticks_off()
ax.set_xticks(range(3), [label for _, label in BUCKETS])
ax.set_xlim(-0.55, 2.55)
ax.set_ylabel("Variants found in ≥2 patients\n(count per cohort, log scale)")
ax.yaxis.grid(True, color="0.88", linewidth=0.5); ax.set_axisbelow(True)
ax.legend(handles=legend_handles(["AA", "EA", "IN"], ["AA", "EA", "IN"]), ncol=3, frameon=False,
          loc="lower center", bbox_to_anchor=(0.5, 1.0), handlelength=1.6, handleheight=1.0, columnspacing=1.2)
fig.tight_layout(pad=0.2)
fig.savefig(f"{out}/fig_burden.png", dpi=600)
plt.close(fig)
print(f"wrote {out}/fig_burden.png")

# ---- Fig. 3a: within-population difference in % of patients with a variant (high - low grade)
GENES = ["PTEN", "TP53", "PIK3CA", "SMAD4", "GNAS", "APC"]
N_HIGH, N_LOW = 2, 3                         # whole-genome-sequenced patients per population
cross = {r["gene"]: r for r in read("tableIVE_cross_pop.csv")}


def lean(hi_lo):
    nh, nl = (int(x) for x in hi_lo.split("/"))
    return round(100 * nh / N_HIGH - 100 * nl / N_LOW)


gap = {k: [lean(cross[g][f"{k}_hi/lo"]) for g in GENES] for k in ("AA", "EA")}
fig, ax = plt.subplots(figsize=(8.3 * CM, 4.6 * CM))
h = 0.36
for j, k in enumerate(["AA", "EA"]):           # AA above EA within each gene
    ys = [len(GENES) - 1 - g + (0.5 - j) * h for g in range(len(GENES))]
    ax.barh(ys, gap[k], h * 0.9, color=STYLE[k]["color"], hatch=STYLE[k]["hatch"], edgecolor="black", linewidth=0.5)
    for y, v in zip(ys, gap[k]):
        ax.text(v - 3 if v < 0 else v + 3, y, f"{v:d}".replace("-", "−"),
                ha="right" if v < 0 else "left", va="center", fontsize=6)
ax.set_yticks(range(len(GENES)), GENES[::-1], fontstyle="italic")
ax.set_xlim(-150, 115); ax.set_ylim(-0.6, len(GENES) - 0.4)
ax.axvline(0, color="0.3", linewidth=0.6)
ax.xaxis.grid(True, color="0.88", linewidth=0.5); ax.set_axisbelow(True)
ax.set_xlabel("% of high-grade patients with a variant\nminus % of low-grade patients with a variant")
ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}".replace("-", "−")))
ax.legend(handles=legend_handles(["AA", "EA"], ["African American", "European American"]), frameon=False,
          loc="upper left", fontsize=6, handlelength=1.8, handleheight=1.0)
fig.tight_layout(pad=0.2)
fig.savefig(f"{out}/fig_lean.png", dpi=600)
plt.close(fig)
print(f"wrote {out}/fig_lean.png")
