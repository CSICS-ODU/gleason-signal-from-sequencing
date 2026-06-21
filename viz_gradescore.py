# -*- coding: utf-8 -*-
"""viz_gradescore.py -- proposed Fig 8 (Option A): a one-axis 'grade-association score'.

For every sample we compute, in the cohort-normalized (per-gene-z) space, how much more
it resembles the high-grade reference group than the low-grade group:

    score(x) = dist(x, low-mean) - dist(x, high-mean)      (higher => more high-like)

Reference samples are scored leave-one-out (self excluded from its own class mean) so the
two reference distributions are honest; the 10 AA+EA query patients are scored against the
full reference class means. This is the quantity the comparison actually rests on, so it
surfaces the signal that PCA buries in a central blob. No classifier / no p-value here.
"""
import types, functools
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import distance_knn as dk
print = functools.partial(print, flush=True)

METRIC = "spearman"           # the headline recovery metric
OUR = ["PTEN", "TP53", "PIK3CA", "SMAD4", "GNAS", "APC", "ATM", "CDH1"]
HCOL, LCOL = "#c0392b", "#2980b9"

# ---- load (per make_final_fig) --------------------------------------------
ref_mat = dk.load_reference(types.SimpleNamespace(ref_rcc_dir="GSE244631", ref_matrix=None))
ref_labels = dk.load_ref_labels("GSE244631_ref_labels.csv")
ref_samples = [s for s in ref_mat.columns if s in ref_labels]
ref_mat = ref_mat[ref_samples]
ref_lab = np.array([ref_labels[s] for s in ref_samples])
qprof = dk.load_query_profiles("cohort_matched")
truth = dk.load_query_truth("genefilelists/metadata_orig.csv", "ID", "Gleason")
qpids = [p for p in qprof if p in truth and not str(p).isdigit()]
qlab = np.array([truth[p] for p in qpids])
listB = set()
for p in qprof:
    if str(p).isdigit():
        listB |= {g for g in qprof[p].index if qprof[p].get(g, 0) > 0}
genes = sorted((set(OUR) | listB) & set(ref_mat.index))
sub = ref_mat.reindex(genes).fillna(0.0)

# ---- cohort-normalized (per-gene-z) space ---------------------------------
ref_g = dk.zscore_per_gene(sub).values.T                              # 51 x 23
Q = pd.DataFrame({p: qprof[p].reindex(genes).fillna(0.0) for p in qpids})
q_g = dk.zscore_per_gene(Q).values.T                                  # 10 x 23
hi = ref_lab == "high"
lo = ref_lab == "low"
H_all, L_all = ref_g[hi].mean(0), ref_g[lo].mean(0)


def score(x, hm, lm):
    return dk.compute_dist(x, lm, METRIC) - dk.compute_dist(x, hm, METRIC)


# reference scores, leave-one-out
ref_score = np.empty(len(ref_lab))
idx = np.arange(len(ref_lab))
for i in range(len(ref_lab)):
    if ref_lab[i] == "high":
        hm = ref_g[hi & (idx != i)].mean(0); lm = L_all
    else:
        hm = H_all; lm = ref_g[lo & (idx != i)].mean(0)
    ref_score[i] = score(ref_g[i], hm, lm)
q_score = np.array([score(q_g[i], H_all, L_all) for i in range(len(qpids))])

for nm, m in (("ref high", ref_score[hi]), ("ref low", ref_score[lo]),
              ("query high", q_score[qlab == "high"]), ("query low", q_score[qlab == "low"])):
    print(f"{nm:<11} median {np.median(m):+.3f}  mean {m.mean():+.3f}  n={len(m)}")

# ---- figure ---------------------------------------------------------------
plt.rcParams.update({"font.size": 12})
rng = np.random.default_rng(0)
fig, ax = plt.subplots(figsize=(7.4, 4.2))
rows = [("Reference\nhigh grade", ref_score[hi], HCOL, 3, False),
        ("Reference\nlow grade", ref_score[lo], LCOL, 2, False),
        ("Query AA+EA\nhigh grade", q_score[qlab == "high"], HCOL, 1, True),
        ("Query AA+EA\nlow grade", q_score[qlab == "low"], LCOL, 0, True)]

ax.axvline(0, color="grey", ls="--", lw=1, zorder=1)
ax.text(-0.03, 3.78, "more low-grade-like  $\\leftarrow$", ha="right", va="center",
        fontsize=10, color="dimgrey")
ax.text(0.03, 3.78, "$\\rightarrow$  more high-grade-like", ha="left", va="center",
        fontsize=10, color="dimgrey")

for label, vals, col, y, is_query in rows:
    if not is_query and len(vals) > 4:
        vp = ax.violinplot(vals, positions=[y], vert=False, widths=0.75,
                           showmeans=False, showextrema=False)
        for b in vp["bodies"]:
            b.set_facecolor(col); b.set_alpha(0.16); b.set_edgecolor("none")
            v = b.get_paths()[0].vertices          # clip the density tail to the data
            v[:, 0] = np.clip(v[:, 0], vals.min(), vals.max())
    jit = rng.uniform(-0.12, 0.12, len(vals))
    ax.scatter(vals, y + jit, c=col, s=(170 if is_query else 42),
               edgecolor="black", linewidth=(1.5 if is_query else 0.4),
               alpha=(0.95 if is_query else 0.6), zorder=4)
    ax.scatter(np.median(vals), y, marker="D", c="white", edgecolor=col,
               linewidth=2.0, s=110, zorder=5)

ax.axhline(1.5, color="black", lw=0.7, alpha=0.3)        # divider: reference | query
ax.set_yticks([r[3] for r in rows])
ax.set_yticklabels([r[0] for r in rows])
ax.set_ylim(-0.6, 4.0)
ax.set_xlabel("Grade-association score\n(resemblance to high-grade $-$ to low-grade reference group)")
fig.tight_layout()
fig.savefig("fig_gradescore.png", dpi=160)
print("\nwrote fig_gradescore.png")
