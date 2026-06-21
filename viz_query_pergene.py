# -*- coding: utf-8 -*-
"""viz_query_pergene.py -- project the AA+EA query patients into the PER-GENE-Z
reference space (the RIGHT subplot of fig_ref_pca_panel23.png) and see where they land.

A real classifier CANNOT do this: per-gene z-scoring the query needs the other query
patients to set each gene's mean/SD, which is test-time leakage. But per-gene z is also
the only thing that makes two different modalities comparable -- it turns both expression
counts and variant counts into dimensionless z-scores (each gene centred within its own
modality). So as a DIAGNOSTIC we grant ourselves the leak and ask the real question:

   if the query were allowed into per-gene-z space, would high/low queries fall near the
   matching reference grade region, and would a nearest-class-mean call get them right?

If YES  -> the cross-modality signal exists; the only thing killing the classifier was the
           leakage-forced self-normalization (a bigger/held-out query cohort could fix it).
If NO   -> the modality gap is fundamental and per-gene z does not rescue it.

We z-score the query per gene ACROSS THE 10 QUERY PATIENTS (the leak), z-score the
reference per gene across its 51 samples, then project the query through the reference's
per-gene-z PCA and also run a nearest-class-mean call in that space.

Output: fig_query_pergene_panel23.png + printed "cheating" balanced accuracy.
"""
import types
import functools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import distance_knn as dk
print = functools.partial(print, flush=True)

OUR = ["PTEN", "TP53", "PIK3CA", "SMAD4", "GNAS", "APC", "ATM", "CDH1"]
HCOL, LCOL = "#c0392b", "#2980b9"

# ---- load (same as before) ------------------------------------------------
ref_mat = dk.load_reference(types.SimpleNamespace(ref_rcc_dir="GSE244631", ref_matrix=None))
ref_labels = dk.load_ref_labels("GSE244631_ref_labels.csv")
ref_samples = [s for s in ref_mat.columns if s in ref_labels]
ref_mat = ref_mat[ref_samples]
ref_lab = np.array([ref_labels[s] for s in ref_samples])

qprof = dk.load_query_profiles("cohort_matched")
truth = dk.load_query_truth("genefilelists/metadata_orig.csv", "ID", "Gleason")
qpids = [p for p in qprof if p in truth and not str(p).isdigit()]   # AA+EA only
qlab = np.array([truth[p] for p in qpids])

listB = set()
for p in qprof:
    if str(p).isdigit():
        listB |= {g for g in qprof[p].index if qprof[p].get(g, 0) > 0}
panel23 = sorted((set(OUR) | listB) & set(ref_mat.index))
print(f"panel: {len(panel23)} genes | reference 51 ({(ref_lab=='high').sum()}H/"
      f"{(ref_lab=='low').sum()}L) | query {len(qpids)} ({(qlab=='high').sum()}H/"
      f"{(qlab=='low').sum()}L, AA+EA)\n")

# ---- per-gene z on each side (the leak is on the query side) ---------------
refz = dk.zscore_per_gene(ref_mat.reindex(panel23).fillna(0.0))     # genes x 51
Xref = refz.values.T                                               # 51 x 23

Q = pd.DataFrame({p: qprof[p].reindex(panel23).fillna(0.0) for p in qpids})  # genes x 10
Qz = dk.zscore_per_gene(Q)                                          # genes x 10 (LEAK)
Xq = Qz.values.T                                                   # 10 x 23

# ---- "cheating" nearest-class-mean in per-gene-z space --------------------
def balacc(lab, preds):
    recs = [(preds[lab == c] == c).mean() for c in ("high", "low") if (lab == c).any()]
    return float(np.mean(recs))


print(f"{'metric':<10}{'leaky bal_acc':>15}{'perm_p':>10}   "
      f"(vs legit self-z 0.62/0.67, panel-restricted)")
mh, ml = Xref[ref_lab == "high"].mean(0), Xref[ref_lab == "low"].mean(0)
for metric in ("cosine", "spearman"):
    preds = np.array(["high" if dk.compute_dist(Xq[i], mh, metric)
                      < dk.compute_dist(Xq[i], ml, metric) else "low"
                      for i in range(len(qpids))])
    ba = balacc(qlab, preds)
    rng = np.random.default_rng(12345)
    ge = sum(balacc(rng.permutation(qlab), preds) >= ba - 1e-9 for _ in range(1000))
    pval = (ge + 1) / 1001
    hits = " ".join(f"{p}:{t}->{pr}{'*' if t==pr else ''}"
                    for p, t, pr in zip(qpids, qlab, preds))
    print(f"{metric:<10}{ba:>15.3f}{pval:>10.4f}")
    print(f"            {hits}")
print()

# ---- PCA on reference per-gene-z, project query into it -------------------
mu = Xref.mean(0)
Xc = Xref - mu
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
PCs = Vt[:2]
evr = (S ** 2) / (S ** 2).sum()
ref_scores = Xc @ PCs.T
q_scores = (Xq - mu) @ PCs.T

fig, ax = plt.subplots(figsize=(8.2, 6.6))
for cls, col in (("high", HCOL), ("low", LCOL)):
    m = ref_lab == cls
    ax.scatter(ref_scores[m, 0], ref_scores[m, 1], c=col, s=55, alpha=0.55,
               edgecolor="white", linewidth=0.5, label=f"ref {cls} (n={m.sum()})", zorder=2)
    ax.scatter(*ref_scores[m].mean(0), c=col, marker="*", s=620, edgecolor="black",
               linewidth=1.3, zorder=6)
for cls, mk in (("high", "^"), ("low", "v")):
    m = qlab == cls
    col = HCOL if cls == "high" else LCOL
    ax.scatter(q_scores[m, 0], q_scores[m, 1], c=col, marker=mk, s=150,
               edgecolor="black", linewidth=1.6, zorder=7,
               label=f"query {cls} (n={m.sum()})")
    ax.scatter(*q_scores[m].mean(0), c=col, marker="X", s=340, edgecolor="black",
               linewidth=1.6, zorder=8)
for i, p in enumerate(qpids):
    ax.annotate(p.replace("PT", ""), q_scores[i], fontsize=7,
                xytext=(3, 3), textcoords="offset points", zorder=9)
ax.set_xlabel(f"PC1 ({evr[0]*100:.0f}%)")
ax.set_ylabel(f"PC2 ({evr[1]*100:.0f}%)")
ax.set_title("Query projected into PER-GENE-Z reference space (23-gene panel)\n"
             "stars = reference grade prototypes; X = query grade centroids\n"
             "(leaky diagnostic: a real classifier cannot per-gene-z the query)")
ax.legend(fontsize=8, framealpha=0.9, loc="best")
ax.axhline(0, color="grey", lw=0.5, alpha=0.4)
ax.axvline(0, color="grey", lw=0.5, alpha=0.4)
fig.tight_layout()
fig.savefig("fig_query_pergene_panel23.png", dpi=140)
print("wrote fig_query_pergene_panel23.png")
