# -*- coding: utf-8 -*-
"""viz_reference.py -- is the 'high prototype vs low prototype' idea actually present
in the Brazilian reference (GSE244631) space?

The cross-modality classifier collapses onto one class. That can happen for two very
different reasons, and this script separates them:

  (A) the reference high/low classes are NOT separable even among themselves
      (no grade structure in the expression space) -> the premise fails outright; or
  (B) they ARE separable within the reference, but our DNA-variant query vectors all
      land in one region of that space (modality gap) -> the premise is fine, the
      cross-modality projection is what breaks.

We test both:
  * within-reference separability: leave-one-out nearest-class-mean balanced accuracy
    (expression-vs-expression, the friendly within-modality case) + class-prototype
    cosine separation, on the full panel and on the 23-gene candidate panel, under
    both normalizations (the classifier's self-signature z and the proper per-gene z).
  * geometry: PCA of the 51 reference samples coloured by grade, with the two class
    prototypes marked, and -- on the 23-gene panel -- our 10 AA+EA query patients
    overlaid in the SAME projection, so we can literally see whether they all fall on
    the 'low' side.

Outputs: fig_ref_pca_full.png, fig_ref_pca_panel23.png, and a printed diagnostic table.
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
HCOL, LCOL = "#c0392b", "#2980b9"   # high = red, low = blue

# ---- load reference + labels + query (same as the classifier) -------------
ref_mat = dk.load_reference(types.SimpleNamespace(ref_rcc_dir="GSE244631", ref_matrix=None))
ref_labels = dk.load_ref_labels("GSE244631_ref_labels.csv")
ref_samples = [s for s in ref_mat.columns if s in ref_labels]
ref_mat = ref_mat[ref_samples]
ref_lab = np.array([ref_labels[s] for s in ref_samples])

qprof = dk.load_query_profiles("cohort_matched")
truth = dk.load_query_truth("genefilelists/metadata_orig.csv", "ID", "Gleason")
qpids = [p for p in qprof if p in truth and not str(p).isdigit()]   # AA+EA only
qlab = np.array([truth[p] for p in qpids])

# 23-gene candidate panel = our candidates U Indian-detected, intersect reference
listB = set()
for p in qprof:
    if str(p).isdigit():
        listB |= {g for g in qprof[p].index if qprof[p].get(g, 0) > 0}
panel23 = sorted((set(OUR) | listB) & set(ref_mat.index))
full_panel = sorted(ref_mat.index)
print(f"reference: {ref_mat.shape[0]} genes x {len(ref_samples)} samples "
      f"({(ref_lab=='high').sum()} high / {(ref_lab=='low').sum()} low)")
print(f"23-gene candidate panel present in reference: {len(panel23)}")
print(f"AA+EA query patients: {len(qpids)} "
      f"({(qlab=='high').sum()} high / {(qlab=='low').sum()} low)\n")


# ---- separability diagnostics ---------------------------------------------
def loo_ncm(X, labels, metric):
    """Leave-one-out nearest-class-mean balanced accuracy. X: samples x genes."""
    n = len(labels)
    preds = []
    for i in range(n):
        keep = np.arange(n) != i
        Xt, yt = X[keep], labels[keep]
        mh, ml = Xt[yt == "high"].mean(0), Xt[yt == "low"].mean(0)
        dh, dl = dk.compute_dist(X[i], mh, metric), dk.compute_dist(X[i], ml, metric)
        preds.append("high" if dh < dl else "low")
    preds = np.array(preds)
    recs = [(preds[labels == c] == c).mean() for c in ("high", "low")]
    return float(np.mean(recs))


def proto_sep(X, labels, metric):
    """Distance between the high and low class prototypes (bigger = more separated)."""
    mh, ml = X[labels == "high"].mean(0), X[labels == "low"].mean(0)
    return dk.compute_dist(mh, ml, metric)


def norm_samples(sub, kind):
    """sub: genes x samples -> samples x genes, normalized. kind: 'selfz'|'genez'."""
    if kind == "selfz":
        return dk.within_sample_normalize_matrix(sub, "z").values.T
    return dk.zscore_per_gene(sub).values.T


print(f"{'panel':<10}{'norm':<8}{'metric':<10}{'LOO bal_acc':>12}{'proto_sep':>11}")
for pname, gl in (("full", full_panel), ("panel23", panel23)):
    sub = ref_mat.reindex(gl).fillna(0.0)
    for kind in ("selfz", "genez"):
        X = norm_samples(sub, kind)
        for metric in ("cosine", "spearman"):
            acc = loo_ncm(X, ref_lab, metric)
            sep = proto_sep(X, ref_lab, metric)
            print(f"{pname:<10}{kind:<8}{metric:<10}{acc:>12.3f}{sep:>11.3f}")
print("\n(LOO bal_acc ~0.5 => grade is NOT separable even within the reference;"
      "\n >0.5 => structure exists in expression space and the failure is the modality gap)\n")


# ---- PCA helpers -----------------------------------------------------------
def pca_fit(X):
    mu = X.mean(0)
    Xc = X - mu
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    evr = (S ** 2) / (S ** 2).sum()
    return mu, Vt[:2], Xc @ Vt[:2].T, evr[:2]


def pca_apply(mu, PCs, Y):
    return (Y - mu) @ PCs.T


def scatter_ref(ax, scores, evr, title, query_scores=None):
    for cls, col, mk in (("high", HCOL, "o"), ("low", LCOL, "o")):
        m = ref_lab == cls
        ax.scatter(scores[m, 0], scores[m, 1], c=col, marker=mk, s=55,
                   edgecolor="white", linewidth=0.6, alpha=0.9,
                   label=f"ref {cls} (n={m.sum()})", zorder=3)
        cen = scores[m].mean(0)
        ax.scatter(*cen, c=col, marker="*", s=520, edgecolor="black",
                   linewidth=1.2, zorder=5)
    if query_scores is not None:
        for cls, mk in (("high", "^"), ("low", "v")):
            m = qlab == cls
            ax.scatter(query_scores[m, 0], query_scores[m, 1], c="black",
                       marker=mk, s=90, linewidth=1.4, zorder=6,
                       label=f"query {cls} (n={m.sum()})")
    ax.set_xlabel(f"PC1 ({evr[0]*100:.0f}%)")
    ax.set_ylabel(f"PC2 ({evr[1]*100:.0f}%)")
    ax.set_title(title)
    ax.legend(fontsize=8, framealpha=0.9)
    ax.axhline(0, color="grey", lw=0.5, alpha=0.4)
    ax.axvline(0, color="grey", lw=0.5, alpha=0.4)


# ---- Figure 1: full panel, reference only (both normalizations) -----------
fig, axes = plt.subplots(1, 2, figsize=(13, 5.4))
for ax, kind, ttl in zip(axes, ("selfz", "genez"),
                         ("self-signature z (classifier space)", "per-gene z (proper)")):
    sub = ref_mat.reindex(full_panel).fillna(0.0)
    X = norm_samples(sub, kind)
    mu, PCs, scores, evr = pca_fit(X)
    scatter_ref(ax, scores, evr, f"Reference, full panel\n{ttl}")
fig.suptitle("Brazilian reference (GSE244631): is grade separable in expression space?",
             fontsize=12)
fig.tight_layout()
fig.savefig("fig_ref_pca_full.png", dpi=140)
print("wrote fig_ref_pca_full.png")

# ---- Figure 2: 23-gene panel, reference + query overlay -------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5.4))
for ax, kind, ttl in zip(axes, ("selfz", "genez"),
                         ("self-signature z (classifier space)", "per-gene z (proper)")):
    sub = ref_mat.reindex(panel23).fillna(0.0)
    X = norm_samples(sub, kind)
    mu, PCs, scores, evr = pca_fit(X)
    # project the AA+EA query patients into the SAME space (classifier-space only)
    qscores = None
    if kind == "selfz":
        Y = np.array([dk.normalize_query_vector(
            qprof[p].reindex(panel23).fillna(0.0).values, "within_patient_z")
            for p in qpids])
        qscores = pca_apply(mu, PCs, Y)
    scatter_ref(ax, scores, evr, f"23-gene candidate panel\n{ttl}", qscores)
fig.suptitle("23-gene panel: reference grade structure + where our AA+EA queries land",
             fontsize=12)
fig.tight_layout()
fig.savefig("fig_ref_pca_panel23.png", dpi=140)
print("wrote fig_ref_pca_panel23.png")
