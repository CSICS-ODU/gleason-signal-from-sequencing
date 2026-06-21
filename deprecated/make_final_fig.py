# -*- coding: utf-8 -*-
"""make_final_fig.py -- the new closing figure + numbers for the manuscript.

ONE classifier (the panel-restricted, subsample-balanced nearest-class-mean used for
Table 'distance_knn_bm'), ONE query (the 10 broadly sequenced AA+EA patients, 4H/6L),
ONE gene panel (the fixed 23-gene candidate panel). The ONLY thing that changes between
the two panels is the normalization:

  LEFT  self-signature z  -- leakage-free; each patient normalized against itself only.
                            This is what a real classifier must use. Signal hidden.
  RIGHT per-gene z         -- leaky diagnostic; each gene centred across the query cohort,
                            the axis on which grade actually lives. Signal recovered.

Isolating the normalization shows the genome-wide failure is methodological (the
leakage-forced self-normalization at n=10), not the DNA-vs-mRNA modality gap.

Writes manuscript/fig_recovery_pergene.png and prints the balanced accuracies + perm p.
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

ROUNDS, SEED, NPERM = 200, 0, 1000
OUR = ["PTEN", "TP53", "PIK3CA", "SMAD4", "GNAS", "APC", "ATM", "CDH1"]
HCOL, LCOL = "#c0392b", "#2980b9"

# ---- load -----------------------------------------------------------------
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
genes = sorted((set(OUR) | listB) & set(ref_mat.index))
sub = ref_mat.reindex(genes).fillna(0.0)                     # genes x 51


def balacc(lab, preds):
    recs = [(preds[lab == c] == c).mean() for c in ("high", "low") if (lab == c).any()]
    return float(np.mean(recs))


def classify(ref_arr, qvecs):
    """Same classifier as Table distance_knn_bm: subsample-balanced nearest-class-mean
    over the individual reference samples."""
    rng = np.random.default_rng(SEED)
    preds = []
    for qv in qvecs:
        dists = np.array([dk.compute_dist(qv, ref_arr[j], m) for j in range(len(ref_arr))])
        p, _ = dk.predict_balanced(dists, ref_lab, "classmean", None, "subsample", ROUNDS, rng)
        preds.append(p)
    return np.array(preds)


def perm_p(preds, obs):
    rng = np.random.default_rng(12345)
    ge = sum(balacc(rng.permutation(qlab), preds) >= obs - 1e-9 for _ in range(NPERM))
    return (ge + 1) / (NPERM + 1)


def pca_fit(X):
    mu = X.mean(0)
    U, S, Vt = np.linalg.svd(X - mu, full_matrices=False)
    return mu, Vt[:2], (X - mu) @ Vt[:2].T, (S ** 2) / (S ** 2).sum()


# ---- build both spaces -----------------------------------------------------
# self-signature z (leakage-free): ref within-sample z, query within-patient z
ref_self = dk.within_sample_normalize_matrix(sub, "z").values.T          # 51 x 23
q_self = np.array([dk.normalize_query_vector(
    qprof[p].reindex(genes).fillna(0.0).values, "within_patient_z") for p in qpids])

# per-gene z (leaky diagnostic): ref per-gene z; query per-gene z across the 10 queries
ref_gene = dk.zscore_per_gene(sub).values.T                              # 51 x 23
Q = pd.DataFrame({p: qprof[p].reindex(genes).fillna(0.0) for p in qpids})
q_gene = dk.zscore_per_gene(Q).values.T                                  # 10 x 23

print(f"panel {len(genes)} genes | AA+EA query {len(qpids)} ({(qlab=='high').sum()}H/"
      f"{(qlab=='low').sum()}L) | reference 51 ({(ref_lab=='high').sum()}H/{(ref_lab=='low').sum()}L)\n")
print(f"{'space':<14}{'metric':<10}{'bal_acc':>9}{'perm_p':>9}")
ann = {}
for space, R, Qv in (("self-signature z", ref_self, q_self),
                     ("per-gene z (leak)", ref_gene, q_gene)):
    for metric in ("cosine", "spearman"):
        # rebuild classify with the chosen metric
        rng = np.random.default_rng(SEED)
        preds = []
        for qv in Qv:
            dists = np.array([dk.compute_dist(qv, R[j], metric) for j in range(len(R))])
            p, _ = dk.predict_balanced(dists, ref_lab, "classmean", None, "subsample", ROUNDS, rng)
            preds.append(p)
        preds = np.array(preds)
        obs = balacc(qlab, preds)
        pv = perm_p(preds, obs)
        ann[(space, metric)] = (obs, pv)
        print(f"{space:<14}{metric:<10}{obs:>9.3f}{pv:>9.4f}")

# ---- figure ----------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5.8))
for ax, (title, R, Qv, kind) in zip(axes, (
        ("Leakage-free space (self-signature z)\nwhat a real classifier must use",
         ref_self, q_self, "self"),
        ("Per-gene-z diagnostic\nthe axis on which grade lives",
         ref_gene, q_gene, "gene"))):
    mu, PCs, rs, evr = pca_fit(R)
    qs = (Qv - mu) @ PCs.T
    for cls, col in (("high", HCOL), ("low", LCOL)):
        m = ref_lab == cls
        ax.scatter(rs[m, 0], rs[m, 1], c=col, s=46, alpha=0.5, edgecolor="white",
                   linewidth=0.5, label=f"ref {cls} (n={m.sum()})", zorder=2)
        ax.scatter(*rs[m].mean(0), c=col, marker="*", s=560, edgecolor="black",
                   linewidth=1.2, zorder=6)
    for cls, mk in (("high", "^"), ("low", "v")):
        m = qlab == cls
        col = HCOL if cls == "high" else LCOL
        ax.scatter(qs[m, 0], qs[m, 1], c=col, marker=mk, s=140, edgecolor="black",
                   linewidth=1.5, zorder=7, label=f"query {cls} (n={m.sum()})")
        ax.scatter(*qs[m].mean(0), c=col, marker="X", s=300, edgecolor="black",
                   linewidth=1.5, zorder=8)
    co, cp = ann[("self-signature z" if kind == "self" else "per-gene z (leak)", "cosine")]
    so, sp = ann[("self-signature z" if kind == "self" else "per-gene z (leak)", "spearman")]
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(f"PC1 ({evr[0]*100:.0f}%)")
    ax.set_ylabel(f"PC2 ({evr[1]*100:.0f}%)")
    ax.text(0.02, 0.98,
            f"balanced acc.\ncosine {co:.2f} (p={cp:.2f})\nspearman {so:.2f} (p={sp:.2f})",
            transform=ax.transAxes, va="top", ha="left", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", ec="grey", alpha=0.9))
    ax.legend(fontsize=8, framealpha=0.9, loc="lower right")
    ax.axhline(0, color="grey", lw=0.5, alpha=0.4)
    ax.axvline(0, color="grey", lw=0.5, alpha=0.4)
fig.tight_layout()
fig.savefig("manuscript/fig_recovery_pergene.png", dpi=150)
print("\nwrote manuscript/fig_recovery_pergene.png")
