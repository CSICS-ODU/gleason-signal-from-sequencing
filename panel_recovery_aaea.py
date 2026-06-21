# -*- coding: utf-8 -*-
"""panel_recovery_aaea.py -- the FINAL panel-restricted classification, re-run on the
broadly sequenced AA+EA stratum ONLY (the Indian targeted-panel cohort is excluded
from the query).

What is the SAME as panel_recovery.py:
  * the gene panel. It is still the FIXED 23-gene union
        { our named AA+EA candidates }  UNION  { every gene detected in the Indian-panel
          patients }   intersected with the GSE244631 reference panel.
    The Indian patients still BUILD the panel (listB); they are only dropped from the
    set of patients we classify.
  * the reference (GSE244631, 19 high / 32 low), the cross-modality cosine/Spearman
    nearest-class-mean classifier, the class-balanced subsampling (majority class
    subsampled to the 19-high minority -> 19 high / 19 low per round), and the
    self-normalization on both sides.

What CHANGES:
  * the query is the 10 broadly sequenced AA+EA patients (4 high / 6 low) instead of
    all 20. The Indian targeted-panel patients are not classified.
  * because none of the AA+EA query patients contribute to the panel through listB,
    the panel is genuinely fixed across every query patient -- there is no per-fold
    rebuild to do, so we apply the single 23-gene panel unchanged to all 10.

The permutation null is a straight label shuffle on the fixed predictions over the
10 AA+EA patients.
"""
import types
import functools
import numpy as np
import pandas as pd
import distance_knn as dk
print = functools.partial(print, flush=True)

ROUNDS = 200
SEED = 0
NPERM = 1000
OUR_CANDIDATES = ["PTEN", "TP53", "PIK3CA", "SMAD4", "GNAS", "APC", "ATM", "CDH1"]

# ---- load exactly as panel_recovery.py ------------------------------------
ref_mat = dk.load_reference(types.SimpleNamespace(ref_rcc_dir="GSE244631", ref_matrix=None))
ref_labels = dk.load_ref_labels("GSE244631_ref_labels.csv")
qprof = dk.load_query_profiles("cohort_matched")
truth = dk.load_query_truth("genefilelists/metadata_orig.csv", "ID", "Gleason")
truth = {p: v for p, v in truth.items() if p in qprof}
pids = [p for p in qprof if p in truth]
ref_samples = [s for s in ref_mat.columns if s in ref_labels]
ref_mat = ref_mat[ref_samples]
ref_lab = np.array([ref_labels[s] for s in ref_samples])
panel = set(ref_mat.index)


def is_indian(p):
    return str(p).isdigit()


# ---- build the FIXED union panel from ALL patients (unchanged) ------------
# The Indian patients still define listB here; only the query set drops them below.
listB = set()
for p in pids:
    if is_indian(p):
        listB |= {g for g in qprof[p].index if qprof[p].get(g, 0) > 0}
listA = set(OUR_CANDIDATES)
genes = sorted((listA | listB) & panel)

# ---- query = AA+EA only ----------------------------------------------------
query_pids = [p for p in pids if not is_indian(p)]
nhi = sum(truth[p] == "high" for p in query_pids)
print(f"reference: {len(ref_lab)} samples "
      f"({(ref_lab=='high').sum()} high / {(ref_lab=='low').sum()} low) "
      f"-> balanced to {min((ref_lab=='high').sum(), (ref_lab=='low').sum())}"
      f"+{min((ref_lab=='high').sum(), (ref_lab=='low').sum())} per round")
print(f"FIXED panel on reference: {len(genes)} genes")
print(f"query (AA+EA only): {len(query_pids)} patients "
      f"({nhi} high / {len(query_pids)-nhi} low)")
print(f"  dropped Indian-panel patients: {sum(is_indian(p) for p in pids)}\n")


def predict_all(metric):
    """Classify each AA+EA query patient against the reference on the FIXED panel."""
    rng = np.random.default_rng(SEED)
    sub = ref_mat.reindex(genes).fillna(0.0)
    ref_arr = dk.within_sample_normalize_matrix(sub, "z").values.T
    preds = {}
    for p in query_pids:
        qv = dk.normalize_query_vector(
            qprof[p].reindex(genes).fillna(0.0).values, "within_patient_z")
        dists = np.array([dk.compute_dist(qv, ref_arr[j], metric)
                          for j in range(len(ref_arr))])
        preds[p], _ = dk.predict_balanced(dists, ref_lab, "classmean", None,
                                          "subsample", ROUNDS, rng)
    return preds


def bal_acc(tmap, preds):
    t = np.array([tmap[p] for p in query_pids])
    pr = np.array([preds[p] for p in query_pids])
    recs = [(pr[t == c] == c).mean() for c in ("high", "low") if (t == c).any()]
    return float(np.mean(recs))


labels = np.array([truth[p] for p in query_pids])
print(f"{'metric':<10}{'bal_acc':>9}{'perm_p':>10}")
for metric in ("cosine", "spearman"):
    preds = predict_all(metric)
    obs = bal_acc(truth, preds)
    rng = np.random.default_rng(12345)
    ge = 1
    for i in range(NPERM):
        perm = rng.permutation(labels)
        if bal_acc({p: perm[j] for j, p in enumerate(query_pids)}, preds) >= obs - 1e-9:
            ge += 1
    pval = ge / (NPERM + 1)
    print(f"{metric:<10}{obs:>9.2f}{pval:>10.4f}")
