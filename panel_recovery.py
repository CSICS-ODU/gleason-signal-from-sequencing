# -*- coding: utf-8 -*-
"""panel_recovery.py -- rebuild the candidate set as the honest UNION, then re-run the
cross-modality recovery on that FIXED set.

  set = { our named AA+EA candidates }  UNION  { every gene detected in Indian-panel
          patients (clinically supported; no statistical test required) }
        intersected with the GSE244631 reference panel (only ref-measurable genes usable).

The set is defined by EXTERNAL warrants (our analysis + the clinical panel), not by the
test labels, so there is no feature-selection leakage to control: we evaluate the fixed
set under leave-one-patient-out. Predictions come only from the reference and the query
profile (never the query's own label), so the permutation null is a straight label
shuffle on fixed predictions.
"""
import types
import functools
import os
import shutil
import numpy as np
import pandas as pd
import distance_knn as dk
print = functools.partial(print, flush=True)

ROUNDS = 200
SEED = 0
NPERM = 1000
OUR_CANDIDATES = ["PTEN", "TP53", "PIK3CA", "SMAD4", "GNAS", "APC", "ATM", "CDH1"]

# ---- load exactly as nested_loo.py ----------------------------------------
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


# ---- build the union set --------------------------------------------------
listB = set()
for p in pids:
    if is_indian(p):
        listB |= {g for g in qprof[p].index if qprof[p].get(g, 0) > 0}
listA = set(OUR_CANDIDATES)
unionset = (listA | listB) & panel
genes = sorted(unionset)

nin = sum(is_indian(p) for p in pids)
print(f"cohort: {len(pids)} patients ({nin} Indian-panel)")
print(f"List A (our named candidates): {len(listA)}  {sorted(listA)}")
print(f"  on reference panel: {sorted(listA & panel)}")
print(f"  dropped, not on ref panel: {sorted(listA - panel)}")
print(f"List B (genes detected in Indian patients): {len(listB)} "
      f"(on reference panel: {len(listB & panel)})")
print(f"A-only on panel (ours, not in Indian set): {sorted((listA & panel) - listB)}")
print(f"UNION on reference panel: {len(genes)} genes\n")

if os.path.exists("biomarker_subset.csv"):
    shutil.copy("biomarker_subset.csv", "biomarker_subset_OLD.csv")
pd.DataFrame({"gene": genes}).to_csv("biomarker_subset.csv", index=False)
print(f"wrote biomarker_subset.csv ({len(genes)} genes); old -> biomarker_subset_OLD.csv\n")


# ---- proper LOO: rebuild List B from the OTHER Indian patients per fold ----
# The feature set never sees the held-out patient's own variants. (List A is a fixed
# external list, so it does not change across folds.) Predictions still depend only on
# the reference + the held-out patient's own profile, never on its Gleason label, so the
# permutation null is a straight label shuffle on the fixed LOO predictions.
def listB_excluding(held):
    s = set()
    for p in pids:
        if p != held and is_indian(p):
            s |= {g for g in qprof[p].index if qprof[p].get(g, 0) > 0}
    return s


fold_set = {p: sorted((listA | listB_excluding(p)) & panel) for p in pids}
sizes = [len(fold_set[p]) for p in pids]
print(f"per-fold LOO set size: min {min(sizes)} max {max(sizes)} (full union {len(genes)})")


def predict_all_loo(metric):
    rng = np.random.default_rng(SEED)
    preds = {}
    for p in pids:
        gl = fold_set[p]
        sub = ref_mat.reindex(gl).fillna(0.0)
        ref_arr = dk.within_sample_normalize_matrix(sub, "z").values.T
        qv = dk.normalize_query_vector(
            qprof[p].reindex(gl).fillna(0.0).values, "within_patient_z")
        dists = np.array([dk.compute_dist(qv, ref_arr[j], metric)
                          for j in range(len(ref_arr))])
        preds[p], _ = dk.predict_balanced(dists, ref_lab, "classmean", None,
                                          "subsample", ROUNDS, rng)
    return preds


def bal_acc(tmap, preds):
    t = np.array([tmap[p] for p in pids])
    pr = np.array([preds[p] for p in pids])
    recs = [(pr[t == c] == c).mean() for c in ("high", "low") if (t == c).any()]
    return float(np.mean(recs))


plist = list(pids)
labels = np.array([truth[p] for p in plist])
print(f"\n{'metric':<10}{'bal_acc':>9}{'perm_p':>10}")
for metric in ("cosine", "spearman"):
    preds = predict_all_loo(metric)
    obs = bal_acc(truth, preds)
    rng = np.random.default_rng(12345)
    ge = 1
    for i in range(NPERM):
        perm = rng.permutation(labels)
        if bal_acc({p: perm[j] for j, p in enumerate(plist)}, preds) >= obs - 1e-9:
            ge += 1
    pval = ge / (NPERM + 1)
    print(f"{metric:<10}{obs:>9.2f}{pval:>10.4f}")
