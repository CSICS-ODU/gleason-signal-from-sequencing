# -*- coding: utf-8 -*-
"""Diagnostic: why is the permutation p so high for the AA+EA-only run?
Reuses panel_recovery_aaea's setup, prints per-patient predictions, the confusion,
and the FULL discrete permutation null (all C(10,4)=210 distinct label assignments)."""
import types, functools, itertools
import numpy as np
import pandas as pd
import distance_knn as dk
print = functools.partial(print, flush=True)

ROUNDS, SEED = 200, 0
OUR = ["PTEN", "TP53", "PIK3CA", "SMAD4", "GNAS", "APC", "ATM", "CDH1"]

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
is_indian = lambda p: str(p).isdigit()

listB = set()
for p in pids:
    if is_indian(p):
        listB |= {g for g in qprof[p].index if qprof[p].get(g, 0) > 0}
genes = sorted((set(OUR) | listB) & panel)
qpids = [p for p in pids if not is_indian(p)]


def predict_all(metric):
    rng = np.random.default_rng(SEED)
    ref_arr = dk.within_sample_normalize_matrix(ref_mat.reindex(genes).fillna(0.0), "z").values.T
    preds = {}
    for p in qpids:
        qv = dk.normalize_query_vector(qprof[p].reindex(genes).fillna(0.0).values, "within_patient_z")
        dists = np.array([dk.compute_dist(qv, ref_arr[j], metric) for j in range(len(ref_arr))])
        preds[p], _ = dk.predict_balanced(dists, ref_lab, "classmean", None, "subsample", ROUNDS, rng)
    return preds


def balacc(tarr, parr):
    recs = [(parr[tarr == c] == c).mean() for c in ("high", "low") if (tarr == c).any()]
    return float(np.mean(recs))


truth_arr = np.array([truth[p] for p in qpids])
for metric in ("cosine", "spearman"):
    preds = predict_all(metric)
    parr = np.array([preds[p] for p in qpids])
    obs = balacc(truth_arr, parr)
    print(f"\n=== {metric} : observed balanced acc = {obs:.3f} ===")
    print(f"{'patient':<9}{'truth':<7}{'pred':<7}{'ok'}")
    for p in qpids:
        print(f"{p:<9}{truth[p]:<7}{preds[p]:<7}{'Y' if preds[p]==truth[p] else '.'}")
    # EXACT discrete null: every way to call 4 of the 10 patients 'high'
    n = len(qpids)
    null = []
    for hi_idx in itertools.combinations(range(n), 4):
        lab = np.array(["low"] * n, dtype=object)  # object dtype: no string truncation
        lab[list(hi_idx)] = "high"
        null.append(balacc(lab, parr))
    null = np.array(null)
    vals, cnts = np.unique(np.round(null, 4), return_counts=True)
    print(f"\nexact null over all C(10,4)={len(null)} label assignments:")
    print(f"{'balacc':>8}{'count':>8}{'frac':>8}{'cum>=':>9}")
    cum = 0
    for v, c in sorted(zip(vals, cnts), reverse=True):
        cum += c
        print(f"{v:>8.3f}{c:>8}{c/len(null):>8.3f}{cum/len(null):>9.3f}")
    p_exact = (null >= obs - 1e-9).mean()
    print(f"exact permutation p (>= {obs:.3f}) = {p_exact:.3f}")
