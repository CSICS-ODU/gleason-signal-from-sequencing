# -*- coding: utf-8 -*-
"""quick_norm_sweep.py -- compare normalization schemes for the panel comparison.

Three normalization schemes for the panel-restricted cross-modality classifier
(23-gene panel, subsample-balanced nearest-class-mean), each on the full 20-patient
cohort and on the AA+EA-only stratum:

  A  ref self-signature z , query self-signature z   <- the manuscript's comparison
                                                         (panel_recovery*.py)
  B  ref per-gene z        , query self-signature z   <- mixed alternative, not reported
  C  ref per-gene z        , query per-gene z         <- cohort-standardized (transductive)

A is the manuscript's fixed-panel comparison (Table 2b rows 1-2); C is its
cohort-standardized row (Table 2b row 4, Fig. 8). p-values here use 1,000 random
shuffles; table2_exact.py gives the exact p-values the manuscript reports.
"""
import types, functools
import numpy as np, pandas as pd
import distance_knn as dk
print = functools.partial(print, flush=True)

ROUNDS, SEED, NPERM = 200, 0, 1000
OUR = ["PTEN", "TP53", "PIK3CA", "SMAD4", "GNAS", "APC", "ATM", "CDH1"]

ref_mat = dk.load_reference(types.SimpleNamespace(ref_rcc_dir="GSE244631", ref_matrix=None))
ref_labels = dk.load_ref_labels("GSE244631_ref_labels.csv")
ref_samples = [s for s in ref_mat.columns if s in ref_labels]
ref_mat = ref_mat[ref_samples]
ref_lab = np.array([ref_labels[s] for s in ref_samples])
qprof = dk.load_query_profiles("cohort_matched")
truth = dk.load_query_truth("genefilelists/metadata_orig.csv", "ID", "Gleason")
truth = {p: v for p, v in truth.items() if p in qprof}
allp = [p for p in qprof if p in truth]
listB = set()
for p in allp:
    if str(p).isdigit():
        listB |= {g for g in qprof[p].index if qprof[p].get(g, 0) > 0}
genes = sorted((set(OUR) | listB) & set(ref_mat.index))
sub = ref_mat.reindex(genes).fillna(0.0)
print(f"panel {len(genes)} genes\n")


def balacc(lab, preds):
    return float(np.mean([(preds[lab == c] == c).mean()
                          for c in ("high", "low") if (lab == c).any()]))


def ref_arr(scheme):
    if scheme == "self":
        return dk.within_sample_normalize_matrix(sub, "z").values.T
    return dk.zscore_per_gene(sub).values.T          # per-gene


def query_vecs(qpids, scheme):
    if scheme == "self":
        return np.array([dk.normalize_query_vector(
            qprof[p].reindex(genes).fillna(0.0).values, "within_patient_z") for p in qpids])
    Q = pd.DataFrame({p: qprof[p].reindex(genes).fillna(0.0) for p in qpids})
    return dk.zscore_per_gene(Q).values.T            # per-gene (leak)


def run(qpids, refscheme, qscheme, metric):
    R = ref_arr(refscheme)
    Qv = query_vecs(qpids, qscheme)
    lab = np.array([truth[p] for p in qpids])
    rng = np.random.default_rng(SEED)
    preds = []
    for qv in Qv:
        dists = np.array([dk.compute_dist(qv, R[j], metric) for j in range(len(R))])
        p, _ = dk.predict_balanced(dists, ref_lab, "classmean", None, "subsample", ROUNDS, rng)
        preds.append(p)
    preds = np.array(preds)
    obs = balacc(lab, preds)
    rng2 = np.random.default_rng(12345)
    ge = sum(balacc(rng2.permutation(lab), preds) >= obs - 1e-9 for _ in range(NPERM))
    return obs, (ge + 1) / (NPERM + 1)


aaea = [p for p in allp if not str(p).isdigit()]
schemes = [("A self/self", "self", "self"),
           ("B pergene/self", "gene", "self"),
           ("C pergene/pergene", "gene", "gene")]
for cohort_name, qp in (("FULL 20 (13H/7L)", allp), ("AA+EA 10 (4H/6L)", aaea)):
    print(f"=== {cohort_name} ===")
    print(f"{'scheme':<24}{'cos bal(p)':>16}{'spear bal(p)':>16}")
    for name, rs, qs in schemes:
        co, cp = run(qp, rs, qs, "cosine")
        so, sp = run(qp, rs, qs, "spearman")
        print(f"{name:<24}{f'{co:.2f} ({cp:.3f})':>16}{f'{so:.2f} ({sp:.3f})':>16}")
    print()
