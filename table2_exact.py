# -*- coding: utf-8 -*-
"""table2_exact.py -- Table 2b of the manuscript (23-gene candidate panel), every row,
with the EXACT permutation p-values the table reports.

Rows (balanced accuracy, cosine and Spearman):
  1  Full cohort (13H/7L), fixed panel      -- as panel_recovery.py (per-fold List B)
  2  AA+EA (4H/6L), fixed panel             -- as panel_recovery_aaea.py
  3  AA+EA, nested selection                -- as nested_loo.py AA,EA (11 genes per tail,
                                               re-selected inside every leave-one-out fold)
  4  AA+EA, cohort-standardized             -- quick_norm_sweep.py scheme C (each gene
                                               z-scored across samples, reference and query
                                               separately; transductive)

Every row uses the same scorer: nearest-class-mean distance, class-balanced by
subsampling the 32 low-grade reference samples to 19 over 200 rounds, majority vote.

Exact nulls:
  * full cohort: all C(20,7) = 77,520 allocations of 7 low-grade labels;
  * AA+EA: the 100 population-preserving allocations (2 high / 3 low within each of
    AA and EA, C(5,2)^2), so no result can arise from population alone.
p = fraction of allocations whose balanced accuracy is >= the observed value.
For rows 1, 2 and 4 the predictions never use a query label, so they are computed once
and only the labels are permuted. For row 3 gene selection uses the training labels, so
the whole nested pipeline is re-run for each of the 100 allocations.
"""
import types
import itertools
import functools
import numpy as np
import pandas as pd
import distance_knn as dk
print = functools.partial(print, flush=True)

ROUNDS = 200
SEED = 0
N_EACH = 11
OUR_CANDIDATES = ["PTEN", "TP53", "PIK3CA", "SMAD4", "GNAS", "APC", "ATM", "CDH1"]
METRICS = ("cosine", "spearman")

# ---- load (as panel_recovery.py / nested_loo.py) ----------------------------
ref_mat = dk.load_reference(types.SimpleNamespace(ref_rcc_dir="GSE244631", ref_matrix=None))
ref_labels = dk.load_ref_labels("GSE244631_ref_labels.csv")
ref_samples = [s for s in ref_mat.columns if s in ref_labels]
ref_mat = ref_mat[ref_samples]
ref_lab = np.array([ref_labels[s] for s in ref_samples])
panel = set(ref_mat.index)
qprof = dk.load_query_profiles("cohort_matched")
truth = dk.load_query_truth("genefilelists/metadata_orig.csv", "ID", "Gleason")
races = dk.load_query_races("genefilelists/metadata_orig.csv", "ID", "race")
truth = {p: v for p, v in truth.items() if p in qprof}
pids = [p for p in qprof if p in truth]
aaea = [p for p in pids if races.get(p) in ("AA", "EA")]


def is_indian(p):
    return races.get(p) == "IN"


def listB(exclude=None):
    """Genes with >=1 variant in the Indian clinical-panel patients."""
    s = set()
    for p in pids:
        if is_indian(p) and p != exclude:
            s |= {g for g in qprof[p].index if qprof[p].get(g, 0) > 0}
    return s


fixed_genes = sorted((set(OUR_CANDIDATES) | listB()) & panel)


def predict(qvec, ref_arr, metric, rng):
    dists = np.array([dk.compute_dist(qvec, r, metric) for r in ref_arr])
    return dk.predict_balanced(dists, ref_lab, "classmean", None, "subsample", ROUNDS, rng)[0]


def self_z_preds(query, genes_for, metric):
    """Each profile z-scored across its own genes (per-tumor), reference likewise."""
    rng = np.random.default_rng(SEED)
    preds = {}
    for p in query:
        g = genes_for(p)
        ref_arr = dk.within_sample_normalize_matrix(ref_mat.reindex(g).fillna(0.0), "z").values.T
        qv = dk.normalize_query_vector(qprof[p].reindex(g).fillna(0.0).values, "within_patient_z")
        preds[p] = predict(qv, ref_arr, metric, rng)
    return preds


def cohort_z_preds(query, metric):
    """Each gene z-scored across samples: reference across its 51, query across its 10."""
    rng = np.random.default_rng(SEED)
    sub = ref_mat.reindex(fixed_genes).fillna(0.0)
    ref_arr = dk.zscore_per_gene(sub).values.T
    Q = pd.DataFrame({p: qprof[p].reindex(fixed_genes).fillna(0.0) for p in query})
    qarr = dk.zscore_per_gene(Q).values.T
    return {p: predict(qarr[i], ref_arr, metric, rng) for i, p in enumerate(query)}


def bal_acc(tmap, preds, query):
    recs = [np.mean([preds[p] == c for p in query if tmap[p] == c]) for c in ("high", "low")]
    return float(np.mean(recs))


def allocations_full():
    n_low = sum(truth[p] == "low" for p in pids)
    for lows in itertools.combinations(pids, n_low):
        lows = set(lows)
        yield {p: ("low" if p in lows else "high") for p in pids}


def allocations_aaea():
    aa = [p for p in aaea if races[p] == "AA"]
    ea = [p for p in aaea if races[p] == "EA"]
    na = sum(truth[p] == "high" for p in aa)
    ne = sum(truth[p] == "high" for p in ea)
    for ha in itertools.combinations(aa, na):
        for he in itertools.combinations(ea, ne):
            highs = set(ha) | set(he)
            yield {p: ("high" if p in highs else "low") for p in aaea}


def exact_p(obs, null):
    null = np.asarray(null)
    return float(np.mean(null >= obs - 1e-9)), len(null)


# ---- nested selection (nested_loo.py, AA+EA) --------------------------------
universe = sorted(panel & set().union(*[set(qprof[p].index) for p in aaea]))
carrier = {p: {g: (qprof[p].get(g, 0) > 0) for g in universe} for p in aaea}


def select(train, tmap):
    """11 most low-enriched + 11 most high-enriched genes by carrier-rate gap.
    Genes are sorted by gap with ties in gene-symbol order, and the two ends taken."""
    hi = [p for p in train if tmap[p] == "high"]
    lo = [p for p in train if tmap[p] == "low"]
    gap = {g: (sum(carrier[p][g] for p in hi) / len(hi)
               - sum(carrier[p][g] for p in lo) / len(lo)) for g in universe}
    order = sorted(universe, key=lambda g: gap[g])
    return order[:N_EACH] + order[-N_EACH:]


_dcache = {}


def nested_balacc(tmap, metric):
    rng = np.random.default_rng(SEED)
    preds = {}
    for p in aaea:
        g = select([q for q in aaea if q != p], tmap)
        key = (p, tuple(g))
        if key not in _dcache:
            ref_arr = dk.within_sample_normalize_matrix(ref_mat.reindex(g).fillna(0.0), "z").values.T
            qv = dk.normalize_query_vector(qprof[p].reindex(g).fillna(0.0).values, "within_patient_z")
            _dcache[key] = {m: np.array([dk.compute_dist(qv, r, m) for r in ref_arr]) for m in METRICS}
        preds[p] = dk.predict_balanced(_dcache[key][metric], ref_lab, "classmean", None,
                                       "subsample", ROUNDS, rng)[0]
    return bal_acc(tmap, preds, aaea)


# ---- run ------------------------------------------------------------------------
nh = lambda q: sum(truth[p] == "high" for p in q)
print(f"reference: {len(ref_lab)} samples ({(ref_lab == 'high').sum()} high / "
      f"{(ref_lab == 'low').sum()} low); fixed panel: {len(fixed_genes)} genes")
print(f"full cohort: {len(pids)} ({nh(pids)}H/{len(pids) - nh(pids)}L); "
      f"AA+EA: {len(aaea)} ({nh(aaea)}H/{len(aaea) - nh(aaea)}L); "
      f"nested gene universe: {len(universe)} of {len(panel)} reference genes\n")

full_allocs = list(allocations_full())
aaea_allocs = list(allocations_aaea())
rows = []

res = {}
for m in METRICS:
    preds = self_z_preds(pids, lambda p: sorted((set(OUR_CANDIDATES) | listB(exclude=p)) & panel), m)
    obs = bal_acc(truth, preds, pids)
    res[m] = (obs,) + exact_p(obs, [bal_acc(a, preds, pids) for a in full_allocs])
rows.append(("Full cohort, fixed panel", res))

res = {}
for m in METRICS:
    preds = self_z_preds(aaea, lambda p: fixed_genes, m)
    obs = bal_acc(truth, preds, aaea)
    res[m] = (obs,) + exact_p(obs, [bal_acc(a, preds, aaea) for a in aaea_allocs])
rows.append(("AA+EA, fixed panel", res))

res = {}
for m in METRICS:
    obs = nested_balacc(truth, m)
    res[m] = (obs,) + exact_p(obs, [nested_balacc(a, m) for a in aaea_allocs])
rows.append(("AA+EA, nested selection", res))

res = {}
for m in METRICS:
    preds = cohort_z_preds(aaea, m)
    obs = bal_acc(truth, preds, aaea)
    res[m] = (obs,) + exact_p(obs, [bal_acc(a, preds, aaea) for a in aaea_allocs])
rows.append(("AA+EA, cohort-standardized", res))

print(f"{'Table 2b row':<30}{'cosine bal (p)':>20}{'Spearman bal (p)':>20}{'allocations':>13}")
for name, r in rows:
    c, s = r["cosine"], r["spearman"]
    print(f"{name:<30}{f'{c[0]:.3f} ({c[1]:.3f})':>20}{f'{s[0]:.3f} ({s[1]:.3f})':>20}{c[2]:>13,}")
