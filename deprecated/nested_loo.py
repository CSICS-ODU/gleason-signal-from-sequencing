# -*- coding: utf-8 -*-
"""
nested_loo.py -- the honest version of "filter the reference by our biomarkers".

The biomarkers were selected using the cohort's own high/low labels, so scoring
the same cohort with them leaks. Here we redo selection INSIDE each leave-one-out
fold (never touching the held-out patient), so the gene choice can't peek at the
answer. We compare three things head to head:

  leaky    : select biomarkers on ALL patients (incl. the test one)  <- the 0.83
  nested   : select biomarkers on the n-1 training patients only      <- honest
  k-curve  : balanced accuracy at k=1,3,5 + classmean, to dissect the k=1 worry

Selection mirrors the study's disproportion logic: a gene is "carried" if the
patient has >=1 variant in it (variant_count>0); rate_gap = carrier-rate(high) -
carrier-rate(low); take the top N high-enriched + N low-enriched genes that are
ALSO on the reference panel (only those are measurable in the reference).
"""
import sys
import types
import numpy as np
import pandas as pd
import distance_knn as dk
import functools
print = functools.partial(print, flush=True)   # every line written immediately

N_EACH = 11          # genes per side -> ~22 features, matching biomarker_subset
ROUNDS = 200
SEED = 0
METRICS = ["cosine", "spearman", "euclidean"]
KS = [1, 3, 5]

races_keep = None
_arg1 = sys.argv[1] if len(sys.argv) > 1 else ""
if _arg1 and not _arg1.startswith("perm"):  # e.g. python nested_loo.py AA,EA
    races_keep = {r.strip().upper() for r in _arg1.split(",")}

# ---- load -----------------------------------------------------------------
ref_mat = dk.load_reference(types.SimpleNamespace(ref_rcc_dir="GSE244631", ref_matrix=None))
ref_labels = dk.load_ref_labels("GSE244631_ref_labels.csv")
qprof = dk.load_query_profiles("cohort_matched")   # matched-layer = biomarker representation
truth = dk.load_query_truth("genefilelists/metadata_orig.csv", "ID", "Gleason")
if races_keep:
    races = dk.load_query_races("genefilelists/metadata_orig.csv", "ID", "race")
    qprof = {p: v for p, v in qprof.items() if races.get(p) in races_keep}
truth = {p: v for p, v in truth.items() if p in qprof}

pids = [p for p in qprof if p in truth]
ref_samples = [s for s in ref_mat.columns if s in ref_labels]
ref_mat = ref_mat[ref_samples]
ref_lab = np.array([ref_labels[s] for s in ref_samples])

panel = set(ref_mat.index)
qgenes = set().union(*[set(qprof[p].index) for p in pids])
panel_genes = sorted(panel & qgenes)                 # candidate features
carrier = {p: {g: (qprof[p].get(g, 0) > 0) for g in panel_genes} for p in pids}
print(f"cohort: {len(pids)} patients "
      f"({sum(truth[p]=='high' for p in pids)}h/{sum(truth[p]=='low' for p in pids)}l) "
      f"| candidate panel genes: {len(panel_genes)}")


def select_with(train_ids, n_each, tmap):
    hi = [p for p in train_ids if tmap[p] == "high"]
    lo = [p for p in train_ids if tmap[p] == "low"]
    if not hi or not lo:
        return []
    gap = {g: (sum(carrier[p][g] for p in hi) / len(hi)
               - sum(carrier[p][g] for p in lo) / len(lo)) for g in panel_genes}
    order = sorted(panel_genes, key=lambda g: gap[g])
    return order[:n_each] + order[-n_each:]           # low-enriched + high-enriched


def select(train_ids, n_each):
    return select_with(train_ids, n_each, truth)


def classify(test_pid, genes, metric, classifier, k, rng):
    sub = ref_mat.reindex(genes).fillna(0.0)
    ref_norm = dk.within_sample_normalize_matrix(sub, "z")
    ref_arr = ref_norm.values.T
    qv = dk.normalize_query_vector(
        qprof[test_pid].reindex(genes).fillna(0.0).values, "within_patient_z")
    dists = np.array([dk.compute_dist(qv, ref_arr[j], metric)
                      for j in range(len(ref_arr))])
    pred, _ = dk.predict_balanced(dists, ref_lab, classifier, k,
                                  "subsample", ROUNDS, rng)
    return pred


def bal_acc(rows):
    t = np.array([r[0] for r in rows]); p = np.array([r[1] for r in rows])
    recs = [(p[t == c] == c).mean() for c in ("high", "low") if (t == c).any()]
    acc = (p == t).mean()
    return acc, float(np.mean(recs))


# ---- nested vs leaky -------------------------------------------------------
leaky_genes = select(pids, N_EACH)                    # uses ALL patients (leaks)
fold_genes = {p: select([q for q in pids if q != p], N_EACH) for p in pids}

# gene-set stability across folds (speaks to overfitting / k=1 fragility)
from collections import Counter
cnt = Counter(g for gs in fold_genes.values() for g in gs)
stable = [g for g, c in cnt.items() if c == len(pids)]
print(f"selected genes per fold: {N_EACH*2} | distinct across folds: {len(cnt)} | "
      f"in EVERY fold: {len(stable)}\n")

# class-mean only, cosine + spearman only -- the configs the tables actually use.
print(f"{'mode':<7}{'metric':<10}{'acc':>6}{'bal_acc':>9}")
for mode, getg in (("leaky", lambda p: leaky_genes),
                   ("nested", lambda p: fold_genes[p])):
    for metric in ("cosine", "spearman"):
        rng = np.random.default_rng(SEED)
        rows = [(truth[p], classify(p, getg(p), metric, "classmean", None, rng))
                for p in pids]
        acc, bal = bal_acc(rows)
        print(f"{mode:<7}{metric:<10}{acc:>6.2f}{bal:>9.2f}")
    print()


# ---- all-genes (genome-wide) baseline, same pipeline, no gene selection ----
_rng = np.random.default_rng(SEED)
_rows = [(truth[p], classify(p, panel_genes, "spearman", "classmean", None, _rng))
         for p in pids]
_acc0, _bal0 = bal_acc(_rows)
print(f"all-genes ({len(panel_genes)}) spearman classmean: "
      f"acc {_acc0:.2f}  bal_acc {_bal0:.2f}\n")


# ---- nested permutation significance test (spearman, nearest-class-mean) ----
# Shuffle the query labels, RESELECT genes per fold on the shuffled labels, score
# against the shuffle. This wraps the selection step, so the null reflects the whole
# leaky-prone pipeline -- the honest p for the nested result.
def nested_balacc(tmap, metric="spearman", seed=SEED):
    fg = {p: select_with([q for q in pids if q != p], N_EACH, tmap) for p in pids}
    rng = np.random.default_rng(seed)
    rows = [(tmap[p],
             classify(p, fg[p], metric, "classmean", None, rng) if fg[p] else "high")
            for p in pids]
    return bal_acc(rows)[1]


if any(a.startswith("perm") for a in sys.argv[1:]):
    nperm = next((int(a.split(":", 1)[1]) for a in sys.argv[1:] if a.startswith("perm:")),
                 1000)
    plist = list(pids)
    labels = np.array([truth[p] for p in plist])
    for metric in ("cosine", "spearman"):
        obs = nested_balacc(truth, metric)
        print(f"[{metric}] observed nested bal_acc = {obs:.3f}; running {nperm} shuffles...")
        rng = np.random.default_rng(12345)
        null = np.empty(nperm)
        ge = 0                                   # count of null >= obs (running tally)
        for i in range(nperm):
            perm = rng.permutation(labels)
            null[i] = nested_balacc({p: perm[j] for j, p in enumerate(plist)}, metric, seed=i)
            if null[i] >= obs - 1e-9:
                ge += 1
            if (i + 1) % 25 == 0:                # progress line every 25 shuffles
                print(f"  [{metric}] {i+1}/{nperm} | null>=obs={ge} | "
                      f"running p~{(ge + 1) / (i + 2):.4f}")
        pval = (ge + 1) / (nperm + 1)
        print(f"NESTED PERM ({metric}, classmean), {nperm} shuffles: "
              f"obs {obs:.3f}, null {null.mean():.3f} "
              f"[{np.quantile(null, 0.025):.3f},{np.quantile(null, 0.975):.3f}], p={pval:.4f}")
