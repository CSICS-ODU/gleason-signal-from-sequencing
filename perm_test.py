# -*- coding: utf-8 -*-
"""Permutation test for the suspicious nested k=3 spike.

Shuffle the query high/low labels, rerun the FULL nested-LOO pipeline (biomarker
selection inside each fold + balanced kNN classification) for one config, and
compare the observed balanced accuracy to the null distribution. A lone-k peak
that is real should beat the null; if p is large it was an artifact of the
vote mechanics at that k.
"""
import sys
import types
import numpy as np
import distance_knn as dk

CFG = dict(metric="cosine", classifier="knn", k=3)
ROUNDS = 100
N_PERM = 400
N_EACH = 11

ref_mat = dk.load_reference(types.SimpleNamespace(ref_rcc_dir="GSE244631", ref_matrix=None))
ref_labels = dk.load_ref_labels("GSE244631_ref_labels.csv")
qprof = dk.load_query_profiles("cohort_samples")
truth0 = dk.load_query_truth("genefilelists/metadata_orig.csv", "ID", "Gleason")
if len(sys.argv) > 1:
    keep = {r.strip().upper() for r in sys.argv[1].split(",")}
    races = dk.load_query_races("genefilelists/metadata_orig.csv", "ID", "race")
    qprof = {p: v for p, v in qprof.items() if races.get(p) in keep}
truth0 = {p: v for p, v in truth0.items() if p in qprof}

pids = [p for p in qprof if p in truth0]
ref_samples = [s for s in ref_mat.columns if s in ref_labels]
ref_mat = ref_mat[ref_samples]
ref_lab = np.array([ref_labels[s] for s in ref_samples])
panel_genes = sorted(set(ref_mat.index) & set().union(*[set(qprof[p].index) for p in pids]))
carrier = {p: {g: (qprof[p].get(g, 0) > 0) for g in panel_genes} for p in pids}
# precompute normalized query signature is gene-set dependent, so done per fold


def run(labels, rng):
    """nested-LOO balanced accuracy under a given label dict."""
    def select(train):
        hi = [p for p in train if labels[p] == "high"]
        lo = [p for p in train if labels[p] == "low"]
        if not hi or not lo:
            return panel_genes[:2 * N_EACH]
        gap = {g: (sum(carrier[p][g] for p in hi) / len(hi)
                   - sum(carrier[p][g] for p in lo) / len(lo)) for g in panel_genes}
        order = sorted(panel_genes, key=lambda g: gap[g])
        return order[:N_EACH] + order[-N_EACH:]

    preds, tr = [], []
    for p in pids:
        genes = select([q for q in pids if q != p])
        ref_norm = dk.within_sample_normalize_matrix(ref_mat.reindex(genes).fillna(0.0), "z")
        ref_arr = ref_norm.values.T
        qv = dk.normalize_query_vector(qprof[p].reindex(genes).fillna(0.0).values,
                                       "within_patient_z")
        dists = np.array([dk.compute_dist(qv, ref_arr[j], CFG["metric"])
                          for j in range(len(ref_arr))])
        pred, _ = dk.predict_balanced(dists, ref_lab, CFG["classifier"], CFG["k"],
                                      "subsample", ROUNDS, rng)
        preds.append(pred); tr.append(labels[p])
    preds, tr = np.array(preds), np.array(tr)
    recs = [(preds[tr == c] == c).mean() for c in ("high", "low") if (tr == c).any()]
    return float(np.mean(recs))


rng = np.random.default_rng(0)
obs = run(truth0, rng)
labs = np.array([truth0[p] for p in pids])
null = []
for _ in range(N_PERM):
    perm = {p: l for p, l in zip(pids, rng.permutation(labs))}
    null.append(run(perm, rng))
null = np.array(null)
p = (null >= obs).mean()
print(f"config: nested {CFG} rounds={ROUNDS}")
print(f"observed balanced_acc = {obs:.3f}")
print(f"null: mean={null.mean():.3f} sd={null.std():.3f} "
      f"95th={np.percentile(null,95):.3f} max={null.max():.3f}")
print(f"permutation p-value = {p:.3f}  (N={N_PERM})")
