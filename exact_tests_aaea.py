# -*- coding: utf-8 -*-
"""exact_tests_aaea.py -- the exact population-preserving tests on the whole-genome-
sequenced AA+EA patients (4 high / 6 low), plus the gene-level correlation check.

  1  Variant burden (Results, "Variants and genes found in only one grade"; Discussion):
     per patient, the catalogued, non-benign variants found in >=2 patients summed over
     genes (cohort_matched/ variant_count), and the number of genes with >=1 such
     variant. Statistic: mean(low) - mean(high). Also reported: distinct variants per
     patient (a variant in overlapping genes counted once), from the SNV master table.
  2  Fig. 3b p and q: for each candidate gene, does the high-minus-low difference in
     the fraction of patients with a variant differ between AA and EA?
     Statistic |(h_EA/2 - l_EA/3) - (h_AA/2 - l_AA/3)|; Benjamini-Hochberg q over the
     eight candidate genes.
  3  Gene-level Spearman correlation between mean reference expression (GSE244631)
     and mean variant count (Results, "Predicting Gleason grade from an external
     dataset": "near zero").

Null for 1 and 2: all 100 allocations that keep 2 high / 3 low within each of AA and
EA (C(5,2)^2), so no result can arise from population alone; p = fraction of
allocations at least as extreme as observed.
"""
import csv
import glob
import types
import itertools
import functools
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import distance_knn as dk
print = functools.partial(print, flush=True)

CANDIDATES = ["PTEN", "TP53", "PIK3CA", "SMAD4", "GNAS", "APC", "ATM", "CDH1"]
MASTER = "snv_rebuild/PT510_low/matched/ALL/**/ALL_master_matched_genes.csv"

qprof = dk.load_query_profiles("cohort_matched")
truth = dk.load_query_truth("genefilelists/metadata_orig.csv", "ID", "Gleason")
races = dk.load_query_races("genefilelists/metadata_orig.csv", "ID", "race")
AA = [p for p in qprof if races.get(p) == "AA"]
EA = [p for p in qprof if races.get(p) == "EA"]
pts = AA + EA
high = {p for p in pts if truth[p] == "high"}


def allocations():
    na = sum(p in high for p in AA)
    ne = sum(p in high for p in EA)
    for ha in itertools.combinations(AA, na):
        for he in itertools.combinations(EA, ne):
            yield set(ha) | set(he)


allocs = list(allocations())
print(f"AA+EA: {len(pts)} patients ({len(high)} high / {len(pts) - len(high)} low); "
      f"{len(allocs)} population-preserving allocations\n")

# ---- 1  variant burden ------------------------------------------------------------
summed = {p: float(qprof[p].sum()) for p in pts}
ngenes = {p: int((qprof[p] > 0).sum()) for p in pts}
distinct = {p: 0 for p in pts}
with open(glob.glob(MASTER, recursive=True)[0], newline="") as fh:
    for r in csv.DictReader(fh):
        for p in set(r["high_patients"].split(";") + r["low_patients"].split(";")) & set(pts):
            distinct[p] += 1


def low_minus_high(x, H):
    return (np.mean([x[p] for p in pts if p not in H]) - np.mean([x[p] for p in H]))


print("1  Variant burden, high vs low grade (exact one-sided p, low > high)")
for name, x in (("variants summed over genes", summed), ("genes with a variant", ngenes),
                ("distinct variants", distinct)):
    obs = low_minus_high(x, high)
    p = np.mean([low_minus_high(x, H) >= obs - 1e-9 for H in allocs])
    mh = np.mean([x[q] for q in high]); ml = np.mean([x[q] for q in pts if q not in high])
    print(f"   {name:<28} high {mh:>9,.1f}   low {ml:>9,.1f}   p = {p:.3f}")

# ---- 2  Fig. 3b: AA-versus-EA difference in the high-minus-low lean ----------------
carrier = {p: set(qprof[p].index[qprof[p] > 0]) for p in pts}


def lean_diff(g, H):
    def lean(grp):
        h = [p for p in grp if p in H]
        l = [p for p in grp if p not in H]
        return np.mean([g in carrier[p] for p in h]) - np.mean([g in carrier[p] for p in l])
    return lean(EA) - lean(AA)


def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    q = np.empty_like(p)
    prev = 1.0
    for rank in range(len(p) - 1, -1, -1):
        prev = min(prev, p[order[rank]] * len(p) / (rank + 1))
        q[order[rank]] = prev
    return q


pvals = []
for g in CANDIDATES:
    obs = abs(lean_diff(g, high))
    pvals.append(np.mean([abs(lean_diff(g, H)) >= obs - 1e-9 for H in allocs]))
qvals = bh(pvals)
print("\n2  Fig. 3b: patients with a variant (H/L) and the AA-vs-EA test")
print(f"   {'gene':<8}{'AA H/L':>8}{'EA H/L':>8}{'p':>7}{'q':>7}")
for g, p, q in zip(CANDIDATES, pvals, qvals):
    hl = lambda grp: (f"{sum(g in carrier[x] for x in grp if x in high)}/"
                      f"{sum(g in carrier[x] for x in grp if x not in high)}")
    print(f"   {g:<8}{hl(AA):>8}{hl(EA):>8}{p:>7.2f}{q:>7.2f}")

# ---- 3  gene-level correlation between the two data types --------------------------
ref = dk.load_reference(types.SimpleNamespace(ref_rcc_dir="GSE244631", ref_matrix=None))
ref_labels = dk.load_ref_labels("GSE244631_ref_labels.csv")
ref = ref[[s for s in ref.columns if s in ref_labels]]
print("\n3  Spearman correlation across genes: mean reference expression vs mean variant count")
for name, group in (("full cohort", list(qprof)), ("AA+EA", pts)):
    genes = sorted(set(ref.index) & set().union(*[set(qprof[p].index) for p in group]))
    Q = pd.DataFrame({p: qprof[p].reindex(genes).fillna(0.0) for p in group})
    rho = spearmanr(ref.loc[genes].mean(axis=1), Q.mean(axis=1))[0]
    print(f"   {name:<12} {len(genes)} genes   rho = {rho:+.2f}")
