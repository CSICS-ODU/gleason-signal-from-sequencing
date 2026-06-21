"""
compare_signatures.py - in-house, leave-one-out severity classification of the 20
patients, comparing two SIGNATURE representations on the same harness.

Same concept as the prior distance/classification test (genefilelists/distance.py):
represent each patient as a vector over a gene subspace, compare by distance, assign
to the nearer severity reference. TWO differences, both deliberate:

  1. NO AVERAGING. The old high/low references were averaged over an external
     55-patient cohort -- averaging collapses molecular heterogeneity into a phantom
     mean patient, so high and low centroids regress together. Here the reference is
     the SET of real in-house patients: each held-out patient is scored against its
     NEAREST real high patient vs nearest real low patient (leave-one-out, k=1).
  2. CANDIDATE SUBSPACE. Restrict to Layer-1 disproportion genes, not the ~17k-gene
     universe, so the discriminative signal isn't diluted by passenger/coverage noise.

Two representations compared (the head-to-head of invariants):
  * relief    -- per-gene variant COUNT, normalized to sum 1 (the magnitude
                 topography hypothesis).
  * relational-- BINARY presence (which candidate genes are mutated; units discarded,
                 only co-occurrence structure of the gene set retained).
Both scored with the same four metrics (cosine/canberra/euclidean/spearman) so the
only thing that varies is the vector. Cohort = PT510=low (13 high / 7 low).

Usage:
  python compare_signatures.py                         # PT510_low matched ALL master
  python compare_signatures.py <master_genes.csv> --top-k 40 --min-total 4
  python compare_signatures.py <f> --genes APC,ATM,TP53,...   # explicit subspace
  python compare_signatures.py <f> --all-genes         # reproduce the diluted version
"""

import argparse
import csv
import glob
from collections import defaultdict

import numpy as np
from scipy.spatial.distance import canberra, cosine, euclidean
from scipy.stats import spearmanr

METRICS = ["cosine", "canberra", "euclidean", "spearman"]


def split_field(val):
    return [x for x in (val or "").split(";") if x]


def load_master(path):
    """SNV master -> per-patient per-gene variant COUNTS, patient classes, cohort.

    count[patient][gene] = number of that gene's SNVs the patient carries (the relief
    magnitude); presence = count > 0. A patient's class is the list it appears in.
    """
    counts = defaultdict(lambda: defaultdict(int))
    high_all, low_all = set(), set()
    high_car, low_car = defaultdict(set), defaultdict(set)
    saw_genes = False
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if "genes" in r:
                saw_genes = True
            hp, lp = split_field(r.get("high_patients")), split_field(r.get("low_patients"))
            high_all.update(hp)
            low_all.update(lp)
            for g in split_field(r.get("genes")):
                for p in hp:
                    counts[p][g] += 1
                for p in lp:
                    counts[p][g] += 1
                if hp:
                    high_car[g].update(hp)
                if lp:
                    low_car[g].update(lp)
    if not saw_genes:
        raise SystemExit(f"{path} has no `genes` column")
    cls = {p: "high" for p in high_all}
    cls.update({p: "low" for p in low_all - high_all})
    return {"counts": counts, "cls": cls,
            "H": len(high_all), "L": len(low_all),
            "high_car": high_car, "low_car": low_car}


def pick_candidates(data, top_k, min_total, explicit, all_genes):
    genes = set(data["high_car"]) | set(data["low_car"])
    if all_genes:
        return sorted(genes)
    if explicit:
        return [g for g in explicit if g in genes]
    H, L = data["H"], data["L"]
    scored = []
    for g in genes:
        nh, nl = len(data["high_car"].get(g, ())), len(data["low_car"].get(g, ()))
        if nh + nl <= min_total:
            continue
        scored.append((g, nh / H - nl / L))
    highs = sorted(scored, key=lambda s: -s[1])[:top_k]
    lows = sorted(scored, key=lambda s: s[1])[:top_k]
    return sorted({s[0] for s in highs + lows})


def vectors(data, cands, representation):
    """Patient -> normalized vector over the candidate subspace."""
    out = {}
    for p in data["cls"]:
        raw = np.array([data["counts"][p].get(g, 0) for g in cands], dtype=float)
        if representation == "relational":
            raw = (raw > 0).astype(float)
        s = raw.sum()
        out[p] = raw / s if s > 0 else raw  # sum-1 normalize (relief is a distribution)
    return out


def dist(a, b, metric):
    if a.sum() == 0 or b.sum() == 0:
        return np.inf                       # carries none of the subspace -> no signal
    if metric == "cosine":
        return cosine(a, b)
    if metric == "canberra":
        return canberra(a, b)
    if metric == "euclidean":
        return euclidean(a, b)
    coef, _ = spearmanr(a, b)
    return 1 - coef if not np.isnan(coef) else np.inf


def classify(data, vecs, metric):
    """Leave-one-out, nearest REAL high vs nearest REAL low (no averaging)."""
    cls = data["cls"]
    patients = list(cls)
    rows = []
    for p in patients:
        others = [q for q in patients if q != p]
        d_high = min((dist(vecs[p], vecs[q], metric) for q in others if cls[q] == "high"),
                     default=np.inf)
        d_low = min((dist(vecs[p], vecs[q], metric) for q in others if cls[q] == "low"),
                    default=np.inf)
        delta = d_low - d_high                       # >0 -> closer to a high patient
        pred = "high" if delta > 0 else "low"
        rows.append({"patient": p, "true": cls[p], "pred": pred,
                     "d_high": d_high, "d_low": d_low, "correct": pred == cls[p]})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("master", nargs="?")
    ap.add_argument("--top-k", type=int, default=40,
                    help="top-K HIGH and top-K LOW Layer-1 genes (effect size)")
    ap.add_argument("--min-total", type=int, default=4)
    ap.add_argument("--genes", help="explicit candidate subspace, comma-separated")
    ap.add_argument("--all-genes", action="store_true",
                    help="use full gene universe (reproduces the diluted version)")
    a = ap.parse_args()

    path = a.master or next(iter(sorted(glob.glob(
        "snv_rebuild/PT510_low/matched/ALL/**/ALL_master_matched*_genes.csv",
        recursive=True))), None)
    if not path:
        raise SystemExit("no PT510_low master found; pass one explicitly")

    data = load_master(path)
    explicit = [g.strip() for g in a.genes.split(",")] if a.genes else None
    cands = pick_candidates(data, a.top_k, a.min_total, explicit, a.all_genes)
    if len(cands) < 2:
        raise SystemExit("need >=2 candidate genes")

    print("=" * 64)
    print(f"{path}\n  cohort: {data['H']} high / {data['L']} low | "
          f"subspace: {len(cands)} genes | majority baseline: "
          f"{max(data['H'], data['L'])}/{data['H'] + data['L']} = "
          f"{max(data['H'], data['L']) / (data['H'] + data['L']):.2f}")
    print(f"\n  {'representation':<14}{'metric':<11}{'correct':>9}{'accuracy':>10}")
    results = {}
    for rep in ("relief", "relational"):
        vecs = vectors(data, cands, rep)
        for metric in METRICS:
            rows = classify(data, vecs, metric)
            nc = sum(r["correct"] for r in rows)
            acc = nc / len(rows)
            results[(rep, metric)] = rows
            print(f"  {rep:<14}{metric:<11}{nc:>6}/{len(rows):<2}{acc:>10.2f}")

    # save per-patient predictions for the best cell of each representation
    for rep in ("relief", "relational"):
        best = max(METRICS, key=lambda m: sum(r["correct"] for r in results[(rep, m)]))
        with open(f"compare_{rep}_{best}_preds.csv", "w", newline="",
                  encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["patient", "true", "pred", "d_high", "d_low", "correct"])
            for r in results[(rep, best)]:
                w.writerow([r["patient"], r["true"], r["pred"],
                            f"{r['d_high']:.4f}", f"{r['d_low']:.4f}", r["correct"]])
        print(f"  wrote compare_{rep}_{best}_preds.csv (best {rep} metric: {best})")


if __name__ == "__main__":
    main()
