"""
cohort_coverage_analysis.py - characterize WHY the Indian (IN) cohort contributes
no discriminative (severity / co-occurrence) signal. Not a salvage; a documented
analysis of the data-generating limitation, for the paper's data-characterization /
limitations section and the cohort-expansion motivation.

Mechanism under test: discriminative analysis needs VARIATION in variant presence
across patients. A targeted panel with no matched-normal subtraction produces a
near-constant, germline-dominated readout -- the same variants in everyone -- which
carries no individual signal regardless of method.

Metrics per ancestry (matched bucket, PT510=low):
  * coverage      -- distinct genes / SNVs (panel breadth)
  * %universal    -- SNVs present in ALL patients of the group (germline/constant share)
  * %private      -- SNVs in exactly ONE patient (individual-specific variation;
                     somatic-rich data has many, a constant readout has ~none)
  * mean prevalence -- average fraction of the group carrying an SNV (cross-n
                     comparable; ->1 means uniform/constant, low means variable)
Then: the IN variable SNVs (the only non-constant content) and whether they can
support any high/low contrast given the 9-high/1-low IN split.
"""

import csv
import glob
from collections import Counter, defaultdict


def load(anc, bucket="matched"):
    snv_car = defaultdict(set)
    gene_car = defaultdict(set)
    rows = []
    for p in glob.glob(
            f"snv_rebuild/PT510_low/{bucket}/{anc}/**/{anc}_master_{bucket}*_genes.csv",
            recursive=True):
        with open(p, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                hp = [x for x in (r["high_patients"] or "").split(";") if x]
                lp = [x for x in (r["low_patients"] or "").split(";") if x]
                ps = hp + lp
                snv_car[r["SNV"]].update(ps)
                genes = [g for g in (r["genes"] or "").split(";") if g]
                for g in genes:
                    gene_car[g].update(ps)
                rows.append((r["SNV"], genes, len(hp), len(lp)))
    pats = set().union(*snv_car.values()) if snv_car else set()
    return snv_car, gene_car, pats, rows


def main():
    print("Per-ancestry coverage & invariance (matched, PT510=low)")
    print(f"  {'anc':<4}{'pts':>5}{'genes':>8}{'SNVs':>8}{'%univ':>7}"
          f"{'%private':>10}{'mean_prev':>11}")
    store = {}
    for anc in ("IN", "AA", "EA"):
        snv_car, gene_car, pats, rows = load(anc)
        store[anc] = (snv_car, gene_car, pats, rows)
        n, m = len(pats), len(snv_car)
        univ = sum(1 for v in snv_car.values() if len(v) == n)
        priv = sum(1 for v in snv_car.values() if len(v) == 1)
        meanprev = sum(len(v) / n for v in snv_car.values()) / m if m else 0
        print(f"  {anc:<4}{n:>5}{len(gene_car):>8}{m:>8}{100*univ/m:>6.0f}%"
              f"{100*priv/m:>9.0f}%{meanprev:>11.2f}")

    # IN variable content: the only non-constant SNVs, and can they discriminate?
    snv_car, gene_car, pats, rows = store["IN"]
    n = len(pats)
    var = [(s, g, nh, nl) for (s, g, nh, nl) in rows if len(snv_car[s]) < n]
    print(f"\nIN: of {len(snv_car)} SNVs, {len(var)} are variable (not in all {n}).")
    genes = Counter(g for _, gs, _, _ in var for g in gs)
    print(f"  variable SNVs sit in {len(genes)} genes; top: "
          f"{', '.join(f'{g}({c})' for g, c in genes.most_common(10))}")
    split = Counter((nh, nl) for _, _, nh, nl in var)
    print(f"  their (n_high, n_low) within IN's 9-high/1-low split:")
    for (nh, nl), c in sorted(split.items(), key=lambda x: -x[1]):
        print(f"    high={nh} low={nl}: {c} SNVs")
    disc = sum(c for (nh, nl), c in split.items() if nl >= 1 and nh == 0)
    print(f"  SNVs that are low-only (the 1 low patient, no highs): {disc} "
          f"-> no severity contrast possible from a single low patient.")

    print("\nConclusion: IN is a near-constant, germline-dominated panel readout "
          "(~all variants in all patients, ~0 private). AA/EA carry the variation; "
          "IN cannot feed severity or co-occurrence analysis by construction.")


if __name__ == "__main__":
    main()
