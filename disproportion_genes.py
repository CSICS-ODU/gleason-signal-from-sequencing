"""
disproportion_genes.py - the GENE-level analogue of disproportion.py.

disproportion.py tests individual SNVs for high-vs-low skew. This rolls the same
test up to genes: a patient "carries" a gene if they carry >=1 SNV in it
(presence/absence union -- NOT a count, so a gene with 5 SNVs in one patient still
counts as one carrier). For each recurrent gene we compare per-grade RATES
(n_high/H vs n_low/L) and test the 2x2 carrier table
  [[high carriers, high non-carriers], [low carriers, low non-carriers]]
with Fisher's exact test (odds ratio + p-value).

Why rebuild from the SNV master instead of the gene matrices:
  The gene-results CSVs (high_only/low_only/high_low_common) are partitioned
  PER SNV, so one gene can appear in several of them with per-category carrier
  counts -- concatenating double-counts. The SNV master carries, per SNV, its
  `genes`, `high_patients`, `low_patients`, so a clean union per gene is exact.

Usage:
  python disproportion_genes.py                       # default: PT510_low ALL
                                                      #   annotated masters,
                                                      #   benign + matched + unmatched
  python disproportion_genes.py <ALL_master_*_genes.csv> ...
  python disproportion_genes.py <f> --min-total 10 --top 25
"""

import argparse
import csv
import glob
from collections import defaultdict

from scipy.stats import fisher_exact


def split_field(val):
    """Master list-fields are ;-delimited; tolerate stray blanks."""
    return [x for x in (val or "").split(";") if x]


def keep(pid, no_indian):
    """Indian targeted-panel patients have all-numeric IDs; AA/EA are PT-prefixed.
    --no-indian restricts to the AA+EA broad-panel stratum (the IN panel is a
    near-constant germline readout; see cohort_coverage_analysis.py)."""
    return not (no_indian and pid.isdigit())


def aggregate_genes(path, no_indian=False):
    """SNV master -> per-gene carrier sets. Returns (H, L, {gene: (high_set, low_set)}).

    H, L = distinct high / low patients anywhere in the file (the cohort the rates
    are taken over). A patient is a gene carrier if they carry any SNV in the gene.
    """
    high_all, low_all = set(), set()
    high_car = defaultdict(set)
    low_car = defaultdict(set)
    saw_genes_col = False
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if "genes" in r:
                saw_genes_col = True
            hp = [x for x in split_field(r.get("high_patients")) if keep(x, no_indian)]
            lp = [x for x in split_field(r.get("low_patients")) if keep(x, no_indian)]
            high_all.update(hp)
            low_all.update(lp)
            for g in split_field(r.get("genes")):
                if hp:
                    high_car[g].update(hp)
                if lp:
                    low_car[g].update(lp)
    if not saw_genes_col:
        return None  # plain (unannotated) master -- no gene column to roll up
    genes = {g: (high_car.get(g, set()), low_car.get(g, set()))
             for g in set(high_car) | set(low_car)}
    return len(high_all), len(low_all), genes


def add_qvalues(items):
    """Benjamini-Hochberg FDR q-values across ALL tested items (both directions
    together -- every Fisher test we ran is one hypothesis). q is monotone in p."""
    m = len(items)
    if not m:
        return
    order = sorted(range(m), key=lambda i: items[i]["p"])
    running = 1.0
    for rank in range(m, 0, -1):           # largest p first, cumulative-min downward
        i = order[rank - 1]
        running = min(running, items[i]["p"] * m / rank)
        items[i]["q"] = min(running, 1.0)


def analyze(path, min_total, no_indian=False):
    agg = aggregate_genes(path, no_indian)
    if agg is None:
        return None
    H, L, genes = agg
    out = []
    for gene, (hset, lset) in genes.items():
        nh, nl = len(hset), len(lset)
        nt = nh + nl
        if nt <= min_total:
            continue
        odds, p = fisher_exact([[nh, H - nh], [nl, L - nl]])
        hr, lr = (nh / H if H else 0.0), (nl / L if L else 0.0)
        out.append({
            "gene": gene, "nh": nh, "nl": nl, "hr": hr, "lr": lr, "diff": hr - lr,
            "or": odds, "p": p, "dir": "HIGH" if hr > lr else ("LOW" if lr > hr else "="),
        })
    add_qvalues(out)
    return H, L, out


def show(title, items, n):
    print(f"\n  {title}  (top {n}, ranked by effect size)")
    print(f"    {'gene':<16}{'hi':>6}{'lo':>6}{'high%':>7}{'low%':>7}"
          f"{'OR':>8}{'p':>9}{'q':>9}")
    for d in items[:n]:
        orr = "inf" if d["or"] == float("inf") else f"{d['or']:.2f}"
        print(f"    {d['gene'][:15]:<16}"
              f"{d['nh']:>6}{d['nl']:>6}{d['hr']*100:>6.0f}%{d['lr']*100:>6.0f}%"
              f"{orr:>8}{d['p']:>9.2g}{d['q']:>9.2g}")


def default_files():
    """The annotated ALL master per layer (PT510_low). Only *_genes.csv masters
    carry a `genes` column, so the plain unmatched master is skipped here."""
    files = []
    for layer in ("benign", "matched", "unmatched"):
        files.extend(glob.glob(
            f"snv_rebuild/PT510_low/{layer}/ALL/**/ALL_master_{layer}*_genes.csv",
            recursive=True))
    return sorted(files)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="*",
                    help="annotated ALL_master_*_genes.csv (default: PT510_low ALL)")
    ap.add_argument("--min-total", type=int, default=10,
                    help="keep genes carried by strictly more than this many patients")
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--no-indian", action="store_true",
                    help="exclude Indian (numeric-ID) patients -> AA+EA broad-panel stratum")
    ap.add_argument("--exclusive", action="store_true",
                    help="keep only high-only (n_low=0) and low-only (n_high=0) genes; "
                         "drop anything carried by BOTH high and low (likely shared/germline)")
    a = ap.parse_args()

    files = a.files or default_files()
    if not files:
        raise SystemExit("no annotated ALL master files found; pass paths explicitly")

    for path in files:
        # unmatched tops out low (Indians contribute ~no unmatched variants),
        # so a >10 cutoff finds nothing -- use 5 there, as in disproportion.py.
        mt = 5 if "unmatched" in path else a.min_total
        res = analyze(path, mt, a.no_indian)
        print("=" * 70)
        if res is None:
            print(f"{path}\n  (no `genes` column -- skipped)")
            continue
        H, L, items = res
        if a.exclusive:
            items = [d for d in items if d["nh"] == 0 or d["nl"] == 0]
            add_qvalues(items)   # recompute FDR over the reduced (exclusive) set
        print(f"{path}\n  cohort: {H} high / {L} low | "
              f"genes carried by >{mt} patients: {len(items)}"
              f"{' [exclusive: high-only + low-only]' if a.exclusive else ''}")
        # rank by effect size (rate gap), p as tiebreak -- in a small cohort the
        # magnitude of separation matters more than the nominal p.
        high = sorted([d for d in items if d["dir"] == "HIGH"],
                      key=lambda d: (-d["diff"], d["p"]))
        low = sorted([d for d in items if d["dir"] == "LOW"],
                     key=lambda d: (d["diff"], d["p"]))
        sig_h = sum(1 for d in high if d["p"] < a.alpha)
        sig_l = sum(1 for d in low if d["p"] < a.alpha)
        q_h = sum(1 for d in high if d["q"] < a.alpha)
        q_l = sum(1 for d in low if d["q"] < a.alpha)
        print(f"  high-skewed: {len(high)} (p<{a.alpha}: {sig_h}, q<{a.alpha}: {q_h}) | "
              f"low-skewed: {len(low)} (p<{a.alpha}: {sig_l}, q<{a.alpha}: {q_l})")
        show("HIGH-skewed (over-represented in high-severity)", high, a.top)
        show("LOW-skewed (over-represented in low-severity)", low, a.top)


if __name__ == "__main__":
    main()
