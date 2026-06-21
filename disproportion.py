"""
disproportion.py - among recurrent SNVs (n_total > threshold) in an annotated
master, find those carried DISPROPORTIONATELY by high- vs low-severity patients.

Why rates, not counts: the cohort is imbalanced (e.g. ALL = 14 high / 6 low). A SNV
in 10 high + 4 low looks "high-heavy" but is ~71% of high and ~67% of low -- NOT
skewed. So we compare per-grade RATES (n_high/H vs n_low/L) and test the 2x2 table
  [[high carriers, high non-carriers], [low carriers, low non-carriers]]
with Fisher's exact test (odds ratio + p-value).

Usage:
  python disproportion.py                          # default: PT510_low ALL masters,
                                                   #   benign + matched + unmatched
  python disproportion.py <file_genes.csv> ...     # specific annotated master(s)
  python disproportion.py <f> --min-total 10 --top 25
"""

import argparse
import csv
import glob

from scipy.stats import fisher_exact


def keep(pid, no_indian):
    """Indian targeted-panel patients have all-numeric IDs; AA/EA are PT-prefixed.
    The IN panel is a near-constant germline readout (see cohort_coverage_analysis.py),
    so --no-indian restricts to the AA+EA broad-panel stratum that actually varies."""
    return not (no_indian and pid.isdigit())


def cohort_size(rows, no_indian=False):
    """H, L = distinct high / low patients appearing in the file."""
    high, low = set(), set()
    for r in rows:
        high |= {p for p in (r.get("high_patients") or "").split(";") if p and keep(p, no_indian)}
        low |= {p for p in (r.get("low_patients") or "").split(";") if p and keep(p, no_indian)}
    return len(high), len(low)


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
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    H, L = cohort_size(rows, no_indian)
    out = []
    for r in rows:
        # recompute counts from the patient lists so --no-indian is exact (the stored
        # n_high/n_low columns include Indian patients).
        hp = [x for x in (r.get("high_patients") or "").split(";") if x and keep(x, no_indian)]
        lp = [x for x in (r.get("low_patients") or "").split(";") if x and keep(x, no_indian)]
        nh, nl = len(hp), len(lp)
        nt = nh + nl
        if nt <= min_total:
            continue
        odds, p = fisher_exact([[nh, H - nh], [nl, L - nl]])
        hr, lr = (nh / H if H else 0.0), (nl / L if L else 0.0)
        out.append({
            "snv": r["SNV"], "gene": r.get("genes", ""), "cat": r.get("category", ""),
            "nh": nh, "nl": nl, "hr": hr, "lr": lr, "diff": hr - lr,
            "or": odds, "p": p, "dir": "HIGH" if hr > lr else ("LOW" if lr > hr else "="),
        })
    add_qvalues(out)
    return H, L, out


def show(title, items, n):
    print(f"\n  {title}  (top {n}, ranked by effect size)")
    print(f"    {'SNV':<22}{'gene':<12}{'hi':>6}{'lo':>6}{'high%':>7}{'low%':>7}"
          f"{'OR':>8}{'p':>9}{'q':>9}")
    for d in items[:n]:
        orr = "inf" if d["or"] == float("inf") else f"{d['or']:.2f}"
        print(f"    {d['snv']:<22}{(d['gene'] or '-')[:11]:<12}"
              f"{d['nh']:>6}{d['nl']:>6}{d['hr']*100:>6.0f}%{d['lr']*100:>6.0f}%"
              f"{orr:>8}{d['p']:>9.2g}{d['q']:>9.2g}")


def default_files():
    """The ALL master per layer (PT510_low): prefer the annotated *_genes.csv,
    fall back to the plain master (e.g. unmatched, which has no gene annotation)."""
    files = []
    for layer in ("benign", "matched", "unmatched"):
        found = glob.glob(
            f"snv_rebuild/PT510_low/{layer}/ALL/**/ALL_master_{layer}*.csv",
            recursive=True)
        genes = [f for f in found if f.endswith("_genes.csv")]
        plain = [f for f in found if not f.endswith("_genes.csv")]
        files.extend(genes if genes else plain)
    return sorted(files)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="*",
                    help="annotated *_genes.csv (default: PT510_low ALL masters)")
    ap.add_argument("--min-total", type=int, default=10,
                    help="keep SNVs with n_total strictly greater than this")
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--no-indian", action="store_true",
                    help="exclude Indian (numeric-ID) patients -> AA+EA broad-panel stratum")
    ap.add_argument("--exclusive", action="store_true",
                    help="keep only high-only (n_low=0) and low-only (n_high=0) SNVs; "
                         "drop anything carried by BOTH high and low (likely shared/germline)")
    a = ap.parse_args()

    files = a.files or default_files()
    if not files:
        raise SystemExit("no ALL master files found; pass paths explicitly")

    for path in files:
        # unmatched tops out near 10 (Indians contribute ~no unmatched variants),
        # so a >10 cutoff finds nothing -- use 5 there.
        mt = 5 if "unmatched" in path else a.min_total
        H, L, res = analyze(path, mt, a.no_indian)
        if a.exclusive:
            res = [d for d in res if d["nh"] == 0 or d["nl"] == 0]
            add_qvalues(res)   # recompute FDR over the reduced (exclusive) set
        print("=" * 78)
        print(f"{path}\n  cohort: {H} high / {L} low | "
              f"SNVs with n_total>{mt}: {len(res)}"
              f"{' [exclusive: high-only + low-only]' if a.exclusive else ''}")
        # rank by effect size (rate gap), p as tiebreak -- in a small cohort the
        # magnitude of separation matters more than the nominal p.
        high = sorted([d for d in res if d["dir"] == "HIGH"],
                      key=lambda d: (-d["diff"], d["p"]))
        low = sorted([d for d in res if d["dir"] == "LOW"],
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
