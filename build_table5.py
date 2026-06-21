"""Regenerate Table 5 (tab:gene_map): GENE-level exclusive high-only / low-only
recurrent genes in the pathogenic (matched) bucket, AA+EA stratum, >=3 carriers,
ranked by effect size. Mirrors disproportion_genes.py --no-indian --min-total 2
--exclusive. A patient carries a gene if it carries >=1 matched SNV in it.
"""
import csv
import glob
from collections import defaultdict
from scipy.stats import fisher_exact

MASTER = glob.glob("snv_rebuild/PT510_low/matched/ALL/**/ALL_master_matched_genes.csv",
                   recursive=True)[0]
MIN_CARRIERS = 3


def split(v):
    return [x for x in (v or "").split(";") if x]


def aa_ea(ids):
    return [p for p in ids if not p.isdigit()]


hc, lc = defaultdict(set), defaultdict(set)
H, L = set(), set()
with open(MASTER, newline="", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        hp, lp = aa_ea(split(r.get("high_patients"))), aa_ea(split(r.get("low_patients")))
        H.update(hp); L.update(lp)
        for g in split(r.get("genes")):
            if hp:
                hc[g].update(hp)
            if lp:
                lc[g].update(lp)
H, L = len(H), len(L)

out = []
for g in set(hc) | set(lc):
    nh, nl = len(hc.get(g, set())), len(lc.get(g, set()))
    if nh + nl < MIN_CARRIERS:
        continue
    if nh and nl:                       # exclusive only
        continue
    odds, p = fisher_exact([[nh, H - nh], [nl, L - nl]])
    out.append({"gene": g, "nh": nh, "nl": nl, "diff": nh / H - nl / L,
                "or": odds, "p": p, "dir": "high-only" if nl == 0 else "low-only"})

m = len(out)
order = sorted(range(m), key=lambda i: out[i]["p"])
run = 1.0
for rank in range(m, 0, -1):
    i = order[rank - 1]
    run = min(run, out[i]["p"] * m / rank)
    out[i]["q"] = min(run, 1.0)

hi = sorted([d for d in out if d["dir"] == "high-only"], key=lambda d: (-d["diff"], d["p"]))
lo = sorted([d for d in out if d["dir"] == "low-only"], key=lambda d: (d["diff"], d["p"]))

with open("paper_results/table5_exclusive_aaea.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["set", "gene", "n_high", "n_low", "rate_gap", "OR", "p", "q"])
    for d in hi + lo:
        orr = "inf" if d["or"] == float("inf") else f"{d['or']:.2f}"
        w.writerow([d["dir"], d["gene"], d["nh"], d["nl"], f"{d['diff']:.2f}",
                    orr, f"{d['p']:.4f}", f"{d['q']:.3f}"])

print(f"cohort AA+EA: {H} high / {L} low | exclusive recurrent (>={MIN_CARRIERS}) "
      f"matched GENES: {len(out)}  (high-only {len(hi)}, low-only {len(lo)})")
print(f"min q: {min(d['q'] for d in out):.3f}  (q<0.05: {sum(d['q'] < 0.05 for d in out)})\n")
for title, items in (("HIGH-only genes (top 18)", hi[:18]), ("LOW-only genes (top 18)", lo[:18])):
    print(title)
    for d in items:
        print(f"  {d['gene'][:20]:<22}{d['nh']}/{d['nl']:<3} gap={d['diff']:+.2f} "
              f"p={d['p']:.4f} q={d['q']:.3f}")
    print()
