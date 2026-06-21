"""Regenerate Table 4 (tab:snv_sets): exclusive high-only / low-only recurrent
variants in the pathogenic (matched) bucket, AA+EA stratum, >=3 carriers,
ranked by effect size. Writes paper_results/table4_exclusive_aaea.csv and prints
a representative top-N for the manuscript.
"""
import csv
import glob
from scipy.stats import fisher_exact

MASTER = glob.glob("snv_rebuild/PT510_low/matched/ALL/**/ALL_master_matched_genes.csv",
                   recursive=True)[0]
MIN_CARRIERS = 3          # >=3 carriers (the --min-total 2 run)


def split(v):
    return [x for x in (v or "").split(";") if x]


def aa_ea(ids):                      # drop Indian (numeric-ID) patients
    return [p for p in ids if not p.isdigit()]


rows, H, L = [], set(), set()
with open(MASTER, newline="", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        hp = aa_ea(split(r.get("high_patients")))
        lp = aa_ea(split(r.get("low_patients")))
        H.update(hp); L.update(lp)
        rows.append((r["SNV"], split(r.get("genes")), hp, lp))
H, L = len(H), len(L)

out = []
for snv, genes, hp, lp in rows:
    nh, nl = len(hp), len(lp)
    nt = nh + nl
    if nt < MIN_CARRIERS:
        continue
    if nh and nl:                    # exclusive: drop anything in BOTH
        continue
    odds, p = fisher_exact([[nh, H - nh], [nl, L - nl]])
    out.append({"snv": snv, "gene": ";".join(genes) or "-", "nh": nh, "nl": nl,
                "diff": nh / H - nl / L, "or": odds, "p": p,
                "dir": "high-only" if nl == 0 else "low-only"})

# Benjamini-Hochberg q over the exclusive set
m = len(out)
order = sorted(range(m), key=lambda i: out[i]["p"])
run = 1.0
for rank in range(m, 0, -1):
    i = order[rank - 1]
    run = min(run, out[i]["p"] * m / rank)
    out[i]["q"] = min(run, 1.0)

hi = sorted([d for d in out if d["dir"] == "high-only"], key=lambda d: (-d["diff"], d["p"]))
lo = sorted([d for d in out if d["dir"] == "low-only"], key=lambda d: (d["diff"], d["p"]))

with open("paper_results/table4_exclusive_aaea.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["set", "variant", "gene", "n_high", "n_low", "rate_gap", "OR", "p", "q"])
    for d in hi + lo:
        orr = "inf" if d["or"] == float("inf") else f"{d['or']:.2f}"
        w.writerow([d["dir"], d["snv"], d["gene"], d["nh"], d["nl"],
                    f"{d['diff']:.2f}", orr, f"{d['p']:.4f}", f"{d['q']:.3f}"])

print(f"cohort AA+EA: {H} high / {L} low | exclusive recurrent (>={MIN_CARRIERS}) "
      f"matched SNVs: {len(out)}  (high-only {len(hi)}, low-only {len(lo)})")
print(f"min q over the set: {min(d['q'] for d in out):.3f}  (FDR survivors q<0.05: "
      f"{sum(d['q'] < 0.05 for d in out)})\n")
for title, items in (("HIGH-only (top 15)", hi[:15]), ("LOW-only (top 15)", lo[:15])):
    print(title)
    print(f"  {'variant':<22}{'gene':<18}{'h/l':>6}{'gap':>6}{'p':>9}{'q':>7}")
    for d in items:
        hl = f"{d['nh']}/{d['nl']}"
        print(f"  {d['snv']:<22}{d['gene'][:17]:<18}{hl:>6}"
              f"{d['diff']:>6.2f}{d['p']:>9.4f}{d['q']:>7.3f}")
    print()
