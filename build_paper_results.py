"""
build_paper_results.py - generate the complete Results table set for the manuscript,
at the locked PT510=low setting, in one reproducible pass.

This replaces the discarded manuscript_results/*.tex draft. It reads the rebuilt SNV
masters directly (snv_rebuild/PT510_low/<layer>/<group>/.../<group>_master_<layer>*.csv)
and emits one CSV per manuscript result block into paper_results/:

  tableIVB_bucket_counts.csv   IV.B  SNV counts per annotation bucket x group x severity
  tableIVC_snv_sets.csv        IV.C  common / high-enriched / low-enriched / unknown SNVs
  tableIVD_gene_mapping.csv    IV.D  high/low-enriched genes + per-group recurrent genes
  tableIVE_cross_pop.csv       IV.E  AA-vs-EA cross-population consistency / divergence
  tableIV_candidates.csv       Tbl IV  computational status of the 7 qPCR genes

Design choices (consistent with the locked analysis):
  * PT510 = low everywhere -> ALL cohort is 13 high / 7 low.
  * v_j = presence/absence at the set/overlap level (a patient carries a gene/SNV if it
    appears in their variant profile); the IV.A distance probe separately uses normalized
    variant counts -- this split is stated in the manuscript.
  * Severity skew = per-grade RATE comparison (n_high/H vs n_low/L) with Fisher's exact
    test + Benjamini-Hochberg q across all tests; ranked by effect size (rate gap), since
    at n=20 nothing survives FDR and magnitude is more informative than nominal p.
  * IN = clinical-confirmation arm: a targeted hotspot panel run on suspected cancer
    patients to verify/rule out known cancer genes. Its near-universal recurrent calls are
    reported as confirmation of known genes, NOT as a within-IN high/low contrast (1 low
    IN patient makes that vacuous). Cross-population work is the AA+EA broad stratum.

Run:  python build_paper_results.py
Manuscript mapping: tableIVB -> Table 1b and Fig. 2; tableIVE -> Fig. 3; tableIV_candidates
-> Table 4a, DNA-variant column. (The IV.x labels are the section numbers of an earlier
draft.)
"""

import argparse
import csv
import glob
import os
from collections import defaultdict

from scipy.stats import fisher_exact

BASE = "snv_rebuild/PT510_low"
GROUPS = ("ALL", "IN", "EA", "AA")
LAYERS = ("matched", "benign", "unmatched")
OUTDIR = "paper_results"

# Folder name -> manuscript annotation bucket (the terminology the paper uses).
BUCKET_LABEL = {"matched": "pathogenic/known", "benign": "benign/likely-benign",
                "unmatched": "unknown/unannotated"}

# qPCR / cell-line genes that Table IV must report a computational status for.
QPCR_GENES = ["APC", "ATM", "CDH1", "TP53", "E2F1", "SIRT1", "PRKAB2"]

# Known cancer drivers used for the cross-population (IV.E) read; restricting to a
# biologically-motivated set is what makes AA-vs-EA consistency/divergence interpretable
# rather than chance at 2-3 patients per cell.
KNOWN = ["APC", "ATM", "CDH1", "TP53", "E2F1", "SIRT1", "PRKAB2", "PIK3CD", "MAPK12",
         "NPM1", "KRAS", "IDH1", "MET", "NOTCH1", "SMAD4", "EGFR", "KDR", "STK11",
         "PTEN", "RB1", "FBXW7", "BRAF", "FGFR3", "RET", "PIK3CA", "AKT1", "CTNNB1",
         "GNAS", "ALK", "ERBB2", "ERBB4", "FLT3", "KIT", "CDKN2A", "HRAS", "JAK3",
         "PDGFRA", "SRC", "VHL", "MLH1", "MSH2", "NRAS", "SMO", "SMARCB1"]


# --------------------------------------------------------------------------- io

def split_field(val):
    return [x for x in (val or "").split(";") if x]


def keep(pid, no_indian):
    """Indian targeted-panel patients have all-numeric IDs; AA/EA are PT-prefixed."""
    return not (no_indian and pid.isdigit())


def master_path(group, layer):
    """Consolidated master for a group+layer; prefer the gene-annotated one."""
    found = glob.glob(f"{BASE}/{layer}/{group}/**/{group}_master_{layer}*.csv",
                      recursive=True)
    genes = [f for f in found if f.endswith("_genes.csv")]
    return (genes or found or [None])[0]


def load(path, no_indian=False):
    """Read a master -> (rows, H, L). Each row keeps recomputed high/low carrier lists
    (so --no-indian is exact) plus n_total, category, genes (if present)."""
    rows, H, L = [], set(), set()
    has_genes = False
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if "genes" in r:
                has_genes = True
            hp = [x for x in split_field(r.get("high_patients")) if keep(x, no_indian)]
            lp = [x for x in split_field(r.get("low_patients")) if keep(x, no_indian)]
            H.update(hp)
            L.update(lp)
            rows.append({"snv": r["SNV"], "cat": r.get("category", ""),
                         "genes": split_field(r.get("genes")), "hp": hp, "lp": lp})
    return rows, len(H), len(L), has_genes


# ------------------------------------------------------------------ statistics

def bh_qvalues(items):
    """Benjamini-Hochberg q across all Fisher tests (both directions = one family)."""
    m = len(items)
    if not m:
        return
    order = sorted(range(m), key=lambda i: items[i]["p"])
    running = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, items[i]["p"] * m / rank)
        items[i]["q"] = min(running, 1.0)


def fisher_rows(units, H, L, min_total):
    """units: iterable of (name, gene_or_blank, high_set, low_set). Per-unit rate test."""
    out = []
    for name, gene, hset, lset in units:
        nh, nl = len(hset), len(lset)
        if nh + nl <= min_total:
            continue
        odds, p = fisher_exact([[nh, H - nh], [nl, L - nl]])
        hr, lr = (nh / H if H else 0.0), (nl / L if L else 0.0)
        out.append({"name": name, "gene": gene, "nh": nh, "nl": nl,
                    "hr": hr, "lr": lr, "diff": hr - lr, "or": odds, "p": p,
                    "dir": "HIGH" if hr > lr else ("LOW" if lr > hr else "=")})
    bh_qvalues(out)
    return out


def gene_carriers(rows):
    """Per-gene presence union: {gene: (high_set, low_set)} (>=1 SNV => carrier)."""
    hc, lc = defaultdict(set), defaultdict(set)
    for r in rows:
        for g in r["genes"]:
            if r["hp"]:
                hc[g].update(r["hp"])
            if r["lp"]:
                lc[g].update(r["lp"])
    return {g: (hc.get(g, set()), lc.get(g, set())) for g in set(hc) | set(lc)}


def write_csv(name, header, rows):
    os.makedirs(OUTDIR, exist_ok=True)
    path = os.path.join(OUTDIR, name)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print(f"  wrote {path}  ({len(rows)} rows)")
    return path


# ----------------------------------------------------------------- IV.B counts

def table_IVB():
    """SNV counts per annotation bucket x group x severity category."""
    out = []
    print("\n[IV.B] SNV bucket counts (PT510_low)")
    print(f"  {'bucket':<20}{'group':<5}{'H':>3}{'L':>4}{'nSNV':>9}"
          f"{'hi_only':>9}{'lo_only':>9}{'common':>9}")
    for layer in LAYERS:
        for group in GROUPS:
            path = master_path(group, layer)
            if not path:
                continue
            rows, H, L, _ = load(path)
            cat = defaultdict(int)
            for r in rows:
                cat[r["cat"]] += 1
            n = len(rows)
            ho, lo, co = cat["high_only"], cat["low_only"], cat["high_low_common"]
            out.append([BUCKET_LABEL[layer], group, H, L, n, ho, lo, co])
            print(f"  {BUCKET_LABEL[layer]:<20}{group:<5}{H:>3}{L:>4}{n:>9}"
                  f"{ho:>9}{lo:>9}{co:>9}")
    write_csv("tableIVB_bucket_counts.csv",
              ["bucket", "group", "n_high", "n_low", "n_snv",
               "high_only", "low_only", "high_low_common"], out)


# --------------------------------------------------------------- IV.C SNV sets

def table_IVC(top):
    """Common / high-enriched / low-enriched (matched ALL) + unknown candidates."""
    rows, H, L, _ = load(master_path("ALL", "matched"))
    units = [(r["snv"], ";".join(r["genes"]), set(r["hp"]), set(r["lp"]))
             for r in rows if len(r["hp"]) + len(r["lp"]) > 10]
    res = fisher_rows(units, H, L, min_total=-1)  # already prefiltered to recurrent
    high = sorted([d for d in res if d["dir"] == "HIGH"], key=lambda d: (-d["diff"], d["p"]))
    low = sorted([d for d in res if d["dir"] == "LOW"], key=lambda d: (d["diff"], d["p"]))
    n_common = sum(1 for r in rows if r["cat"] == "high_low_common")

    # unknown candidates: recurrent novel variants (no annotation) -- report by burden.
    urows, _, _, _ = load(master_path("ALL", "unmatched"))
    unk = sorted(urows, key=lambda r: -(len(r["hp"]) + len(r["lp"])))

    out = []
    for label, items in (("high_enriched", high[:top]), ("low_enriched", low[:top])):
        for d in items:
            out.append([label, d["name"], d["gene"], d["nh"], d["nl"],
                        f"{d['hr']:.2f}", f"{d['lr']:.2f}", f"{d['diff']:.2f}",
                        ("inf" if d["or"] == float("inf") else f"{d['or']:.2f}"),
                        f"{d['p']:.2g}", f"{d['q']:.2g}"])
    for r in unk[:top]:
        out.append(["unknown_candidate", r["snv"], "", len(r["hp"]), len(r["lp"]),
                    "", "", "", "", "", ""])
    write_csv("tableIVC_snv_sets.csv",
              ["set", "snv", "genes", "n_high", "n_low", "high_rate", "low_rate",
               "rate_gap", "odds_ratio", "p", "q"], out)
    sig_q = sum(1 for d in res if d["q"] < 0.05)
    print(f"\n[IV.C] matched ALL: cohort {H}h/{L}l | recurrent SNVs tested {len(res)} | "
          f"common {n_common} | high-enriched {len(high)} | low-enriched {len(low)} | "
          f"q<0.05: {sig_q}")
    print(f"  unknown(unmatched) ALL recurrent variants: {len(urows)} "
          f"(IN contributes ~none -- targeted panel)")


# ----------------------------------------------------------- IV.D gene mapping

def table_IVD(top):
    """High/low-enriched genes (matched ALL) + per-group recurrent genes."""
    rows, H, L, has_genes = load(master_path("ALL", "matched"))
    if not has_genes:
        raise SystemExit("matched ALL master has no genes column")
    carr = gene_carriers(rows)
    units = [(g, g, hs, ls) for g, (hs, ls) in carr.items() if len(hs) + len(ls) > 6]
    res = fisher_rows(units, H, L, min_total=-1)
    high = sorted([d for d in res if d["dir"] == "HIGH"], key=lambda d: (-d["diff"], d["p"]))
    low = sorted([d for d in res if d["dir"] == "LOW"], key=lambda d: (d["diff"], d["p"]))

    # per-group recurrent genes (>=2 carriers) -> population-enriched view.
    group_genes = {}
    for group in ("IN", "EA", "AA"):
        grows, _, _, hg = load(master_path(group, "matched"))
        gc = gene_carriers(grows) if hg else {}
        rec = {g: len(hs | ls) for g, (hs, ls) in gc.items() if len(hs | ls) >= 2}
        group_genes[group] = rec

    out = []
    for label, items in (("high_enriched", high[:top]), ("low_enriched", low[:top])):
        for d in items:
            out.append([label, d["name"], d["nh"], d["nl"],
                        f"{d['diff']:.2f}",
                        ("inf" if d["or"] == float("inf") else f"{d['or']:.2f}"),
                        f"{d['p']:.2g}", f"{d['q']:.2g}"])
    for group in ("IN", "EA", "AA"):
        for g, n in sorted(group_genes[group].items(), key=lambda kv: -kv[1])[:top]:
            out.append([f"recurrent_{group}", g, n, "", "", "", "", ""])
    write_csv("tableIVD_gene_mapping.csv",
              ["set", "gene", "n_high", "n_low", "rate_gap", "odds_ratio", "p", "q"], out)
    print(f"\n[IV.D] matched ALL genes tested {len(res)} | "
          f"high-enriched {len(high)} | low-enriched {len(low)} | "
          f"recurrent genes: IN {len(group_genes['IN'])}, "
          f"EA {len(group_genes['EA'])}, AA {len(group_genes['AA'])}")


# -------------------------------------------------------- IV.E cross-population

def direction(nh, nl, H, L):
    hr, lr = (nh / H if H else 0), (nl / L if L else 0)
    return "HIGH" if hr > lr else ("LOW" if lr > hr else "=")


def verdict(da, de):
    if da == de and da in ("HIGH", "LOW"):
        return f"CONSISTENT-{da}"
    if {da, de} == {"HIGH", "LOW"}:
        return "DIVERGENT"
    if da in ("HIGH", "LOW") and de == "=":
        return f"AA-only {da}"
    if de in ("HIGH", "LOW") and da == "=":
        return f"EA-only {de}"
    return "flat"


def table_IVE():
    """AA-vs-EA consistency/divergence over known genes (matched layer)."""
    aa, Ha, La, _ = load(master_path("AA", "matched"))
    ea, He, Le, _ = load(master_path("EA", "matched"))
    ca, ce = gene_carriers(aa), gene_carriers(ea)
    # cross-population recurrent presence: genes recurrent (>=2) in IN, AA, AND EA.
    inr, _, _, _ = load(master_path("IN", "matched"))
    cin = gene_carriers(inr)

    def rec(carr, g):
        hs, ls = carr.get(g, (set(), set()))
        return len(hs | ls)

    rows = []
    for g in KNOWN:
        nha, nla = len(ca.get(g, (set(), set()))[0]), len(ca.get(g, (set(), set()))[1])
        nhe, nle = len(ce.get(g, (set(), set()))[0]), len(ce.get(g, (set(), set()))[1])
        if nha + nla + nhe + nle == 0:
            continue
        da, de = direction(nha, nla, Ha, La), direction(nhe, nle, He, Le)
        xpop = "yes" if (rec(cin, g) >= 2 and rec(ca, g) >= 2 and rec(ce, g) >= 2) else ""
        rows.append([g, f"{nha}/{nla}", da, f"{nhe}/{nle}", de, verdict(da, de), xpop])
    order = {"CONSISTENT-HIGH": 0, "CONSISTENT-LOW": 0, "DIVERGENT": 1}
    rows.sort(key=lambda r: (order.get(r[5], 2 if "only" in r[5] else 3), r[0]))
    write_csv("tableIVE_cross_pop.csv",
              ["gene", "AA_hi/lo", "AA_dir", "EA_hi/lo", "EA_dir", "verdict",
               "recurrent_all3_pops"], rows)
    print(f"\n[IV.E] AA cohort {Ha}h/{La}l  EA cohort {He}h/{Le}l  "
          f"(broad stratum) | {len(rows)} known genes present")
    for r in rows:
        print(f"  {r[0]:<9} AA {r[1]:>4} {r[2]:<5} EA {r[3]:>4} {r[4]:<5} {r[5]}"
              f"{'  [recurrent in all 3]' if r[6] else ''}")


# ----------------------------------------------------------- Table IV lookup

def table_IV_candidates():
    """Computational status of the 7 qPCR genes for Table IV's empty column."""
    allrows, H, L, _ = load(master_path("ALL", "matched"))
    callr = gene_carriers(allrows)
    per_group = {}
    for group in ("IN", "EA", "AA"):
        grows, _, _, hg = load(master_path(group, "matched"))
        per_group[group] = gene_carriers(grows) if hg else {}
    aa, Ha, La, _ = load(master_path("AA", "matched"))
    ea, He, Le, _ = load(master_path("EA", "matched"))
    ca, ce = gene_carriers(aa), gene_carriers(ea)

    out = []
    for g in QPCR_GENES:
        hs, ls = callr.get(g, (set(), set()))
        nh, nl = len(hs), len(ls)
        present = nh + nl > 0
        groups = [grp for grp in ("IN", "EA", "AA")
                  if sum(len(s) for s in per_group[grp].get(g, (set(), set()))) > 0]
        if present:
            d = direction(nh, nl, H, L)
            nha = sum(len(s) for s in ca.get(g, (set(), set())))
            nhe = sum(len(s) for s in ce.get(g, (set(), set())))
            da = direction(len(ca.get(g, (set(), set()))[0]),
                           len(ca.get(g, (set(), set()))[1]), Ha, La) if nha else "absent"
            de = direction(len(ce.get(g, (set(), set()))[0]),
                           len(ce.get(g, (set(), set()))[1]), He, Le) if nhe else "absent"
            context = (f"recurrent (ALL {nh}h/{nl}l, {d}); groups: {','.join(groups)}; "
                       f"AA {da}/EA {de}")
        else:
            context = "not detected in cohort variant data (experimental/pathway context only)"
        out.append([g, "yes" if present else "no", nh, nl, ";".join(groups), context])
    write_csv("tableIV_candidates.csv",
              ["gene", "in_variant_data", "n_high", "n_low", "groups",
               "computational_signal_context"], out)
    print("\n[Table IV] qPCR-gene computational status:")
    for r in out:
        print(f"  {r[0]:<8} {r[5]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=25, help="rows per ranked set")
    a = ap.parse_args()
    print(f"Building manuscript result set from {BASE}  ->  {OUTDIR}/")
    table_IVB()
    table_IVC(a.top)
    table_IVD(a.top)
    table_IVE()
    table_IV_candidates()
    print("\nDone. Five CSVs in paper_results/ ready to populate IV.B-E and Table IV.")


if __name__ == "__main__":
    main()
