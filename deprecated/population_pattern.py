"""
population_pattern.py - per-population (AA vs EA) high/low pattern for a focused gene
set, so population-specific severity signal can be read directly.

Why a focused set: at 2 high / 3 low per population the perfect-split p-floor is 0.1
and thousands of genes hit it by chance, so genome-wide AA-vs-EA differences are
indistinguishable from noise. Restricting to biologically-motivated genes (known
cancer drivers) shrinks the hypothesis space so that cross-population CONSISTENCY
(same direction in both = robust) and DIVERGENCE (opposite = population-specific) are
interpretable rather than chance.

Usage:
  python population_pattern.py                      # known cancer genes, PT510_low matched
  python population_pattern.py --genes APC,TP53,ATM,...
"""

import argparse
import csv
import glob
from collections import defaultdict

KNOWN = ["APC", "ATM", "CDH1", "TP53", "E2F1", "SIRT1", "PRKAB2", "PIK3CD", "MAPK12",
         "NPM1", "KRAS", "IDH1", "MET", "NOTCH1", "SMAD4", "EGFR", "KDR", "STK11",
         "PTEN", "RB1", "FBXW7", "BRAF", "FGFR3", "RET", "PIK3CA", "AKT1", "CTNNB1",
         "GNAS", "ALK", "ERBB2", "ERBB4", "FLT3", "KIT", "CDKN2A", "HRAS", "JAK3",
         "PDGFRA", "SRC", "VHL", "MLH1", "MSH2", "NRAS", "SMO", "SMARCB1"]


def load(paths, somatic=None):
    """Union carriers across one or more layer masters. If `somatic` is given
    (patient -> set of somatic-tagged SNV keys), a patient counts as a gene-carrier
    only via SNVs DRAGEN tagged Somatic -- so the rates reflect somatic, not germline,
    variants. Cohort sizes H/L stay the FULL cohort (denominator), unfiltered."""
    hc, lc = defaultdict(set), defaultdict(set)
    H, L = set(), set()
    for path in paths:
        with open(path, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                snv = r["SNV"]
                hp = [x for x in (r["high_patients"] or "").split(";") if x]
                lp = [x for x in (r["low_patients"] or "").split(";") if x]
                H.update(hp)
                L.update(lp)
                if somatic is not None:
                    hp = [x for x in hp if snv in somatic.get(x, ())]
                    lp = [x for x in lp if snv in somatic.get(x, ())]
                for g in (r["genes"] or "").split(";"):
                    if g:
                        hc[g].update(hp)
                        lc[g].update(lp)
    return len(H), len(L), hc, lc


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", help="comma-separated gene set (default: known drivers)")
    ap.add_argument("--layers", default="matched",
                    help="comma list of layers to union: e.g. matched,benign,unmatched")
    ap.add_argument("--somatic", help="somatic_keys.csv: restrict carriers to "
                    "DRAGEN Somatic-tagged variants (germline stripped)")
    a = ap.parse_args()

    layers = [x.strip() for x in a.layers.split(",")]

    somatic = None
    if a.somatic:
        somatic = {}
        with open(a.somatic, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                somatic.setdefault(r["patient"], set()).add(r["snv"])

    def paths_for(anc):
        ps = []
        for layer in layers:
            ps += glob.glob(
                f"snv_rebuild/PT510_low/{layer}/{anc}/**/{anc}_master_{layer}*_genes.csv",
                recursive=True)
        return ps

    aa_paths, ea_paths = paths_for("AA"), paths_for("EA")
    genes = [g.strip() for g in a.genes.split(",")] if a.genes else KNOWN

    Ha, La, hca, lca = load(aa_paths, somatic)
    He, Le, hce, lce = load(ea_paths, somatic)
    tag = "  [SOMATIC-only]" if somatic else ""
    print(f"layers: {','.join(layers)}{tag}  ({len(aa_paths)} AA files, {len(ea_paths)} EA files)")
    print(f"AA cohort {Ha}h/{La}l   EA cohort {He}h/{Le}l   |  {len(genes)} candidate genes")
    print(f"\n  {'gene':<10}{'AA hi/lo':>9}{'AA dir':>8}{'EA hi/lo':>9}{'EA dir':>8}   verdict")
    rows = []
    for g in genes:
        nha, nla = len(hca.get(g, ())), len(lca.get(g, ()))
        nhe, nle = len(hce.get(g, ())), len(lce.get(g, ()))
        if nha + nla + nhe + nle == 0:
            continue                                   # absent from both broad panels
        da = direction(nha, nla, Ha, La)
        de = direction(nhe, nle, He, Le)
        rows.append((g, nha, nla, da, nhe, nle, de, verdict(da, de)))

    # order: consistent first, then divergent, then population-specific, then flat
    rank = {"CONSISTENT-HIGH": 0, "CONSISTENT-LOW": 0, "DIVERGENT": 1}
    rows.sort(key=lambda r: (rank.get(r[7], 2 if "only" in r[7] else 3), r[0]))
    for g, nha, nla, da, nhe, nle, de, v in rows:
        print(f"  {g:<10}{f'{nha}/{nla}':>9}{da:>8}{f'{nhe}/{nle}':>9}{de:>8}   {v}")
    if not rows:
        print("  (none of the candidate genes appear in the AA/EA broad panels)")


if __name__ == "__main__":
    main()
