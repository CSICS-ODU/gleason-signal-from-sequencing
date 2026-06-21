"""
cooccurrence.py - Layer 3: gene-gene co-occurrence patterns for BIOLOGICAL
VERIFICATION. Produces co-mutated gene pairs/modules with a severity direction,
anchored on the paper's already-validated candidate genes, so biologists can verify
them experimentally and extend the existing qPCR / cell-line component.

Same concept as the prior classification test (genefilelists/distance.py): each gene
is compared to a HIGH reference and a LOW reference and assigned to the nearer one --
here the "reference" is simply the in-house high vs low patient group (NO averaging),
so a gene's direction is its carrier rate in high (n_high/H) vs low (n_low/L). Genes
that CO-OCCUR (carried by the same patients, above chance) propagate that label:
"this co-occurs with the highs => it is high too."

In-house, PT510=low (13 high / 7 low). Two traps handled:
  * burden/ubiquity -- score co-occurrence ABOVE chance (lift + Fisher on the 2x2
    patient table), never raw overlap.
  * multi-gene-SNV artifact -- genes co-annotated to the SAME SNV co-occur
    mechanically; each pair carries a shared-SNV count and such edges are excluded
    from "real" co-occurrence and from propagation.

Usage:
  python cooccurrence.py                       # PT510_low matched ALL master
  python cooccurrence.py <master_genes.csv> --top-k 25 --min-total 4 --alpha 0.1
  python cooccurrence.py <f> --genes APC,ATM,TP53,CDH1,NPM1,KRAS   # explicit subspace
"""

import argparse
import csv
import glob
from collections import defaultdict
from itertools import combinations

from scipy.stats import fisher_exact

# The paper's experimentally-grounded genes (Section V qPCR / cell-line) plus the
# strongest disproportion drivers -- co-occurrences touching these are the
# verification-ready ones. Override with --known.
KNOWN = ["APC", "ATM", "CDH1", "TP53", "E2F1", "SIRT1", "PRKAB2", "PIK3CD", "MAPK12",
         "NPM1", "KRAS", "IDH1", "MET", "NOTCH1", "SMAD4", "EGFR", "KDR", "STK11",
         "PTEN", "RB1", "FBXW7", "BRAF", "FGFR3", "RET"]


def split_field(val):
    return [x for x in (val or "").split(";") if x]


def keep(pid, no_indian):
    """Indian targeted-panel patients have all-numeric IDs; AA/EA are PT-prefixed."""
    return not (no_indian and pid.isdigit())


def load_master(path, no_indian=False):
    """SNV master -> per-gene high/low carrier sets + per-gene SNV sets + cohort."""
    high_all, low_all = set(), set()
    high_car, low_car = defaultdict(set), defaultdict(set)
    gene_snvs = defaultdict(set)
    saw_genes = False
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if "genes" in r:
                saw_genes = True
            snv = r.get("SNV", "")
            hp = [x for x in split_field(r.get("high_patients")) if keep(x, no_indian)]
            lp = [x for x in split_field(r.get("low_patients")) if keep(x, no_indian)]
            high_all.update(hp)
            low_all.update(lp)
            for g in split_field(r.get("genes")):
                gene_snvs[g].add(snv)
                if hp:
                    high_car[g].update(hp)
                if lp:
                    low_car[g].update(lp)
    if not saw_genes:
        raise SystemExit(f"{path} has no `genes` column")
    genes = {}
    for g in set(high_car) | set(low_car):
        hi, lo = high_car.get(g, set()), low_car.get(g, set())
        genes[g] = {"high": hi, "low": lo, "carriers": hi | lo, "snvs": gene_snvs[g]}
    return {"H": len(high_all), "L": len(low_all), "N": len(high_all | low_all),
            "genes": genes}


def direction(rec, H, L):
    nh, nl = len(rec["high"]), len(rec["low"])
    hr, lr = (nh / H if H else 0.0), (nl / L if L else 0.0)
    d = "HIGH" if hr > lr else ("LOW" if lr > hr else "=")
    return d, hr, lr, hr - lr


def pick_subspace(data, top_k, min_total, explicit, known):
    genes = data["genes"]
    if explicit:
        return [g for g in explicit if g in genes]
    H, L = data["H"], data["L"]
    scored = []
    for g, rec in genes.items():
        if len(rec["carriers"]) <= min_total:
            continue
        _, _, _, diff = direction(rec, H, L)
        scored.append((g, diff))
    highs = [g for g, _ in sorted(scored, key=lambda s: -s[1])[:top_k]]
    lows = [g for g, _ in sorted(scored, key=lambda s: s[1])[:top_k]]
    anchored = [g for g in known if g in genes]          # always include known genes
    return sorted(set(highs) | set(lows) | set(anchored))


def build(data, cands, known, alpha):
    H, L, N = data["H"], data["L"], data["N"]
    genes = data["genes"]
    sign = {"HIGH": 1, "LOW": -1, "=": 0}
    kset = set(known)

    nodes = {}
    for g in cands:
        d, hr, lr, diff = direction(genes[g], H, L)
        nodes[g] = {"gene": g, "dir": d, "hr": hr, "lr": lr, "diff": diff,
                    "n": len(genes[g]["carriers"]), "known": g in kset}

    edges = []
    for a, b in combinations(cands, 2):
        ca, cb = genes[a]["carriers"], genes[b]["carriers"]
        na, nb = len(ca), len(cb)
        both = len(ca & cb)
        if na == 0 or nb == 0:
            continue
        a_only, b_only = na - both, nb - both
        neither = N - both - a_only - b_only
        odds, p = fisher_exact([[both, a_only], [b_only, neither]])
        lift = both * N / (na * nb)
        shared = len(genes[a]["snvs"] & genes[b]["snvs"])
        edges.append({
            "a": a, "b": b, "both": both,
            "jaccard": both / (na + nb - both) if (na + nb - both) else 0.0,
            "lift": lift, "or": odds, "p": p, "shared_snvs": shared,
            "same_dir": nodes[a]["dir"] == nodes[b]["dir"],
            "n_known": (a in kset) + (b in kset),
            "real": lift > 1 and p < alpha and shared == 0,
        })

    # propagation: implied direction = lift-weighted vote of REAL neighbors
    adj = defaultdict(list)
    for e in edges:
        if e["real"]:
            adj[e["a"]].append((e["b"], e["lift"]))
            adj[e["b"]].append((e["a"], e["lift"]))
    for g in cands:
        score = sum(w * sign[nodes[o]["dir"]] for o, w in adj[g])
        deg = len(adj[g])
        implied = "HIGH" if score > 0 else ("LOW" if score < 0 else "-")
        nodes[g].update(degree=deg, implied=implied,
                        agree=(deg > 0 and implied == nodes[g]["dir"]),
                        rescue=(deg >= 2 and abs(nodes[g]["diff"]) < 0.2 and implied != "-"))

    # modules: connected components of the real-edge graph
    seen, modules = set(), []
    for g in cands:
        if g in seen or not adj[g]:
            continue
        stack, comp = [g], []
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            comp.append(x)
            stack.extend(o for o, _ in adj[x] if o not in seen)
        if len(comp) >= 2:
            s = sum(sign[nodes[x]["dir"]] for x in comp)
            modules.append({"genes": sorted(comp),
                            "label": "HIGH" if s > 0 else ("LOW" if s < 0 else "mixed"),
                            "n_known": sum(nodes[x]["known"] for x in comp)})
    return nodes, edges, modules


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("master", nargs="?")
    ap.add_argument("--top-k", type=int, default=25)
    ap.add_argument("--min-total", type=int, default=4)
    ap.add_argument("--alpha", type=float, default=0.1)
    ap.add_argument("--genes", help="explicit candidate subspace, comma-separated")
    ap.add_argument("--known", help="override the verification anchor-gene list")
    ap.add_argument("--no-indian", action="store_true",
                    help="exclude Indian (numeric-ID) patients -> AA+EA broad-panel stratum")
    ap.add_argument("--drop-ubiquitous", action="store_true",
                    help="drop genes carried by ALL stratum patients (no co-occurrence info)")
    ap.add_argument("--out", default="cooccurrence_results")
    a = ap.parse_args()

    path = a.master or next(iter(sorted(glob.glob(
        "snv_rebuild/PT510_low/matched/ALL/**/ALL_master_matched*_genes.csv",
        recursive=True))), None)
    if not path:
        raise SystemExit("no PT510_low master found; pass one explicitly")

    data = load_master(path, a.no_indian)
    if a.drop_ubiquitous:
        N = data["N"]
        data["genes"] = {g: rec for g, rec in data["genes"].items()
                         if 0 < len(rec["carriers"]) < N}
    known = [g.strip() for g in a.known.split(",")] if a.known else KNOWN
    explicit = [g.strip() for g in a.genes.split(",")] if a.genes else None
    cands = pick_subspace(data, a.top_k, a.min_total, explicit, known)
    if len(cands) < 2:
        raise SystemExit("need >=2 candidate genes")

    nodes, edges, modules = build(data, cands, known, a.alpha)
    real = [e for e in edges if e["real"]]
    same = sum(1 for e in real if e["same_dir"])
    nk = sum(nodes[g]["known"] for g in cands)

    print("=" * 70)
    print(f"{path}\n  cohort: {data['H']} high / {data['L']} low (N={data['N']}) | "
          f"subspace: {len(cands)} genes ({nk} known) | alpha={a.alpha}")
    print(f"  real co-occurrence edges (lift>1, p<{a.alpha}, no shared SNV): {len(real)}"
          + (f" | same-direction: {same}/{len(real)} ({100*same/len(real):.0f}%)" if real else ""))

    print("\n  VERIFICATION-READY PAIRS (>=1 known gene, real, same-direction)")
    print(f"    {'gene A':<13}{'gene B':<13}{'dir':>5}{'both':>6}{'lift':>7}{'p':>9}{'known':>7}")
    vp = [e for e in real if e["same_dir"] and e["n_known"] >= 1]
    for e in sorted(vp, key=lambda e: (-e["n_known"], -e["lift"]))[:30]:
        dd = nodes[e["a"]]["dir"][0] + nodes[e["b"]]["dir"][0]
        print(f"    {e['a'][:12]:<13}{e['b'][:12]:<13}{dd:>5}{e['both']:>6}"
              f"{e['lift']:>7.2f}{e['p']:>9.2g}{e['n_known']:>7}")

    print("\n  CO-OCCURRENCE MODULES (connected co-mutated gene sets)")
    for m in sorted(modules, key=lambda m: (-m["n_known"], -len(m["genes"]))):
        print(f"    [{m['label']:<5} | {m['n_known']} known] {', '.join(m['genes'])}")

    rescued = [g for g in cands if nodes[g]["rescue"]]
    if rescued:
        print("\n  RESCUED (weak own signal, labeled by co-occurrence neighbors)")
        for g in rescued:
            n = nodes[g]
            print(f"    {g:<13} own={n['dir']:<5} implied={n['implied']:<5} "
                  f"deg={n['degree']} {'(known)' if n['known'] else ''}")

    import os
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "cooc_edges.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["gene_a", "gene_b", "both", "jaccard", "lift", "odds_ratio", "p",
                    "shared_snvs", "same_dir", "n_known", "real"])
        for e in sorted(edges, key=lambda e: -e["lift"]):
            w.writerow([e["a"], e["b"], e["both"], f"{e['jaccard']:.4f}", f"{e['lift']:.4f}",
                        e["or"], e["p"], e["shared_snvs"], e["same_dir"], e["n_known"], e["real"]])
    with open(os.path.join(a.out, "cooc_nodes.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["gene", "direction", "high_rate", "low_rate", "diff", "n_carriers",
                    "known", "degree", "implied", "agree", "rescue"])
        for g in cands:
            n = nodes[g]
            w.writerow([g, n["dir"], f"{n['hr']:.3f}", f"{n['lr']:.3f}", f"{n['diff']:.3f}",
                        n["n"], n["known"], n["degree"], n["implied"], n["agree"], n["rescue"]])
    print(f"\n  wrote {a.out}/cooc_edges.csv and cooc_nodes.csv")


if __name__ == "__main__":
    main()
