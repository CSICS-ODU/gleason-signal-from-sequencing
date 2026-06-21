"""
query_genes.py - query the rebuilt gene x patient matrices.

A "SKU" is one ancestry x category slice, within a layer (default matched)
and PT510 setting (default high).

Short aliases (case-insensitive)
--------------------------------
  ancestry:  A = ALL    AA = AA    EA = EA    IN = IN
  category:  h = high_only   l = low_only   c = high_low_common   a = all (pool of the 3)
  layer:     b = benign      m = matched    u = unmatched
  pt510:     H = high        L = low        (CAPITAL, to avoid clashing with category h/l)

A SKU is written ANCESTRY:CATEGORY, e.g.  AA:h   EA:l   A:a   IN:c

Commands
--------
  list                              every SKU and its gene count
  top --sku AA:h                    top genes in one slice
  top --all                         top genes cohort-wide  (= --sku A:a)
  top --all --category h            cohort-wide, one category
      [--metric total|patients|high_pat|low_pat] [--n 20]
  absent --sku AA:h                 genes in the ancestry's 'all' but NOT in this slice
  gene -g TP53 [-a AA]              per-patient breakdown for one gene (default ancestry ALL)

Common options:  --layer b|m|u   --pt510 h|l

Examples
--------
  python query_genes.py top --sku A:a --n 10
  python query_genes.py top --sku EA:h --metric total
  python query_genes.py absent --sku AA:h
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path("snv_rebuild")
ANCESTRIES = ["AA", "EA", "IN", "ALL"]
CATEGORIES = ["high_only", "low_only", "high_low_common"]
METRICS = ["total", "patients", "high_pat", "low_pat"]

ANC_ALIAS = {"a": "ALL", "all": "ALL", "aa": "AA",
             "e": "EA", "ea": "EA", "i": "IN", "in": "IN"}
CAT_ALIAS = {"a": "all", "all": "all",
             "h": "high_only", "ho": "high_only", "high_only": "high_only",
             "l": "low_only", "lo": "low_only", "low_only": "low_only",
             "c": "high_low_common", "hlc": "high_low_common",
             "common": "high_low_common", "high_low_common": "high_low_common"}
LAYER_ALIAS = {"b": "benign", "benign": "benign", "m": "matched",
               "matched": "matched", "u": "unmatched", "unmatched": "unmatched"}
# PT510 is CASE-SENSITIVE: capital H/L only, so it never clashes with category h/l.
PT510_ALIAS = {"H": "high", "L": "low", "high": "high", "low": "low",
               "HIGH": "high", "LOW": "low"}
METRIC_ALIAS = {"t": "total", "total": "total", "p": "patients",
                "patients": "patients", "high_pat": "high_pat",
                "low_pat": "low_pat", "hp": "high_pat", "lp": "low_pat"}


def resolve(token, alias, what, case_sensitive=False):
    key = str(token) if case_sensitive else str(token).lower()
    out = alias.get(key)
    if out is None:
        raise SystemExit(f"unknown {what}: {token!r}  (valid: "
                         f"{sorted(set(alias.values()))})")
    return out


def gene_file(ancestry, category, layer, pt510):
    return (ROOT / f"PT510_{pt510}" / layer / ancestry / "gene results"
            / f"{ancestry}_{category}_genes_{layer}.csv")


def load_sku(ancestry, category, layer, pt510):
    """gene -> {'total':int, 'high':set(patients), 'low':set(patients)}.

    Breadth is patient SETS so pooling across categories unions patients
    instead of double-counting them.
    """
    path = gene_file(ancestry, category, layer, pt510)
    out = {}
    if not path.exists():
        return out
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        gi, ti = header.index("Gene"), header.index("Total")
        patient_cols = [(i, header[i].rsplit(":", 1)) for i in range(gi + 1, ti)]
        for row in reader:
            high, low, total = set(), set(), 0
            for i, (pid, grade) in patient_cols:
                v = int(row[i])
                total += v
                if v > 0:
                    (high if grade == "high" else low).add(pid)
            out[row[gi]] = {"total": total, "high": high, "low": low}
    return out


def merge(dicts):
    """Combine gene records across SKUs: SUM magnitude, UNION patient sets."""
    pooled = defaultdict(lambda: {"total": 0, "high": set(), "low": set()})
    for d in dicts:
        for g, rec in d.items():
            pooled[g]["total"] += rec["total"]
            pooled[g]["high"] |= rec["high"]
            pooled[g]["low"] |= rec["low"]
    return pooled


def load_slice(ancestry, category, layer, pt510):
    """One slice; category 'all' pools the three real categories."""
    if category == "all":
        return merge(load_sku(ancestry, c, layer, pt510) for c in CATEGORIES)
    return load_sku(ancestry, category, layer, pt510)


def patient_list(ancestry, layer, pt510):
    """Ordered [(patient, grade)] for an ancestry (high columns first, then low)."""
    for cat in ("high_low_common", "high_only", "low_only"):
        path = gene_file(ancestry, cat, layer, pt510)
        if path.exists():
            with open(path, newline="", encoding="utf-8") as f:
                header = next(csv.reader(f))
            gi, ti = header.index("Gene"), header.index("Total")
            return [tuple(header[i].rsplit(":", 1)) for i in range(gi + 1, ti)]
    return []


def gene_row(ancestry, category, layer, pt510, gene):
    """{patient: count} for one gene in one slice ({} if the gene is absent)."""
    path = gene_file(ancestry, category, layer, pt510)
    if not path.exists():
        return {}
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        gi, ti = header.index("Gene"), header.index("Total")
        pcols = [(i, header[i].rsplit(":", 1)[0]) for i in range(gi + 1, ti)]
        for row in reader:
            if row[gi] == gene:
                return {pid: int(row[i]) for i, pid in pcols}
    return {}


def metric_value(rec, metric):
    if metric == "total":
        return rec["total"]
    if metric == "patients":
        return len(rec["high"]) + len(rec["low"])
    if metric == "high_pat":
        return len(rec["high"])
    if metric == "low_pat":
        return len(rec["low"])
    raise ValueError(metric)


def parse_sku(s):
    try:
        anc, cat = s.split(":", 1)
    except ValueError:
        raise SystemExit(f"--sku must look like AA:h, got {s!r}")
    return resolve(anc, ANC_ALIAS, "ancestry"), resolve(cat, CAT_ALIAS, "category")


def print_ranked(title, genes, metric, n):
    ranked = sorted(genes.items(), key=lambda kv: metric_value(kv[1], metric),
                    reverse=True)
    print(title)
    if not ranked:
        print("  (no genes - empty slice)")
        return
    print(f"  {'#':>3}  {'gene':<16}{metric:>10}   (total/high_pat/low_pat)")
    for i, (g, rec) in enumerate(ranked[:n], 1):
        print(f"  {i:>3}  {g:<16}{metric_value(rec, metric):>10}   "
              f"({rec['total']}/{len(rec['high'])}/{len(rec['low'])})")


def cmd_list(args):
    print(f"SKUs for layer={args.layer} PT510={args.pt510}  (genes per slice)")
    print(f"  {'ancestry':<6} " + " ".join(f"{c:>16}" for c in CATEGORIES)
          + f"{'all':>16}")
    for anc in ANCESTRIES:
        cells = [str(len(load_sku(anc, c, args.layer, args.pt510))) for c in CATEGORIES]
        cells.append(str(len(load_slice(anc, "all", args.layer, args.pt510))))
        print(f"  {anc:<6} " + " ".join(f"{c:>16}" for c in cells))


def cmd_top(args):
    if args.sku:
        anc, cat = parse_sku(args.sku)
    else:  # --all
        anc, cat = "ALL", (resolve(args.category, CAT_ALIAS, "category")
                           if args.category else "all")
    genes = load_slice(anc, cat, args.layer, args.pt510)
    metric = resolve(args.metric, METRIC_ALIAS, "metric")
    print_ranked(
        f"Top {args.n} genes | {anc}:{cat} | layer={args.layer} "
        f"PT510={args.pt510} | by {metric}", genes, metric, args.n)


def cmd_absent(args):
    anc, cat = parse_sku(args.sku)
    present = set(load_slice(anc, cat, args.layer, args.pt510))
    allgenes = set(load_slice(anc, "all", args.layer, args.pt510))  # the relevant 'all'
    absent = sorted(allgenes - present)
    print(f"Genes ABSENT from {anc}:{cat}  (vs {anc}:all) | "
          f"layer={args.layer} PT510={args.pt510}")
    print(f"  {anc}:all genes: {len(allgenes)} | present in {cat}: {len(present)} | "
          f"absent: {len(absent)}")
    for g in absent:
        print(f"  {g}")


def cmd_gene(args):
    anc = resolve(args.ancestry, ANC_ALIAS, "ancestry") if args.ancestry else "ALL"
    g = args.gene
    pats = patient_list(anc, args.layer, args.pt510)
    if not pats:
        raise SystemExit(f"no gene matrices for ancestry {anc}")
    cat_rows = {c: gene_row(anc, c, args.layer, args.pt510, g) for c in CATEGORIES}
    if not any(cat_rows.values()):
        print(f"Gene {g!r} not found in {anc} (layer={args.layer}, PT510={args.pt510}).")
        return

    print(f"Gene {g} | ancestry={anc} layer={args.layer} PT510={args.pt510}")
    print(f"  {'patient':<14}{'grade':<6}{'high_only':>10}{'low_only':>10}"
          f"{'common':>9}{'total':>8}")
    tot = [0, 0, 0]
    n_high = n_low = 0
    for pid, grade in pats:
        vals = [cat_rows[c].get(pid, 0) for c in CATEGORIES]
        s = sum(vals)
        for i in range(3):
            tot[i] += vals[i]
        if s > 0:
            n_high += grade == "high"
            n_low += grade == "low"
        flag = "" if s else "   (absent)"
        print(f"  {pid:<14}{grade:<6}{vals[0]:>10}{vals[1]:>10}{vals[2]:>9}{s:>8}{flag}")
    print(f"  {'TOTAL':<14}{'':<6}{tot[0]:>10}{tot[1]:>10}{tot[2]:>9}{sum(tot):>8}")
    print(f"  carried by {n_high} high + {n_low} low patient(s)")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--layer", default="matched")
    p.add_argument("--pt510", default="high")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    t = sub.add_parser("top")
    t.add_argument("--sku")
    t.add_argument("--all", action="store_true")
    t.add_argument("--category")
    t.add_argument("--metric", default="total")
    t.add_argument("--n", type=int, default=20)
    a = sub.add_parser("absent")
    a.add_argument("--sku", required=True)
    gn = sub.add_parser("gene")
    gn.add_argument("-g", "--gene", required=True)
    gn.add_argument("-a", "--ancestry", default=None,
                    help="default ALL (whole cohort)")

    args = p.parse_args()
    args.layer = resolve(args.layer, LAYER_ALIAS, "layer")
    args.pt510 = resolve(args.pt510, PT510_ALIAS, "pt510", case_sensitive=True)
    if args.cmd == "top" and not args.sku and not args.all:
        raise SystemExit("top needs --sku ANC:CAT (e.g. AA:h) or --all")
    {"list": cmd_list, "top": cmd_top, "absent": cmd_absent,
     "gene": cmd_gene}[args.cmd](args)


if __name__ == "__main__":
    main()
