"""
manuscript_results.py - generate paste-ready result tables for the manuscript
("Population-Aware Genomic Analysis ...") from the snv_rebuild/ outputs.

It does NOT touch the paper. It writes LaTeX snippets + a readable summary into
manuscript_results/, aligned to the draft's placeholders:

  IV.B/C  SNV sets by population x severity, per annotation bucket (layer)
  IV.D    gene-mapping counts by population x severity, per layer
  IV.E    cross-population overlap (shared) vs population-enriched genes
  Tab.IV  computational-signal context for the candidate genes

Terminology bridge (paper <-> pipeline):
  pathogenic = matched   benign = benign   unknown = unmatched
  high-severity-enriched = high_only   low-severity-enriched = low_only
  shared/common = high_low_common       PT510 = high (matches draft Table II)

Run:  python manuscript_results.py            (PT510=high, all layers)
"""

import csv
import os
from pathlib import Path

from query_genes import (ANCESTRIES, CATEGORIES, load_sku, load_slice,
                         gene_file, metric_value)

PT510 = "high"
LAYERS = [("matched", "pathogenic / likely pathogenic"),
          ("benign", "benign / likely benign"),
          ("unmatched", "unknown / unannotated")]
# present populations in the order the paper uses
POP_ORDER = ["ALL", "IN", "EA", "AA"]
CANDIDATES = ["APC", "ATM", "CDH1", "TP53", "E2F1", "SIRT1", "PRKAB2",
              "EGFR", "KDR", "STK11"]
OUT = Path("manuscript_results")
ROOT = Path("snv_rebuild")


def read_summary(layer):
    """{group: {...}} from snv_rebuild/PT510_x/<layer>/summary_<layer>.csv."""
    path = ROOT / f"PT510_{PT510}" / layer / f"summary_{layer}.csv"
    rows = {}
    if path.exists():
        with open(path, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                rows[r["Group"]] = r
    return rows


def tex(s):
    """Escape LaTeX-special chars in a gene name / free string."""
    s = str(s)
    for a, b in (("\\", ""), ("_", r"\_"), ("%", r"\%"), ("&", r"\&"),
                 ("#", r"\#"), ("$", r"\$")):
        s = s.replace(a, b)
    return s


def latex_table(caption, label, header, body_rows):
    cols = "l" + " c" * (len(header) - 1)
    lines = [r"\begin{table}[t]", r"\centering", rf"\caption{{{caption}}}",
             rf"\label{{{label}}}", r"\small",
             rf"\begin{{tabular}}{{{cols}}}", r"\hline",
             " & ".join(header) + r" \\", r"\hline"]
    lines += [" & ".join(str(c) for c in row) + r" \\" for row in body_rows]
    lines += [r"\hline", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lines)


def snv_set_tables():
    """IV.B / IV.C - SNV sets by population x severity, one table per layer."""
    out = ["% ===== IV.B / IV.C : SNV sets by population x severity =====\n"]
    for layer, pretty in LAYERS:
        s = read_summary(layer)
        body = []
        for g in POP_ORDER:
            if g not in s:
                continue
            r = s[g]
            total = int(r["HighOnly"]) + int(r["LowOnly"]) + int(r["HighLowCommon"])
            body.append([g, r["HighPatients"], r["LowPatients"],
                         r["HighOnly"], r["LowOnly"], r["HighLowCommon"], total])
        out.append(latex_table(
            f"Recurrent SNV sets by population and severity "
            f"({pretty} bucket).",
            f"tab:snv_{layer}",
            ["Population", "High pts", "Low pts", "High-enriched",
             "Low-enriched", "Shared", "Total"],
            body))
    return "\n".join(out)


def gene_counts(layer):
    """{(pop,cat): n_genes} for one layer."""
    return {(p, c): len(load_sku(p, c, layer, PT510))
            for p in ANCESTRIES for c in CATEGORIES}


def gene_map_tables():
    """IV.D - mapped-gene counts by population x severity, per layer."""
    out = ["% ===== IV.D : mapped gene counts by population x severity =====\n"]
    for layer, pretty in LAYERS:
        gc = gene_counts(layer)
        body = []
        for g in POP_ORDER:
            row = [g] + [gc[(g, c)] for c in CATEGORIES]
            row.append(len(load_slice(g, "all", layer, PT510)))
            body.append(row)
        out.append(latex_table(
            f"Mapped-gene counts by population and severity ({pretty} bucket).",
            f"tab:gene_{layer}",
            ["Population", "High-enriched", "Low-enriched", "Shared", "All"],
            body))
    return "\n".join(out)


def topn(pop, cat, layer, n):
    genes = load_sku(pop, cat, layer, PT510)
    return [g for g, _ in sorted(genes.items(),
            key=lambda kv: metric_value(kv[1], "total"), reverse=True)[:n]]


def cross_group_tables(layer="matched", n=20):
    """IV.E - shared (recurrent across populations) vs population-enriched."""
    out = [f"% ===== IV.E : cross-population (layer={layer}) =====\n"]

    # shared/recurrent: genes in high_low_common of every sub-population
    common_sets = [set(load_sku(p, "high_low_common", layer, PT510))
                   for p in ["AA", "EA", "IN"]]
    shared = sorted(set.intersection(*common_sets)) if common_sets else []
    out.append(f"% Recurrent genes shared across AA, EA, IN (high_low_common): "
               f"{len(shared)}\n% " + ", ".join(shared[:50]) + "\n")

    # population-enriched: top high_only genes unique to one population
    rows = []
    for p in ["AA", "EA", "IN"]:
        mine = set(load_sku(p, "high_only", layer, PT510))
        others = set()
        for q in ["AA", "EA", "IN"]:
            if q != p:
                others |= set(load_sku(q, "high_only", layer, PT510))
        uniq = topn_filtered(p, "high_only", layer, n, mine - others)
        rows.append([p, ", ".join(tex(x) for x in uniq[:10])])
    out.append(latex_table(
        f"Population-enriched high-severity genes (top, high\\_only unique to "
        f"each population; {layer}).",
        "tab:pop_enriched",
        ["Population", "Top population-unique high-enriched genes"], rows))

    # high_only ALL top-N with per-population high-frequency (the reconciled view)
    out.append(highonly_population_table(layer, n))
    return "\n".join(out)


def topn_filtered(pop, cat, layer, n, keep):
    genes = load_sku(pop, cat, layer, PT510)
    genes = {g: r for g, r in genes.items() if g in keep}
    return [g for g, _ in sorted(genes.items(),
            key=lambda kv: metric_value(kv[1], "total"), reverse=True)[:n]]


def highonly_population_table(layer, n):
    """Top-N ALL high_only genes; per-population HIGH-patient frequency (additive,
    gap-free reconciliation discussed for IV.E)."""
    all_high = load_sku("ALL", "high_only", layer, PT510)
    top = [g for g, _ in sorted(all_high.items(),
           key=lambda kv: kv[1]["total"], reverse=True)[:n]]
    # high-patient counts per population for each gene
    pop_high = {p: load_slice(p, "all", layer, PT510) for p in ["IN", "AA", "EA"]}
    pop_npat = {p: _n_high_patients(p, layer) for p in ["IN", "AA", "EA"]}
    body = []
    for i, g in enumerate(top, 1):
        cells = [i, tex(g), all_high[g]["total"]]
        for p in ["IN", "AA", "EA"]:
            rec = pop_high[p].get(g)
            if rec and len(rec["high"]) > 0:
                cells.append(f"{len(rec['high'])}/{pop_npat[p]}")
            else:
                cells.append("--")
        body.append(cells)
    return latex_table(
        f"Top {n} cohort high-severity-enriched genes with per-population "
        f"high-patient support (carriers/high patients; {layer}).",
        "tab:highonly_pop",
        ["\\#", "Gene", "ALL total", "IN", "AA", "EA"], body)


def _n_high_patients(pop, layer):
    s = read_summary(layer)
    return int(s[pop]["HighPatients"]) if pop in s else 0


def candidate_table(layer="matched"):
    """Table IV - computational-signal context for the candidate genes."""
    out = ["% ===== Table IV : candidate-gene computational signal context =====\n"]
    body = []
    for gene in CANDIDATES:
        present = []
        for c in CATEGORIES:
            rec = load_sku("ALL", c, layer, PT510).get(gene)
            if rec:
                tag = {"high_only": "high-enriched", "low_only": "low-enriched",
                       "high_low_common": "shared"}[c]
                present.append(f"{tag} ({rec['total']})")
        pops = [p for p in ["AA", "EA", "IN"]
                if gene in load_slice(p, "all", layer, PT510)]
        ctx = "; ".join(present) if present else "not recurrent"
        if pops:
            ctx += "; in " + "/".join(pops)
        body.append([tex(gene), ctx])
    out.append(latex_table(
        f"Candidate-gene computational signal context ({layer} bucket, PT510=high).",
        "tab:candidate_context",
        ["Candidate gene", "Computational signal context"], body))
    return "\n".join(out)


def main():
    OUT.mkdir(exist_ok=True)
    sections = {
        "IV_B_C_snv_sets.tex": snv_set_tables(),
        "IV_D_gene_map.tex": gene_map_tables(),
        "IV_E_cross_group.tex": cross_group_tables(),
        "TableIV_candidates.tex": candidate_table(),
    }
    for name, text in sections.items():
        (OUT / name).write_text(text, encoding="utf-8")
        print(f"wrote {OUT/name}")
    combined = "\n\n".join(sections.values())
    (OUT / "all_results.tex").write_text(combined, encoding="utf-8")
    print(f"wrote {OUT/'all_results.tex'}")
    print("\n--- preview ---\n")
    print(combined)


if __name__ == "__main__":
    main()
