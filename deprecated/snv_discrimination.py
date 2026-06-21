"""
snv_discrimination.py - high/low severity discrimination at the SNV level, LABEL-FREE.

The matched/benign/unmatched layers come from a ClinVar/COSMIC LOOKUP. Restricting the
signal to "matched" lets that external annotation decide what matters (and imports its
known-gene / population biases). This does the discrimination on ALL SNVs WITHOUT using
the lookup: each SNV (a specific chr:pos:ref>alt) is a label-free feature; we just ask
which SNVs are carried disproportionately by high vs low patients. The ClinVar layer is
kept only as a POST-HOC TAG, to check whether the discriminating SNVs are pathogenic,
benign, or novel -- i.e. whether the signal depends on the lookup at all.

SNV level (not gene level) on purpose: unioning layers at the GENE level makes big genes
universal and washes the signal out; SNVs stay specific, so all-variant discrimination
is preserved.

Per ancestry (PT510_low): union matched+benign+unmatched, report top discriminating SNVs
(gene + source layer) and the layer breakdown of the discriminators.
"""

import argparse
import glob
import csv
from collections import Counter

from scipy.stats import fisher_exact


def layer_file(anc, layer):
    """One master per layer: prefer the annotated *_genes.csv, else the plain master
    (unmatched AA/EA has no gene annotation)."""
    fs = glob.glob(f"snv_rebuild/PT510_low/{layer}/{anc}/**/{anc}_master_{layer}*.csv",
                   recursive=True)
    genes = [f for f in fs if f.endswith("_genes.csv")]
    return (genes or fs or [None])[0]


def load_all(anc):
    """Union all SNVs across layers; tag each with its ClinVar layer (post-hoc only)."""
    H, L, snvs = set(), set(), {}
    for layer in ("matched", "benign", "unmatched"):
        f = layer_file(anc, layer)
        if not f:
            continue
        with open(f, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                hp = {x for x in (r.get("high_patients") or "").split(";") if x}
                lp = {x for x in (r.get("low_patients") or "").split(";") if x}
                H |= hp
                L |= lp
                gs = [g for g in (r.get("genes") or "").split(";") if g]
                snvs[r["SNV"]] = {"hp": hp, "lp": lp, "layer": layer,
                                  "gene": gs[0] if gs else "-"}
    return len(H), len(L), snvs


def analyze(anc, min_total):
    H, L, snvs = load_all(anc)
    out = []
    for s, d in snvs.items():
        nh, nl = len(d["hp"]), len(d["lp"])
        if nh + nl <= min_total:
            continue
        hr, lr = (nh / H if H else 0), (nl / L if L else 0)
        odds, p = fisher_exact([[nh, H - nh], [nl, L - nl]])
        out.append({"snv": s, "gene": d["gene"], "layer": d["layer"],
                    "nh": nh, "nl": nl, "hr": hr, "lr": lr, "diff": hr - lr, "p": p,
                    "dir": "HIGH" if hr > lr else ("LOW" if lr > hr else "=")})
    return H, L, out


def show(anc, H, L, out, top):
    disc = [d for d in out if d["dir"] != "="]
    print("=" * 66)
    print(f"{anc}: {H} high / {L} low | {len(out)} SNVs tested (>min-total), "
          f"{len(disc)} discriminating")
    # the key question: do discriminators depend on the ClinVar lookup?
    by_layer = Counter(d["layer"] for d in disc)
    tot_layer = Counter(d["layer"] for d in out)
    print("  discriminating SNVs by source layer (vs all tested):")
    for layer in ("matched", "benign", "unmatched"):
        print(f"    {layer:<10} {by_layer.get(layer,0):>6} / {tot_layer.get(layer,0):<6}")
    for label, key in (("HIGH-skewed", lambda d: (-d["diff"], d["p"])),
                       ("LOW-skewed", lambda d: (d["diff"], d["p"]))):
        side = [d for d in disc if d["dir"] == ("HIGH" if "HIGH" in label else "LOW")]
        print(f"\n  {label} (top {top}, effect size)")
        print(f"    {'SNV':<22}{'gene':<12}{'layer':<10}{'hi':>3}{'lo':>3}{'p':>8}")
        for d in sorted(side, key=key)[:top]:
            print(f"    {d['snv']:<22}{d['gene'][:11]:<12}{d['layer']:<10}"
                  f"{d['nh']:>3}{d['nl']:>3}{d['p']:>8.2g}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-total", type=int, default=1,
                    help="keep SNVs carried by strictly more than this many patients")
    ap.add_argument("--top", type=int, default=12)
    ap.add_argument("--anc", default="AA,EA", help="comma list of ancestries")
    a = ap.parse_args()
    for anc in [x.strip() for x in a.anc.split(",")]:
        H, L, out = analyze(anc, a.min_total)
        show(anc, H, L, out, a.top)


if __name__ == "__main__":
    main()
