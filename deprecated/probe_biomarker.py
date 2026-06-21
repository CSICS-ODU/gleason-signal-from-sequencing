# -*- coding: utf-8 -*-
"""probe_biomarker.py -- where do the 22 genes in biomarker_subset.csv come from?

Reproduces the EXACT selection the leakage-free pipeline uses (nested_loo.py):
candidate pool = genes shared between the GSE244631 reference panel and the cohort;
rank by carrier-rate gap = carrier%(high) - carrier%(low); take top-11 high-enriched
+ top-11 low-enriched = 22. Then compares that computed set against the genes actually
listed in biomarker_subset.csv, to see if the file is reproducible or hand-curated.
"""
import types
import numpy as np
import pandas as pd
import distance_knn as dk

# ---- load exactly as nested_loo.py does -----------------------------------
ref_mat = dk.load_reference(types.SimpleNamespace(ref_rcc_dir="GSE244631", ref_matrix=None))
qprof = dk.load_query_profiles("cohort_matched")
truth = dk.load_query_truth("genefilelists/metadata_orig.csv", "ID", "Gleason")
truth = {p: v for p, v in truth.items() if p in qprof}
pids = [p for p in qprof if p in truth]

panel = set(ref_mat.index)                       # genes measurable in the reference
qgenes = set().union(*[set(qprof[p].index) for p in pids])
panel_genes = sorted(panel & qgenes)             # the candidate pool
carrier = {p: {g: (qprof[p].get(g, 0) > 0) for g in panel_genes} for p in pids}

hi = [p for p in pids if truth[p] == "high"]
lo = [p for p in pids if truth[p] == "low"]
print(f"cohort: {len(pids)} ({len(hi)} high / {len(lo)} low)")
print(f"reference panel genes: {len(panel)} | cohort genes: {len(qgenes)} | "
      f"shared candidate pool: {len(panel_genes)}\n")

gap = {g: (sum(carrier[p][g] for p in hi) / len(hi)
           - sum(carrier[p][g] for p in lo) / len(lo)) for g in panel_genes}
order = sorted(panel_genes, key=lambda g: gap[g])
low_sel = order[:11]                              # most low-enriched
high_sel = order[-11:]                            # most high-enriched
computed = low_sel + high_sel

print("COMPUTED full-cohort selection (carrier-rate gap, 11 low + 11 high):")
print("  low-enriched :", ", ".join(f"{g}({gap[g]:+.2f})" for g in low_sel))
print("  high-enriched:", ", ".join(f"{g}({gap[g]:+.2f})" for g in high_sel))

# ---- compare against the file ---------------------------------------------
filed = pd.read_csv("biomarker_subset.csv").iloc[:, 0].astype(str).str.strip().tolist()
filed = [g for g in filed if g]
cset, fset = set(computed), set(filed)
print(f"\nbiomarker_subset.csv has {len(filed)} genes")
print(f"computed selection has {len(computed)} genes")
print(f"  overlap (in both)      : {len(cset & fset)}  {sorted(cset & fset)}")
print(f"  in FILE but NOT computed: {len(fset - cset)}  {sorted(fset - cset)}")
print(f"  COMPUTED but not in file: {len(cset - fset)}  {sorted(cset - fset)}")

# are the file's genes even in the candidate pool / reference at all?
not_in_pool = [g for g in filed if g not in set(panel_genes)]
not_in_ref = [g for g in filed if g not in panel]
not_in_cohort = [g for g in filed if g not in qgenes]
print(f"\nfile genes NOT in shared candidate pool: {len(not_in_pool)}  {not_in_pool}")
print(f"file genes NOT on reference panel       : {len(not_in_ref)}  {not_in_ref}")
print(f"file genes NOT in cohort variants       : {len(not_in_cohort)}  {not_in_cohort}")

# ---- the smoking gun: are the "high-enriched" genes just Indian-panel genes? ----
# Indian targeted-panel patients have all-numeric IDs; AA/EA are PT-prefixed.
def origin(p):
    return "IN" if str(p).isdigit() else "AA/EA"

print("\nancestry of carriers for the COMPUTED high-enriched selection:")
print(f"  ({sum(origin(p)=='IN' for p in hi)} of {len(hi)} high patients are Indian-panel)")
for g in high_sel:
    hc = [p for p in hi if carrier[p][g]]
    lc = [p for p in lo if carrier[p][g]]
    nin = sum(origin(p) == "IN" for p in hc)
    print(f"  {g:<8} high-carriers={len(hc):2} (IN={nin:2}, AA/EA={len(hc)-nin}) | "
          f"low-carriers={len(lc)}")

print("\nancestry of carriers for the COMPUTED low-enriched selection:")
for g in low_sel:
    hc = [p for p in hi if carrier[p][g]]
    lc = [p for p in lo if carrier[p][g]]
    nin = sum(origin(p) == "IN" for p in lc)
    print(f"  {g:<8} low-carriers={len(lc):2} (IN={nin:2}, AA/EA={len(lc)-nin}) | "
          f"high-carriers={len(hc)}")

# for each file gene, show gap, rank, and carrier ancestry -> classify the source
print("\nper-file-gene detail (gap, rank, carrier ancestry):")
print(f"  {'gene':<10}{'gap':>7}  {'rank':>8}  high:IN/US  low:IN/US  likely source")
ranked = {g: i for i, g in enumerate(order)}     # 0 = most low-enriched
our_candidates = {"PTEN", "TP53", "PIK3CA", "SMAD4", "GNAS"}   # Table IX (AA+EA)
for g in filed:
    if g not in gap:
        print(f"  {g:<10} NOT in candidate pool")
        continue
    hc = [p for p in hi if carrier[p][g]]
    lc = [p for p in lo if carrier[p][g]]
    h_in = sum(origin(p) == "IN" for p in hc)
    l_in = sum(origin(p) == "IN" for p in lc)
    # heuristic: high-lean + carried only by Indians = clinical-panel; otherwise AA/EA
    if g in our_candidates:
        src = "OUR AA+EA candidate (Table IX)"
    elif gap[g] > 0 and h_in == len(hc) and len(hc) > 0:
        src = "Indian panel (all-IN carriers)"
    elif gap[g] < 0 and l_in == 0:
        src = "AA+EA low-enriched"
    else:
        src = "mixed/other"
    print(f"  {g:<10}{gap[g]:>+7.2f}  {ranked[g]:>4}/{len(order)-1}  "
          f"{len(hc):>5}:{h_in}/{len(hc)-h_in}  {len(lc):>4}:{l_in}/{len(lc)-l_in}  {src}")
