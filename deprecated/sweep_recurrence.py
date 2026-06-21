"""Sweep the recurrence threshold (gene-level, full cohort) and show how the
gene count and FDR survivors change -- testing whether a higher criterion than
the paper's >=2 yields FDR-significant hits (fewer tests -> lighter correction).
"""
import disproportion_genes as dg

print(f"{'layer':<10}{'min_total':>10}{'rec>=':>7}{'n_genes':>9}"
      f"{'hi_p<.05':>10}{'hi_q<.05':>10}{'lo_p<.05':>10}{'lo_q<.05':>10}")
for path in dg.default_files():
    layer = "matched" if "matched" in path else ("benign" if "benign" in path else "unmatched")
    for mt in (1, 2, 3, 4, 5, 6):
        res = dg.analyze(path, mt, no_indian=False)
        if res is None:
            continue
        H, L, out = res
        hi = [d for d in out if d["dir"] == "HIGH"]
        lo = [d for d in out if d["dir"] == "LOW"]
        hp = sum(d["p"] < 0.05 for d in hi); hq = sum(d["q"] < 0.05 for d in hi)
        lp = sum(d["p"] < 0.05 for d in lo); lq = sum(d["q"] < 0.05 for d in lo)
        print(f"{layer:<10}{mt:>10}{mt+1:>7}{len(out):>9}"
              f"{hp:>10}{hq:>10}{lp:>10}{lq:>10}")
    print()
