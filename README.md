# Population-Aware Genomic Analysis of Prostate Cancer Grade — analysis code

Analysis code for the study *"Population-Aware Genomic Analysis of Prostate Cancer Grade:
Exploratory Signals from SNV-Derived Gene-Level Profiles."* The manuscript itself is not
part of this repository.

The study asks whether prostate cancer grade (Gleason grade, high versus low) is
reflected in gene-level profiles of single-nucleotide variants (SNVs) from tumor
sequencing, and whether that structure differs by population, in a 20-patient cohort
(African American, AA; European American, EA; Indian, IN). The profiles are also
compared with an external gene-expression reference cohort, GEO GSE244631 (Brazil;
NanoString Breast Cancer 360 panel).

## Getting started

You need two things: **this repository** (code) and **`data_bundle.zip`** (the
de-identified analysis inputs, available from the corresponding author on reasonable
request; it is not, and must not be, in this repository). Nothing else needs to be
downloaded and no preprocessing is needed.

```bash
# 1. Get the code
git clone https://github.com/CSICS-ODU/Prostate-Cancer-Research.git
cd Prostate-Cancer-Research

# 2. Add the data: extract data_bundle.zip INTO the repository folder, so it creates
#    cohort_matched/, GSE244631/, genefilelists/, snv_rebuild/ right here:
unzip -o /path/to/data_bundle.zip                    # Linux / macOS (terminal)
#   macOS (Finder):   double-click the .zip, then move those 4 folders into the repo
#   Windows (PowerShell):  Expand-Archive C:\path\to\data_bundle.zip -DestinationPath . -Force
#   Windows (Explorer):    right-click the .zip -> Extract All... -> into the repo folder

# 3. Set up the Python environment (conda), or: pip install -r requirements.txt
conda env create -f py3-environment.yml
conda activate py3

# 4. Run any analysis (full list below), e.g.:
python table2_exact.py
```

After step 2 the folder looks like this:

```
Prostate-Cancer-Research/
├── distance_knn.py  table2_exact.py  ...   (code, from the repository)
├── cohort_matched/          20 per-patient gene-count CSVs   ┐
├── GSE244631/               NanoString .RCC reference files  │  from
├── genefilelists/           metadata_orig.csv                │  data_bundle.zip
├── snv_rebuild/             SNV master tables                ┘
└── ...
```

## Reproducing the manuscript

Run from the repository folder, with the data unzipped in place. Every script reads its
inputs from the unzipped folders and writes its outputs (printed tables, CSVs, PNGs) to
the current folder. Each runs in seconds; `table2_exact.py` takes under a minute.

| Manuscript item | Command | Output |
|---|---|---|
| Table 1 (patients; variants by category and group) | `python build_paper_results.py` | `paper_results/tableIVB_bucket_counts.csv` |
| Fig. 2 and Fig. 3a | `python make_paper_figures.py` (after `build_paper_results.py`) | `fig_burden.png`, `fig_lean.png` |
| Table 2a (genome-wide comparison) | `python distance_knn.py --ref-rcc-dir GSE244631 --ref-labels GSE244631_ref_labels.csv --query-dir cohort_matched --query-meta genefilelists/metadata_orig.csv --metrics cosine,spearman` | printed summary: cosine 0.50 (p 1.00); Spearman best k=3, 0.60 (p 0.27) |
| Table 2b (23-gene panel), all four rows with their exact *p* | `python table2_exact.py` | printed table |
| Fig. 7 (the 23-gene candidate panel) | `python panel_recovery.py` | `biomarker_subset.csv` |
| Table 3a (variants found in only one grade) | `python build_table4.py` | `paper_results/table4_exclusive_aaea.csv` |
| Table 3b (genes found in only one grade) | `python build_table5.py` | `paper_results/table5_exclusive_aaea.csv` |
| Fig. 3b counts | `python build_paper_results.py` | `paper_results/tableIVE_cross_pop.csv` |
| Fig. 3b *p* and *q*; variant-burden test (*p* = 0.030); gene-level correlation | `python exact_tests_aaea.py` | printed |
| Table 4a, DNA-variant column | `python build_paper_results.py` | `paper_results/tableIV_candidates.csv` |
| Fig. 8 (grade-association score) | `python viz_gradescore.py` | `fig_gradescore.png` |
| No gene survives FDR correction (floor *q* ≈ 0.08) | `python build_paper_results.py` | `paper_results/tableIVD_gene_mapping.csv` (lowest *q* 0.081) |
| Per-variant and per-gene grade tests at other carrier thresholds | `python disproportion.py`, `python disproportion_genes.py` (add `--no-indian` for AA+EA) | printed |

Expected Table 2b output (`table2_exact.py`):

```
Table 2b row                        cosine bal (p)    Spearman bal (p)  allocations
Full cohort, fixed panel             0.813 (0.012)       0.780 (0.022)       77,520
AA+EA, fixed panel                   0.625 (0.400)       0.667 (0.340)          100
AA+EA, nested selection              0.250 (0.910)       0.708 (0.200)          100
AA+EA, cohort-standardized           0.708 (0.270)       0.833 (0.090)          100
```

Notes:

- **Exact *p*-values.** `table2_exact.py` and `exact_tests_aaea.py` enumerate every label
  allocation: all 77,520 for the full cohort (13 high / 7 low), and for the 10 AA+EA
  patients the 100 allocations that keep 2 high and 3 low within each of AA and EA, so no
  result can arise from population alone. `panel_recovery.py`, `panel_recovery_aaea.py`
  and `quick_norm_sweep.py` print the same balanced accuracies with *p*-values from
  1,000 random label shuffles, which differ slightly from the exact values.
- **Nested selection.** `python nested_loo.py AA,EA` is the nested-selection analysis on
  its own (cosine 0.25, Spearman 0.71). Within each leave-one-out fold, genes are ranked
  by the difference in the fraction of high- and low-grade training patients with a
  variant, and the 11 lowest and 11 highest are kept from a list sorted by that
  difference with ties in gene-symbol order.
- **Output names.** `build_paper_results.py` names its CSVs after the section numbers of
  an earlier draft (`tableIVB`, `tableIVE`, ...). Its annotation categories map to the
  manuscript as: `pathogenic/known` = catalogued, not benign; `benign/likely-benign` =
  catalogued benign/likely benign; `unknown/unannotated` = not catalogued.
- **Supporting diagnostics** (not in the manuscript's tables): `panel_recovery_aaea.py`,
  `quick_norm_sweep.py`, `diag_aaea.py`, `viz_reference.py`, `viz_query_pergene.py`.

## What is in `data_bundle.zip`

Everything the code needs, all de-identified:

| Created on unzip | Contents | Read by |
|---|---|---|
| `cohort_matched/` | 20 per-patient `gene,variant_count` CSVs: catalogued, non-benign SNVs found in at least two patients, counted per gene | the comparison scripts, `exact_tests_aaea.py` |
| `snv_rebuild/PT510_low/…` | SNV master tables per annotation category and group | `build_paper_results.py`, `build_table4/5.py`, `disproportion*.py`, `exact_tests_aaea.py` |
| `genefilelists/metadata_orig.csv` | cohort `ID, race, Gleason` | all scripts |
| `GSE244631/` | the public NanoString `.RCC` reference files | the comparison scripts |

`GSE244631_ref_labels.csv` (reference high/low labels; Gleason 3+4 counts as low, as in
the cohort) and `biomarker_subset.csv` (the 23-gene candidate panel) are part of the
repository. `panel_recovery.py` re-writes the panel file, keeping the previous copy as
`biomarker_subset_OLD.csv`.

## Data, privacy, and provenance

- **Code only in git.** Patient-level data (VCFs, sequencing reads, identifiers) are
  never committed, and `.gitignore` enforces this. The empty data folders carry a
  `PLACEHOLDER.md` until the bundle is unzipped over them.
- **Raw-data rebuild (optional).** The bundle is the ready-to-run product of the
  upstream processing. How it is built from raw data (genome build GRCh37; ClinVar,
  COSMIC and GENCODE versions; download locations) is described in **`UPSTREAM.md`**. It
  is not needed for the analyses above.
- **History.** Superseded or exploratory scripts are kept under **`deprecated/`** (see
  `deprecated/README.md`); none are needed to reproduce the manuscript.

## Environment

Python with `numpy`, `pandas`, `scipy` and `matplotlib` (`py3-environment.yml` or
`requirements.txt`). The results were produced under Python 3.13 and re-verified under
Python 3.14 (numpy 2.3.5, pandas 2.3.3, scipy 1.16.3, matplotlib 3.10.7). The raw-VCF
parsing stage in `UPSTREAM.md` additionally uses `pysam`.

## Scripts

- **`distance_knn.py`**: shared loaders and the cross-modality comparison
  (reference/query loading, normalization, distances, class-balanced voting); run
  directly, it gives the genome-wide comparison (Table 2a).
- **`table2_exact.py`**: Table 2b, every row, with exact permutation *p*-values.
- **`panel_recovery.py`, `panel_recovery_aaea.py`**: the panel comparison on the full
  cohort and on the AA+EA patients; `panel_recovery.py` also writes the panel.
- **`nested_loo.py`**: nested (leave-one-out) gene selection.
- **`quick_norm_sweep.py`, `diag_aaea.py`**: normalization schemes and the permutation
  null behind the AA+EA results.
- **`exact_tests_aaea.py`**: variant-burden test, Fig. 3b *p*/*q*, gene-level correlation.
- **`build_paper_results.py`, `build_table4.py`, `build_table5.py`**: variant and gene
  tables.
- **`disproportion.py`, `disproportion_genes.py`**: Fisher's exact tests of grade
  differences per variant and per gene, with Benjamini–Hochberg FDR.
- **`make_paper_figures.py`, `viz_gradescore.py`**: Figs. 2, 3a and 8.
- **`viz_reference.py`, `viz_query_pergene.py`**: separability and projection diagnostics.

## Citation

Manuscript submitted. Funding: AIM-AHEAD Program, Agreement No. 1OT2OD032581.
