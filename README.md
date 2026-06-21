# Population-Aware Genomic Analysis of Prostate Cancer Severity — analysis code

Analysis code for the study *"Population-Aware Genomic Analysis of Prostate Cancer
Severity: Exploratory Signals from SNV-Derived Gene-Level Profiles."* The manuscript is
maintained separately (Overleaf) and is not part of this repository.

The study asks whether prostate-cancer severity (Gleason grade) is encoded in
sequencing-derived, gene-level single-nucleotide-variant (SNV) profiles, and whether
that structure differs by ancestry, in a 20-patient multi-ancestry cohort (African
American, European American, Indian) compared against the external NanoString BC360
reference cohort GSE244631 (Brazil).

## Getting started

You need two things: **this repository** (code) and **`data_bundle.zip`** (the
de-identified data, sent to you separately — it is not, and must not be, in the repo).
There is nothing else to download and no preprocessing to run.

```bash
# 1. Get the code
git clone <repo-url>
cd <repo>

# 2. Add the data: extract data_bundle.zip INTO the repo folder, so it creates
#    cohort_matched/, GSE244631/, genefilelists/, snv_rebuild/ right here:
unzip /path/to/data_bundle.zip                       # Linux / macOS (terminal)
#   macOS (Finder):   double-click the .zip, then move those 4 folders into the repo
#   Windows (PowerShell):  Expand-Archive C:\path\to\data_bundle.zip -DestinationPath .
#   Windows (Explorer):    right-click the .zip -> Extract All... -> into the repo folder

# 3. Set up the Python environment
conda env create -f py3-environment.yml
conda activate py3

# 4. Run any analysis (full list under "Reproducing the main results"), e.g.:
python panel_recovery_aaea.py
```

After step 2 the folder looks like this — repo code and bundle data side by side:

```
<repo>/
├── distance_knn.py  panel_recovery.py  viz_gradescore.py  ...   (code, from the repo)
├── cohort_matched/          20 gene-count CSVs            ┐
├── GSE244631/               NanoString .RCC reference     │  from
├── genefilelists/           metadata_orig.csv             │  data_bundle.zip
├── snv_rebuild/             SNV master tables             ┘
└── ...
```

## Reproducing the main results

With the data unzipped in place and `conda activate py3` done, run any of these. Each
reads its inputs from the unzipped folders and writes outputs (CSVs / PNGs) into the
current directory.

```bash
# Cross-modality severity analysis (the panel-restricted comparison)
python panel_recovery.py        # full 20-patient cohort       -> cosine 0.81 / Spearman 0.78
python panel_recovery_aaea.py   # broadly sequenced AA+EA only -> chance (0.62 / 0.67)
python quick_norm_sweep.py      # self- vs per-gene normalization sweep

# Final figure: cohort-normalized grade-association score
python viz_gradescore.py        # writes fig_gradescore.png

# Supporting diagnostics
python viz_reference.py         # within-reference separability + PCA
python viz_query_pergene.py     # per-gene-z query projection
python diag_aaea.py             # exact permutation-null granularity

# Within-cohort variant / gene analyses
python disproportion.py         # Fisher's exact test, per-variant severity skew
python disproportion_genes.py   # same, aggregated to genes
python build_paper_results.py   # result tables
python build_paper_figures.py   # variant-burden / population-divergence figures
python build_table4.py          # Table: severity-exclusive variants
python build_table5.py          # Table: severity-exclusive genes
```

## What's in `data_bundle.zip`

Everything the code needs, all de-identified (keyed by code, no patient identifiers):

| Created on unzip | Contents | Read by |
|---|---|---|
| `cohort_matched/` | 20 per-patient `gene,variant_count` CSVs | the cross-modality scripts |
| `snv_rebuild/PT510_low/…` | SNV master tables | `disproportion*`, `build_paper_results`, `build_table4/5` |
| `genefilelists/metadata_orig.csv` | cohort `ID, race, Gleason` (truth labels) | all scripts |
| `GSE244631/` | public NanoString `.RCC` reference | the cross-modality scripts |

`GSE244631_ref_labels.csv` (reference high/low labels) and `biomarker_subset.csv` (the
23-gene candidate panel) ship **in the repo**; the panel file is also re-written by
`panel_recovery.py`.

## Data, privacy, and provenance

- **Code only in git.** Raw/patient-identifiable data (per-patient VCFs, FASTQ, names)
  are never committed; they live only in the separately-shared bundle, and a
  `.gitignore` enforces this. The empty data folders carry a `PLACEHOLDER.md` until the
  bundle is unzipped over them.
- **Raw-data rebuild (optional).** The bundle is the light, ready-to-run product. To
  regenerate it from raw sequencing data — including where to download GEO GSE244631,
  ClinVar, COSMIC, and GENCODE — see **`UPSTREAM.md`**. Not needed for the analyses above.
- **History.** Superseded/scratch scripts are preserved under **`deprecated/`** (see
  `deprecated/README.md`); none are needed to reproduce the results.

## Environment

Python via conda (`py3`): `numpy`, `pandas`, `scipy`, `matplotlib`. The raw-VCF parsing
stage additionally uses `pysam` (WSL/Linux), but that is only for the `UPSTREAM.md`
rebuild, not for the analyses above.

## Scripts

- **`distance_knn.py`** — shared loaders and the cross-modality classifier
  (reference/query loading, normalization, distances, class-balanced voting).
- **`panel_recovery.py` / `panel_recovery_aaea.py`** — panel-restricted comparison on
  the full cohort vs the broadly sequenced AA+EA stratum.
- **`quick_norm_sweep.py`, `diag_aaea.py`** — normalization sweep and permutation-null
  diagnostics behind the AA+EA result.
- **`viz_gradescore.py`** — the grade-association-score figure.
- **`viz_reference.py`, `viz_query_pergene.py`** — separability / projection diagnostics.
- **`disproportion.py`, `disproportion_genes.py`** — Fisher's-exact severity-skew tests
  (variant- and gene-level), FDR-corrected.
- **`build_paper_results.py`, `build_paper_figures.py`, `build_table4.py`,
  `build_table5.py`, `perm_test.py`** — result tables, figures, permutation testing.

## Citation

Manuscript in preparation. Funding: AIM-AHEAD Program, Agreement No. 1OT2OD032581.
