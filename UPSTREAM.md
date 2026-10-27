# Upstream data processing (how the analysis inputs are built)

The scripts in this repository start from **processed gene-level inputs**
(`cohort_matched/`, `GSE244631_ref_labels.csv`). This document describes how those
inputs are produced from raw data. The raw-processing scripts themselves are **not
included** here: they require `pysam`, the raw sequencing data, and large annotation
databases, and some reference patient identifiers. They are retained in the authors'
working environment and can be shared on request, subject to the data constraints
below.

> **Genome build:** GRCh37 (hg19). All variant coordinates and annotation databases must
> use the same build.

---

## 1. Data sources

| Data | Source / where to download | Notes |
|------|----------------------------|-------|
| **GSE244631 reference** (NanoString BC360 expression) | NCBI GEO, accession **GSE244631** — `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244631` (download the supplementary `.RCC` files) | 51 labeled prostate samples used here (19 high / 32 low). Place in `GSE244631/`. |
| **Cohort: African American + European American** (5 + 5) | Moffitt Cancer Center — tumor-only whole-genome sequencing called with DRAGEN 4.4.4, `*.hard-filtered.vcf.gz` (GRCh37) | Not public. Institutional data-use agreement; de-identified. |
| **Cohort: Indian** (10) | Thakurpukur Cancer Center / Saroj Gupta Cancer Centre & Research Institute — CHPV2 targeted clinical panel, per-sample `*-full.tsv` | Not public. **Source files carry patient names; de-identify before any use.** |
| **ClinVar** (annotation) | release of 31 August 2025: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2025/clinvar_20250831.vcf.gz` | Used for benign / pathogenic / VUS classification. |
| **COSMIC** (annotation) | `https://cancer.sanger.ac.uk/cosmic/download` (v102, GRCh37: `COSMIC_MUTANTS_all_v102_GRCh37.vcf.gz`; free registration) | Used for the "database-matched" bucket. |
| **GENCODE** (gene model) | `https://www.gencodegenes.org/human/release_19.html` (release 19, the GRCh37 gene model: `gencode.v19.annotation.gtf.gz`) | For SNV -> gene mapping. |

---

## 2. Environment

Variant parsing/annotation uses `pysam`, which runs under **WSL/Linux**, not native
Windows:

```bash
conda activate py3        # same env as the analysis (py3-environment.yml) + pysam
```

---

## 3. Rebuild steps

The original (non-shipped) scripts are named in brackets so they can be matched to the
authors' working copy.

1. **Acquire and place raw data** as in the table above:
   `*.hard-filtered.vcf.gz` (AA/EA), `*-full.tsv` (IN, de-identified), `GSE244631/`
   (RCC), and the ClinVar / COSMIC / GENCODE databases.

2. **Parse variants to per-patient SNV lists.** Read the DRAGEN VCFs and the CHPV2
   TSVs into a common `(chrom, pos, ref, alt)` SNV representation per patient, keeping
   high-confidence calls. *(`read_vcf.py` / `extract_somatic.py`, plus the CHPV2 TSV
   parser.)*

3. **Annotate each SNV** against ClinVar and COSMIC and assign the three mutually
   exclusive buckets used throughout the paper — *benign / likely-benign* (ClinVar
   Benign/Likely-benign), *database-matched non-benign* (any remaining SNV catalogued
   in COSMIC or ClinVar), and *unmatched* (in neither). *(`annotate_master.py`, driven
   over all samples by `annotate_all.py`.)*

4. **Build the SNV master matrices** (patient x SNV, per bucket and per
   population/severity group). *(`rebuild_snv_matrix.py` -> `snv_rebuild/`.)*

5. **Map SNV -> gene (GENCODE) and aggregate to per-patient gene vectors**
   (`gene, variant_count`), one CSV per patient. *(`build_matched_query.py` ->
   `cohort_matched/`.)* This directory is the input the analysis scripts read.

6. **Reference labels.** Apply the study's Gleason high/low rule to the GSE244631
   sample metadata to produce `GSE244631_ref_labels.csv` (`sample_id,class`). The RCC
   files are parsed directly by `distance_knn.py` at analysis time.

7. **Cohort labels.** `genefilelists/metadata_orig.csv` holds `ID,Gleason` (+ ancestry)
   for the 20 patients. Supply a de-identified copy alongside `cohort_matched/`.

---

## 4. Result

After step 5 you have `cohort_matched/` and (step 6) `GSE244631_ref_labels.csv`. With
those in place, every script in the parent folder runs as described in `README.md`.
