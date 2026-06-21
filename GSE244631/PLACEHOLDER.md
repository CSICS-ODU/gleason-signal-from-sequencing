# GSE244631/ — data slot (not in repository)

Drop the **NanoString BC360 `.RCC` reference files** here.

- **Source (public):** NCBI GEO accession **GSE244631** —
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244631 (download the
  supplementary `.RCC` files). 51 labeled prostate samples are used (19 high / 32 low).
- **Read by:** `distance_knn.load_reference(ref_rcc_dir="GSE244631")` (used by all the
  cross-modality scripts). The high/low labels are in `../GSE244631_ref_labels.csv`
  (included in the repo).

Everything in this folder except this file is git-ignored.
