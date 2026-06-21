# genefilelists/ — data slot (not in repository)

Drop the **de-identified cohort label file** `metadata_orig.csv` here, with columns
`ID,Gleason` (and an ancestry/`race` column).

- **Source:** the cohort metadata; **shared separately** (private — contains cohort
  identifiers, not committed to this repo).
- **Read by:** `distance_knn.load_query_truth()` and `load_query_races()` (used by
  `panel_recovery*.py`, `viz_*.py`, `quick_norm_sweep.py`, `diag_aaea.py`,
  `perm_test.py`).

Everything in this folder except this file is git-ignored.
