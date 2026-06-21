# cohort_matched/ — data slot (not in repository)

Drop the **20 per-patient gene-count CSVs** here, one file per patient named `<ID>.csv`
with columns `gene,variant_count` (e.g. `PT028.csv`, `50702111163.csv`).

- **Source:** derived from the cohort's tumor sequencing; **shared separately under a
  data-use agreement** (private — not committed to this repo).
- **How it is produced:** see `../UPSTREAM.md`.
- **Read by:** `distance_knn.load_query_profiles()` (used by `panel_recovery*.py`,
  `viz_*.py`, `quick_norm_sweep.py`, `diag_aaea.py`, `perm_test.py`).

Everything in this folder except this file is git-ignored.
