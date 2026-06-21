# Deprecated / superseded code

These scripts were used during the project but are **not** part of the current
analysis. They are kept for provenance and history. None of them are needed to
reproduce the manuscript results; the current equivalents live in the parent folder.

Do not build on these without checking against the current scripts first.

| Script | What it did | Why deprecated |
|--------|-------------|----------------|
| `nested_loo.py` | Nested leave-one-out severity classifier with per-fold biomarker reselection; produced the early **Spearman 0.89** recovery. | That result was an assay/ancestry confound (the reselected genes tracked the Indian targeted panel, not severity). Replaced by `panel_recovery.py` + `panel_recovery_aaea.py`. **Dead — do not resurrect.** |
| `make_final_fig.py` | Two-panel self-signature-z vs per-gene-z recovery figure (`fig_recovery_pergene`). | Replaced by the single grade-association-score figure from `viz_gradescore.py`. |
| `make_recovery_fig.py` | Grouped-bar genome-wide-vs-panel recovery figure (`fig_recovery`). | Same; superseded by `viz_gradescore.py`. |
| `probe_biomarker.py` | One-off diagnostic: provenance and per-gene carrier ancestry of the candidate panel (showed the panel ≈ the Indian clinical-panel genes). | Diagnostic only; conclusions folded into the manuscript. |
| `compare_signatures.py` | Early in-house leave-one-out severity classifier on gene profiles. | Superseded by `distance_knn.py` / `panel_recovery*.py`. |
| `snv_discrimination.py` | Early SNV-level high/low discrimination. | Superseded by `disproportion.py`. |
| `group_analysis.py`, `group_gene_analysis.py` | Early per-group (population x severity) variant / gene tallies. | Superseded by `build_paper_results.py` and `disproportion*.py`. |
| `manuscript_results.py` | Early paste-ready result tables. | Superseded by `build_paper_results.py`. |
| `population_pattern.py` | Per-population AA-vs-EA high/low lean for focus genes. | Exploratory support for the population-divergence result; not in the final pipeline. |
| `cohort_coverage_analysis.py` | Why the Indian cohort contributes far fewer variants (targeted panel vs broad sequencing). | Background analysis; conclusion stated in the manuscript, script not rerun. |
| `cooccurrence.py` | Exploratory gene-gene co-occurrence ("Layer 3"). | Never used in the manuscript. |
| `sweep_recurrence.py` | Recurrence-threshold sweep (gene-level). | Exploratory. |
| `query_genes.py` | Utility to query the rebuilt gene x patient matrices. | Ad hoc utility. |
| `distance.py`, `copy_of_distance.py` | Original prototype distance / classification (HIGH/LOW prototype averaging), from `genefilelists/`. | Replaced by `distance_knn.py` (which drops the prototype averaging). |
| `b.py` | Scratch comparison script. | Throwaway. |
