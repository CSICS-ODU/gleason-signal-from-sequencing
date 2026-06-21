# -*- coding: utf-8 -*-
"""
distance_knn.py  --  classifier-viability experiment, averaging removed.

WHAT THIS IS
------------
A rebuild of distance.py that DROPS the HIGH/LOW prototype averaging and instead
maps each of our cohort patients directly against EVERY individual reference
sample (GSE244631, 51 NanoString BC360 expression samples), then classifies.

The question is whether the "invariant" idea survives a rigorous method:
is there a signal shared across two different MEASUREMENTS
    reference = mRNA abundance (NanoString panel expression)
    query     = DNA variant burden (variant_count per gene, from VCFs)
such that a high-grade / high-risk tumor "looks alike" across both.

TRAIN / TEST DISCIPLINE (this is a classifier-viability test, so it matters)
---------------------------------------------------------------------------
  * Reference = TRAINING set. Per-gene-across-sample statistics are allowed
    (z-score each gene over the 51 reference samples). Fit ONCE on the reference.
  * Query (our cohort) = TEST set. Each patient is normalized using ONLY itself
    (+ fixed annotations like gene length). NEVER the other query patients --
    that would leak information you would not have when classifying one new
    patient in deployment.

The asymmetry (reference normalized per-gene, query normalized per-patient) is
the honest consequence of "no test-time leakage", not an oversight. Because the
two sides end up standardized on different axes, prefer scale-tolerant metrics
(cosine, spearman).

This script SWEEPS the open choices so you can explore:
    query normalization : within_patient_z | within_patient_rank | none
    length correction   : on | off        (variant_count -> variants per kb)
    reference norm       : zscore_per_gene | none
    metric              : cosine | spearman | euclidean | canberra
    classifier          : knn (k=1,3,5...) | classmean
and reports accuracy against the Gleason high/low truth, alongside a
majority-class baseline (so "no better than base rate" is impossible to miss).

Optional within-modality positive controls (leave-one-out) live behind --control.
They are OFF by default so the main run is just the cross-modality test.

INPUTS (point the CONFIG / CLI flags at these once the data is pulled)
---------------------------------------------------------------------
  --ref-rcc-dir DIR      directory of GSE244631 *.RCC files (parsed here), OR
  --ref-matrix FILE      a genes x samples CSV (auto-detects orientation)
  --ref-labels FILE      CSV with columns: sample_id,class  (class in {high,low};
                         map recurrent / CAPRA-S high-risk -> high when you build it)
  --query-dir DIR        per-patient CSVs (gene,variant_count), filename = patient id
  --query-meta FILE      metadata CSV with an ID column + a Gleason column
  --gene-lengths FILE    optional CSV: gene,length_bp  (required if --length-correct on)

Nothing here is run for you; build/edit only -- you run and test it.
"""

import argparse
import glob
import os
import re

import numpy as np
import pandas as pd
from scipy.spatial.distance import canberra, cosine, euclidean
from scipy.stats import spearmanr, rankdata


# ============================================================================
# ID + label helpers
# ============================================================================

def normalize_id(x):
    """Canonicalize an ID or filename to 'PT<n>' or a bare digit string."""
    if x is None:
        return None
    s = os.path.basename(str(x).strip())
    s = re.sub(r"\.csv$", "", s, flags=re.IGNORECASE)
    if re.fullmatch(r"\d+", s):
        return s
    if (m := re.fullmatch(r"PT[_-]?(\d+)", s, flags=re.IGNORECASE)):
        return f"PT{m.group(1)}"
    if (m := re.search(r"PT[_-]?(\d+)", s, flags=re.IGNORECASE)):
        return f"PT{m.group(1)}"  # matches embedded forms e.g. S01_PT_028_RNA
    if (m := re.search(r"sample_PT[_-]?(\d+)(?:_|$)", s, flags=re.IGNORECASE)):
        return f"PT{m.group(1)}"
    if (m := re.search(r"sample_(\d+)(?:_|$)", s, flags=re.IGNORECASE)):
        return m.group(1)
    return None


def classify_gleason(gleason_value):
    """Raw Gleason -> 'high' / 'low'. 4+3 high, 3+4 low, >7 high else low."""
    s = str(gleason_value).strip()
    if "4+3" in s:
        return "high"
    if "3+4" in s:
        return "low"
    try:
        return "high" if int(s) > 7 else "low"
    except ValueError:
        raise ValueError(f"Unrecognized Gleason value: {gleason_value!r}")


# ============================================================================
# Reference loading  (GSE244631 NanoString expression)
# ============================================================================

def parse_rcc_file(path):
    """Parse one NanoString .RCC file -> Series {gene_name: count}.

    RCC files are plain text with a <Code_Summary> section whose rows are
    'CodeClass,Name,Accession,Count'. We keep Endogenous + Housekeeping genes.
    """
    counts = {}
    in_summary = False
    with open(path, "r", encoding="latin1") as fh:
        for line in fh:
            line = line.strip()
            if line == "<Code_Summary>":
                in_summary = True
                continue
            if line == "</Code_Summary>":
                break
            if not in_summary or not line:
                continue
            parts = line.split(",")
            if len(parts) < 4 or parts[0] == "CodeClass":
                continue
            code_class, name = parts[0], parts[1]
            if code_class not in ("Endogenous", "Housekeeping"):
                continue  # drop Positive/Negative spike-in controls
            try:
                counts[name] = float(parts[3])
            except ValueError:
                continue
    return pd.Series(counts, name=normalize_id(path) or os.path.basename(path))


def load_reference_from_rcc_dir(rcc_dir):
    """Build a genes x samples DataFrame from a directory of .RCC files."""
    paths = sorted(glob.glob(os.path.join(rcc_dir, "*.RCC")) +
                   glob.glob(os.path.join(rcc_dir, "*.rcc")))
    if not paths:
        raise FileNotFoundError(f"No .RCC files found in {rcc_dir}")
    cols = {}
    for p in paths:
        s = parse_rcc_file(p)
        base = os.path.basename(p)
        m = re.search(r"GSM\d+", base)          # key by GSM accession so it
        key = m.group(0) if m else os.path.splitext(base)[0]  # joins to ref-labels
        cols[key] = s
    mat = pd.DataFrame(cols)
    return mat


def load_reference_matrix(path):
    """Load a genes x samples CSV. Auto-detect orientation by assuming the
    axis with more entries is genes (panels have ~770 genes, ~51 samples)."""
    df = pd.read_csv(path, index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")
    if df.shape[0] < df.shape[1]:
        # more columns than rows -> samples are rows; transpose to genes x samples
        df = df.T
    return df


def load_reference(args):
    if args.ref_rcc_dir:
        mat = load_reference_from_rcc_dir(args.ref_rcc_dir)
    elif args.ref_matrix:
        mat = load_reference_matrix(args.ref_matrix)
    else:
        raise ValueError("Provide --ref-rcc-dir or --ref-matrix")
    # collapse duplicate genes by mean
    mat = mat.groupby(mat.index).mean()
    return mat  # genes x samples


def load_ref_labels(path):
    """CSV with columns sample_id,class (class in {high,low})."""
    df = pd.read_csv(path)
    df.columns = [c.strip().lower() for c in df.columns]
    if "sample_id" not in df.columns or "class" not in df.columns:
        raise ValueError("--ref-labels must have columns: sample_id,class")
    labels = {}
    for _, r in df.iterrows():
        cls = str(r["class"]).strip().lower()
        if cls not in ("high", "low"):
            raise ValueError(f"ref label class must be high/low, got {cls!r}")
        labels[str(r["sample_id"]).strip()] = cls
    return labels


# ============================================================================
# Query loading  (our cohort variant burden)
# ============================================================================

def load_query_profiles(query_dir):
    """dir of per-patient CSVs (gene,variant_count) -> {patient_id: Series}."""
    profiles = {}
    for path in sorted(glob.glob(os.path.join(query_dir, "*.csv"))):
        pid = normalize_id(path)
        if pid is None:
            continue
        df = pd.read_csv(path)
        df.columns = [c.strip() for c in df.columns]
        if "gene" not in df.columns:
            continue
        val_col = ("variant_count" if "variant_count" in df.columns
                   else "average_expression" if "average_expression" in df.columns
                   else df.columns[1])
        s = pd.Series(pd.to_numeric(df[val_col], errors="coerce").values,
                      index=df["gene"].astype(str).values).dropna()
        profiles[pid] = s.groupby(level=0).mean()
    return profiles


def load_query_truth(meta_path, id_col="ID", gleason_col="Gleason"):
    df = pd.read_csv(meta_path, skipinitialspace=True)
    df.columns = df.columns.str.strip()
    truth = {}
    for _, r in df.iterrows():
        pid = normalize_id(r[id_col])
        if pid is None:
            continue
        truth[pid] = classify_gleason(r[gleason_col])
    return truth


def load_query_races(meta_path, id_col="ID", race_col="race"):
    """{patient_id: race} (e.g. IN/AA/EA), for ancestry filtering."""
    df = pd.read_csv(meta_path, skipinitialspace=True)
    df.columns = df.columns.str.strip()
    if race_col not in df.columns:
        return {}
    races = {}
    for _, r in df.iterrows():
        pid = normalize_id(r[id_col])
        if pid is not None:
            races[pid] = str(r[race_col]).strip().upper()
    return races


def load_gene_lengths(path):
    """CSV gene,length_bp -> Series gene->length_kb."""
    df = pd.read_csv(path)
    df.columns = [c.strip().lower() for c in df.columns]
    gcol = "gene" if "gene" in df.columns else df.columns[0]
    lcol = ("length_bp" if "length_bp" in df.columns
            else "length" if "length" in df.columns else df.columns[1])
    s = pd.Series(pd.to_numeric(df[lcol], errors="coerce").values,
                  index=df[gcol].astype(str).values).dropna()
    s = s.groupby(level=0).mean() / 1000.0  # bp -> kb
    return s[s > 0]


# ============================================================================
# Normalization
# ============================================================================

def zscore_per_gene(mat):
    """Per-gene z-score ACROSS samples (rows=genes, cols=samples).
    Allowed for the reference/training set only. Genes with zero variance -> 0."""
    mu = mat.mean(axis=1)
    sd = mat.std(axis=1, ddof=0).replace(0, np.nan)
    z = mat.sub(mu, axis=0).div(sd, axis=0)
    return z.fillna(0.0)


def normalize_query_vector(vec, method):
    """Leakage-free per-PATIENT normalization. Uses only this one vector.

    method:
      within_patient_z    : (x - mean_genes) / std_genes  for this patient
      within_patient_rank : rank of each gene within this patient
      none                : raw values
    """
    v = np.asarray(vec, dtype=float)
    if method == "none":
        return v
    if method == "within_patient_z":
        sd = v.std(ddof=0)
        if sd == 0:
            return np.zeros_like(v)
        return (v - v.mean()) / sd
    if method == "within_patient_rank":
        return rankdata(v)
    raise ValueError(f"Unknown query-norm method {method!r}")


def within_sample_normalize_matrix(mat, kind):
    """Normalize each SAMPLE (column) across its own genes -- the SAME self-signature
    transform applied to the query, so reference and query live on the same axis.
    mat is genes x samples; kind is 'z' or 'rank'. Vectorized over columns (numerically
    identical to per-column normalize_query_vector, but ~30x faster in tight loops)."""
    arr = mat.to_numpy(dtype=float)
    if kind == "z":
        mu = arr.mean(axis=0, keepdims=True)
        sd = arr.std(axis=0, ddof=0, keepdims=True)
        out = np.where(sd == 0.0, 0.0, (arr - mu) / np.where(sd == 0.0, 1.0, sd))
    else:  # rank
        out = np.apply_along_axis(rankdata, 0, arr)
    return pd.DataFrame(out, index=mat.index, columns=mat.columns)


# ============================================================================
# Distance + classification
# ============================================================================

def compute_dist(a, b, metric):
    if metric == "canberra":
        return canberra(a, b)
    if metric == "cosine":
        return cosine(a, b)
    if metric == "euclidean":
        return euclidean(a, b)
    if metric == "spearman":
        coef, _ = spearmanr(a, b)
        return 1.0 - (coef if np.isfinite(coef) else 0.0)
    raise ValueError(f"Metric {metric!r} not recognized")


def knn_predict(dists, ref_labels_ordered, k):
    """dists: 1D distances to each reference sample (same order as labels).
    Majority label of the k nearest. Ties -> 'high' is not assumed; falls back
    to the single nearest on a tie."""
    order = np.argsort(dists)
    top = order[:k]
    votes = [ref_labels_ordered[i] for i in top]
    n_high = votes.count("high")
    n_low = votes.count("low")
    if n_high == n_low:
        return ref_labels_ordered[order[0]]  # nearest breaks the tie
    return "high" if n_high > n_low else "low"


def classmean_predict(dists, ref_labels_ordered):
    """Mean distance to all HIGH refs vs all LOW refs; closer class wins.
    This is the old prototype idea WITHOUT pre-averaging the profiles."""
    d = np.asarray(dists, dtype=float)
    lab = np.asarray(ref_labels_ordered)
    dh = d[lab == "high"].mean() if (lab == "high").any() else np.inf
    dl = d[lab == "low"].mean() if (lab == "low").any() else np.inf
    return "high" if dh < dl else "low", dh, dl


# --- class-balanced voting -------------------------------------------------
# Raw kNN votes regress to the REFERENCE prior (19 high / 32 low = low-majority),
# which is inverted vs the query (high-majority) -> below-random accuracy.
# Two fixes: (a) inverse-frequency 'weight' the votes; (b) 'subsample' the
# majority class down to the minority count, repeated over rounds so every
# low sample is used across rounds, then majority-vote the rounds.

def _knn_vote(dists, lab, idx, k, weighted):
    sub = idx[np.argsort(dists[idx])[:k]]
    votes = lab[sub]
    if weighted:
        ch = max((lab[idx] == "high").sum(), 1)
        cl = max((lab[idx] == "low").sum(), 1)
        wh, wl = (votes == "high").sum() / ch, (votes == "low").sum() / cl
    else:
        wh, wl = (votes == "high").sum(), (votes == "low").sum()
    if wh == wl:
        return lab[sub[0]]
    return "high" if wh > wl else "low"


def _classmean_vote(dists, lab, idx):
    dh = dists[idx][lab[idx] == "high"]
    dl = dists[idx][lab[idx] == "low"]
    mh = dh.mean() if len(dh) else np.inf
    ml = dl.mean() if len(dl) else np.inf
    return ("high" if mh < ml else "low"), mh, ml


def predict_balanced(dists, lab, classifier, k, mode, rounds, rng):
    """Return (prediction, (d_high_mean, d_low_mean)) under a balancing mode."""
    full = np.arange(len(lab))

    def one(idx, weighted):
        if classifier == "knn":
            return _knn_vote(dists, lab, idx, k, weighted), (np.nan, np.nan)
        p, mh, ml = _classmean_vote(dists, lab, idx)
        return p, (mh, ml)

    if mode == "none":
        return one(full, False)
    if mode == "weight":
        return one(full, True)
    # subsample: match majority down to minority, repeat, majority-vote rounds
    hi = np.where(lab == "high")[0]
    lo = np.where(lab == "low")[0]
    m = min(len(hi), len(lo))
    preds, mhs, mls = [], [], []
    for _ in range(rounds):
        idx = np.concatenate([rng.choice(hi, m, replace=False),
                              rng.choice(lo, m, replace=False)])
        p, (mh, ml) = one(idx, False)
        preds.append(p); mhs.append(mh); mls.append(ml)
    nh = preds.count("high")
    pred = "high" if nh >= len(preds) - nh else "low"
    return pred, (float(np.nanmean(mhs)), float(np.nanmean(mls)))


# ============================================================================
# Build aligned, normalized matrices
# ============================================================================

def build_matrices(ref_mat, ref_labels, query_profiles, args, gene_lengths=None):
    """Return reference matrix (genes x samples, normalized), ordered ref labels,
    and a {patient: normalized vector} dict -- all on a shared gene universe."""
    # keep only labeled reference samples
    ref_samples = [s for s in ref_mat.columns if s in ref_labels]
    if not ref_samples:
        raise ValueError("No reference samples matched --ref-labels sample_id values. "
                          f"ref cols e.g. {list(ref_mat.columns)[:3]}, "
                          f"label ids e.g. {list(ref_labels)[:3]}")
    ref_mat = ref_mat[ref_samples]

    # shared gene universe = reference genes  INTERSECT  union of query genes
    query_genes = set().union(*[set(s.index) for s in query_profiles.values()])
    universe = sorted(set(ref_mat.index) & query_genes)
    if not universe:
        raise ValueError("Empty gene universe: reference and query share no gene names. "
                          "Check gene symbol nomenclature on both sides.")

    # optional restriction to a hypothesis-driven gene subset (e.g. our biomarkers).
    # NOTE: if the subset was selected using the QUERY high/low labels, classifying
    # the same cohort with it is feature-selection leakage -> optimistic. Read it as
    # a biology probe, not a clean held-out viability number.
    subset = getattr(args, "gene_subset", None)
    if subset:
        want = set(pd.read_csv(subset).iloc[:, 0].astype(str))
        universe = sorted(set(universe) & want)
        if not universe:
            raise ValueError(f"Gene subset {subset} shares no genes with the universe.")
        print(f"  gene-subset {subset}: restricted universe -> {len(universe)} genes")

    ref_mat = ref_mat.reindex(universe).fillna(0.0)

    # reference normalization.
    #   within_sample_z/rank : each reference sample self-normalized across its own
    #       genes -- the SAME transform as the query (the intended "signature" method).
    #   zscore_per_gene      : per-gene across samples (training stats; different axis).
    if args.ref_norm == "within_sample_z":
        ref_norm = within_sample_normalize_matrix(ref_mat, "z")
    elif args.ref_norm == "within_sample_rank":
        ref_norm = within_sample_normalize_matrix(ref_mat, "rank")
    elif args.ref_norm == "zscore_per_gene":
        ref_norm = zscore_per_gene(ref_mat)
    elif args.ref_norm == "none":
        ref_norm = ref_mat
    else:
        raise ValueError(f"Unknown ref-norm {args.ref_norm!r}")

    ref_labels_ordered = [ref_labels[s] for s in ref_samples]

    # query: length-correct (stateless) then per-patient normalize (leakage-free)
    if args.length_correct and gene_lengths is None:
        raise ValueError("--length-correct on requires --gene-lengths")

    query_vecs = {}
    for pid, prof in query_profiles.items():
        vec = prof.reindex(universe).fillna(0.0)
        if args.length_correct:
            lens = gene_lengths.reindex(universe)
            vec = vec / lens.replace(0, np.nan)
            vec = vec.fillna(0.0)  # genes without length -> 0 contribution
        nv = normalize_query_vector(vec.values, args.query_norm)
        query_vecs[pid] = nv

    return ref_norm, ref_labels_ordered, query_vecs, universe


# ============================================================================
# Main experiment sweep
# ============================================================================

def run_once(ref_norm, ref_labels_ordered, query_vecs, query_truth,
             metric, classifier, k, balance="none", rounds=50, seed=0):
    """One (metric, classifier, k) configuration -> per-patient rows."""
    ref_arr = ref_norm.values.T  # samples x genes
    lab = np.asarray(ref_labels_ordered)
    rng = np.random.default_rng(seed)
    rows = []
    for pid, qv in query_vecs.items():
        truth = query_truth.get(pid)
        dists = np.array([compute_dist(qv, ref_arr[j], metric)
                          for j in range(ref_arr.shape[0])])
        pred, (dh, dl) = predict_balanced(dists, lab, classifier, k,
                                          balance, rounds, rng)
        extra = {"k": k} if classifier == "knn" else {"d_high_mean": dh,
                                                       "d_low_mean": dl}
        rows.append({
            "sample_id": pid, "metric": metric, "classifier": classifier,
            "balance": balance,
            "true_class": truth, "predicted_class": pred,
            "correct": (None if truth is None else pred == truth),
            **extra,
        })
    return rows


def majority_baseline(query_truth):
    vals = [v for v in query_truth.values() if v in ("high", "low")]
    if not vals:
        return None, None
    n_high = vals.count("high")
    maj = "high" if n_high >= len(vals) - n_high else "low"
    return maj, max(n_high, len(vals) - n_high) / len(vals)


def balanced_accuracy(sub):
    """Mean of per-class recall -- neutralises the reference/query prior mismatch.
    0.5 == chance regardless of class imbalance."""
    recs = []
    for cls in ("high", "low"):
        m = sub["true_class"] == cls
        if m.any():
            recs.append((sub.loc[m, "predicted_class"] == cls).mean())
    return float(np.mean(recs)) if recs else float("nan")


def permutation_balacc(sub, n_perm=5000, seed=0):
    """Permutation null for balanced accuracy: keep the model's predictions fixed,
    shuffle the TRUE labels n_perm times, recompute balanced accuracy each time.
    Returns (observed, p_value, null_2.5pct, null_97.5pct), where p is the upper-
    tail fraction of permutations >= observed (i.e. P(chance does this well or
    better)). This is the defensible test for a negative classification result:
    a value is 'no better than chance' iff it falls inside this null."""
    s = sub[sub["true_class"].isin(["high", "low"])]
    t = s["true_class"].to_numpy()
    p = s["predicted_class"].to_numpy()
    if t.size == 0:
        return float("nan"), float("nan"), float("nan"), float("nan")

    def bal(tt):
        recs = []
        for cls in ("high", "low"):
            m = tt == cls
            if m.any():
                recs.append((p[m] == cls).mean())
        return float(np.mean(recs)) if recs else float("nan")

    obs = bal(t)
    rng = np.random.default_rng(seed)
    null = np.array([bal(rng.permutation(t)) for _ in range(n_perm)])
    p_val = (np.sum(null >= obs - 1e-12) + 1) / (n_perm + 1)
    lo, hi = float(np.percentile(null, 2.5)), float(np.percentile(null, 97.5))
    return obs, float(p_val), lo, hi


def summarize(all_rows, query_truth, n_perm=5000, seed=0):
    df = pd.DataFrame(all_rows)
    df["correct"] = df["correct"].astype("boolean")
    grp_cols = [c for c in ["balance", "metric", "classifier", "k"] if c in df.columns]
    rows = []
    for key, sub in df.groupby(grp_cols, dropna=False):
        rec = dict(zip(grp_cols, key if isinstance(key, tuple) else (key,)))
        rec["n"] = int(sub["correct"].notna().sum())
        rec["accuracy"] = float(sub["correct"].fillna(False).mean())
        rec["balanced_acc"] = balanced_accuracy(sub)
        _, rec["perm_p"], rec["null_lo"], rec["null_hi"] = permutation_balacc(
            sub, n_perm=n_perm, seed=seed)
        rows.append(rec)
    summary = pd.DataFrame(rows).sort_values("balanced_acc", ascending=False)
    maj, base = majority_baseline(query_truth)
    summary["majority_baseline"] = base
    summary["beats_balanced_chance"] = summary["balanced_acc"] > 0.5
    return df, summary


# --------------------------- LOO within-modality controls --------------------

def loo_control(mat, labels_ordered, sample_names, metric, k, norm="zscore_per_gene"):
    """Leave-one-out kNN within a single modality (positive control).
    mat: genes x samples (raw). Refits per-gene z on the n-1 remaining each fold."""
    cols = list(sample_names)
    rows = []
    for i, held in enumerate(cols):
        keep = [c for c in cols if c != held]
        sub = mat[keep]
        if norm == "zscore_per_gene":
            mu = sub.mean(axis=1)
            sd = sub.std(axis=1, ddof=0).replace(0, np.nan)
            sub_n = sub.sub(mu, axis=0).div(sd, axis=0).fillna(0.0)
            held_v = ((mat[held] - mu) / sd).fillna(0.0).values
        else:
            sub_n = sub
            held_v = mat[held].values
        keep_labels = [labels_ordered[cols.index(c)] for c in keep]
        dists = np.array([compute_dist(held_v, sub_n[c].values, metric) for c in keep])
        pred = knn_predict(dists, keep_labels, k)
        rows.append({"held": held, "true": labels_ordered[i],
                     "pred": pred, "correct": pred == labels_ordered[i]})
    acc = np.mean([r["correct"] for r in rows]) if rows else float("nan")
    return acc, rows


# ============================================================================
# Manuscript tables: all-genes (failure) vs biomarker-subset (recovery)
# ============================================================================

def write_paper_summary(ref_mat, ref_labels, query_profiles, query_truth, args,
                        gene_lengths):
    """Run the metric x classifier sweep twice -- once on the full gene universe
    and once restricted to the --paper-tables biomarker gene list -- and write a
    single tidy CSV (knn_paper_summary.csv) tagged by gene_set. Feeds the
    manuscript's Table 1 (gene_set='all', the failure) and the Section VI.B
    mirror table (gene_set='biomarker', the recovery) from one consistent run."""
    metrics = [m.strip() for m in args.metrics.split(",") if m.strip()]
    ks = [int(x) for x in args.k.split(",") if x.strip()]
    classifiers = (["knn", "classmean"] if args.classifier == "both"
                   else [args.classifier])
    parts = []
    for gene_set, subset in (("all", None), ("biomarker", args.paper_tables)):
        args.gene_subset = subset
        ref_norm, ref_labels_ordered, query_vecs, universe = build_matrices(
            ref_mat, ref_labels, query_profiles, args, gene_lengths)
        rows = []
        for metric in metrics:
            for clf in classifiers:
                for k in (ks if clf == "knn" else [None]):
                    rows.extend(run_once(ref_norm, ref_labels_ordered, query_vecs,
                                         query_truth, metric, clf, k,
                                         balance=args.balance_ref,
                                         rounds=args.balance_rounds, seed=args.seed))
        _, summary = summarize(rows, query_truth)
        summary.insert(0, "gene_set", gene_set)
        summary.insert(1, "n_genes", len(universe))
        parts.append(summary)
    combined = pd.concat(parts, ignore_index=True)
    combined.to_csv("knn_paper_summary.csv", index=False)
    print(f"\nwrote knn_paper_summary.csv "
          f"({len(combined)} rows: all-genes vs biomarker subset)")
    return combined


# ============================================================================
# CLI
# ============================================================================

def get_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    # reference
    p.add_argument("--ref-rcc-dir", help="dir of GSE244631 *.RCC files")
    p.add_argument("--ref-matrix", help="genes x samples CSV (alternative to --ref-rcc-dir)")
    p.add_argument("--ref-labels", required=True, help="CSV: sample_id,class (high/low)")
    # query
    p.add_argument("--query-dir", required=True, help="dir of per-patient gene,variant_count CSVs")
    p.add_argument("--query-meta", required=True, help="metadata CSV with ID + Gleason")
    p.add_argument("--id-col", default="ID")
    p.add_argument("--gleason-col", default="Gleason")
    p.add_argument("--race-col", default="race")
    p.add_argument("--query-races", default="",
                   help="restrict query to these ancestries, e.g. 'AA,EA' "
                        "(drops the IN targeted-panel cohort; empty = all)")
    # normalization knobs
    p.add_argument("--ref-norm", default="within_sample_z",
                   choices=["within_sample_z", "within_sample_rank",
                            "zscore_per_gene", "none"])
    p.add_argument("--query-norm", default="within_patient_z",
                   choices=["within_patient_z", "within_patient_rank", "none"])
    p.add_argument("--length-correct", action="store_true",
                   help="convert variant_count -> variants/kb (needs --gene-lengths)")
    p.add_argument("--gene-lengths", help="CSV: gene,length_bp")
    p.add_argument("--gene-subset", help="CSV whose first column is a gene list; "
                   "restrict the distance to these genes (e.g. biomarker_subset.csv)")
    # sweep
    p.add_argument("--metrics", default="cosine,spearman,euclidean,canberra")
    p.add_argument("--k", default="1,3,5", help="comma list of k for kNN")
    p.add_argument("--classifier", default="both", choices=["knn", "classmean", "both"])
    p.add_argument("--balance-ref", default="subsample",
                   choices=["none", "subsample", "weight"],
                   help="balance the imbalanced reference vote (19 high/32 low): "
                        "'subsample' matches majority->minority over rounds; "
                        "'weight' uses inverse-frequency votes")
    p.add_argument("--balance-rounds", type=int, default=200,
                   help="subsampling rounds (subsample mode)")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--sweep-norms", action="store_true",
                   help="also sweep every query-norm x length-correct combo")
    # controls
    p.add_argument("--control", action="store_true",
                   help="run leave-one-out within-modality positive controls")
    # output
    p.add_argument("--out", default="knn_results.csv")
    p.add_argument("--summary-out", default="knn_summary.csv")
    p.add_argument("--paper-tables",
                   help="biomarker gene-list CSV; run the sweep on all genes and "
                        "on this subset and write knn_paper_summary.csv (tagged by "
                        "gene_set) for the manuscript's failure and recovery tables")
    return p.parse_args()


def main():
    args = get_args()

    ref_mat = load_reference(args)
    ref_labels = load_ref_labels(args.ref_labels)
    query_profiles = load_query_profiles(args.query_dir)
    query_truth = load_query_truth(args.query_meta, args.id_col, args.gleason_col)
    gene_lengths = load_gene_lengths(args.gene_lengths) if args.gene_lengths else None

    # optional ancestry restriction (e.g. AA,EA -> drop the IN targeted-panel cohort)
    if args.query_races.strip():
        keep = {r.strip().upper() for r in args.query_races.split(",") if r.strip()}
        races = load_query_races(args.query_meta, args.id_col, args.race_col)
        query_profiles = {p: v for p, v in query_profiles.items()
                          if races.get(p) in keep}
        query_truth = {p: v for p, v in query_truth.items() if races.get(p) in keep}
        print(f"ancestry filter {sorted(keep)}: kept {len(query_profiles)} query patients")

    print(f"reference: {ref_mat.shape[0]} genes x {ref_mat.shape[1]} samples "
          f"({len(ref_labels)} labeled)")
    n_hi = sum(v == "high" for v in query_truth.values())
    print(f"query: {len(query_profiles)} patients, {len(query_truth)} with truth "
          f"({n_hi} high / {len(query_truth)-n_hi} low)")

    metrics = [m.strip() for m in args.metrics.split(",") if m.strip()]
    ks = [int(x) for x in args.k.split(",") if x.strip()]
    classifiers = (["knn", "classmean"] if args.classifier == "both"
                   else [args.classifier])

    # which normalization combos to explore
    if args.sweep_norms:
        norm_combos = [(qn, lc)
                       for qn in ("within_patient_z", "within_patient_rank", "none")
                       for lc in (False, True)]
    else:
        norm_combos = [(args.query_norm, args.length_correct)]

    all_rows = []
    for qn, lc in norm_combos:
        if lc and gene_lengths is None:
            print(f"  [skip] query-norm={qn} length-correct=on : no --gene-lengths")
            continue
        args.query_norm, args.length_correct = qn, lc
        ref_norm, ref_labels_ordered, query_vecs, universe = build_matrices(
            ref_mat, ref_labels, query_profiles, args, gene_lengths)
        tag = f"{qn}|len={'on' if lc else 'off'}"
        print(f"  config {tag}: universe = {len(universe)} shared genes")
        for metric in metrics:
            for clf in classifiers:
                kk = ks if clf == "knn" else [None]
                for k in kk:
                    rows = run_once(ref_norm, ref_labels_ordered, query_vecs,
                                    query_truth, metric, clf, k,
                                    balance=args.balance_ref,
                                    rounds=args.balance_rounds, seed=args.seed)
                    for r in rows:
                        r["norm"] = tag
                    all_rows.extend(rows)

    results_df = pd.DataFrame(all_rows)
    full_df, summary = summarize(all_rows, query_truth)

    # group summary additionally by norm config when sweeping
    if args.sweep_norms:
        grp = [c for c in ["norm", "balance", "metric", "classifier", "k"]
               if c in full_df.columns]
        rows = []
        for key, sub in full_df.groupby(grp, dropna=False):
            rec = dict(zip(grp, key if isinstance(key, tuple) else (key,)))
            rec["n"] = int(sub["correct"].notna().sum())
            rec["accuracy"] = float(sub["correct"].fillna(False).mean())
            rec["balanced_acc"] = balanced_accuracy(sub)
            rows.append(rec)
        summary = pd.DataFrame(rows).sort_values("balanced_acc", ascending=False)
        maj, base = majority_baseline(query_truth)
        summary["majority_baseline"] = base
        summary["beats_balanced_chance"] = summary["balanced_acc"] > 0.5

    results_df.to_csv(args.out, index=False)
    summary.to_csv(args.summary_out, index=False)

    if getattr(args, "paper_tables", None):
        write_paper_summary(ref_mat, ref_labels, query_profiles, query_truth,
                            args, gene_lengths)

    maj, base = majority_baseline(query_truth)
    print(f"\nmajority-class baseline: predict '{maj}' -> {base:.0%}")
    print("\n=== accuracy summary (sorted) ===")
    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(summary.to_string(index=False))

    # ---------------- optional positive controls ----------------
    if args.control:
        print("\n=== within-modality positive controls (leave-one-out kNN) ===")
        # reference-only: expression vs expression
        ref_samples = [s for s in ref_mat.columns if s in ref_labels]
        rmat = ref_mat[ref_samples]
        rlab = [ref_labels[s] for s in ref_samples]
        for metric in metrics:
            for k in ks:
                acc, _ = loo_control(rmat, rlab, ref_samples, metric, k)
                print(f"  reference (expr) LOO  metric={metric:<9} k={k}: acc={acc:.2f}")
        # query-only: variant vs variant (uses our cohort as its own reference;
        # this DOES use other query samples -- it is a control, not the main test)
        q_ids = [p for p in query_profiles if p in query_truth]
        if len(q_ids) >= 3:
            universe_q = sorted(set().union(*[set(query_profiles[p].index) for p in q_ids]))
            qmat = pd.DataFrame({p: query_profiles[p].reindex(universe_q).fillna(0.0)
                                 for p in q_ids})
            qlab = [query_truth[p] for p in q_ids]
            for metric in metrics:
                for k in ks:
                    acc, _ = loo_control(qmat, qlab, q_ids, metric, k)
                    print(f"  query (variant) LOO   metric={metric:<9} k={k}: acc={acc:.2f}")

    return results_df, summary


if __name__ == "__main__":
    main()
