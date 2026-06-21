import csv
import glob
import os
import re
from pathlib import Path
from collections import defaultdict, Counter
import pysam

INPUT_ROOT = Path("processed_filtered")

MATCHED_DIR = "matched"
COSMIC_DIR = "cosmic"
CLINVAR_DIR = "clinvar"

# --- gene-matrix rebuild (operates on snv_rebuild/ master files) ---
SNV_ROOT = Path("snv_rebuild")          # output of rebuild_snv_matrix.py
PROCESSED = Path("processed")           # per-patient VCFs (for Indian inline genes)
METADATA = "metadata.csv"
LAYERS = ["benign", "matched", "unmatched"]
CATEGORIES = ["high_only", "low_only", "high_low_common"]

# Indian VCFs carry the gene inline, e.g. INFO ...FUNC=[{'gene':'MET',...}]
GENE_INLINE_RE = re.compile(r"'gene'\s*:\s*'([^']+)'")


def find_indexed_vcfs():
    """
    Find bgzip-compressed VCFs that have tabix indexes.
    Requires both:
        file.vcf.gz
        file.vcf.gz.tbi
    """
    vcf_files = []

    for folder in [MATCHED_DIR, COSMIC_DIR, CLINVAR_DIR]:
        vcf_files.extend(glob.glob(os.path.join(folder, "*.vcf.gz")))

    indexed = []

    for path in sorted(vcf_files):
        if os.path.exists(path + ".tbi"):
            indexed.append(Path(path))
        else:
            print("Skipping unindexed VCF:", path)

    return indexed


def parse_info_field(info):
    """
    Convert VCF INFO string into a dictionary.
    Example:
        GENE=TP53;AF=0.1
    becomes:
        {"GENE": "TP53", "AF": "0.1"}
    """
    info_dict = {}

    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info_dict[key] = value
        else:
            info_dict[item] = True

    return info_dict


def clean_gene_name(gene):
    gene = gene.strip()

    if not gene:
        return None

    bad_values = {".", "-", "na", "n/a", "none", "null", "unknown"}
    if gene.lower() in bad_values:
        return None

    return gene


def extract_genes_from_info(info, ann_gene_idx=3, csq_gene_idx=3):
    """
    Tries several common VCF annotation styles.

    This is intentionally tolerant because COSMIC, ClinVar, ANNOVAR,
    VEP, and SnpEff may encode genes differently.

    ann_gene_idx / csq_gene_idx are the sub-field positions of the gene
    symbol within the SnpEff ANN and VEP CSQ pipe-delimited values. They
    default to 3 but are normally derived from the VCF header's Format line
    (see extract_annotation_format).
    """
    genes = set()
    info_dict = parse_info_field(info)

    # Simple direct gene fields
    direct_gene_keys = [
        "GENE",
        "Gene",
        "SYMBOL",
        "Gene.refGene",
        "GENEINFO",
        "GENENAME",
        "HGNC",
    ]

    for key in direct_gene_keys:
        if key not in info_dict:
            continue

        value = str(info_dict[key])

        # ClinVar-style GENEINFO can look like:
        # TP53:7157|BRCA1:672
        if key == "GENEINFO":
            chunks = re.split(r"[|,]", value)
            for chunk in chunks:
                gene = chunk.split(":")[0]
                gene = clean_gene_name(gene)
                if gene:
                    genes.add(gene)
        else:
            chunks = re.split(r"[,|]", value)
            for chunk in chunks:
                gene = clean_gene_name(chunk)
                if gene:
                    genes.add(gene)

    # SnpEff ANN format:
    # ANN=Allele|Annotation|Impact|Gene_Name|Gene_ID|...
    if "ANN" in info_dict:
        annotations = str(info_dict["ANN"]).split(",")

        for ann in annotations:
            fields = ann.split("|")
            if len(fields) > ann_gene_idx:
                gene = clean_gene_name(fields[ann_gene_idx])
                if gene:
                    genes.add(gene)

    # VEP CSQ format:
    # CSQ=Allele|Consequence|IMPACT|SYMBOL|Gene|...
    #
    # csq_gene_idx is taken from the header's "Format:" line when available
    # (extract_annotation_format), so we no longer hard-assume SYMBOL is field 3.
    if "CSQ" in info_dict:
        annotations = str(info_dict["CSQ"]).split(",")

        for csq in annotations:
            fields = csq.split("|")
            if len(fields) > csq_gene_idx:
                gene = clean_gene_name(fields[csq_gene_idx])
                if gene:
                    genes.add(gene)

    return genes


def chrom_candidates(chrom):
    chrom = chrom.strip()

    if chrom.startswith("chr"):
        return [chrom, chrom[3:]]

    return [chrom, "chr" + chrom]


def parse_variant_key(variant_key):
    """
    Parse a variant key into (chrom, pos, ref, alt).

    Accepts:
        chr5:132198287:C>G
        5:132198287:C>G

    Returns None for anything that doesn't match this layout (so callers
    can skip+warn instead of crashing the whole run on one bad row).
    """
    try:
        chrom, pos, rest = variant_key.split(":", 2)
        ref, alt = rest.split(">")
        return chrom, int(pos), ref.upper(), alt.upper()
    except (ValueError, AttributeError):
        return None


def parse_format_fields(header_line):
    """
    Pull the pipe-delimited field list out of a VEP/SnpEff INFO header line.

    VEP:    ...Description="...Format: Allele|Consequence|IMPACT|SYMBOL|Gene|...">
    SnpEff: ...Description="Functional annotations: 'Allele | Annotation | ... | Gene_Name | ...'">
    """
    match = re.search(r"Format:\s*([^\">]+)", header_line)
    if not match:
        match = re.search(r"'\s*([^']*\|[^']*)'", header_line)
    if not match:
        return None
    return [f.strip() for f in match.group(1).split("|")]


def gene_field_index(fields, preferred):
    """Index of the first preferred field name (case-insensitive), else None."""
    if not fields:
        return None
    lowered = [f.lower() for f in fields]
    for name in preferred:
        if name.lower() in lowered:
            return lowered.index(name.lower())
    return None


def extract_annotation_format(tbx):
    """
    Inspect the VCF header to find which CSQ/ANN sub-field holds the gene
    symbol, rather than hard-coding index 3. Falls back to 3 when the header
    is absent or unparseable.
    """
    csq_idx = None
    ann_idx = None

    try:
        header_lines = list(tbx.header)
    except Exception:
        header_lines = []

    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ"):
            csq_idx = gene_field_index(parse_format_fields(line),
                                       ["SYMBOL", "Gene_Name", "Gene"])
        elif line.startswith("##INFO=<ID=ANN"):
            ann_idx = gene_field_index(parse_format_fields(line),
                                       ["Gene_Name", "SYMBOL", "Gene"])

    return {
        "csq_gene_idx": csq_idx if csq_idx is not None else 3,
        "ann_gene_idx": ann_idx if ann_idx is not None else 3,
    }


def lookup_variant_in_open_tabix(handle, parsed_variant):
    """
    handle: dict with keys {"tbx", "csq_gene_idx", "ann_gene_idx"}.
    parsed_variant: (chrom, pos, ref, alt) from parse_variant_key().

    Returns the union of gene symbols annotated at this exact variant.
    """
    chrom, pos, ref, alt = parsed_variant
    tbx = handle["tbx"]

    genes = set()

    for query_chrom in chrom_candidates(chrom):
        try:
            records = tbx.fetch(query_chrom, pos - 1, pos)
        except ValueError:
            continue

        for line in records:
            parts = line.rstrip("\n").split("\t")

            if len(parts) < 8:
                continue

            rec_pos = int(parts[1])
            rec_ref = parts[3].upper()
            rec_alts = parts[4].upper().split(",")
            rec_info = parts[7]

            if rec_pos == pos and rec_ref == ref and alt in rec_alts:
                # Accumulate (don't early-return): a position can have several
                # matching records, each carrying different gene annotations.
                genes.update(extract_genes_from_info(
                    rec_info,
                    ann_gene_idx=handle["ann_gene_idx"],
                    csq_gene_idx=handle["csq_gene_idx"],
                ))

    return genes


def lookup_variant_across_vcfs(handles, parsed_variant):
    genes = set()

    for handle in handles:
        genes.update(lookup_variant_in_open_tabix(handle, parsed_variant))

    return genes


def load_ancestry():
    """patient_id -> ancestry (AA / EA / IN)."""
    anc = {}
    with open(METADATA, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            pid = (row.get("ID") or "").strip()
            if pid:
                anc[pid] = row["race"].strip()
    return anc


def build_indian_gene_map(indian_ids):
    """
    SNV key -> set(genes), read INLINE from the Indian per-patient VCFs.
    Indian variants annotate the gene in INFO, e.g. FUNC=[{'gene':'MET',...}].
    """
    gene_map = defaultdict(set)
    for layer in LAYERS:
        for pid in indian_ids:
            path = PROCESSED / layer / f"{pid}.{layer}.vcf"
            if not path.exists():
                continue
            with open(path, encoding="utf-8", errors="replace") as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 8:
                        continue
                    key = f"{p[0]}:{p[1]}:{p[3]}>{p[4]}"
                    for g in GENE_INLINE_RE.findall(p[7]):
                        g = clean_gene_name(g)
                        if g:
                            gene_map[key].add(g)
    return gene_map


def genes_for_snv(snv, handles, variant_cache, indian_gene_map):
    """Genes for an SNV: tabix-DB lookup UNION the Indian inline gene (if any)."""
    if snv not in variant_cache:
        parsed = parse_variant_key(snv)
        variant_cache[snv] = (lookup_variant_across_vcfs(handles, parsed)
                              if parsed is not None else set())
    return set(variant_cache[snv]) | indian_gene_map.get(snv, set())


def split_patients(field):
    return [p.strip() for p in (field or "").split(";") if p.strip()]


def master_to_gene_matrices(master_file, handles, variant_cache, indian_gene_map):
    """
    Read one *_master_*.csv and emit a gene x patient COUNT matrix per category.
    Cell = how many of that gene's SNVs the patient carries (magnitude, not a set).
    """
    counts = {cat: defaultdict(Counter) for cat in CATEGORIES}  # cat -> gene -> pt -> n
    high_pts, low_pts = set(), set()
    no_gene = 0

    with open(master_file, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None or "SNV" not in reader.fieldnames:
            print("  Skipping (no 'SNV'):", master_file)
            return
        for row in reader:
            snv = (row.get("SNV") or "").strip()
            cat = (row.get("category") or "").strip()
            if cat not in counts:
                continue
            hp = split_patients(row.get("high_patients"))
            lp = split_patients(row.get("low_patients"))
            high_pts.update(hp)
            low_pts.update(lp)
            genes = genes_for_snv(snv, handles, variant_cache, indian_gene_map)
            if not genes:
                no_gene += 1
                continue
            for gene in genes:
                for pt in hp + lp:
                    counts[cat][gene][pt] += 1

    # columns: high patients first, then low, so the high/low contrast reads left-to-right
    columns = sorted(high_pts) + sorted(low_pts)
    grade = {p: "high" for p in high_pts}
    grade.update({p: "low" for p in low_pts})

    gene_dir = master_file.parent.parent / "gene results"
    gene_dir.mkdir(parents=True, exist_ok=True)
    stem = master_file.name.replace("_master_", "_{cat}_genes_").rsplit(".", 1)[0]
    for cat in CATEGORIES:
        gp = counts[cat]
        out = gene_dir / (stem.format(cat=cat) + ".csv")
        with open(out, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["Gene"] + [f"{p}:{grade[p]}" for p in columns]
                       + ["Total", "n_high_pat", "n_low_pat"])
            for gene in sorted(gp):
                pc = gp[gene]
                vec = [pc.get(p, 0) for p in columns]
                n_high = sum(1 for p in high_pts if pc.get(p, 0) > 0)
                n_low = sum(1 for p in low_pts if pc.get(p, 0) > 0)
                w.writerow([gene] + vec + [sum(vec), n_high, n_low])

    summary = ", ".join(f"{c}={len(counts[c])}" for c in CATEGORIES)
    tail = f"  (no-gene SNVs: {no_gene})" if no_gene else ""
    print(f"  {master_file.parent.parent.name}/{master_file.name}: {summary}{tail}")


def main():
    ancestry = load_ancestry()
    indian_ids = [p for p, a in ancestry.items() if a == "IN"]
    print(f"patients: {len(ancestry)}  (Indian: {len(indian_ids)})")

    indexed_vcfs = find_indexed_vcfs()
    print("Indexed VCF files found:", len(indexed_vcfs))
    if not indexed_vcfs:
        raise RuntimeError("No indexed .vcf.gz files with .tbi indexes were found.")
    handles = []
    for path in indexed_vcfs:
        tbx = pysam.TabixFile(str(path), encoding="utf-8")
        handles.append({"tbx": tbx, **extract_annotation_format(tbx)})

    print("Building Indian inline gene map ...")
    indian_gene_map = build_indian_gene_map(indian_ids)
    print(f"  inline-annotated SNVs: {len(indian_gene_map)}")

    master_files = sorted(SNV_ROOT.rglob("*_master_*.csv"))
    print(f"master files: {len(master_files)}")

    variant_cache = {}
    try:
        for mf in master_files:
            master_to_gene_matrices(mf, handles, variant_cache, indian_gene_map)
    finally:
        for handle in handles:
            handle["tbx"].close()
    print("Unique SNVs gene-mapped:", len(variant_cache))


if __name__ == "__main__":
    main()