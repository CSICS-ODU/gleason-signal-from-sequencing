import pandas as pd
from pathlib import Path
from tqdm import tqdm
from collections import defaultdict
import csv
import re
import pdb

# ==============================
# CONFIG
# ==============================
metadata_path = "metadata.csv"
processed_root = Path("processed")   # IMPORTANT FIX
output_root = Path("processed_filtered")

modes = ["unmatched", "matched", "benign"]

# ==============================
# LOAD METADATA
# ==============================
metadata = pd.read_csv(metadata_path)

# ==============================
# GLEASON PARSING
# ==============================
def parse_gleason(val):
    try:
        val = str(val).strip()
        if "+" in val:
            a, b = val.split("+")
            return int(a), int(b)
        else:
            return int(val), 0
    except:
        return None, None

# ==============================
# INDIAN EXTRACTION (B or A2)
# ==============================
def extract_indian_pattern(extra):
    try:
        extra = str(extra)

        match_b = re.search(r'B\s*\((\d+\+\d+)\)', extra)
        if match_b:
            return match_b.group(1)

        match_a2 = re.search(r'A2\s*\((\d+\+\d+)\)', extra)
        if match_a2:
            return match_a2.group(1)

        return None
    except:
        return None

# ==============================
# FINAL GLEASON GROUP
# ==============================
def gleason_group(row):
    print (row)
    # pdb.set_trace()
    race = row["race"]

    if race == "IN":
        pattern = extract_indian_pattern(row["Extra"])
        if pattern:
            a, b = parse_gleason(pattern)
        else:
            a, b = parse_gleason(row["Gleason"])
    else:
        a, b = parse_gleason(row["Gleason"])

    if a is None:
        return None

    # SPECIAL RULE
    if a == 3 and b == 4:
        return "low"
    elif a == 4 and b == 3:
        return "high"

    score = a + b
    return "low" if score <= 6 else "high"

metadata["gleason_group"] = metadata.apply(gleason_group, axis=1)

# ==============================
# SNV FUNCTIONS
# ==============================
def is_snv(line):
    parts = line.split("\t")
    return len(parts[3]) == 1 and all(len(a) == 1 for a in parts[4].split(","))

def snv_key(line):
    parts = line.strip().split("\t")
    return f"{parts[0]}:{parts[1]}:{parts[3]}>{parts[4]}"

# ==============================
# GROUP DEFINITIONS
# ==============================
races = ["ALL", "IN", "EA", "AA"]
gleasons = ["high", "low", "ALL"]

# ==============================
# MAIN LOOP
# ==============================
for mode in modes:

    print(f"\n\n========== MODE: {mode.upper()} ==========")

    table_counts = defaultdict(lambda: {"high": 0, "low": 0, "ALL": 0})
    snv_counts = defaultdict(lambda: {"high": 0, "low": 0, "ALL": 0})

    # ✅ USE CORRECT INPUT FOLDER
    mode_input_folder = processed_root / mode
    files = list(mode_input_folder.rglob("*.vcf"))

    print(f"Reading from: {mode_input_folder}")
    print("VCF files found:", len(files))

    for race in races:
        for gleason in gleasons:

            group_name = f"{race}_{gleason}_{mode}"
            print(f"\nProcessing: {group_name}")

            df = metadata.copy()

            if race != "ALL":
                df = df[df["race"] == race]

            if gleason != "ALL":
                df = df[df["gleason_group"] == gleason]

            valid_ids = set(df["ID"].astype(str))

            table_counts[race][gleason] = len(valid_ids)

            group_dir = output_root / mode / group_name
            group_dir.mkdir(parents=True, exist_ok=True)

            snv_dict = defaultdict(set)

            # ==============================
            # PROCESS FILES
            # ==============================
            for file_path in tqdm(files):

                sample_id = file_path.stem.split(".")[0]

                if sample_id not in valid_ids:
                    continue

                with open(file_path) as f:
                    lines = f.readlines()

                headers = []
                records = []

                for line in lines:
                    if line.startswith("#"):
                        headers.append(line)
                    else:
                        if is_snv(line):
                            records.append(line)
                            snv_dict[snv_key(line)].add(sample_id)

                # SAVE VCF
                out_vcf = group_dir / f"{sample_id}.{mode}.vcf"
                with open(out_vcf, "w") as out:
                    out.writelines(headers + records)

            # ==============================
            # COMMON SNV COUNT
            # ==============================
            common_snvs = [snv for snv, p in snv_dict.items() if len(p) >= 2]
            snv_counts[race][gleason] = len(common_snvs)

            # ==============================
            # SAVE COMMON SNVs
            # ==============================
            out_csv = group_dir / "common_snvs.csv"
            with open(out_csv, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["SNV", "Count", "Patients"])

                for snv in common_snvs:
                    patients = snv_dict[snv]
                    writer.writerow([snv, len(patients), ",".join(sorted(patients))])

        # ✅ FIX ALL COLUMN (NO DOUBLE COUNTING)
        table_counts[race]["ALL"] = (
            table_counts[race]["high"] + table_counts[race]["low"]
        )

        snv_counts[race]["ALL"] = (
            snv_counts[race]["high"] + snv_counts[race]["low"]
        )

    # ==============================
    # SAVE TABLES
    # ==============================
    mode_dir = output_root / mode

    # PATIENT TABLE
    table_path = mode_dir / "group_summary_table.csv"
    with open(table_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Group", "High", "Low", "ALL"])

        for race in races:
            writer.writerow([
                race,
                table_counts[race]["high"],
                table_counts[race]["low"],
                table_counts[race]["ALL"]
            ])

    # SNV TABLE
    snv_table_path = mode_dir / "group_snv_summary_table.csv"
    with open(snv_table_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Group", "High", "Low", "ALL"])

        for race in races:
            writer.writerow([
                race,
                snv_counts[race]["high"],
                snv_counts[race]["low"],
                snv_counts[race]["ALL"]
            ])

    print(f"\nSaved tables for {mode}")

print("\nAll modes processed successfully.")