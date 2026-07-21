import sys
import re
from collections import defaultdict

import pandas as pd


def read_barcode_tsv(path, lineage_include=None):
    sep = "," if path.endswith(".csv") else "\t"
    try:
        df = pd.read_csv(path, sep=sep)
    except Exception as e:
        sys.exit(f"[ERROR] Failed to read barcode table '{path}': {e}")

    if df.shape[1] < 2:
        sys.exit("[ERROR] Barcode table must have at least 2 columns (lineage + one marker).")

    first_col = df.columns[0]
    if first_col.lower() not in {"lineage", "lin", "name"}:
        df = df.rename(columns={first_col: "lineage"})
    else:
        df = df.rename(columns={first_col: "lineage"})

    all_lineages = []
    long_records = []
    pat = re.compile(r"^([ACGT])(\d+)([ACGT])$")

    for _, row in df.iterrows():
        lin = row["lineage"]
        if lineage_include and lin not in lineage_include:
            continue
        all_lineages.append(lin)

        for col in df.columns[1:]:
            m = pat.match(col)
            if not m:
                continue
            ref, pos_str, alt = m.group(1), m.group(2), m.group(3)
            val = row[col]

            try:
                present = int(val) == 1
            except Exception:
                present = str(val).strip() == "1"

            if not present:
                continue

            try:
                pos = int(pos_str)
            except ValueError:
                continue

            long_records.append((lin, pos, ref, alt))

    if not all_lineages:
        sys.exit("[ERROR] No lineages loaded (check --lineage-include or barcode CSV/TSV).")
    if not long_records:
        sys.exit(
            "[ERROR] No marker mutations found across any lineage in the barcode "
            "table. At least one lineage must have a defining mutation to anchor "
            "the reference base for each position."
        )

    lineage_sites = {lin: {} for lin in all_lineages}
    pos_to_ref = {}
    all_positions = set()

    long_records.sort(key=lambda x: (x[0], x[1], x[3]))

    for lin, pos, ref, alt in long_records:
        ref = str(ref).upper()
        alt = str(alt).upper()
        if len(ref) != 1 or len(alt) != 1:
            continue

        if pos in pos_to_ref and pos_to_ref[pos] != ref:
            sys.exit(
                f"[ERROR] Conflicting REF bases at pos {pos} in barcode table: "
                f"{pos_to_ref[pos]} vs {ref}"
            )
        pos_to_ref[pos] = ref

        lineage_sites[lin][pos] = (ref, alt)
        all_positions.add(pos)

    empty_lineages = [lin for lin in all_lineages if not lineage_sites[lin]]
    if empty_lineages:
        print(
            f"[INFO] {len(empty_lineages)} lineage(s) have no defining mutations "
            f"in the barcode (reference-equivalent). These are detectable in EM "
            f"mode only.",
            file=sys.stderr,
        )

    site_to_lineage_alt = defaultdict(dict)
    for lin, sites in lineage_sites.items():
        for pos, (ref, alt) in sites.items():
            site_to_lineage_alt[pos][lin] = alt

    return lineage_sites, site_to_lineage_alt, sorted(all_positions), pos_to_ref
