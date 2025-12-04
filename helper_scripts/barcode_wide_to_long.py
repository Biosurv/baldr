import re, sys, pandas as pd

if len(sys.argv) < 3:
    sys.exit("usage: barcodes_wide_to_long.py <wide.tsv/csv> <out.tsv>")

inp, outp = sys.argv[1], sys.argv[2]

sep = "," if inp.endswith(".csv") else "\t"

df = pd.read_csv(inp, sep=sep)

if df.columns[0].lower() not in {"lineage","lin","name"}:
    df = df.rename(columns={df.columns[0]: "lineage"})
else:
    df = df.rename(columns={df.columns[0]: "lineage"})

long = []
pat = re.compile(r"^([ACGT])(\d+)([ACGT])$")
for _, row in df.iterrows():
    lin = row["lineage"]
    for col in df.columns[1:]:
        m = pat.match(col)
        if not m:
            continue
        ref, pos, alt = m.group(1), int(m.group(2)), m.group(3)
        val = row[col]
        try:
            present = int(val) == 1
        except:
            present = str(val).strip() == "1"
        if present:
            long.append((lin, pos, ref, alt))

out = pd.DataFrame(long, columns=["lineage","pos","ref","alt"])
out.sort_values(["lineage","pos","alt"], inplace=True)
out.to_csv(outp, sep="\t", index=False)
print(f"Wrote {len(out)} lineage-site records to {outp}")
