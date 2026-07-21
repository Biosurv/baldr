#!/usr/bin/env python3
"""Generate the static test fixtures committed under tests/data/.

Run this by hand when you want to add or change a fixture, then commit the
regenerated files:

    python tests/data/make_fixtures.py
"""

import json
import os
import random

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))

REF_NAME = "AMPLICON_REF"
REF_LEN = 600

MARKER_POSITIONS = [30, 70, 110, 150, 190, 230, 270, 310, 350, 385]
MARKER_FREE_START = 420
MARKER_FREE_END = 580


LINEAGE_MARKERS = {
    "A.1": [30, 70, 110, 150, 190],
    "B.1": [150, 190, 230, 270, 310],
    "C.1": [270, 310, 350, 385],
}


UNMODELLED_MARKERS = {"D.1": [30, 190, 350]}

BASE_QUALITY = 30
ERROR_RATE = 0.005

SEED = 42


def make_reference(rng):
    return "".join(rng.choice("ACGT") for _ in range(REF_LEN))


def make_alt_bases(reference, rng):
    alts = {}
    for pos in MARKER_POSITIONS:
        ref_base = reference[pos - 1]
        alts[pos] = rng.choice([b for b in "ACGT" if b != ref_base])
    return alts


def marker_column(reference, alts, pos):
    return f"{reference[pos - 1]}{pos}{alts[pos]}"


def write_reference(reference, path):
    with open(path, "w") as fh:
        fh.write(f">{REF_NAME}\n")
        for i in range(0, len(reference), 60):
            fh.write(reference[i : i + 60] + "\n")


def write_barcode(reference, alts, lineages, path):


    columns = [marker_column(reference, alts, p) for p in MARKER_POSITIONS]
    with open(path, "w") as fh:
        fh.write("lineage," + ",".join(columns) + "\n")
        for lin, positions in lineages.items():
            carried = set(positions)
            cells = ["1" if p in carried else "0" for p in MARKER_POSITIONS]
            fh.write(lin + "," + ",".join(cells) + "\n")


def build_haplotype(reference, alts, positions):
    seq = list(reference)
    for pos in positions:
        seq[pos - 1] = alts[pos]
    return "".join(seq)


def sequence_with_errors(haplotype, start, end, rng, error_rate=ERROR_RATE):
    seq = list(haplotype[start:end])
    for i, base in enumerate(seq):
        if rng.random() < error_rate:
            seq[i] = rng.choice([b for b in "ACGT" if b != base])
    return "".join(seq)


def make_segment(name, seq, start, lineage, flag=0):
    seg = pysam.AlignedSegment()
    seg.query_name = name
    seg.flag = flag
    seg.reference_id = 0
    seg.reference_start = start
    seg.mapping_quality = 60
    seg.cigartuples = [(0, len(seq))]
    seg.query_sequence = seq
    seg.query_qualities = pysam.qualitystring_to_array(
        chr(BASE_QUALITY + 33) * len(seq)
    )
    seg.set_tag("XL", lineage, value_type="Z")
    return seg


def write_bam(segments, path):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": REF_NAME, "LN": REF_LEN}],
    }
    segments = sorted(segments, key=lambda s: s.reference_start)
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for seg in segments:
            out.write(seg)
    pysam.index(path)


def simulate_single(haplotypes, composition, n_reads, rng, error_rate=ERROR_RATE):


    segments = []
    lineages = list(composition)
    counts = {lin: int(round(composition[lin] * n_reads)) for lin in lineages}
    drift = n_reads - sum(counts.values())
    if drift:
        counts[max(counts, key=counts.get)] += drift

    read_id = 0
    for lin in lineages:
        for _ in range(counts[lin]):
            start = rng.randint(0, 20)
            end = rng.randint(390, 400)
            seq = sequence_with_errors(
                haplotypes[lin], start, end, rng, error_rate
            )
            segments.append(
                make_segment(f"r{read_id:05d}", seq, start, lin)
            )
            read_id += 1
    return segments, counts


def simulate_no_sites(reference, n_reads, rng):
    segments = []
    for read_id in range(n_reads):
        start = rng.randint(MARKER_FREE_START, MARKER_FREE_START + 20)
        end = rng.randint(MARKER_FREE_END - 20, MARKER_FREE_END)
        seq = sequence_with_errors(reference, start, end, rng)
        segments.append(
            make_segment(f"nosite{read_id:05d}", seq, start, "none")
        )
    return segments


def simulate_paired(haplotypes, composition, n_fragments, rng):
    segments = []
    lineages = list(composition)
    counts = {lin: int(round(composition[lin] * n_fragments)) for lin in lineages}
    drift = n_fragments - sum(counts.values())
    if drift:
        counts[max(counts, key=counts.get)] += drift

    frag_id = 0
    for lin in lineages:
        for _ in range(counts[lin]):
            name = f"frag{frag_id:05d}"
            r1_start, r1_end = rng.randint(0, 10), rng.randint(200, 210)
            r2_start, r2_end = rng.randint(215, 225), rng.randint(390, 400)

            r1_seq = sequence_with_errors(haplotypes[lin], r1_start, r1_end, rng)
            r2_seq = sequence_with_errors(haplotypes[lin], r2_start, r2_end, rng)


            r1 = make_segment(name, r1_seq, r1_start, lin, flag=99)
            r2 = make_segment(name, r2_seq, r2_start, lin, flag=147)

            tlen = r2_end - r1_start
            r1.next_reference_id = 0
            r1.next_reference_start = r2_start
            r1.template_length = tlen
            r2.next_reference_id = 0
            r2.next_reference_start = r1_start
            r2.template_length = -tlen

            segments.extend([r1, r2])
            frag_id += 1
    return segments, counts


def simulate_short(haplotypes, composition, n_reads, rng, read_length=60):
    segments = []
    lineages = list(composition)
    counts = {lin: int(round(composition[lin] * n_reads)) for lin in lineages}
    drift = n_reads - sum(counts.values())
    if drift:
        counts[max(counts, key=counts.get)] += drift

    read_id = 0
    for lin in lineages:
        for _ in range(counts[lin]):
            start = rng.randint(0, 400 - read_length)
            seq = sequence_with_errors(
                haplotypes[lin], start, start + read_length, rng
            )
            segments.append(make_segment(f"short{read_id:05d}", seq, start, lin))
            read_id += 1
    return segments, counts


def simulate_paired_uninformative(reference, haplotypes, composition,
                                  n_fragments, rng):
    segments = []
    lineages = list(composition)
    counts = {lin: int(round(composition[lin] * n_fragments)) for lin in lineages}
    drift = n_fragments - sum(counts.values())
    if drift:
        counts[max(counts, key=counts.get)] += drift

    frag_id = 0
    for lin in lineages:
        for _ in range(counts[lin]):
            name = f"splitfrag{frag_id:05d}"
            r1_start, r1_end = rng.randint(0, 10), rng.randint(390, 400)
            r2_start = rng.randint(MARKER_FREE_START, MARKER_FREE_START + 20)
            r2_end = rng.randint(MARKER_FREE_END - 20, MARKER_FREE_END)

            r1_seq = sequence_with_errors(haplotypes[lin], r1_start, r1_end, rng)
            r2_seq = sequence_with_errors(reference, r2_start, r2_end, rng)

            r1 = make_segment(name, r1_seq, r1_start, lin, flag=99)
            r2 = make_segment(name, r2_seq, r2_start, lin, flag=147)

            tlen = r2_end - r1_start
            r1.next_reference_id = 0
            r1.next_reference_start = r2_start
            r1.template_length = tlen
            r2.next_reference_id = 0
            r2.next_reference_start = r1_start
            r2.template_length = -tlen

            segments.extend([r1, r2])
            frag_id += 1
    return segments, counts


def rng_for(name):
    return random.Random(f"{SEED}:{name}")


def main():
    reference = make_reference(rng_for("reference"))
    alts = make_alt_bases(reference, rng_for("alts"))

    write_reference(reference, os.path.join(HERE, "reference.fasta"))
    write_barcode(reference, alts, LINEAGE_MARKERS, os.path.join(HERE, "barcode.csv"))

    edge_lineages = dict(LINEAGE_MARKERS)
    edge_lineages["B.1.dup"] = list(LINEAGE_MARKERS["B.1"])
    edge_lineages["REF.1"] = []
    write_barcode(
        reference, alts, edge_lineages, os.path.join(HERE, "barcode_edge.csv")
    )

    haplotypes = {
        lin: build_haplotype(reference, alts, positions)
        for lin, positions in LINEAGE_MARKERS.items()
    }
    haplotypes.update(
        {
            lin: build_haplotype(reference, alts, positions)
            for lin, positions in UNMODELLED_MARKERS.items()
        }
    )

    samples = []

    single_specs = [
        ("pure_A.bam", {"A.1": 1.0}, 300, ERROR_RATE, True,
         "Single lineage, no mixture. Sanity anchor."),
        ("mix_75_25.bam", {"A.1": 0.75, "B.1": 0.25}, 400, ERROR_RATE, True,
         "Two-lineage mixture at 75/25."),
        ("mix_50_50.bam", {"A.1": 0.50, "B.1": 0.50}, 400, ERROR_RATE, True,
         "Two-lineage mixture at 50/50."),
        ("mix_95_05.bam", {"A.1": 0.95, "B.1": 0.05}, 1200, ERROR_RATE, True,
         "Skewed mixture. The minor lineage is where EM's prior feedback "
         "would push reads toward the dominant lineage."),
        ("mix_three.bam", {"A.1": 0.50, "B.1": 0.30, "C.1": 0.20}, 1000,
         ERROR_RATE, True, "Three-lineage mixture."),
        ("mix_noisy.bam", {"A.1": 0.60, "B.1": 0.40}, 800, 0.08, True,
         "Same shape as the clean mixtures but at a 16x higher error rate, "
         "so assignment is no longer trivially perfect. Gives accuracy and "
         "ambiguity rates that can actually regress."),
        ("unmodelled.bam", {"A.1": 0.70, "D.1": 0.30}, 500, ERROR_RATE, False,
         "30% of reads come from D.1, which is absent from barcode.csv. "
         "Measures how much an unmodelled lineage inflates listed lineages."),
    ]

    for (filename, composition, n_reads, error_rate, recoverable,
         description) in single_specs:
        segments, counts = simulate_single(
            haplotypes, composition, n_reads, rng_for(filename), error_rate
        )
        write_bam(segments, os.path.join(HERE, filename))
        samples.append(
            {
                "bam": filename,
                "barcode": "barcode.csv",
                "n_reads": sum(counts.values()),
                "read_counts": counts,
                "truth": {
                    lin: n / sum(counts.values()) for lin, n in counts.items()
                },
                "error_rate": error_rate,
                "paired": False,
                "recoverable": recoverable,
                "description": description,
            }
        )

    short_segments, short_counts = simulate_short(
        haplotypes, {"A.1": 0.50, "B.1": 0.50}, 600, rng_for("short_reads.bam")
    )
    write_bam(short_segments, os.path.join(HERE, "short_reads.bam"))
    samples.append(
        {
            "bam": "short_reads.bam",
            "barcode": "barcode.csv",
            "n_reads": sum(short_counts.values()),
            "read_counts": short_counts,
            "truth": {
                lin: n / sum(short_counts.values())
                for lin, n in short_counts.items()
            },
            "error_rate": ERROR_RATE,
            "read_length": 60,
            "paired": False,
            "recoverable": True,
            "description": "60 bp reads scattered across the marker region. "
            "Most cover only one or two markers, many of which are shared "
            "between lineages, so assignment is genuinely uncertain.",
        }
    )

    no_site_segments = simulate_no_sites(reference, 100, rng_for("no_sites.bam"))
    write_bam(no_site_segments, os.path.join(HERE, "no_sites.bam"))
    samples.append(
        {
            "bam": "no_sites.bam",
            "barcode": "barcode.csv",
            "n_reads": len(no_site_segments),
            "read_counts": {},
            "truth": {},
            "paired": False,
            "recoverable": False,
            "description": "Reads aligned outside the marker region entirely. "
            "Every read should be unassignable for lack of evidence.",
        }
    )

    paired_segments, paired_counts = simulate_paired(
        haplotypes, {"A.1": 0.60, "C.1": 0.40}, 200, rng_for("paired.bam")
    )
    write_bam(paired_segments, os.path.join(HERE, "paired.bam"))
    samples.append(
        {
            "bam": "paired.bam",
            "barcode": "barcode.csv",
            "n_reads": len(paired_segments),
            "n_fragments": sum(paired_counts.values()),
            "read_counts": paired_counts,
            "truth": {
                lin: n / sum(paired_counts.values())
                for lin, n in paired_counts.items()
            },
            "paired": True,
            "recoverable": True,
            "description": "Paired fragments whose mates cover disjoint marker "
            "sets. Both mates share a query name.",
        }
    )

    split_segments, split_counts = simulate_paired_uninformative(
        reference, haplotypes, {"A.1": 0.50, "C.1": 0.50}, 100,
        rng_for("paired_split.bam")
    )
    write_bam(split_segments, os.path.join(HERE, "paired_split.bam"))
    samples.append(
        {
            "bam": "paired_split.bam",
            "barcode": "barcode.csv",
            "n_reads": len(split_segments),
            "n_fragments": sum(split_counts.values()),
            "read_counts": split_counts,
            "truth": {
                lin: n / sum(split_counts.values())
                for lin, n in split_counts.items()
            },
            "paired": True,
            "recoverable": True,
            "description": "Paired fragments where read 1 spans every marker "
            "and read 2 lands in the marker-free tail. Every fragment is "
            "identifiable, but the two mates carry very different evidence "
            "under a shared query name.",
        }
    )

    manifest = {
        "seed": SEED,
        "reference": "reference.fasta",
        "reference_name": REF_NAME,
        "reference_length": REF_LEN,
        "default_barcode": "barcode.csv",
        "edge_barcode": "barcode_edge.csv",
        "base_quality": BASE_QUALITY,
        "error_rate": ERROR_RATE,
        "marker_positions": MARKER_POSITIONS,
        "lineage_markers": LINEAGE_MARKERS,
        "unmodelled_markers": UNMODELLED_MARKERS,
        "truth_tag": "XL",
        "samples": samples,
    }
    with open(os.path.join(HERE, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")

    print(f"Wrote {len(samples)} samples to {HERE}")
    for sample in samples:
        print(f"  {sample['bam']:20s} {sample['n_reads']:5d} reads")


if __name__ == "__main__":
    main()
