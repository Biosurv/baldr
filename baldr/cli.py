import sys
import os
import argparse
import glob

import pysam

from ._version import __version__
from .barcode import read_barcode_tsv
from .assign import assign_reads_on_bam


def parse_args():
    p = argparse.ArgumentParser(
        description=f"Assign reads to lineages using barcode SNVs (v{__version__})."
    )
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--bam", help="Input BAM (indexed).")
    mode.add_argument("--bam-dir", help="Directory containing BAM files (process all *.bam).")

    p.add_argument(
        "--barcode-csv",
        required=True,
        help=(
            "Barcode table in wide format (CSV or TSV). First column is lineage name "
            "(or lin/name). Remaining columns are markers like A123G, C456T, etc.; "
            "cells with value 1 indicate that lineage has that SNV."
        ),
    )
    p.add_argument(
        "--lineage-include",
        nargs="*",
        default=None,
        help="Optional explicit list of lineage names to consider. If omitted, all in table are used.",
    )

    p.add_argument("--mapq-min", type=int, default=0, help="Minimum MAPQ to consider a read [0].")

    p.add_argument(
        "--baseq-min",
        type=int,
        default=7,
        help="Minimum base quality to count a site [7].",
    )
    p.add_argument(
        "--min-sites",
        type=int,
        default=0,
        help=(
            "Minimum informative sites covered for assignment [0]. In EM mode this "
            "is a per-read floor (number of barcode positions covered). In voting "
            "mode it is per-lineage (number of the winning lineage's defining "
            "sites covered)."
        ),
    )
    p.add_argument(
        "--min-margin",
        type=int,
        default=1,
        help="(Voting mode only) Minimum score margin over second-best lineage to assign [1].",
    )
    p.add_argument(
        "--ignore-ref-mismatch",
        action="store_true",
        help=(
            "If set, do not subtract a vote when the read shows the reference base "
            "(treat as 0, not -1)."
        ),
    )
    p.add_argument(
        "--allow-secondary",
        action="store_true",
        help="Include secondary/supplementary alignments in assignment (default: primary only).",
    )

    p.add_argument(
        "--voting",
        action="store_true",
        help=(
            "Use voting-based assignment instead of EM. Voting cannot detect "
            "reference-equivalent lineages (rows that are all 0 in the barcode)."
        ),
    )
    p.add_argument(
        "--post-min",
        type=float,
        default=0.90,
        help="Posterior threshold to assign a read to a lineage in EM mode [0.90].",
    )
    p.add_argument(
        "--em-max-iters",
        type=int,
        default=50,
        help="Max EM iterations [50].",
    )
    p.add_argument(
        "--em-tol",
        type=float,
        default=1e-4,
        help="EM convergence tolerance on mixture weights [1e-4].",
    )

    p.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Threads for BAM reading (pysam).",
    )
    p.add_argument(
        "--outdir",
        required=True,
        help="Output directory (parent for folder mode).",
    )
    p.add_argument(
        "--prefix",
        default=None,
        help="Prefix for outputs. In --bam-dir mode, defaults to each BAM basename if not set.",
    )
    p.add_argument(
        "--write-bams",
        action="store_true",
        help="Write per-lineage BAMs.",
    )
    p.add_argument(
        "--write-fastq",
        action="store_true",
        help="Write per-lineage FASTQs (implies iterating reads twice).",
    )
    p.add_argument(
        "--write-unassigned",
        action="store_true",
        help="Also write BAM/FASTQ for ambiguous/unassigned.",
    )
    p.add_argument(
        "--bam-glob",
        default="*.bam",
        help="Only used with --bam-dir: glob pattern to match BAMs [*.bam].",
    )
    p.add_argument(
        "--write-mix-summary",
        action="store_true",
        help="Write per-sample mix summary TSV with read counts and fractions per lineage.",
    )

    return p.parse_args()


def main():
    args = parse_args()
    args.use_em = not args.voting

    lineage_sites, site_to_lineage_alt, all_positions, pos_to_ref = read_barcode_tsv(
        args.barcode_csv, args.lineage_include
    )

    if args.voting:
        for lin, sites in sorted(lineage_sites.items()):
            n_sites = len(sites)
            if 0 < n_sites < args.min_sites:
                print(
                    f"[WARN] Lineage '{lin}' has only {n_sites} defining site(s), "
                    f"fewer than --min-sites {args.min_sites}. "
                    f"Effective minimum capped at {n_sites} for this lineage.",
                    file=sys.stderr,
                )
            elif n_sites == 0:
                print(
                    f"[WARN] Lineage '{lin}' has no defining mutations and is "
                    f"undetectable in voting mode. Use EM (default) to detect "
                    f"reference-equivalent lineages.",
                    file=sys.stderr,
                )

    if args.bam:
        bam_path = args.bam
        if not os.path.exists(f"{bam_path}.bai"):
            try:
                pysam.index(bam_path)
            except Exception as e:
                sys.exit(f"[ERROR] Failed to index BAM: {bam_path}\n{e}")
        outdir = args.outdir
        prefix = args.prefix if args.prefix else os.path.splitext(
            os.path.basename(bam_path)
        )[0]
        rc = assign_reads_on_bam(
            bam_path,
            outdir,
            prefix,
            args,
            lineage_sites,
            site_to_lineage_alt,
            all_positions,
            pos_to_ref,
        )
        sys.exit(rc)

    os.makedirs(args.outdir, exist_ok=True)
    bam_paths = sorted(glob.glob(os.path.join(args.bam_dir, args.bam_glob)))
    if not bam_paths:
        sys.exit(
            f"[ERROR] No BAMs matched pattern {args.bam_glob} in {args.bam_dir}"
        )

    overall_rc = 0
    for bam_path in bam_paths:
        sample = os.path.splitext(os.path.basename(bam_path))[0]
        sample_outdir = os.path.join(args.outdir, sample)
        prefix = args.prefix if args.prefix else sample
        if not os.path.exists(f"{bam_path}.bai"):
            try:
                pysam.index(bam_path)
            except Exception as e:
                print(
                    f"[WARN] Failed to index {bam_path}: {e}", file=sys.stderr
                )
                overall_rc = 1
                continue

        print(f">> Processing {sample}")
        rc = assign_reads_on_bam(
            bam_path,
            sample_outdir,
            prefix,
            args,
            lineage_sites,
            site_to_lineage_alt,
            all_positions,
            pos_to_ref,
        )
        if rc != 0:
            overall_rc = rc

    sys.exit(overall_rc)


if __name__ == "__main__":
    main()
