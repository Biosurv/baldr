import sys
import os
import argparse
from collections import defaultdict, Counter
import glob
import math
import pysam
import re
import pandas as pd

from ._version import __version__

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None


# Argument parsing
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
        default=2,
        help="Minimum informative sites covered for assignment [2].",
    )
    p.add_argument(
        "--min-margin",
        type=int,
        default=1,
        help="(Non-EM mode) Minimum score margin over second-best lineage to assign [1].",
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
        "--use-em",
        action="store_true",
        help="Enable EM to estimate lineage mixture and soft-assign reads.",
    )
    p.add_argument(
        "--post-min",
        type=float,
        default=0.90,
        help="Posterior threshold to assign a read to a lineage when --use-em is on [0.90].",
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


# Barcode reading
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

    long_records = []
    pat = re.compile(r"^([ACGT])(\d+)([ACGT])$")

    for _, row in df.iterrows():
        lin = row["lineage"]
        if lineage_include and lin not in lineage_include:
            continue

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

    if not long_records:
        sys.exit("[ERROR] No lineage sites loaded (check --lineage-include or barcode CSV/TSV).")

    lineage_sites = defaultdict(dict)
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

    if not lineage_sites:
        sys.exit("[ERROR] No lineage sites loaded after processing barcode CSV/TSV.")

    site_to_lineage_alt = defaultdict(dict)
    for lin, sites in lineage_sites.items():
        for pos, (ref, alt) in sites.items():
            site_to_lineage_alt[pos][lin] = alt

    return lineage_sites, site_to_lineage_alt, sorted(all_positions), pos_to_ref


def iter_read_aligned_pairs(read, want_positions_set):
    for qpos, rpos in read.get_aligned_pairs(matches_only=False, with_seq=False):
        if rpos is None or qpos is None:
            continue
        pos1 = rpos + 1
        if pos1 in want_positions_set:
            yield pos1, qpos


def phred_to_err(q):
    return min(0.1, max(10 ** (-q / 10.0), 1e-6))


def logsumexp(vals):
    m = max(vals)
    return m + math.log(sum(math.exp(v - m) for v in vals))


def read_loglik_for_lineage(covered_obs, qual_by_pos, lineage_sites, lin, pos_to_ref):
    ll = 0.0
    lin_sites = lineage_sites[lin]
    for pos, b in covered_obs.items():
        ref = pos_to_ref.get(pos)
        if ref is None:
            continue

        if pos in lin_sites:
            exp = lin_sites[pos][1]
        else:
            exp = ref

        q = qual_by_pos.get(pos, 30)
        e = phred_to_err(q)

        if b == exp:
            ll += math.log(1 - e)
        elif b in "ACGT":
            ll += math.log(e)
        else:
            ll += math.log(0.5)
    return ll


# assignment per BAM
def assign_reads_on_bam(
    bam_path,
    outdir,
    prefix,
    args,
    lineage_sites,
    site_to_lineage_alt,
    all_positions,
    pos_to_ref,
):
    os.makedirs(outdir, exist_ok=True)
    want_positions = set(all_positions)

    bam = pysam.AlignmentFile(bam_path, "rb", threads=args.threads)
    ref_lengths = dict(zip(bam.references, bam.lengths))
    contigs = list(ref_lengths.keys())

    # Init counters/storage
    read_names = []
    per_read_cov_sites = []
    per_read_ll = []
    names_per_lineage = {lin: set() for lin in lineage_sites.keys()}
    ambiguous_names = set()

    read_count = 0
    dropped_mapq = 0
    sec_count = 0

    baseq_min = args.baseq_min
    ignore_ref_mismatch = args.ignore_ref_mismatch
    fetch_targets = contigs

    # prog bar
    total_mapped = bam.mapped if bam.mapped is not None else None
    pbar = None
    if tqdm is not None:
        pbar = tqdm(total=total_mapped, desc=f"{prefix}", unit="reads", leave=False)

    # compute votes or log-likelihoods
    for chrom in fetch_targets:
        for read in bam.fetch(chrom):
            read_count += 1
            if pbar is not None:
                pbar.update(1)

            if read.is_unmapped:
                continue
            if not args.allow_secondary and (read.is_secondary or read.is_supplementary):
                sec_count += 1
                continue
            if read.mapping_quality < args.mapq_min:
                dropped_mapq += 1
                continue

            seq = read.query_sequence
            quals = read.query_qualities
            covered_obs = {}
            qual_by_pos = {}
            seen_pos = set()
            for pos1, qpos in iter_read_aligned_pairs(read, want_positions):
                if pos1 in seen_pos:
                    continue
                seen_pos.add(pos1)
                if qpos < 0 or qpos >= len(seq):
                    continue
                qbase = seq[qpos].upper()
                qqual = quals[qpos] if quals is not None else 30
                if qbase == "N" or qqual < baseq_min:
                    continue
                covered_obs[pos1] = qbase
                qual_by_pos[pos1] = qqual

            # ambiguous
            if not covered_obs:
                ambiguous_names.add(read.query_name)
                continue

            read_names.append(read.query_name)

            if args.use_em:
                # log-likelihood per lineage
                ll_dict = {}
                for lin in lineage_sites.keys():
                    ll = read_loglik_for_lineage(
                        covered_obs, qual_by_pos, lineage_sites, lin, pos_to_ref
                    )
                    ll_dict[lin] = ll
                per_read_ll.append(ll_dict)
                per_read_cov_sites.append(len(covered_obs))
            else:
                # voting + margin approach
                votes = Counter()
                covered_sites_per_lin = Counter()

                for pos1, b in covered_obs.items():
                    ref_base = pos_to_ref.get(pos1, None)
                    lin2alt = site_to_lineage_alt.get(pos1, {})
                    for lin, alt in lin2alt.items():
                        covered_sites_per_lin[lin] += 1
                        if b == alt:
                            votes[lin] += 1
                        elif (
                            not ignore_ref_mismatch
                            and ref_base
                            and b == ref_base
                        ):
                            votes[lin] -= 1
                        else:
                            pass

                if not votes:
                    ambiguous_names.add(read.query_name)
                    continue

                best_lin, best_score = None, -10**9
                second_best = -10**9
                for lin in lineage_sites.keys():
                    sc = votes.get(lin, 0)
                    if sc > best_score:
                        second_best = best_score
                        best_score = sc
                        best_lin = lin
                    elif sc > second_best:
                        second_best = sc

                cov_best = covered_sites_per_lin.get(best_lin, 0)
                per_read_cov_sites.append(cov_best)

                if (
                    cov_best >= args.min_sites
                    and best_score >= second_best + args.min_margin
                ):
                    names_per_lineage[best_lin].add(read.query_name)
                else:
                    ambiguous_names.add(read.query_name)

    bam.close()
    if pbar is not None:
        pbar.close()

    # EM mixture & soft assignment
    gamma = None
    pi = None
    if args.use_em and read_names:
        L = list(lineage_sites.keys())
        R = len(per_read_ll)

        # uniform mixture
        pi = {lin: 1.0 / len(L) for lin in L}
        gamma = [{lin: 0.0 for lin in L} for _ in range(R)]

        # EM loop
        for _ in range(args.em_max_iters):
            for r in range(R):
                ll = per_read_ll[r]
                terms = [math.log(pi[lin]) + ll.get(lin, -1e9) for lin in L]
                denom = logsumexp(terms)
                for i, lin in enumerate(L):
                    gamma[r][lin] = math.exp(terms[i] - denom)

            pi_new = {lin: 0.0 for lin in L}
            for r in range(R):
                for lin in L:
                    pi_new[lin] += gamma[r][lin]
            for lin in L:
                pi_new[lin] /= R

            delta = sum(abs(pi_new[lin] - pi[lin]) for lin in L)
            pi = pi_new
            if delta < args.em_tol:
                break

        # hard membership sets using posterior cutoff + min-sites
        name_to_idx = {n: i for i, n in enumerate(read_names)}
        names_per_lineage = {lin: set() for lin in lineage_sites.keys()}
        ambiguous_names = set()
        for n, idx in name_to_idx.items():
            assigned = False
            if per_read_cov_sites[idx] < args.min_sites:
                ambiguous_names.add(n)
                continue
            for lin in L:
                if gamma[idx][lin] >= args.post_min:
                    names_per_lineage[lin].add(n)
                    assigned = True
                    break
            if not assigned:
                ambiguous_names.add(n)

    for lin, names in names_per_lineage.items():
        if names:
            with open(
                os.path.join(outdir, f"{prefix}.{lin}.names.txt"), "w"
            ) as out:
                for n in sorted(names):
                    out.write(n + "\n")
    if ambiguous_names:
        with open(
            os.path.join(outdir, f"{prefix}.ambiguous.names.txt"), "w"
        ) as out:
            for n in sorted(ambiguous_names):
                out.write(n + "\n")

    if args.write_bams or args.write_fastq or args.write_unassigned:
        src = pysam.AlignmentFile(bam_path, "rb", threads=args.threads)
        writers_bam = {}
        writers_fastq = {}
        created_bam_paths = set()

        def get_bam_writer(tag):
            if tag not in writers_bam:
                path = os.path.join(outdir, f"{prefix}.{tag}.bam")
                writers_bam[tag] = pysam.AlignmentFile(path, "wb", template=src)
                created_bam_paths.add(path)
            return writers_bam[tag]

        def get_fastq_writer(tag):
            if tag not in writers_fastq:
                path = os.path.join(outdir, f"{prefix}.{tag}.fastq")
                writers_fastq[tag] = open(path, "w")
            return writers_fastq[tag]

        mem = {lin: s for lin, s in names_per_lineage.items() if s}
        include_ambig = args.write_unassigned and len(ambiguous_names) > 0
        if include_ambig:
            mem["ambiguous"] = ambiguous_names

        for chrom in src.references:
            for read in src.fetch(chrom):
                if read.is_unmapped:
                    continue
                if not args.allow_secondary and (
                    read.is_secondary or read.is_supplementary
                ):
                    continue

                tag = None
                for lin in mem.keys():
                    if lin == "ambiguous":
                        continue
                    if read.query_name in mem[lin]:
                        tag = lin
                        break
                if tag is None and include_ambig and read.query_name in ambiguous_names:
                    tag = "ambiguous"
                if tag is None:
                    continue

                if args.write_bams:
                    get_bam_writer(tag).write(read)
                if args.write_fastq:
                    seq = read.query_sequence
                    qual = read.qual if read.qual else "~" * len(seq)
                    get_fastq_writer(tag).write(
                        f"@{read.query_name}\n{seq}\n+\n{qual}\n"
                    )

        if args.write_bams:
            for w in writers_bam.values():
                w.close()
            for path in created_bam_paths:
                pysam.index(path)
        if args.write_fastq:
            for fh in writers_fastq.values():
                fh.close()
        src.close()

    total_assigned = sum(len(v) for v in names_per_lineage.values())
    outsum = os.path.join(outdir, f"{prefix}.summary.txt")
    with open(outsum, "w") as out:
        out.write(f"reads_processed\t{read_count}\n")
        out.write(f"assigned_reads\t{total_assigned}\n")
        out.write(f"ambiguous_reads\t{len(ambiguous_names)}\n")
        out.write(f"skipped_low_mapq\t{dropped_mapq}\n")
        out.write(f"skipped_secondary\t{sec_count}\n")
        for lin, names in sorted(names_per_lineage.items()):
            out.write(f"lineage_{lin}\t{len(names)}\n")

    if args.write_mix_summary:
        mix_path = os.path.join(outdir, f"{prefix}.mix_summary.tsv")
        with open(mix_path, "w") as fh:
            fh.write("lineage\treads_or_weight\tfrac_assigned\tfrac_total\n")
            if args.use_em and read_names:
                L = list(lineage_sites.keys())
                name_to_idx = {n: i for i, n in enumerate(read_names)}
                weight_per_lin = {lin: 0.0 for lin in L}
                for idx in range(len(read_names)):
                    for lin in L:
                        weight_per_lin[lin] += gamma[idx][lin]
                total_weight = (
                    sum(weight_per_lin.values()) if weight_per_lin else 0.0
                )
                for lin in sorted(weight_per_lin.keys()):
                    w = weight_per_lin[lin]
                    fa = (w / total_weight) if total_weight > 0 else 0.0
                    ft = w / read_count if read_count > 0 else 0.0
                    fh.write(f"{lin}\t{w:.3f}\t{fa:.6f}\t{ft:.6f}\n")
                if len(ambiguous_names) > 0:
                    fh.write(
                        f"ambiguous\t{len(ambiguous_names)}\t0.000000\t"
                        f"{len(ambiguous_names)/read_count if read_count>0 else 0.0:.6f}\n"
                    )
            else:
                for lin, names in sorted(names_per_lineage.items()):
                    n = len(names)
                    if n == 0:
                        continue
                    fa = (n / total_assigned) if total_assigned > 0 else 0.0
                    ft = n / read_count if read_count > 0 else 0.0
                    fh.write(f"{lin}\t{n}\t{fa:.6f}\t{ft:.6f}\n")
                if len(ambiguous_names) > 0:
                    n = len(ambiguous_names)
                    fh.write(
                        f"ambiguous\t{n}\t0.000000\t"
                        f"{n/read_count if read_count>0 else 0.0:.6f}\n"
                    )

    return 0


def main():
    args = parse_args()

    lineage_sites, site_to_lineage_alt, all_positions, pos_to_ref = read_barcode_tsv(
        args.barcode_csv, args.lineage_include
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