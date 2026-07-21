import os
import math
from collections import Counter

import pysam

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None


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
            ll += math.log(e / 3.0)
        else:
            ll += math.log(0.5)
    return ll


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

    total_reads = 0
    read_count = 0
    unmapped_count = 0
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
            total_reads += 1

            if read.is_unmapped:
                unmapped_count += 1
                continue
            if not args.allow_secondary and (read.is_secondary or read.is_supplementary):
                sec_count += 1
                continue
            if read.mapping_quality < args.mapq_min:
                dropped_mapq += 1
                continue

            read_count += 1
            if pbar is not None:
                pbar.update(1)

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

            if args.use_em:
                # per-read evidence floor (uniform across lineages, including
                # reference-equivalent rows)
                if len(covered_obs) < args.min_sites:
                    ambiguous_names.add(read.query_name)
                    continue
                # log-likelihood per lineage
                ll_dict = {}
                for lin in lineage_sites.keys():
                    ll = read_loglik_for_lineage(
                        covered_obs, qual_by_pos, lineage_sites, lin, pos_to_ref
                    )
                    ll_dict[lin] = ll
                per_read_ll.append(ll_dict)
                cov_per_lin = Counter()
                for pos1 in covered_obs:
                    for lin in site_to_lineage_alt.get(pos1, {}):
                        cov_per_lin[lin] += 1
                per_read_cov_sites.append(cov_per_lin)
                read_names.append(read.query_name)
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

                effective_min = min(args.min_sites, len(lineage_sites[best_lin]))
                if (
                    cov_best >= effective_min
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
        log_floor = math.log(1e-300)
        for _ in range(args.em_max_iters):
            for r in range(R):
                ll = per_read_ll[r]
                terms = [
                    (math.log(pi[lin]) if pi[lin] > 0 else log_floor)
                    + ll.get(lin, -1e9)
                    for lin in L
                ]
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

        # hard membership sets using posterior cutoff
        name_to_idx = {n: i for i, n in enumerate(read_names)}
        names_per_lineage = {lin: set() for lin in lineage_sites.keys()}
        for n, idx in name_to_idx.items():
            assigned = False
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
                    quals = read.query_qualities
                    if quals is not None:
                        qual = pysam.qualities_to_qualitystring(quals)
                    else:
                        qual = "~" * len(seq)
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
        out.write(f"reads_total\t{total_reads}\n")
        out.write(f"reads_processed\t{read_count}\n")
        out.write(f"assigned_reads\t{total_assigned}\n")
        out.write(f"ambiguous_reads\t{len(ambiguous_names)}\n")
        out.write(f"skipped_unmapped\t{unmapped_count}\n")
        out.write(f"skipped_low_mapq\t{dropped_mapq}\n")
        out.write(f"skipped_secondary\t{sec_count}\n")
        for lin, names in sorted(names_per_lineage.items()):
            out.write(f"lineage_{lin}\t{len(names)}\n")

    if args.write_mix_summary:
        mix_path = os.path.join(outdir, f"{prefix}.mix_summary.tsv")
        with open(mix_path, "w") as fh:
            fh.write("lineage\treads\tmix_weight\tfrac_assigned\tfrac_total\n")
            for lin, names in sorted(names_per_lineage.items()):
                n = len(names)
                if n == 0:
                    continue
                fa = (n / total_assigned) if total_assigned > 0 else 0.0
                ft = n / read_count if read_count > 0 else 0.0
                mw = f"{pi[lin]:.6f}" if pi is not None and lin in pi else ""
                fh.write(f"{lin}\t{n}\t{mw}\t{fa:.6f}\t{ft:.6f}\n")
            if len(ambiguous_names) > 0:
                n = len(ambiguous_names)
                fh.write(
                    f"ambiguous\t{n}\t\t0.000000\t"
                    f"{n/read_count if read_count>0 else 0.0:.6f}\n"
                )

    return 0
