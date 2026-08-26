import os
import math
from collections import Counter, defaultdict

import numpy as np
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


# per quality log-prob tables
_LOG_MATCH = [math.log(1.0 - phred_to_err(q)) for q in range(256)]
_LOG_MISMATCH = [math.log(phred_to_err(q) / 3.0) for q in range(256)]
_LOG_AMBIG = math.log(0.5)


def build_pos_alt_groups(site_to_lineage_alt, lin_index):
    # Index the barcode by position, then by ALT base.
    # Every lineage carrying the same ALT at the same position takes the same
    # likelihood correction so they can all be updated by one vectorised op.
    
    groups = {}
    for pos, lin2alt in site_to_lineage_alt.items():
        by_alt = defaultdict(list)
        for lin, alt in lin2alt.items():
            i = lin_index.get(lin)
            if i is not None:
                by_alt[alt].append(i)

        entries = []
        for alt, idxs in by_alt.items():
            if len(idxs) == 1:
                entries.append((alt, idxs[0]))
                
            else:
                idxs.sort()
                entries.append((alt, np.array(idxs, dtype=np.intp)))
        groups[pos] = entries
    return groups


def fill_read_loglik(row, covered_obs, qual_by_pos, pos_to_ref, pos_alt_groups):
    row.fill(0.0)
    ll_ref = 0.0

    for pos, b in covered_obs.items():
        ref = pos_to_ref.get(pos)
        if ref is None:
            continue

        q = qual_by_pos.get(pos, 30)
        t_match = _LOG_MATCH[q]
        t_miss = _LOG_MISMATCH[q] if b in "ACGT" else _LOG_AMBIG

        t_ref = t_match if b == ref else t_miss
        ll_ref += t_ref

        for alt, idx in pos_alt_groups.get(pos, ()):
            t_alt = t_match if b == alt else t_miss
            if t_alt != t_ref:
                row[idx] += t_alt - t_ref

    row += ll_ref


class LogLikMatrix:
    # The log-lik matrix, grown in fixed blocks

    def __init__(self, n_lineages, block_bytes=16 << 20):
        self.n_lineages = n_lineages
        self.n_rows = 0
        rows = block_bytes // max(1, n_lineages * 8)
        self._block_rows = int(min(4096, max(256, rows)))
        self._blocks = []
        self._fill = self._block_rows

    def next_row(self):
        if self._fill >= self._block_rows:
            self._blocks.append(
                np.empty((self._block_rows, self.n_lineages), dtype=np.float64)
            )
            self._fill = 0
        row = self._blocks[-1][self._fill]
        self._fill += 1
        self.n_rows += 1
        return row

    def iter_blocks(self):
        start = 0
        last = len(self._blocks) - 1
        for i, block in enumerate(self._blocks):
            n = self._fill if i == last else self._block_rows
            yield start, block[:n]
            start += n


def _log_pi(pi):
    # Mirrors the original scalar guard
    return np.where(pi > 0, np.log(np.where(pi > 0, pi, 1.0)), math.log(1e-300))


def _posterior(block, log_pi):

    terms = block + log_pi
    m = terms.max(axis=1, keepdims=True)
    shifted = terms - m
    np.exp(shifted, out=shifted)
    denom = m + np.log(shifted.sum(axis=1, keepdims=True))
    terms -= denom
    np.exp(terms, out=terms)
    return terms


def run_em(ll_matrix, max_iters, tol, post_min):
    # fit mixture weights, then hard-assign reads on the posterior cutoff.

    n_reads = ll_matrix.n_rows
    n_lin = ll_matrix.n_lineages

    pi = np.full(n_lin, 1.0 / n_lin, dtype=np.float64)
    pi_estep = None  # the pi that the most recent E-step actually used

    for _ in range(max_iters):
        log_pi = _log_pi(pi)
        pi_new = np.zeros(n_lin, dtype=np.float64)
        for _start, block in ll_matrix.iter_blocks():
            pi_new += _posterior(block, log_pi).sum(axis=0)
        pi_new /= n_reads

        pi_estep = pi
        delta = np.abs(pi_new - pi).sum()
        pi = pi_new
        if delta < tol:
            break

    best = np.full(n_reads, -1, dtype=np.int64)
    if pi_estep is not None:

        log_pi = _log_pi(pi_estep)
        for start, block in ll_matrix.iter_blocks():
            hit = _posterior(block, log_pi) >= post_min
            best[start:start + block.shape[0]] = np.where(
                hit.any(axis=1), hit.argmax(axis=1), -1
            )
    elif post_min <= 0.0:

        best[:] = 0

    return pi, best


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
    lineage_order = list(lineage_sites.keys())
    lin_index = {lin: i for i, lin in enumerate(lineage_order)}
    pos_alt_groups = (
        build_pos_alt_groups(site_to_lineage_alt, lin_index) if args.use_em else {}
    )

    read_names = []
    per_read_cov_sites = []
    per_read_ll = LogLikMatrix(len(lineage_order))
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
                # log-likelihood for every lineage, written straight into the
                # matrix so no per-read dict is ever built
                fill_read_loglik(
                    per_read_ll.next_row(),
                    covered_obs,
                    qual_by_pos,
                    pos_to_ref,
                    pos_alt_groups,
                )
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
    pi = None
    if args.use_em and read_names:
        pi_arr, best = run_em(
            per_read_ll, args.em_max_iters, args.em_tol, args.post_min
        )
        pi = {lin: float(pi_arr[i]) for i, lin in enumerate(lineage_order)}

        # hard membership sets using posterior cutoff
        name_to_idx = {n: i for i, n in enumerate(read_names)}
        names_per_lineage = {lin: set() for lin in lineage_sites.keys()}
        for n, idx in name_to_idx.items():
            j = int(best[idx])
            if j >= 0:
                names_per_lineage[lineage_order[j]].add(n)
            else:
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
