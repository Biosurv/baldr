# Shared fixtures and helpers for the BALDR test suite

import json
import os
from types import SimpleNamespace

import pysam
import pytest

from baldr.assign import assign_reads_on_bam
from baldr.barcode import read_barcode_tsv

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


@pytest.fixture(scope="session")
def manifest():
    with open(os.path.join(DATA_DIR, "manifest.json")) as fh:
        return json.load(fh)


@pytest.fixture(scope="session")
def reference_seq():
    lines = []
    with open(os.path.join(DATA_DIR, "reference.fasta")) as fh:
        for line in fh:
            if not line.startswith(">"):
                lines.append(line.strip())
                
    return "".join(lines)


def data_path(name):
    
    return os.path.join(DATA_DIR, name)


def default_args(**overrides):

    args = SimpleNamespace(
        mapq_min=0,
        baseq_min=7,
        min_sites=0,
        min_margin=1,
        ignore_ref_mismatch=False,
        allow_secondary=False,
        use_em=True,
        post_min=0.90,
        em_max_iters=50,
        em_tol=1e-4,
        threads=1,
        write_bams=False,
        write_fastq=False,
        write_unassigned=False,
        write_mix_summary=True,
    )
    
    
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


def run_assign(bam, outdir, barcode="barcode.csv", prefix="sample", **overrides):
    """Run assignment on a fixture BAM and return the parsed outputs."""
    lineage_sites, site_to_lineage_alt, all_positions, pos_to_ref = read_barcode_tsv(
        data_path(barcode)
    )
    args = default_args(**overrides)
    rc = assign_reads_on_bam(
        data_path(bam),
        str(outdir),
        prefix,
        args,
        lineage_sites,
        site_to_lineage_alt,
        all_positions,
        pos_to_ref,
    )
    return Outputs(str(outdir), prefix, rc)



class Outputs:
    """Parsed view of the files a single assignment run writes."""

    def __init__(self, outdir, prefix, rc):
        self.outdir = outdir
        self.prefix = prefix
        self.rc = rc

    def _path(self, suffix):
        return os.path.join(self.outdir, f"{self.prefix}.{suffix}")

    @property
    def summary(self):
        values = {}
        with open(self._path("summary.txt")) as fh:
            for line in fh:
                key, value = line.rstrip("\n").split("\t")
                values[key] = int(value)
        return values

    @property
    def mix(self):
        """lineage -> {reads, mix_weight, frac_assigned, frac_total}."""
        rows = {}
        path = self._path("mix_summary.tsv")
        if not os.path.exists(path):
            return rows
        with open(path) as fh:
            next(fh)
            for line in fh:
                lineage, reads, weight, frac_assigned, frac_total = line.rstrip(
                    "\n"
                ).split("\t")
                rows[lineage] = {
                    "reads": int(reads),
                    "mix_weight": float(weight) if weight else None,
                    "frac_assigned": float(frac_assigned),
                    "frac_total": float(frac_total),
                }
        return rows

    @property
    def assignments(self):
        """lineage -> set of read names assigned to it."""
        result = {}
        for filename in os.listdir(self.outdir):
            if not filename.endswith(".names.txt"):
                continue
            label = filename[len(self.prefix) + 1 : -len(".names.txt")]
            if label == "ambiguous":
                continue
            with open(os.path.join(self.outdir, filename)) as fh:
                result[label] = {line.strip() for line in fh if line.strip()}
        return result

    @property
    def ambiguous(self):
        path = self._path("ambiguous.names.txt")
        if not os.path.exists(path):
            return set()
            
        with open(path) as fh:
            return {line.strip() for line in fh if line.strip()}


def truth_by_read(bam):
    """read name -> true lineage, from the XL tag written by the generator."""
    truth = {}
    with pysam.AlignmentFile(data_path(bam), "rb") as handle:
        for read in handle.fetch(until_eof=True):
            truth[read.query_name] = read.get_tag("XL")
    return truth
