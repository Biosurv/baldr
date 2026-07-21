"""Tests for read assignment and mixture estimation in EM mode."""

import json
import os

import pysam
import pytest

from conftest import DATA_DIR, data_path, run_assign, truth_by_read

with open(os.path.join(DATA_DIR, "manifest.json")) as _fh:
    MANIFEST = json.load(_fh)

SAMPLES = {sample["bam"]: sample for sample in MANIFEST["samples"]}

# Full-length single-end samples whose true composition is fully represented in
# barcode.csv. Paired samples, short reads, and the deliberately unmodelled
# sample are exercised separately.
FULL_LENGTH = [
    sample["bam"]
    for sample in MANIFEST["samples"]
    if sample.get("recoverable")
    and not sample["paired"]
    and "read_length" not in sample
]

MIX_TOLERANCE = 0.02
MIN_ACCURACY = 0.98


@pytest.mark.parametrize("bam", FULL_LENGTH)
def test_mix_weights_match_truth(bam, tmp_path):
    truth = SAMPLES[bam]["truth"]
    mix = run_assign(bam, tmp_path).mix
    for lineage, expected in truth.items():
        observed = mix.get(lineage, {}).get("mix_weight")
        assert observed is not None, f"{lineage} missing from mix summary"
        assert observed == pytest.approx(expected, abs=MIX_TOLERANCE)


@pytest.mark.parametrize("bam", FULL_LENGTH)
def test_lineages_absent_from_the_sample_stay_near_zero(bam, tmp_path):
    truth = SAMPLES[bam]["truth"]
    mix = run_assign(bam, tmp_path).mix
    for lineage, row in mix.items():
        if lineage in truth or lineage == "ambiguous":
            continue
        assert row["mix_weight"] == pytest.approx(0.0, abs=MIX_TOLERANCE)


@pytest.mark.parametrize("bam", FULL_LENGTH)
def test_assigned_reads_carry_their_true_lineage(bam, tmp_path):
    truth = truth_by_read(bam)
    assignments = run_assign(bam, tmp_path).assignments

    total = sum(len(names) for names in assignments.values())
    correct = sum(
        1
        for lineage, names in assignments.items()
        for name in names
        if truth[name] == lineage
    )
    assert total > 0
    assert correct / total >= MIN_ACCURACY


@pytest.mark.parametrize("bam", FULL_LENGTH)
def test_every_read_is_assigned_or_ambiguous_exactly_once(bam, tmp_path):
    truth = truth_by_read(bam)
    outputs = run_assign(bam, tmp_path)

    seen = list(outputs.ambiguous)
    for names in outputs.assignments.values():
        seen.extend(names)

    assert len(seen) == len(set(seen)), "a read name appears in more than one output"
    assert set(seen) == set(truth), "not every read landed in an output file"


@pytest.mark.parametrize("bam", FULL_LENGTH)
def test_summary_counts_are_internally_consistent(bam, tmp_path):
    outputs = run_assign(bam, tmp_path)
    summary = outputs.summary

    assert summary["assigned_reads"] + summary["ambiguous_reads"] == (
        summary["reads_processed"]
    )
    per_lineage = sum(
        value for key, value in summary.items() if key.startswith("lineage_")
    )
    assert per_lineage == summary["assigned_reads"]


def test_no_sites_sample_assigns_nothing(tmp_path):
    """Reads with no marker coverage carry no evidence and must not be called."""
    outputs = run_assign("no_sites.bam", tmp_path)
    assert outputs.assignments == {}
    assert len(outputs.ambiguous) == SAMPLES["no_sites.bam"]["n_reads"]


def test_reference_equivalent_lineage_absent_from_sample_stays_near_zero(tmp_path):
    """An all-zero barcode row must not absorb weight from a sample without it.
    """
    outputs = run_assign("mix_75_25.bam", tmp_path, barcode="barcode_edge.csv")
    weight = outputs.mix.get("REF.1", {}).get("mix_weight") or 0.0
    assert weight == pytest.approx(0.0, abs=MIX_TOLERANCE)


def test_short_reads_recover_the_mixture_less_precisely(tmp_path):
    """Short reads see few markers, so the estimate is looser but still right."""
    truth = SAMPLES["short_reads.bam"]["truth"]
    outputs = run_assign("short_reads.bam", tmp_path)

    for lineage, expected in truth.items():
        observed = outputs.mix[lineage]["mix_weight"]
        assert observed == pytest.approx(expected, abs=0.05)

    assert len(outputs.ambiguous) > 0


@pytest.mark.parametrize("min_sites", [0, 1, 2, 3])
def test_min_sites_never_increases_the_number_assigned(min_sites, tmp_path):
    outputs = run_assign(
        "short_reads.bam", tmp_path / str(min_sites), min_sites=min_sites
    )
    assigned = sum(len(names) for names in outputs.assignments.values())
    baseline = run_assign(
        "short_reads.bam", tmp_path / "baseline", min_sites=0
    )
    baseline_assigned = sum(
        len(names) for names in baseline.assignments.values()
    )
    assert assigned <= baseline_assigned


def test_min_sites_above_available_coverage_assigns_nothing(tmp_path):
    """No 60 bp read spans three markers, so nothing should survive the floor."""
    outputs = run_assign("short_reads.bam", tmp_path, min_sites=3)
    assert outputs.assignments == {}


def test_repeated_runs_produce_identical_output(tmp_path):
    first = run_assign("mix_three.bam", tmp_path / "first")
    second = run_assign("mix_three.bam", tmp_path / "second")
    assert first.assignments == second.assignments
    assert first.ambiguous == second.ambiguous
    assert first.mix == second.mix


def test_post_min_of_one_rejects_everything(tmp_path):
    """A posterior can never exceed 1, so nothing can clear this threshold."""
    outputs = run_assign("mix_75_25.bam", tmp_path, post_min=1.01)
    assert outputs.assignments == {}


def test_baseq_filter_discards_low_quality_bases(tmp_path):
    """Fixture bases are Q30, so a Q31 floor removes all marker evidence."""
    outputs = run_assign("mix_75_25.bam", tmp_path, baseq_min=31)
    assert outputs.assignments == {}
    assert len(outputs.ambiguous) == SAMPLES["mix_75_25.bam"]["n_reads"]


def test_mapq_filter_discards_reads_below_threshold(tmp_path):
    """Fixture reads are MAPQ 60, so a 61 floor drops every one of them."""
    outputs = run_assign("mix_75_25.bam", tmp_path, mapq_min=61)
    summary = outputs.summary
    assert summary["skipped_low_mapq"] == SAMPLES["mix_75_25.bam"]["n_reads"]
    assert summary["reads_processed"] == 0


def test_write_bams_emits_indexed_bams_matching_the_assignments(tmp_path):
    outputs = run_assign("mix_75_25.bam", tmp_path, write_bams=True)

    for lineage, names in outputs.assignments.items():
        path = os.path.join(str(tmp_path), f"sample.{lineage}.bam")
        assert os.path.exists(path)
        assert os.path.exists(path + ".bai"), f"{lineage} BAM was not indexed"
        with pysam.AlignmentFile(path, "rb") as handle:
            written = {read.query_name for read in handle.fetch(until_eof=True)}
        assert written == names


def test_write_fastq_emits_one_record_per_assigned_read(tmp_path):
    outputs = run_assign("mix_75_25.bam", tmp_path, write_fastq=True)

    for lineage, names in outputs.assignments.items():
        path = os.path.join(str(tmp_path), f"sample.{lineage}.fastq")
        with open(path) as fh:
            lines = fh.read().splitlines()
        assert len(lines) == 4 * len(names)
        written = {line[1:] for line in lines[::4]}
        assert written == names


def test_lineage_include_restricts_the_candidate_set(tmp_path):
    """Dropping B.1 from consideration forces its reads somewhere else."""
    from baldr.barcode import read_barcode_tsv
    from baldr.assign import assign_reads_on_bam
    from conftest import default_args

    lineage_sites, site_to_alt, positions, pos_to_ref = read_barcode_tsv(
        data_path("barcode.csv"), lineage_include=["A.1"]
    )
    assign_reads_on_bam(
        data_path("mix_75_25.bam"),
        str(tmp_path),
        "sample",
        default_args(),
        lineage_sites,
        site_to_alt,
        positions,
        pos_to_ref,
    )
    summary_path = os.path.join(str(tmp_path), "sample.summary.txt")
    with open(summary_path) as fh:
        text = fh.read()
    assert "lineage_B.1\t" not in text
