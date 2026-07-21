"""Tests for barcode table parsing."""

import pytest

from baldr.barcode import read_barcode_tsv
from conftest import data_path


@pytest.fixture(scope="module")
def parsed():
    return read_barcode_tsv(data_path("barcode.csv"))


def test_lineages_match_manifest(parsed, manifest):
    lineage_sites, _, _, _ = parsed
    assert set(lineage_sites) == set(manifest["lineage_markers"])


def test_marker_positions_match_manifest(parsed, manifest):
    lineage_sites, _, all_positions, _ = parsed
    assert all_positions == sorted(manifest["marker_positions"])
    for lineage, positions in manifest["lineage_markers"].items():
        assert sorted(lineage_sites[lineage]) == sorted(positions)


def test_ref_bases_match_the_reference_sequence(parsed, reference_seq):
    """The REF base encoded in each column name must match the actual reference.

    This is the check that catches a barcode built against different coordinates
    from the BAM, which otherwise produces confident but meaningless output.
    """
    _, _, _, pos_to_ref = parsed
    for pos, ref_base in pos_to_ref.items():
        assert reference_seq[pos - 1] == ref_base, (
            f"barcode says position {pos} is {ref_base} but the reference has "
            f"{reference_seq[pos - 1]}"
        )


def test_alt_always_differs_from_ref(parsed):
    lineage_sites, _, _, _ = parsed
    for lineage, sites in lineage_sites.items():
        for pos, (ref, alt) in sites.items():
            assert ref != alt, f"{lineage} at {pos} has ALT equal to REF"


def test_site_to_lineage_alt_inverts_lineage_sites(parsed):
    lineage_sites, site_to_lineage_alt, _, _ = parsed
    for lineage, sites in lineage_sites.items():
        for pos, (_, alt) in sites.items():
            assert site_to_lineage_alt[pos][lineage] == alt

    forward_pairs = {
        (lineage, pos)
        for lineage, sites in lineage_sites.items()
        for pos in sites
    }
    inverse_pairs = {
        (lineage, pos)
        for pos, mapping in site_to_lineage_alt.items()
        for lineage in mapping
    }
    assert forward_pairs == inverse_pairs


def test_lineage_include_filters_lineages():
    lineage_sites, _, _, _ = read_barcode_tsv(
        data_path("barcode.csv"), lineage_include=["A.1"]
    )
    assert set(lineage_sites) == {"A.1"}


def test_reference_equivalent_lineage_parses_with_no_sites():
    """An all-zero barcode row is a valid lineage carrying no defining sites."""
    lineage_sites, _, _, _ = read_barcode_tsv(data_path("barcode_edge.csv"))
    assert "REF.1" in lineage_sites
    assert lineage_sites["REF.1"] == {}


def test_identical_rows_produce_identical_site_maps():
    """B.1 and B.1.dup are indistinguishable from the barcode alone.

    Nothing downstream can separate them, so this pins the property at the
    parsing layer where it originates.
    """
    lineage_sites, _, _, _ = read_barcode_tsv(data_path("barcode_edge.csv"))
    assert lineage_sites["B.1"] == lineage_sites["B.1.dup"]


def test_conflicting_ref_bases_are_rejected(tmp_path):
    """Two columns claiming different REF bases at one position is unresolvable."""
    barcode = tmp_path / "conflict.csv"
    barcode.write_text("lineage,A30C,G30T\nX.1,1,0\nY.1,0,1\n")
    with pytest.raises(SystemExit):
        read_barcode_tsv(str(barcode))


def test_barcode_with_no_markers_is_rejected(tmp_path):
    """Without at least one marker there is no way to anchor any REF base."""
    barcode = tmp_path / "empty.csv"
    barcode.write_text("lineage,notamarker\nX.1,1\n")
    with pytest.raises(SystemExit):
        read_barcode_tsv(str(barcode))
