#!/usr/bin/env python

"""Tests for sccmec_type_classification.py

Unit tests for classifiers, lookup tables, mode switching and find_closest_ccr,
plus integration tests that require BLAST+.
"""

import shutil
import subprocess
import sys

import pytest
from pathlib import Path

from sccmecextractor.blast_utils import BlastResult
from sccmecextractor.sccmec_type_classification import (
    CcrClassifier,
    CcrComplexLookup,
    GeneHit,
    MecClassifier,
    MecComplexLookup,
    SCCmecTypeLookup,
    SCCmecTyper,
    TYPING_HEADER,
    collect_input_files,
)

HAS_BLAST = shutil.which("blastn") is not None


def _make_hit(qseqid, sseqid, pident, length, sstart=1, send=None, bitscore=1000):
    """Helper to create a BlastResult for testing."""
    if send is None:
        send = sstart + length - 1
    return BlastResult(
        qseqid=qseqid,
        sseqid=sseqid,
        pident=pident,
        length=length,
        mismatch=int(length * (100 - pident) / 100),
        gapopen=0,
        qstart=1,
        qend=length,
        sstart=sstart,
        send=send,
        evalue=0.0,
        bitscore=bitscore,
    )


def _make_gene_hit(name, pident=95.0, coverage=98.0, classification="full",
                   strand="+", contig="contig_1", start=1000, end=2000):
    """Helper to create a GeneHit for testing."""
    return GeneHit(
        gene_name=name, pident=pident, coverage=coverage,
        classification=classification, strand=strand,
        contig=contig, start=start, end=end,
    )


# ---------------------------------------------------------------------------
# TestMecClassifier
# ---------------------------------------------------------------------------

class TestMecClassifier:
    """Tests for MecClassifier."""

    @pytest.fixture
    def mec_classifier(self):
        """Create a MecClassifier with bundled reference."""
        from sccmecextractor.blast_utils import get_default_ref

        with get_default_ref("mec_class_reference.fasta") as ref:
            return MecClassifier(str(ref))

    def test_full_hit(self, mec_classifier):
        """High identity + high coverage = full."""
        ref_len = mec_classifier.ref_lengths.get("mecA", 2007)
        hit = _make_hit("mecA", "contig_1", 98.0, int(ref_len * 0.95))

        results = mec_classifier.classify([hit])
        assert len(results) == 1
        assert results[0].classification == "full"
        assert results[0].gene_name == "mecA"

    def test_confirmed_partial(self, mec_classifier):
        """High identity + moderate coverage = confirmed partial."""
        ref_len = mec_classifier.ref_lengths.get("mecA", 2007)
        hit = _make_hit("mecA", "contig_1", 97.0, int(ref_len * 0.80))

        results = mec_classifier.classify([hit])
        assert len(results) == 1
        assert results[0].classification == "partial"

    def test_novel_full(self, mec_classifier):
        """Moderate identity + high coverage = novel full."""
        ref_len = mec_classifier.ref_lengths.get("mecA", 2007)
        hit = _make_hit("mecA", "contig_1", 80.0, int(ref_len * 0.95))

        results = mec_classifier.classify([hit])
        assert len(results) == 1
        assert results[0].classification == "novel_full"

    def test_novel_partial(self, mec_classifier):
        """Moderate identity + moderate coverage = novel partial."""
        ref_len = mec_classifier.ref_lengths.get("mecA", 2007)
        hit = _make_hit("mecA", "contig_1", 80.0, int(ref_len * 0.80))

        results = mec_classifier.classify([hit])
        assert len(results) == 1
        assert results[0].classification == "novel_partial"

    def test_below_threshold_excluded(self, mec_classifier):
        """Hits below minimum thresholds are excluded."""
        hit = _make_hit("mecA", "contig_1", 60.0, 500)

        results = mec_classifier.classify([hit])
        assert len(results) == 0

    def test_no_hits(self, mec_classifier):
        """Empty hits return empty results."""
        results = mec_classifier.classify([])
        assert results == []

    def test_bundled_ref_contains_structural_genes(self, mec_classifier):
        """Bundled mec_class_reference.fasta includes IS and regulatory genes."""
        expected = {"IS1272", "IS431", "mecI", "mecR1", "mecA", "mecC", "mecB", "mecD"}
        assert expected.issubset(set(mec_classifier.ref_lengths.keys()))


# ---------------------------------------------------------------------------
# TestCcrClassifier
# ---------------------------------------------------------------------------

class TestCcrClassifier:
    """Tests for CcrClassifier."""

    @pytest.fixture
    def ccr_classifier(self):
        """Create a CcrClassifier with bundled reference."""
        from sccmecextractor.blast_utils import get_default_ref

        with get_default_ref("ccr_genes.fasta") as ref:
            return CcrClassifier(str(ref))

    def test_confirmed_full(self, ccr_classifier):
        """High identity + high coverage = confirmed full."""
        ref_len = ccr_classifier.ref_lengths.get("ccrA1", 1350)
        hit = _make_hit("ccrA1", "contig_1", 95.0, int(ref_len * 0.95))

        results = ccr_classifier.classify([hit])
        assert len(results) == 1
        assert results[0].classification == "full"

    def test_confirmed_partial(self, ccr_classifier):
        """High identity + moderate coverage = confirmed partial."""
        ref_len = ccr_classifier.ref_lengths.get("ccrA1", 1350)
        hit = _make_hit("ccrA1", "contig_1", 90.0, int(ref_len * 0.80))

        results = ccr_classifier.classify([hit])
        assert len(results) == 1
        assert results[0].classification == "partial"

    def test_novel_full(self, ccr_classifier):
        """Moderate identity + high coverage = novel full."""
        ref_len = ccr_classifier.ref_lengths.get("ccrA1", 1350)
        hit = _make_hit("ccrA1", "contig_1", 75.0, int(ref_len * 0.95))

        results = ccr_classifier.classify([hit])
        assert len(results) == 1
        assert results[0].classification == "novel_full"

    def test_novel_partial(self, ccr_classifier):
        """Moderate identity + moderate coverage = novel partial."""
        ref_len = ccr_classifier.ref_lengths.get("ccrA1", 1350)
        hit = _make_hit("ccrA1", "contig_1", 75.0, int(ref_len * 0.80))

        results = ccr_classifier.classify([hit])
        assert len(results) == 1
        assert results[0].classification == "novel_partial"

    def test_below_threshold_excluded(self, ccr_classifier):
        """Hits below minimum thresholds are excluded."""
        hit = _make_hit("ccrA1", "contig_1", 60.0, 500)

        results = ccr_classifier.classify([hit])
        assert len(results) == 0

    def test_overlapping_hits_resolved(self, ccr_classifier):
        """Overlapping hits on same contig resolved by bitscore."""
        ref_len = ccr_classifier.ref_lengths.get("ccrA1", 1350)
        hit1 = _make_hit(
            "ccrA1", "contig_1", 95.0, int(ref_len * 0.95),
            sstart=1000, send=1000 + int(ref_len * 0.95), bitscore=2000,
        )
        hit2 = _make_hit(
            "ccrA2", "contig_1", 90.0, int(ref_len * 0.95),
            sstart=1050, send=1050 + int(ref_len * 0.95), bitscore=1800,
        )

        results = ccr_classifier.classify([hit1, hit2])
        gene_names = [r.gene_name for r in results]
        assert "ccrA1" in gene_names
        assert "ccrA2" not in gene_names


# ---------------------------------------------------------------------------
# TestMecComplexLookup
# ---------------------------------------------------------------------------

class TestMecComplexLookup:
    """Tests for mec complex class assignment (A-E)."""

    def test_class_a_meca_meci(self):
        """mecA + mecI = Class A."""
        hits = [
            _make_gene_hit("mecA", start=5000, end=7000),
            _make_gene_hit("mecI", start=7500, end=8000),
            _make_gene_hit("IS431", start=3000, end=3800),
        ]
        assert MecComplexLookup.lookup(hits) == "A"

    def test_class_b_meca_is1272(self):
        """mecA + IS1272 (no mecI) = Class B."""
        hits = [
            _make_gene_hit("mecA", start=5000, end=7000),
            _make_gene_hit("IS1272", start=7500, end=9000),
            _make_gene_hit("IS431", start=3000, end=3800),
        ]
        assert MecComplexLookup.lookup(hits) == "B"

    def test_class_c1_same_orientation(self):
        """mecA + flanking IS431 same orientation = Class C1."""
        hits = [
            _make_gene_hit("mecA", start=5000, end=7000),
            _make_gene_hit("IS431", strand="+", start=3000, end=3800),
            _make_gene_hit("IS431", strand="+", start=8000, end=8800),
        ]
        assert MecComplexLookup.lookup(hits) == "C1"

    def test_class_c2_opposite_orientation(self):
        """mecA + flanking IS431 opposite orientation = Class C2."""
        hits = [
            _make_gene_hit("mecA", start=5000, end=7000),
            _make_gene_hit("IS431", strand="+", start=3000, end=3800),
            _make_gene_hit("IS431", strand="-", start=8000, end=8800),
        ]
        assert MecComplexLookup.lookup(hits) == "C2"

    def test_class_d_single_is431(self):
        """mecA + single upstream IS431 (no mecI, no IS1272) = Class D."""
        hits = [
            _make_gene_hit("mecA", start=5000, end=7000),
            _make_gene_hit("IS431", start=3000, end=3800),
            _make_gene_hit("mecR1", start=7100, end=8500),
        ]
        assert MecComplexLookup.lookup(hits) == "D"

    def test_class_e_mecc_meci(self):
        """mecC + mecI = Class E."""
        hits = [
            _make_gene_hit("mecC", start=5000, end=7000),
            _make_gene_hit("mecI", start=7500, end=8000),
        ]
        assert MecComplexLookup.lookup(hits) == "E"

    def test_no_mec_returns_dash(self):
        """No mec genes returns '-'."""
        hits = [_make_gene_hit("IS431", start=1000, end=1800)]
        assert MecComplexLookup.lookup(hits) == "-"

    def test_empty_returns_dash(self):
        """Empty results returns '-'."""
        assert MecComplexLookup.lookup([]) == "-"

    def test_meca_only_not_typeable(self):
        """mecA alone (no IS elements nearby) = not_typeable."""
        hits = [_make_gene_hit("mecA", start=5000, end=7000)]
        assert MecComplexLookup.lookup(hits) == "not_typeable"

    def test_is_distant_from_meca_ignored(self):
        """IS elements >10kb from mecA are not counted for class."""
        hits = [
            _make_gene_hit("mecA", start=5000, end=7000),
            # IS1272 is 20kb away — should be ignored
            _make_gene_hit("IS1272", start=27000, end=28500),
        ]
        assert MecComplexLookup.lookup(hits) == "not_typeable"

    def test_is_on_different_contig_ignored(self):
        """IS elements on different contig from mecA are not counted."""
        hits = [
            _make_gene_hit("mecA", contig="contig_1", start=5000, end=7000),
            _make_gene_hit("IS1272", contig="contig_2", start=5500, end=7000),
            _make_gene_hit("IS431", contig="contig_2", start=3000, end=3800),
        ]
        assert MecComplexLookup.lookup(hits) == "not_typeable"

    def test_mecb_detected_but_not_typeable(self):
        """mecB is detected in output but has no class rule."""
        hits = [_make_gene_hit("mecB", start=5000, end=7000)]
        assert MecComplexLookup.lookup(hits) == "not_typeable"

    def test_mecd_detected_but_not_typeable(self):
        """mecD is detected in output but has no class rule."""
        hits = [_make_gene_hit("mecD", start=5000, end=7000)]
        assert MecComplexLookup.lookup(hits) == "not_typeable"


# ---------------------------------------------------------------------------
# TestSCCmecTypeLookup
# ---------------------------------------------------------------------------

class TestSCCmecTypeLookup:
    """Tests for SCCmec type assignment from (mec_class, ccr_complex) pairs."""

    @pytest.mark.parametrize("mec_class,ccr_type,expected", [
        ("B", "1", "I"),
        ("A", "2", "II"),
        ("A", "3", "III"),
        ("B", "2", "IV"),
        ("C2", "5", "V"),
        ("B", "4", "VI"),
        ("C1", "5", "VII"),
        ("A", "4", "VIII"),
        ("C2", "1", "IX"),
        ("C1", "7", "X"),
        ("E", "8", "XI"),
        ("C2", "9", "XII"),
        ("A", "9", "XIII"),
        ("A", "5", "XIV"),
        ("A", "7", "XV"),
    ])
    def test_known_types(self, mec_class, ccr_type, expected):
        """All 15 known SCCmec types map correctly."""
        assert SCCmecTypeLookup.lookup(mec_class, ccr_type) == expected

    def test_novel_combination(self):
        """Unrecognised (mec_class, ccr_type) pair returns novel_combination."""
        assert SCCmecTypeLookup.lookup("B", "5") == "novel_combination"

    def test_no_mec_returns_not_typeable(self):
        """Missing mec class returns not_typeable."""
        assert SCCmecTypeLookup.lookup("-", "2") == "not_typeable"

    def test_no_ccr_returns_not_typeable(self):
        """Missing ccr complex returns not_typeable."""
        assert SCCmecTypeLookup.lookup("A", "-") == "not_typeable"

    def test_composite_ccr_resolves(self):
        """Composite ccr type (e.g. '2;5') tries each part against mec class."""
        # mec class B + ccr "2;5" → Type IV from (B,2)
        result = SCCmecTypeLookup.lookup("B", "2;5")
        assert "IV" in result

    def test_composite_no_match(self):
        """Composite ccr with no matching part returns novel_combination."""
        assert SCCmecTypeLookup.lookup("E", "2;5") == "novel_combination"


# ---------------------------------------------------------------------------
# TestCcrComplexLookup
# ---------------------------------------------------------------------------

class TestCcrComplexLookup:
    """Tests for CcrComplexLookup."""

    def test_simple_pair_type2(self):
        """ccrA2 + ccrB2 = complex type 2."""
        hits = [
            _make_gene_hit("ccrA2"),
            _make_gene_hit("ccrB2"),
        ]
        assert CcrComplexLookup.lookup(hits) == "2"

    def test_single_ccrC1_type5(self):
        """ccrC1 alone = complex type 5."""
        hits = [_make_gene_hit("ccrC1")]
        assert CcrComplexLookup.lookup(hits) == "5"

    def test_mixed_pair_type8(self):
        """ccrA1 + ccrB3 = complex type 8."""
        hits = [
            _make_gene_hit("ccrA1"),
            _make_gene_hit("ccrB3"),
        ]
        assert CcrComplexLookup.lookup(hits) == "8"

    def test_composite_types_2_and_5(self):
        """ccrA2 + ccrB2 + ccrC1 = composite types 2;5."""
        hits = [
            _make_gene_hit("ccrA2"),
            _make_gene_hit("ccrB2"),
            _make_gene_hit("ccrC1"),
        ]
        assert CcrComplexLookup.lookup(hits) == "2;5"

    def test_triple_composite_types_1_5_9(self):
        """ccrA1 + ccrB1 + ccrC1 + ccrC2 = composite types 1;5;9."""
        hits = [
            _make_gene_hit("ccrA1"),
            _make_gene_hit("ccrB1"),
            _make_gene_hit("ccrC1"),
            _make_gene_hit("ccrC2"),
        ]
        assert CcrComplexLookup.lookup(hits) == "1;5;9"

    def test_novel_genes_flagged(self):
        """Novel classification appends '(novel_genes)' to result."""
        hits = [
            _make_gene_hit("ccrA2", classification="novel_full"),
            _make_gene_hit("ccrB2"),
        ]
        assert CcrComplexLookup.lookup(hits) == "2 (novel_genes)"

    def test_unknown_combination(self):
        """Unrecognised allotype combination returns 'novel_combination'."""
        hits = [_make_gene_hit("ccrB9")]
        assert CcrComplexLookup.lookup(hits) == "novel_combination"

    def test_no_hits(self):
        """No ccr hits returns '-'."""
        assert CcrComplexLookup.lookup([]) == "-"

    # --- New complex types (10-22) ---

    def test_type10_ccrA8B9(self):
        """ccrA8 + ccrB9 = complex type 10 (Xiao et al. 2023)."""
        hits = [_make_gene_hit("ccrA8"), _make_gene_hit("ccrB9")]
        assert CcrComplexLookup.lookup(hits) == "10"

    def test_type11_ccrA9B3(self):
        """ccrA9 + ccrB3 = complex type 11 (Huang et al. 2024)."""
        hits = [_make_gene_hit("ccrA9"), _make_gene_hit("ccrB3")]
        assert CcrComplexLookup.lookup(hits) == "11"

    def test_type12_ccrA10B1(self):
        """ccrA10 + ccrB1 = complex type 12 (Huang et al. 2024)."""
        hits = [_make_gene_hit("ccrA10"), _make_gene_hit("ccrB1")]
        assert CcrComplexLookup.lookup(hits) == "12"

    def test_type13_ccrA10B10(self):
        """ccrA10 + ccrB10 = complex type 13 (Huang et al. 2024)."""
        hits = [_make_gene_hit("ccrA10"), _make_gene_hit("ccrB10")]
        assert CcrComplexLookup.lookup(hits) == "13"

    def test_type20_ccrC3(self):
        """ccrC3 alone = complex type 20 (Huang et al. 2024)."""
        hits = [_make_gene_hit("ccrC3")]
        assert CcrComplexLookup.lookup(hits) == "20"

    def test_type21_ccrC4(self):
        """ccrC4 alone = complex type 21 (Huang et al. 2024)."""
        hits = [_make_gene_hit("ccrC4")]
        assert CcrComplexLookup.lookup(hits) == "21"

    def test_type22_ccrC5(self):
        """ccrC5 alone = complex type 22 (Huang et al. 2024)."""
        hits = [_make_gene_hit("ccrC5")]
        assert CcrComplexLookup.lookup(hits) == "22"

    def test_type13_preferred_over_type12(self):
        """When ccrA10 + ccrB10 + ccrB1 present, matched pair (13) takes priority."""
        hits = [
            _make_gene_hit("ccrA10"),
            _make_gene_hit("ccrB10"),
            _make_gene_hit("ccrB1"),
        ]
        result = CcrComplexLookup.lookup(hits)
        assert result == "13"

    def test_composite_type11_and_type5(self):
        """ccrA9 + ccrB3 + ccrC1 = composite types 11;5."""
        hits = [
            _make_gene_hit("ccrA9"),
            _make_gene_hit("ccrB3"),
            _make_gene_hit("ccrC1"),
        ]
        assert CcrComplexLookup.lookup(hits) == "11;5"


# ---------------------------------------------------------------------------
# TestFindClosestCcr
# ---------------------------------------------------------------------------

class TestFindClosestCcr:
    """Tests for SCCmecTyper.find_closest_ccr()."""

    def test_single_ccr_on_same_contig(self):
        """Single ccr on same contig as mecA is returned."""
        mec = [_make_gene_hit("mecA", contig="c1", start=5000, end=7000)]
        ccr = [_make_gene_hit("ccrA2", contig="c1", start=8000, end=9350)]
        result = SCCmecTyper.find_closest_ccr(mec, ccr)
        assert len(result) == 1
        assert result[0].gene_name == "ccrA2"

    def test_pair_partner_included(self):
        """ccrA/ccrB pair partner within 5kb is included."""
        mec = [_make_gene_hit("mecA", contig="c1", start=5000, end=7000)]
        ccr = [
            _make_gene_hit("ccrA2", contig="c1", start=8000, end=9350),
            _make_gene_hit("ccrB2", contig="c1", start=9400, end=10750),
        ]
        result = SCCmecTyper.find_closest_ccr(mec, ccr)
        names = {r.gene_name for r in result}
        assert names == {"ccrA2", "ccrB2"}

    def test_ccrC_returned_alone(self):
        """ccrC is returned as a single gene (unpaired)."""
        mec = [_make_gene_hit("mecA", contig="c1", start=5000, end=7000)]
        ccr = [_make_gene_hit("ccrC1", contig="c1", start=8000, end=9000)]
        result = SCCmecTyper.find_closest_ccr(mec, ccr)
        assert len(result) == 1
        assert result[0].gene_name == "ccrC1"

    def test_closest_selected_from_multiple(self):
        """When multiple ccr on same contig, closest to mecA is selected."""
        mec = [_make_gene_hit("mecA", contig="c1", start=5000, end=7000)]
        ccr = [
            _make_gene_hit("ccrC1", contig="c1", start=20000, end=21000),
            _make_gene_hit("ccrA2", contig="c1", start=8000, end=9350),
            _make_gene_hit("ccrB2", contig="c1", start=9400, end=10750),
        ]
        result = SCCmecTyper.find_closest_ccr(mec, ccr)
        names = {r.gene_name for r in result}
        assert "ccrA2" in names
        assert "ccrC1" not in names

    def test_fallback_when_no_same_contig(self):
        """Returns all ccr when none on same contig as mecA."""
        mec = [_make_gene_hit("mecA", contig="c1", start=5000, end=7000)]
        ccr = [
            _make_gene_hit("ccrA2", contig="c2", start=8000, end=9350),
            _make_gene_hit("ccrB2", contig="c2", start=9400, end=10750),
        ]
        result = SCCmecTyper.find_closest_ccr(mec, ccr)
        assert len(result) == 2  # All returned as fallback

    def test_no_mec_returns_all_ccr(self):
        """With no mecA, all ccr results are returned."""
        mec = []
        ccr = [_make_gene_hit("ccrA2"), _make_gene_hit("ccrB2")]
        result = SCCmecTyper.find_closest_ccr(mec, ccr)
        assert len(result) == 2

    def test_no_ccr_returns_empty(self):
        """With no ccr, empty list returned."""
        mec = [_make_gene_hit("mecA")]
        result = SCCmecTyper.find_closest_ccr(mec, [])
        assert result == []

    def test_mecC_used_when_no_mecA(self):
        """mecC is used for proximity when mecA is absent."""
        mec = [_make_gene_hit("mecC", contig="c1", start=5000, end=7000)]
        ccr = [
            _make_gene_hit("ccrA1", contig="c1", start=8000, end=9350),
            _make_gene_hit("ccrB1", contig="c1", start=9400, end=10750),
            _make_gene_hit("ccrC1", contig="c1", start=25000, end=26000),
        ]
        result = SCCmecTyper.find_closest_ccr(mec, ccr)
        names = {r.gene_name for r in result}
        assert "ccrA1" in names
        assert "ccrC1" not in names

    def test_mecA_preferred_over_mecC(self):
        """When both mecA and mecC present, mecA is used for proximity."""
        mec = [
            _make_gene_hit("mecA", contig="c1", start=5000, end=7000),
            _make_gene_hit("mecC", contig="c1", start=30000, end=32000),
        ]
        ccr = [
            _make_gene_hit("ccrA2", contig="c1", start=8000, end=9350),
            _make_gene_hit("ccrB2", contig="c1", start=9400, end=10750),
            _make_gene_hit("ccrC1", contig="c1", start=31000, end=32000),
        ]
        result = SCCmecTyper.find_closest_ccr(mec, ccr)
        names = {r.gene_name for r in result}
        # ccrA2/ccrB2 are closest to mecA, not ccrC1 which is closest to mecC
        assert "ccrA2" in names
        assert "ccrC1" not in names


# ---------------------------------------------------------------------------
# TestGeneContentMode
# ---------------------------------------------------------------------------

class TestGeneContentMode:
    """Tests for the gene content mode (custom --mec-ref / --ccr-ref)."""

    def test_default_mode_produces_typing(self):
        """No custom refs → full typing columns populated."""
        typer = SCCmecTyper()
        mec = [_make_gene_hit("mecA", start=5000, end=7000)]
        ccr = [
            _make_gene_hit("ccrA1", start=8000, end=9350),
            _make_gene_hit("ccrB1", start=9400, end=10750),
        ]
        result = typer._format_result("test", mec, ccr)
        assert result["mec_class_type"] != "-"
        assert result["ccr_complex_type"] == "1"
        assert result["SCCmec_Type"] != "-"

    def test_custom_mec_skips_mec_typing(self, tmp_path):
        """Custom --mec-ref → mec_class_type and SCCmec_Type are '-'."""
        # Create a minimal custom ref
        ref = tmp_path / "custom_mec.fasta"
        ref.write_text(">my_gene\nATCGATCG\n")
        typer = SCCmecTyper(mec_ref=str(ref))

        assert typer._custom_mec is True
        assert typer._custom_ccr is False

        mec = [_make_gene_hit("my_gene")]
        ccr = [
            _make_gene_hit("ccrA2", start=8000, end=9350),
            _make_gene_hit("ccrB2", start=9400, end=10750),
        ]
        result = typer._format_result("test", mec, ccr)
        assert result["mec_class_type"] == "-"
        assert result["ccr_complex_type"] == "2"  # ccr still typed
        assert result["SCCmec_Type"] == "-"  # can't type without mec class

    def test_custom_ccr_skips_ccr_typing(self, tmp_path):
        """Custom --ccr-ref → ccr_complex_type and SCCmec_Type are '-'."""
        ref = tmp_path / "custom_ccr.fasta"
        ref.write_text(">my_ccr\nATCGATCG\n")
        typer = SCCmecTyper(ccr_ref=str(ref))

        assert typer._custom_mec is False
        assert typer._custom_ccr is True

        mec = [
            _make_gene_hit("mecA", start=5000, end=7000),
            _make_gene_hit("IS1272", start=7500, end=9000),
            _make_gene_hit("IS431", start=3000, end=3800),
        ]
        ccr = [_make_gene_hit("my_ccr")]
        result = typer._format_result("test", mec, ccr)
        assert result["mec_class_type"] == "B"  # mec still typed
        assert result["ccr_complex_type"] == "-"
        assert result["SCCmec_Type"] == "-"

    def test_both_custom_all_typing_skipped(self, tmp_path):
        """Both custom refs → all typing columns are '-'."""
        mec_ref = tmp_path / "custom_mec.fasta"
        mec_ref.write_text(">my_mec\nATCGATCG\n")
        ccr_ref = tmp_path / "custom_ccr.fasta"
        ccr_ref.write_text(">my_ccr\nATCGATCG\n")
        typer = SCCmecTyper(mec_ref=str(mec_ref), ccr_ref=str(ccr_ref))

        result = typer._format_result(
            "test",
            [_make_gene_hit("my_mec")],
            [_make_gene_hit("my_ccr")],
        )
        assert result["mec_class_type"] == "-"
        assert result["ccr_complex_type"] == "-"
        assert result["SCCmec_Type"] == "-"
        # Gene content still reported
        assert "my_mec" in result["mec_genes"]
        assert "my_ccr" in result["ccr_genes"]


# ---------------------------------------------------------------------------
# TestOutputFormat
# ---------------------------------------------------------------------------

class TestOutputFormat:
    """Tests for SCCmecTyper output formatting."""

    def test_format_with_results(self):
        """Formatting with mec and ccr results produces correct fields."""
        typer = SCCmecTyper()
        result = typer._format_result(
            "test_genome",
            [_make_gene_hit("mecA", start=5000, end=7000)],
            [
                _make_gene_hit("ccrA2", start=8000, end=9350),
                _make_gene_hit("ccrB2", start=9400, end=10750),
            ],
        )

        assert result["Input_File"] == "test_genome"
        assert "mecA(full)" in result["mec_genes"]
        assert "ccrA2(full)" in result["ccr_genes"]
        assert "ccrA2" in result["ccr_allotypes"]
        assert "ccrB2" in result["ccr_allotypes"]
        assert result["ccr_complex_type"] == "2"

    def test_format_no_hits(self):
        """Formatting with no hits produces dashes."""
        typer = SCCmecTyper()
        result = typer._format_result("test_genome", [], [])

        assert result["mec_genes"] == "-"
        assert result["ccr_genes"] == "-"
        assert result["ccr_allotypes"] == "-"
        assert result["ccr_complex_type"] == "-"
        assert result["mec_class_type"] == "-"
        assert result["SCCmec_Type"] == "not_typeable"

    def test_all_typing_header_keys_present(self):
        """Result dict contains all TYPING_HEADER keys."""
        typer = SCCmecTyper()
        result = typer._format_result("test", [], [])
        for key in TYPING_HEADER:
            assert key in result

    def test_sccmec_type_assignment(self):
        """Full typing produces correct SCCmec type."""
        typer = SCCmecTyper()
        # mecA + IS1272 (no mecI) → Class B; ccrA1+ccrB1 → ccr type 1
        # (B, 1) → SCCmec Type I
        result = typer._format_result(
            "test",
            [
                _make_gene_hit("mecA", start=5000, end=7000),
                _make_gene_hit("IS1272", start=7500, end=9000),
                _make_gene_hit("IS431", start=3000, end=3800),
            ],
            [
                _make_gene_hit("ccrA1", start=10000, end=11350),
                _make_gene_hit("ccrB1", start=11400, end=12750),
            ],
        )
        assert result["mec_class_type"] == "B"
        assert result["ccr_complex_type"] == "1"
        assert result["SCCmec_Type"] == "I"

    def test_multiple_ccr_primary_and_secondary(self):
        """Multiple ccr complexes produce primary and secondary SCCmec types."""
        typer = SCCmecTyper()
        # mecA + IS1272 → Class B; ccrA1+ccrB1 closest to mecA → Type I (B,1)
        # ccrC1 further away → secondary ccr type 5 → (B,5) = novel_combination
        result = typer._format_result(
            "test",
            [
                _make_gene_hit("mecA", start=5000, end=7000),
                _make_gene_hit("IS1272", start=7500, end=9000),
                _make_gene_hit("IS431", start=3000, end=3800),
            ],
            [
                _make_gene_hit("ccrA1", start=10000, end=11350),
                _make_gene_hit("ccrB1", start=11400, end=12750),
                _make_gene_hit("ccrC1", start=25000, end=26000),
            ],
        )
        assert result["SCCmec_Type"] == "I"
        assert result["SCCmec_Type_secondary"] != "-"
        assert "(composite)" not in result["SCCmec_Type"]

    def test_single_ccr_no_secondary(self):
        """Single ccr complex produces '-' for secondary."""
        typer = SCCmecTyper()
        result = typer._format_result(
            "test",
            [
                _make_gene_hit("mecA", start=5000, end=7000),
                _make_gene_hit("IS1272", start=7500, end=9000),
                _make_gene_hit("IS431", start=3000, end=3800),
            ],
            [
                _make_gene_hit("ccrA1", start=10000, end=11350),
                _make_gene_hit("ccrB1", start=11400, end=12750),
            ],
        )
        assert result["SCCmec_Type"] == "I"
        assert result["SCCmec_Type_secondary"] == "-"

    def test_secondary_known_type(self):
        """Secondary ccr maps to a known SCCmec type."""
        typer = SCCmecTyper()
        # mec class A; ccrA3+ccrB3 closest → Type III (A,3)
        # ccrC1 further → secondary ccr type 5 → (A,5) = Type XIV
        result = typer._format_result(
            "test",
            [
                _make_gene_hit("mecA", start=5000, end=7000),
                _make_gene_hit("mecI", start=7500, end=8000),
                _make_gene_hit("IS431", start=3000, end=3800),
            ],
            [
                _make_gene_hit("ccrA3", start=10000, end=11350),
                _make_gene_hit("ccrB3", start=11400, end=12750),
                _make_gene_hit("ccrC1", start=25000, end=26000),
            ],
        )
        assert result["SCCmec_Type"] == "III"
        assert result["SCCmec_Type_secondary"] == "XIV"

    def test_secondary_novel_combination(self):
        """Secondary ccr that doesn't map shows novel_combination."""
        typer = SCCmecTyper()
        # mec class B; ccrA2+ccrB2 closest → Type IV (B,2)
        # ccrC1 further → secondary ccr type 5 → (B,5) = novel_combination
        result = typer._format_result(
            "test",
            [
                _make_gene_hit("mecA", start=5000, end=7000),
                _make_gene_hit("IS1272", start=7500, end=9000),
                _make_gene_hit("IS431", start=3000, end=3800),
            ],
            [
                _make_gene_hit("ccrA2", start=10000, end=11350),
                _make_gene_hit("ccrB2", start=11400, end=12750),
                _make_gene_hit("ccrC1", start=25000, end=26000),
            ],
        )
        assert result["SCCmec_Type"] == "IV"
        assert result["SCCmec_Type_secondary"] == "novel_combination"


# ---------------------------------------------------------------------------
# TestCollectInputFiles
# ---------------------------------------------------------------------------

class TestCollectInputFiles:
    """Tests for collect_input_files."""

    def test_single_file(self, tmp_path):
        """Single FASTA file is collected."""
        f = tmp_path / "test.fasta"
        f.write_text(">seq\nATCG\n")

        files = collect_input_files([str(f)])
        assert len(files) == 1
        assert files[0] == str(f)

    def test_directory(self, tmp_path):
        """Directory is scanned for FASTA files."""
        for name in ["a.fasta", "b.fna", "c.fa", "d.txt"]:
            (tmp_path / name).write_text(">seq\nATCG\n")

        files = collect_input_files([str(tmp_path)])
        assert len(files) == 3  # .fasta, .fna, .fa but not .txt

    def test_nonexistent_path(self, tmp_path):
        """Non-existent path is skipped with warning."""
        files = collect_input_files([str(tmp_path / "nonexistent.fasta")])
        assert files == []


# ---------------------------------------------------------------------------
# Integration tests — require BLAST+
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not HAS_BLAST, reason="BLAST+ not installed")
class TestSCCmecTyperIntegration:
    """Integration tests requiring BLAST+."""

    def test_type_test_genome(self, test_genome, temp_output_dir):
        """Type the test genome (which may contain SCCmec content)."""
        typer = SCCmecTyper()
        result = typer.type_file(str(test_genome))

        assert "Input_File" in result
        assert result["Input_File"] == "test_genome"
        for key in TYPING_HEADER:
            assert key in result

    def test_cli_help(self):
        """sccmec-type --help runs without error."""
        result = subprocess.run(
            [sys.executable, "-m", "sccmecextractor.sccmec_type_classification", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "mec-ref" in result.stdout
        assert "ccr-ref" in result.stdout
        assert "mec complex class" in result.stdout.lower()

    def test_cli_typing(self, test_genome, temp_output_dir):
        """CLI produces valid TSV output with all expected columns."""
        output_file = temp_output_dir / "typing.tsv"

        result = subprocess.run(
            [
                sys.executable, "-m", "sccmecextractor.sccmec_type_classification",
                "-f", str(test_genome),
                "-o", str(output_file),
            ],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0, f"Script failed:\n{result.stderr}"
        assert output_file.exists()

        with open(output_file) as f:
            header = f.readline().strip().split("\t")
        for col in TYPING_HEADER:
            assert col in header

    def test_cli_mode_message_default(self, test_genome, temp_output_dir):
        """Default mode prints 'full SCCmec typing' message."""
        output_file = temp_output_dir / "typing.tsv"
        result = subprocess.run(
            [
                sys.executable, "-m", "sccmecextractor.sccmec_type_classification",
                "-f", str(test_genome),
                "-o", str(output_file),
            ],
            capture_output=True,
            text=True,
        )
        assert "full SCCmec typing" in result.stdout


# ---------------------------------------------------------------------------
# Pipeline CLI tests — require BLAST+
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not HAS_BLAST, reason="BLAST+ not installed")
class TestPipelineCustomRefs:
    """Test that pipeline accepts and passes through --mec-ref/--ccr-ref."""

    def test_pipeline_help_shows_ref_args(self):
        """Pipeline --help mentions --mec-ref and --ccr-ref."""
        result = subprocess.run(
            [sys.executable, "-m", "sccmecextractor.pipeline", "--help"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert "--mec-ref" in result.stdout
        assert "--ccr-ref" in result.stdout
