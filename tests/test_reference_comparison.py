#!/usr/bin/env python

"""Tests for reference_comparison.py (hybrid detection module).

Unit tests for interval merging, novel coverage calculation, classification
logic, accession extraction, and typing report parsing.  Integration tests
for type_file and type_batch require BLAST+.
"""

import shutil
import tempfile
import pytest
from pathlib import Path

from sccmecextractor.reference_comparison import (
    HybridTyper,
    HYBRID_SUMMARY_HEADER,
    HYBRID_DETAIL_HEADER,
    _extract_accession,
    classify_with_typing,
    read_typing_report,
)

HAS_BLAST = shutil.which("blastn") is not None


# -----------------------------------------------------------------------
# _extract_accession
# -----------------------------------------------------------------------

class TestExtractAccession:
    def test_sccmec_suffix(self):
        assert _extract_accession("GCF_000009585_SCCmec") == "GCF_000009585"

    def test_scc_suffix(self):
        assert _extract_accession("GCF_000009585_SCC") == "GCF_000009585"

    def test_no_suffix(self):
        assert _extract_accession("GCF_000009585") == "GCF_000009585"

    def test_version_suffix(self):
        assert _extract_accession("GCF_000009585.1_ASM_genomic_SCCmec") == "GCF_000009585.1_ASM_genomic"

    def test_empty_string(self):
        assert _extract_accession("") == ""


# -----------------------------------------------------------------------
# classify_with_typing
# -----------------------------------------------------------------------

class TestClassifyWithTyping:
    """Test the 2x2 matrix of single/multi type x single/multi ccr."""

    def _typing_info(self, n_ccr):
        return {
            "ccr_complex_type": "2" if n_ccr == 1 else "2;5",
            "n_ccr_complexes": n_ccr,
            "Is_Composite": "False",
            "SCCmec_Type": "-",
            "mec_genes": "-",
        }

    def test_canonical(self):
        result = classify_with_typing("canonical", 1, self._typing_info(1))
        assert result == "canonical"

    def test_hybrid(self):
        result = classify_with_typing("multi_type", 2, self._typing_info(1))
        assert result == "hybrid"

    def test_multi_ccr(self):
        result = classify_with_typing("multi_type", 2, self._typing_info(2))
        assert result == "multi_ccr"

    def test_multi_ccr_canonical(self):
        result = classify_with_typing("canonical", 1, self._typing_info(2))
        assert result == "multi_ccr_canonical"

    def test_no_match(self):
        result = classify_with_typing("no_match", 0, self._typing_info(1))
        assert result == "no_match"

    def test_no_typing_info_single(self):
        result = classify_with_typing("canonical", 1, None)
        assert result == "canonical"

    def test_no_typing_info_multi(self):
        result = classify_with_typing("multi_type", 3, None)
        assert result == "multi_type"


# -----------------------------------------------------------------------
# _merge_intervals (static method)
# -----------------------------------------------------------------------

class TestMergeIntervals:
    def setup_method(self):
        self.merge = HybridTyper._merge_intervals

    def test_empty(self):
        assert self.merge([]) == []

    def test_single_interval(self):
        assert self.merge([(100, 500)]) == [[100, 500]]

    def test_non_overlapping(self):
        result = self.merge([(100, 500), (700, 1000)])
        assert result == [[100, 500], [700, 1000]]

    def test_overlapping(self):
        result = self.merge([(100, 500), (400, 800)])
        assert result == [[100, 800]]

    def test_contained(self):
        result = self.merge([(100, 1000), (300, 500)])
        assert result == [[100, 1000]]

    def test_adjacent(self):
        result = self.merge([(100, 500), (500, 800)])
        assert result == [[100, 800]]

    def test_unsorted_input(self):
        result = self.merge([(700, 1000), (100, 500)])
        assert result == [[100, 500], [700, 1000]]

    def test_multiple_overlapping(self):
        result = self.merge([(100, 500), (400, 800), (750, 1200)])
        assert result == [[100, 1200]]

    def test_three_separate(self):
        result = self.merge([(100, 200), (500, 600), (900, 1000)])
        assert result == [[100, 200], [500, 600], [900, 1000]]


# -----------------------------------------------------------------------
# _calculate_novel_coverage (static method)
# -----------------------------------------------------------------------

class TestCalculateNovelCoverage:
    def setup_method(self):
        self.calc = HybridTyper._calculate_novel_coverage

    def test_no_overlap(self):
        """Secondary entirely outside primary."""
        novel, novel_bp = self.calc(
            [[100, 500]], [[700, 1000]]
        )
        assert novel == [[700, 1000]]
        assert novel_bp == 300

    def test_complete_overlap(self):
        """Secondary entirely within primary."""
        novel, novel_bp = self.calc(
            [[100, 1000]], [[300, 500]]
        )
        assert novel == []
        assert novel_bp == 0

    def test_partial_overlap_right(self):
        """Secondary extends beyond primary on the right."""
        novel, novel_bp = self.calc(
            [[100, 500]], [[400, 800]]
        )
        assert novel == [[500, 800]]
        assert novel_bp == 300

    def test_partial_overlap_left(self):
        """Secondary extends beyond primary on the left."""
        novel, novel_bp = self.calc(
            [[500, 1000]], [[200, 600]]
        )
        assert novel == [[200, 500]]
        assert novel_bp == 300

    def test_secondary_spans_primary(self):
        """Secondary extends beyond primary on both sides."""
        novel, novel_bp = self.calc(
            [[400, 600]], [[200, 800]]
        )
        assert novel == [[200, 400], [600, 800]]
        assert novel_bp == 400

    def test_empty_primary(self):
        """No primary — all secondary is novel."""
        novel, novel_bp = self.calc(
            [], [[100, 500]]
        )
        assert novel == [[100, 500]]
        assert novel_bp == 400

    def test_empty_secondary(self):
        """No secondary — nothing is novel."""
        novel, novel_bp = self.calc(
            [[100, 500]], []
        )
        assert novel == []
        assert novel_bp == 0

    def test_multiple_primary_gaps(self):
        """Secondary fills a gap between two primary intervals."""
        novel, novel_bp = self.calc(
            [[100, 400], [600, 900]], [[350, 650]]
        )
        assert novel == [[400, 600]]
        assert novel_bp == 200

    def test_multiple_secondary(self):
        """Multiple secondary intervals, partial overlap."""
        novel, novel_bp = self.calc(
            [[300, 700]], [[100, 400], [600, 900]]
        )
        assert novel == [[100, 300], [700, 900]]
        assert novel_bp == 400


# -----------------------------------------------------------------------
# read_typing_report
# -----------------------------------------------------------------------

class TestReadTypingReport:
    def test_basic_parsing(self, tmp_path):
        tsv = tmp_path / "summary.tsv"
        tsv.write_text(
            "Input_File\tStatus\tccr_complex_type\tIs_Composite\tSCCmec_Type\tmec_genes\n"
            "GCF_001\textracted\t2\tFalse\tII\tmecA(full)\n"
            "GCF_002\textracted\t2;5\tTrue\tV\tmecA(full)\n"
            "GCF_003\tfailed\t3\tFalse\tIII\tmecA(full)\n"
        )
        result = read_typing_report(str(tsv))
        assert "GCF_001" in result
        assert "GCF_002" in result
        assert "GCF_003" not in result  # failed, excluded
        assert result["GCF_001"]["n_ccr_complexes"] == 1
        assert result["GCF_002"]["n_ccr_complexes"] == 2

    def test_composite_extracted_included(self, tmp_path):
        tsv = tmp_path / "summary.tsv"
        tsv.write_text(
            "Input_File\tStatus\tccr_complex_type\tIs_Composite\tSCCmec_Type\tmec_genes\n"
            "GCF_001\tcomposite_extracted\t2;4\tTrue\tIV\tmecA(full)\n"
        )
        result = read_typing_report(str(tsv))
        assert "GCF_001" in result
        assert result["GCF_001"]["n_ccr_complexes"] == 2

    def test_no_ccr(self, tmp_path):
        tsv = tmp_path / "summary.tsv"
        tsv.write_text(
            "Input_File\tStatus\tccr_complex_type\tIs_Composite\tSCCmec_Type\tmec_genes\n"
            "GCF_001\textracted\t-\tFalse\t-\t-\n"
        )
        result = read_typing_report(str(tsv))
        assert result["GCF_001"]["n_ccr_complexes"] == 0


# -----------------------------------------------------------------------
# _classify_hybrid (via _merge_intervals and _calculate_novel_coverage)
# -----------------------------------------------------------------------

class TestClassifyHybrid:
    """Test the full classification pipeline with synthetic type_profiles."""

    def setup_method(self):
        self.typer = HybridTyper.__new__(HybridTyper)

    def _make_profile(self, footprint, covered_bp, best_subtype="Ia",
                      coverage=100.0, ref_length=30000, pident=98.0, hsps=1):
        return {
            "element_footprint": footprint,
            "element_covered_bp": covered_bp,
            "best_subtype": best_subtype,
            "best_coverage": coverage,
            "best_ref_length": ref_length,
            "weighted_pident": pident,
            "num_hsps": hsps,
        }

    def test_single_type_canonical(self):
        profiles = {
            "I": self._make_profile([[100, 25000]], 24900, "Ia"),
        }
        summary, detail = self.typer._classify_hybrid(
            "test_element", 30000, profiles
        )
        assert summary["hybrid_call"] == "canonical"
        assert summary["hybrid_best_match"] == "I"
        assert summary["hybrid_components"] == "I"
        assert len(detail) == 1

    def test_two_types_hybrid(self):
        """Two types covering distinct regions → multi_type (no typing info)."""
        profiles = {
            "IV": self._make_profile([[100, 15000]], 14900, "IVa"),
            "II": self._make_profile([[16000, 28000]], 12000, "IIa"),
        }
        summary, detail = self.typer._classify_hybrid(
            "test_element", 30000, profiles
        )
        assert summary["hybrid_call"] == "multi_type"
        assert "IV" in summary["hybrid_components"]
        assert "II" in summary["hybrid_components"]
        assert len(detail) == 2

    def test_two_types_overlapping_canonical(self):
        """Two types covering the same region → canonical (secondary has no novel territory)."""
        profiles = {
            "IV": self._make_profile([[100, 25000]], 24900, "IVa", pident=98.0),
            "II": self._make_profile([[200, 24000]], 23800, "IIa", pident=95.0),
        }
        summary, detail = self.typer._classify_hybrid(
            "test_element", 30000, profiles
        )
        assert summary["hybrid_call"] == "canonical"
        assert summary["hybrid_best_match"] == "IV"

    def test_no_profiles_no_match(self):
        summary, detail = self.typer._classify_hybrid(
            "test_element", 30000, {}
        )
        assert summary["hybrid_call"] == "no_match"
        assert detail == []

    def test_with_typing_info_hybrid(self):
        """Multi-type + single ccr → hybrid."""
        profiles = {
            "IV": self._make_profile([[100, 15000]], 14900, "IVa"),
            "II": self._make_profile([[16000, 28000]], 12000, "IIa"),
        }
        typing_info = {"n_ccr_complexes": 1, "ccr_complex_type": "2"}
        summary, detail = self.typer._classify_hybrid(
            "test_element", 30000, profiles, typing_info
        )
        assert summary["hybrid_call"] == "hybrid"

    def test_with_typing_info_multi_ccr(self):
        """Multi-type + multiple ccr → multi_ccr."""
        profiles = {
            "IV": self._make_profile([[100, 15000]], 14900, "IVa"),
            "II": self._make_profile([[16000, 28000]], 12000, "IIa"),
        }
        typing_info = {"n_ccr_complexes": 2, "ccr_complex_type": "2;5"}
        summary, detail = self.typer._classify_hybrid(
            "test_element", 30000, profiles, typing_info
        )
        assert summary["hybrid_call"] == "multi_ccr"

    def test_three_components(self):
        """Three types covering distinct regions."""
        profiles = {
            "IV": self._make_profile([[100, 8000]], 7900, "IVa"),
            "II": self._make_profile([[10000, 18000]], 8000, "IIa"),
            "V": self._make_profile([[20000, 28000]], 8000, "Va"),
        }
        summary, detail = self.typer._classify_hybrid(
            "test_element", 30000, profiles
        )
        assert summary["hybrid_call"] == "multi_type"
        assert len(summary["hybrid_components"].split(";")) == 3

    def test_novel_below_threshold_stays_canonical(self):
        """Secondary type has novel territory but below thresholds."""
        profiles = {
            "IV": self._make_profile([[100, 25000]], 24900, "IVa"),
            "II": self._make_profile([[25500, 26000]], 500, "IIa", coverage=10.0, pident=90.0),
        }
        summary, detail = self.typer._classify_hybrid(
            "test_element", 30000, profiles
        )
        # 500 bp novel < 1000 bp threshold
        assert summary["hybrid_call"] == "canonical"

    def test_summary_has_all_headers(self):
        profiles = {
            "IV": self._make_profile([[100, 25000]], 24900, "IVa"),
        }
        summary, detail = self.typer._classify_hybrid(
            "test_element", 30000, profiles
        )
        for col in HYBRID_SUMMARY_HEADER:
            assert col in summary

    def test_detail_has_all_headers(self):
        profiles = {
            "IV": self._make_profile([[100, 25000]], 24900, "IVa"),
        }
        summary, detail = self.typer._classify_hybrid(
            "test_element", 30000, profiles
        )
        assert len(detail) > 0
        for col in HYBRID_DETAIL_HEADER:
            assert col in detail[0]


# -----------------------------------------------------------------------
# Integration tests (require BLAST+)
# -----------------------------------------------------------------------

@pytest.mark.skipif(not HAS_BLAST, reason="BLAST+ not installed")
class TestHybridTyperIntegration:
    """Integration tests using real BLAST against bundled references."""

    @pytest.fixture
    def typer(self):
        return HybridTyper()

    @pytest.fixture
    def simple_element(self, tmp_path):
        """Create a minimal FASTA with a short sequence for testing."""
        fasta = tmp_path / "test_element.fasta"
        # 1000 bp of random-ish sequence — won't match any reference well
        seq = "ATCGATCG" * 125
        fasta.write_text(f">test_element\n{seq}\n")
        return str(fasta)

    def test_type_file_returns_tuple(self, typer, simple_element):
        summary, detail = typer.type_file(simple_element)
        assert isinstance(summary, dict)
        assert isinstance(detail, list)
        for col in HYBRID_SUMMARY_HEADER:
            assert col in summary

    def test_type_file_no_match(self, typer, simple_element):
        """A random sequence should produce no_match."""
        summary, detail = typer.type_file(simple_element)
        assert summary["hybrid_call"] == "no_match"

    def test_type_batch_empty(self, typer):
        result = typer.type_batch([])
        assert result == []

    def test_type_batch_single(self, typer, simple_element):
        results = typer.type_batch([simple_element])
        assert len(results) == 1
        summary, detail = results[0]
        assert summary["hybrid_call"] == "no_match"

    def test_type_batch_with_typing_context(self, typer, simple_element):
        typing_context = {
            "test_element": {
                "n_ccr_complexes": 1,
                "ccr_complex_type": "2",
            }
        }
        results = typer.type_batch(
            [simple_element], typing_context=typing_context
        )
        assert len(results) == 1
        summary, detail = results[0]
        # Still no_match since the sequence doesn't match any reference
        assert summary["hybrid_call"] == "no_match"

    def test_type_file_with_typing_info(self, typer, simple_element):
        typing_info = {
            "n_ccr_complexes": 1,
            "ccr_complex_type": "2",
        }
        summary, detail = typer.type_file(
            simple_element, typing_info=typing_info
        )
        assert summary["hybrid_call"] == "no_match"

    def test_element_id_populated(self, typer, simple_element):
        summary, detail = typer.type_file(simple_element)
        assert summary["element_id"] == "test_element"
