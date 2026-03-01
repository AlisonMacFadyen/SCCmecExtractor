#!/usr/bin/env python

"""Tests for type_sccmec.py

Unit tests for classifiers using mock BlastResult objects, plus integration
tests that require BLAST+.
"""

import shutil
import subprocess

import pytest
from pathlib import Path

from sccmecextractor.blast_utils import BlastResult
from sccmecextractor.type_sccmec import (
    CcrClassifier,
    CcrComplexLookup,
    GeneHit,
    MecClassifier,
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


class TestMecClassifier:
    """Tests for MecClassifier."""

    @pytest.fixture
    def mec_classifier(self):
        """Create a MecClassifier with bundled reference."""
        from sccmecextractor.blast_utils import get_default_ref

        with get_default_ref("mec_genes_allotypes.fasta") as ref:
            return MecClassifier(str(ref))

    def test_full_hit(self, mec_classifier):
        """High identity + high coverage = full."""
        # mecA is ~2007bp in the bundled reference
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
        # ccrA2 should be excluded due to overlap
        assert "ccrA2" not in gene_names


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


class TestOutputFormat:
    """Tests for SCCmecTyper output formatting."""

    def test_format_with_results(self):
        """Formatting with mec and ccr results."""
        result = SCCmecTyper._format_result(
            "test_genome",
            [GeneHit("mecA", 98.5, 99.0, "full")],
            [
                GeneHit("ccrA2", 95.0, 98.0, "full"),
                GeneHit("ccrB2", 93.0, 97.0, "full"),
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
        result = SCCmecTyper._format_result("test_genome", [], [])

        assert result["mec_genes"] == "-"
        assert result["ccr_genes"] == "-"
        assert result["ccr_allotypes"] == "-"
        assert result["ccr_complex_type"] == "-"


@pytest.mark.skipif(not HAS_BLAST, reason="BLAST+ not installed")
class TestSCCmecTyperIntegration:
    """Integration tests requiring BLAST+."""

    def test_type_test_genome(self, test_genome, temp_output_dir):
        """Type the test genome (which may contain SCCmec content)."""
        typer = SCCmecTyper()
        result = typer.type_file(str(test_genome))

        assert "Input_File" in result
        assert result["Input_File"] == "test_genome"
        # Result should have all expected keys
        for key in [
            "mec_genes", "mec_identity", "mec_coverage",
            "ccr_genes", "ccr_allotypes", "ccr_identity",
            "ccr_complex_type",
        ]:
            assert key in result

    def test_cli_help(self):
        """sccmec-type --help runs without error."""
        result = subprocess.run(
            ["python", "-m", "sccmecextractor.type_sccmec", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "mec-ref" in result.stdout
        assert "ccr-ref" in result.stdout

    def test_cli_typing(self, test_genome, temp_output_dir):
        """CLI produces valid TSV output."""
        output_file = temp_output_dir / "typing.tsv"

        result = subprocess.run(
            [
                "python", "-m", "sccmecextractor.type_sccmec",
                "-f", str(test_genome),
                "-o", str(output_file),
            ],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0, f"Script failed:\n{result.stderr}"
        assert output_file.exists()

        # Check TSV format
        with open(output_file) as f:
            header = f.readline().strip().split("\t")
        expected = [
            "Input_File", "mec_genes", "mec_identity", "mec_coverage",
            "ccr_genes", "ccr_allotypes", "ccr_identity", "ccr_complex_type",
        ]
        assert header == expected


class TestCcrComplexLookup:
    """Tests for CcrComplexLookup."""

    def test_simple_pair_type2(self):
        """ccrA2 + ccrB2 = complex type 2."""
        hits = [
            GeneHit("ccrA2", 95.0, 98.0, "full"),
            GeneHit("ccrB2", 93.0, 97.0, "full"),
        ]
        assert CcrComplexLookup.lookup(hits) == "2"

    def test_single_ccrC1_type5(self):
        """ccrC1 alone = complex type 5."""
        hits = [GeneHit("ccrC1", 90.0, 95.0, "full")]
        assert CcrComplexLookup.lookup(hits) == "5"

    def test_mixed_pair_type8(self):
        """ccrA1 + ccrB3 = complex type 8."""
        hits = [
            GeneHit("ccrA1", 92.0, 96.0, "full"),
            GeneHit("ccrB3", 90.0, 94.0, "full"),
        ]
        assert CcrComplexLookup.lookup(hits) == "8"

    def test_composite_types_2_and_5(self):
        """ccrA2 + ccrB2 + ccrC1 = composite types 2;5."""
        hits = [
            GeneHit("ccrA2", 95.0, 98.0, "full"),
            GeneHit("ccrB2", 93.0, 97.0, "full"),
            GeneHit("ccrC1", 90.0, 95.0, "full"),
        ]
        assert CcrComplexLookup.lookup(hits) == "2;5"

    def test_triple_composite_types_1_5_9(self):
        """ccrA1 + ccrB1 + ccrC1 + ccrC2 = composite types 1;5;9."""
        hits = [
            GeneHit("ccrA1", 95.0, 98.0, "full"),
            GeneHit("ccrB1", 93.0, 97.0, "full"),
            GeneHit("ccrC1", 90.0, 95.0, "full"),
            GeneHit("ccrC2", 88.0, 92.0, "full"),
        ]
        assert CcrComplexLookup.lookup(hits) == "1;5;9"

    def test_novel_genes_flagged(self):
        """Novel classification appends '(novel_genes)' to result."""
        hits = [
            GeneHit("ccrA2", 75.0, 92.0, "novel_full"),
            GeneHit("ccrB2", 93.0, 97.0, "full"),
        ]
        assert CcrComplexLookup.lookup(hits) == "2 (novel_genes)"

    def test_unknown_combination(self):
        """Unrecognised allotype combination returns 'novel_combination'."""
        hits = [GeneHit("ccrB9", 90.0, 95.0, "full")]
        assert CcrComplexLookup.lookup(hits) == "novel_combination"

    def test_no_hits(self):
        """No ccr hits returns '-'."""
        assert CcrComplexLookup.lookup([]) == "-"

    # --- New complex types (10-22) ---

    def test_type10_ccrA8B9(self):
        """ccrA8 + ccrB9 = complex type 10 (Xiao et al. 2023)."""
        hits = [
            GeneHit("ccrA8", 92.0, 96.0, "full"),
            GeneHit("ccrB9", 90.0, 94.0, "full"),
        ]
        assert CcrComplexLookup.lookup(hits) == "10"

    def test_type11_ccrA9B3(self):
        """ccrA9 + ccrB3 = complex type 11 (Huang et al. 2024)."""
        hits = [
            GeneHit("ccrA9", 92.0, 96.0, "full"),
            GeneHit("ccrB3", 90.0, 94.0, "full"),
        ]
        assert CcrComplexLookup.lookup(hits) == "11"

    def test_type12_ccrA10B1(self):
        """ccrA10 + ccrB1 = complex type 12 (Huang et al. 2024)."""
        hits = [
            GeneHit("ccrA10", 92.0, 96.0, "full"),
            GeneHit("ccrB1", 90.0, 94.0, "full"),
        ]
        assert CcrComplexLookup.lookup(hits) == "12"

    def test_type13_ccrA10B10(self):
        """ccrA10 + ccrB10 = complex type 13 (Huang et al. 2024)."""
        hits = [
            GeneHit("ccrA10", 92.0, 96.0, "full"),
            GeneHit("ccrB10", 90.0, 94.0, "full"),
        ]
        assert CcrComplexLookup.lookup(hits) == "13"

    def test_type20_ccrC3(self):
        """ccrC3 alone = complex type 20 (Huang et al. 2024)."""
        hits = [GeneHit("ccrC3", 90.0, 95.0, "full")]
        assert CcrComplexLookup.lookup(hits) == "20"

    def test_type21_ccrC4(self):
        """ccrC4 alone = complex type 21 (Huang et al. 2024)."""
        hits = [GeneHit("ccrC4", 90.0, 95.0, "full")]
        assert CcrComplexLookup.lookup(hits) == "21"

    def test_type22_ccrC5(self):
        """ccrC5 alone = complex type 22 (Huang et al. 2024)."""
        hits = [GeneHit("ccrC5", 90.0, 95.0, "full")]
        assert CcrComplexLookup.lookup(hits) == "22"

    def test_type13_preferred_over_type12(self):
        """When ccrA10 + ccrB10 + ccrB1 present, matched pair (13) takes priority."""
        hits = [
            GeneHit("ccrA10", 92.0, 96.0, "full"),
            GeneHit("ccrB10", 90.0, 94.0, "full"),
            GeneHit("ccrB1", 88.0, 92.0, "full"),
        ]
        result = CcrComplexLookup.lookup(hits)
        assert result == "13"

    def test_composite_type11_and_type5(self):
        """ccrA9 + ccrB3 + ccrC1 = composite types 11;5."""
        hits = [
            GeneHit("ccrA9", 92.0, 96.0, "full"),
            GeneHit("ccrB3", 90.0, 94.0, "full"),
            GeneHit("ccrC1", 90.0, 95.0, "full"),
        ]
        assert CcrComplexLookup.lookup(hits) == "11;5"
