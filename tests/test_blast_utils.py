#!/usr/bin/env python

"""Tests for blast_utils.py

Unit tests using pre-written outfmt 6 fixtures (no BLAST+ installation needed
for most tests).
"""

import shutil
import tempfile
import textwrap
from pathlib import Path
from unittest.mock import patch

import pytest

from sccmecextractor.blast_utils import (
    BlastNotFoundError,
    BlastResult,
    BlastRunner,
    filter_hits,
    get_best_non_overlapping_hits,
    get_default_ref,
    parse_blast_output,
)


@pytest.fixture
def blast6_empty(tmp_path):
    """Create an empty BLAST outfmt 6 file."""
    f = tmp_path / "empty.blast6"
    f.write_text("")
    return str(f)


@pytest.fixture
def blast6_single(tmp_path):
    """Create a BLAST outfmt 6 file with a single hit."""
    content = "ccrA1\tcontig_1\t95.5\t1350\t60\t1\t1\t1350\t5000\t6349\t0.0\t2200\n"
    f = tmp_path / "single.blast6"
    f.write_text(content)
    return str(f)


@pytest.fixture
def blast6_multiple(tmp_path):
    """Create a BLAST outfmt 6 file with multiple hits."""
    content = textwrap.dedent(
        """\
        ccrA1\tcontig_1\t95.5\t1350\t60\t1\t1\t1350\t5000\t6349\t0.0\t2200
        ccrB2\tcontig_1\t92.0\t1629\t130\t1\t1\t1629\t7000\t8628\t0.0\t2500
        ccrC1\tcontig_2\t88.5\t1677\t193\t0\t1\t1677\t100\t1776\t1e-50\t1800
        rlmH\tcontig_3\t99.0\t480\t5\t0\t1\t480\t10000\t10479\t1e-100\t900
    """
    )
    f = tmp_path / "multiple.blast6"
    f.write_text(content)
    return str(f)


@pytest.fixture
def blast6_overlapping(tmp_path):
    """Create a BLAST outfmt 6 file with overlapping hits on the same contig."""
    content = textwrap.dedent(
        """\
        ccrA1\tcontig_1\t95.5\t1350\t60\t1\t1\t1350\t5000\t6349\t0.0\t2200
        ccrA2\tcontig_1\t90.0\t1350\t135\t0\t1\t1350\t5100\t6449\t0.0\t2000
        ccrB1\tcontig_1\t93.0\t1629\t114\t1\t1\t1629\t20000\t21628\t0.0\t2400
        ccrC1\tcontig_2\t88.0\t1677\t201\t0\t1\t1677\t100\t1776\t1e-50\t1800
    """
    )
    f = tmp_path / "overlapping.blast6"
    f.write_text(content)
    return str(f)


class TestBlastResultParsing:
    """Tests for parse_blast_output."""

    def test_empty_file(self, blast6_empty):
        """Empty file returns empty list."""
        hits = parse_blast_output(blast6_empty)
        assert hits == []

    def test_nonexistent_file(self, tmp_path):
        """Non-existent file returns empty list."""
        hits = parse_blast_output(str(tmp_path / "nonexistent.blast6"))
        assert hits == []

    def test_single_hit(self, blast6_single):
        """Single hit is parsed correctly."""
        hits = parse_blast_output(blast6_single)
        assert len(hits) == 1

        hit = hits[0]
        assert hit.qseqid == "ccrA1"
        assert hit.sseqid == "contig_1"
        assert hit.pident == 95.5
        assert hit.length == 1350
        assert hit.mismatch == 60
        assert hit.gapopen == 1
        assert hit.qstart == 1
        assert hit.qend == 1350
        assert hit.sstart == 5000
        assert hit.send == 6349
        assert hit.evalue == 0.0
        assert hit.bitscore == 2200

    def test_multiple_hits(self, blast6_multiple):
        """Multiple hits are all parsed."""
        hits = parse_blast_output(blast6_multiple)
        assert len(hits) == 4

        assert hits[0].qseqid == "ccrA1"
        assert hits[1].qseqid == "ccrB2"
        assert hits[2].qseqid == "ccrC1"
        assert hits[3].qseqid == "rlmH"

    def test_malformed_line_skipped(self, tmp_path):
        """Lines with fewer than 12 fields are skipped."""
        content = "too\tfew\tfields\n"
        f = tmp_path / "malformed.blast6"
        f.write_text(content)

        hits = parse_blast_output(str(f))
        assert hits == []


class TestHitFiltering:
    """Tests for filter_hits."""

    def test_filter_by_identity(self, blast6_multiple):
        """Hits below identity threshold are removed."""
        hits = parse_blast_output(blast6_multiple)
        ref_lengths = {"ccrA1": 1350, "ccrB2": 1629, "ccrC1": 1677, "rlmH": 480}

        filtered = filter_hits(hits, min_pident=90.0, min_coverage=0.5, ref_lengths=ref_lengths)
        gene_names = [h.qseqid for h in filtered]

        assert "ccrA1" in gene_names
        assert "ccrB2" in gene_names
        assert "rlmH" in gene_names
        assert "ccrC1" not in gene_names  # 88.5 < 90.0

    def test_filter_by_coverage(self, blast6_multiple):
        """Hits below coverage threshold are removed."""
        hits = parse_blast_output(blast6_multiple)
        ref_lengths = {"ccrA1": 1350, "ccrB2": 1629, "ccrC1": 1677, "rlmH": 480}

        filtered = filter_hits(hits, min_pident=80.0, min_coverage=0.99, ref_lengths=ref_lengths)
        gene_names = [h.qseqid for h in filtered]

        assert "ccrA1" in gene_names  # 1350/1350 = 1.0
        assert "ccrB2" in gene_names  # 1629/1629 = 1.0
        assert "ccrC1" in gene_names  # 1677/1677 = 1.0
        assert "rlmH" in gene_names  # 480/480 = 1.0

    def test_unknown_ref_excluded(self, blast6_single):
        """Hits with unknown reference length are excluded."""
        hits = parse_blast_output(blast6_single)
        ref_lengths = {}  # No reference lengths

        filtered = filter_hits(hits, min_pident=80.0, min_coverage=0.5, ref_lengths=ref_lengths)
        assert filtered == []

    def test_empty_input(self):
        """Empty hit list returns empty."""
        filtered = filter_hits([], min_pident=80.0, min_coverage=0.5, ref_lengths={})
        assert filtered == []


class TestBestNonOverlappingHits:
    """Tests for get_best_non_overlapping_hits."""

    def test_no_overlap(self, blast6_multiple):
        """Non-overlapping hits are all kept."""
        hits = parse_blast_output(blast6_multiple)
        selected = get_best_non_overlapping_hits(hits, overlap_threshold=50)

        assert len(selected) == 4

    def test_overlapping_resolved(self, blast6_overlapping):
        """Overlapping hits on the same contig are resolved by bitscore."""
        hits = parse_blast_output(blast6_overlapping)
        selected = get_best_non_overlapping_hits(hits, overlap_threshold=50)

        # ccrA1 (bitscore 2200) should beat ccrA2 (bitscore 2000) on contig_1
        # ccrB1 is on contig_1 but at position 20000 (no overlap)
        # ccrC1 is on contig_2
        contig1_genes = [h.qseqid for h in selected if h.sseqid == "contig_1"]
        assert "ccrA1" in contig1_genes
        assert "ccrA2" not in contig1_genes
        assert "ccrB1" in contig1_genes

        contig2_genes = [h.qseqid for h in selected if h.sseqid == "contig_2"]
        assert "ccrC1" in contig2_genes

    def test_empty_input(self):
        """Empty hit list returns empty."""
        selected = get_best_non_overlapping_hits([])
        assert selected == []

    def test_single_hit(self, blast6_single):
        """Single hit is always kept."""
        hits = parse_blast_output(blast6_single)
        selected = get_best_non_overlapping_hits(hits)
        assert len(selected) == 1

    def test_reverse_strand_overlap(self, tmp_path):
        """Overlapping hits on reverse strand are detected."""
        content = textwrap.dedent(
            """\
            ccrA1\tcontig_1\t95.0\t1350\t67\t1\t1\t1350\t6349\t5000\t0.0\t2200
            ccrA2\tcontig_1\t90.0\t1350\t135\t0\t1\t1350\t6449\t5100\t0.0\t2000
        """
        )
        f = tmp_path / "reverse.blast6"
        f.write_text(content)

        hits = parse_blast_output(str(f))
        selected = get_best_non_overlapping_hits(hits, overlap_threshold=50)

        assert len(selected) == 1
        assert selected[0].qseqid == "ccrA1"


class TestBlastInstallation:
    """Tests for BLAST+ installation check."""

    def test_blast_not_found_raises(self):
        """BlastNotFoundError raised when BLAST+ not on PATH."""
        with patch("shutil.which", return_value=None):
            with pytest.raises(BlastNotFoundError, match="blastn"):
                BlastRunner()

    def test_makeblastdb_not_found_raises(self):
        """BlastNotFoundError raised when makeblastdb not on PATH."""

        def selective_which(name):
            if name == "blastn":
                return "/usr/bin/blastn"
            return None

        with patch("sccmecextractor.blast_utils.shutil.which", side_effect=selective_which):
            with pytest.raises(BlastNotFoundError, match="makeblastdb"):
                BlastRunner()


class TestGetDefaultRef:
    """Tests for get_default_ref context manager."""

    def test_rlmh_ref_exists(self):
        """Bundled rlmH multi-species reference is accessible."""
        with get_default_ref("rlmH.fasta") as path:
            assert path.exists()
            content = path.read_text()
            assert "_rlmH" in content
            assert ">Staphylococcus_epidermidis_rlmH" in content
            assert ">Staphylococcus_aureus_rlmH" in content

    def test_ccr_ref_exists(self):
        """Bundled ccr reference is accessible."""
        with get_default_ref("ccr_genes.fasta") as path:
            assert path.exists()
            content = path.read_text()
            assert ">ccrA1" in content
            assert ">ccrC1" in content

    def test_mec_ref_exists(self):
        """Bundled mec reference is accessible."""
        with get_default_ref("mec_genes_allotypes.fasta") as path:
            assert path.exists()
            content = path.read_text()
            assert ">mecA" in content
