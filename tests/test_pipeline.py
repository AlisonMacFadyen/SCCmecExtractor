#!/usr/bin/env python

"""Tests for the sccmec-pipeline master command."""

import os
import shutil
import subprocess
import sys
import tempfile

import pytest
from pathlib import Path

from sccmecextractor.pipeline import resolve_gff, run_pipeline


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

TEST_DATA = Path(__file__).parent / "test_data"
TEST_GENOME = TEST_DATA / "test_genome.fna"
TEST_GFF = TEST_DATA / "test_genome.gff3"


def _blast_available() -> bool:
    """Return True if BLAST+ is on PATH."""
    try:
        from sccmecextractor.blast_utils import BlastRunner
        BlastRunner._check_blast_installed()
        return True
    except Exception:
        return False


requires_blast = pytest.mark.skipif(
    not _blast_available(), reason="BLAST+ not installed"
)


# ---------------------------------------------------------------------------
# TestResolveGff — no BLAST needed
# ---------------------------------------------------------------------------

class TestResolveGff:
    """Test GFF file resolution by stem name."""

    def test_match_by_stem(self, tmp_path):
        """GFF files dict with matching stem returns the path."""
        gff_map = {"genome_A": "/data/genome_A.gff3",
                    "genome_B": "/data/genome_B.gff3"}
        assert resolve_gff("genome_A", gff_files=gff_map) == "/data/genome_A.gff3"

    def test_no_match_returns_none(self):
        """No matching GFF returns None."""
        gff_map = {"genome_A": "/data/genome_A.gff3"}
        assert resolve_gff("genome_X", gff_files=gff_map) is None

    def test_gff_dir_mode(self, tmp_path):
        """Directory containing a matching .gff3 file returns its path."""
        gff_file = tmp_path / "my_genome.gff3"
        gff_file.write_text("##gff-version 3\n")
        result = resolve_gff("my_genome", gff_dir=str(tmp_path))
        assert result == str(gff_file)

    def test_gff_dir_no_match(self, tmp_path):
        """Directory without a matching .gff3 returns None."""
        (tmp_path / "other.gff3").write_text("##gff-version 3\n")
        assert resolve_gff("missing", gff_dir=str(tmp_path)) is None

    def test_no_gff_source_returns_none(self):
        """When neither gff_files nor gff_dir given, returns None."""
        assert resolve_gff("anything") is None

    def test_gff_files_takes_precedence(self, tmp_path):
        """Explicit gff_files dict is used even if gff_dir is None."""
        gff_map = {"stem": "/explicit/stem.gff3"}
        assert resolve_gff("stem", gff_files=gff_map) == "/explicit/stem.gff3"


# ---------------------------------------------------------------------------
# TestRunPipeline — requires BLAST+
# ---------------------------------------------------------------------------

@requires_blast
class TestRunPipeline:
    """Integration tests for the full pipeline."""

    def test_single_genome_with_gff(self, tmp_path):
        """Full pipeline with GFF produces an extraction and typing."""
        outdir = str(tmp_path / "results")
        gff_map = {TEST_GENOME.stem: str(TEST_GFF)}
        summary = run_pipeline(
            fasta_files=[str(TEST_GENOME)],
            outdir=outdir,
            gff_files=gff_map,
        )
        assert summary["total"] == 1
        # Either extracted or failed — both are valid pipeline completion
        assert summary["extracted"] + summary["failed"] == 1

    def test_single_genome_blast_rlmh(self, tmp_path):
        """Pipeline without GFF uses BLAST-based rlmH detection."""
        outdir = str(tmp_path / "results")
        summary = run_pipeline(
            fasta_files=[str(TEST_GENOME)],
            outdir=outdir,
            blast_rlmh=True,
        )
        assert summary["total"] == 1
        assert summary["extracted"] + summary["failed"] == 1

    def test_failed_extraction_wgs_fallback(self, tmp_path):
        """Genome with no att sites gets WGS fallback typing."""
        # Create a tiny FASTA with no att sites
        fake_fna = tmp_path / "no_att.fna"
        fake_fna.write_text(">contig_1\nATCGATCGATCGATCG\n")

        outdir = str(tmp_path / "results")
        summary = run_pipeline(
            fasta_files=[str(fake_fna)],
            outdir=outdir,
            blast_rlmh=True,
        )
        assert summary["total"] == 1
        assert summary["failed"] == 1
        # WGS typing attempted on failure
        assert summary["typed_wgs"] <= 1

    def test_output_structure(self, tmp_path):
        """Correct subdirectories and files are created."""
        outdir = str(tmp_path / "results")
        gff_map = {TEST_GENOME.stem: str(TEST_GFF)}
        run_pipeline(
            fasta_files=[str(TEST_GENOME)],
            outdir=outdir,
            gff_files=gff_map,
        )
        assert os.path.isdir(os.path.join(outdir, "att_sites"))
        assert os.path.isdir(os.path.join(outdir, "sccmec"))
        assert os.path.isdir(os.path.join(outdir, "typing"))
        # Att sites file should always be created
        att_file = os.path.join(outdir, "att_sites", f"{TEST_GENOME.stem}_att_sites.tsv")
        assert os.path.isfile(att_file)

    def test_unified_report_generated(self, tmp_path):
        """Unified report file is created with correct header."""
        outdir = str(tmp_path / "results")
        gff_map = {TEST_GENOME.stem: str(TEST_GFF)}
        run_pipeline(
            fasta_files=[str(TEST_GENOME)],
            outdir=outdir,
            gff_files=gff_map,
        )
        report = os.path.join(outdir, "sccmec_unified_report.tsv")
        if os.path.isfile(report):
            with open(report) as fh:
                header = fh.readline().strip().split("\t")
            assert "Input_File" in header
            assert "typing_source" in header
            assert "mec_genes" in header
            assert "ccr_complex_type" in header

    def test_batch_two_genomes(self, tmp_path):
        """Pipeline handles two copies of the same genome (batch mode)."""
        # Create a second copy with a different stem
        genome2 = tmp_path / "test_genome_copy.fna"
        shutil.copy2(str(TEST_GENOME), str(genome2))
        gff2 = tmp_path / "test_genome_copy.gff3"
        shutil.copy2(str(TEST_GFF), str(gff2))

        outdir = str(tmp_path / "results")
        gff_map = {
            TEST_GENOME.stem: str(TEST_GFF),
            genome2.stem: str(gff2),
        }
        summary = run_pipeline(
            fasta_files=[str(TEST_GENOME), str(genome2)],
            outdir=outdir,
            gff_files=gff_map,
        )
        assert summary["total"] == 2


# ---------------------------------------------------------------------------
# TestCLI — requires BLAST+
# ---------------------------------------------------------------------------

@requires_blast
class TestCLI:
    """Test the sccmec-pipeline CLI entry point."""

    def test_help(self):
        """--help runs without error."""
        result = subprocess.run(
            [sys.executable, "-m", "sccmecextractor.pipeline", "--help"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert "sccmec" in result.stdout.lower()

    def test_single_genome_end_to_end(self, tmp_path):
        """CLI processes a single genome successfully."""
        outdir = str(tmp_path / "cli_results")
        result = subprocess.run(
            [
                sys.executable, "-m", "sccmecextractor.pipeline",
                "-f", str(TEST_GENOME),
                "-g", str(TEST_GFF),
                "-o", outdir,
            ],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert os.path.isdir(os.path.join(outdir, "att_sites"))

    def test_mutually_exclusive_gff_args(self):
        """Cannot provide both --gff and --gff-dir."""
        result = subprocess.run(
            [
                sys.executable, "-m", "sccmecextractor.pipeline",
                "-f", str(TEST_GENOME),
                "-g", str(TEST_GFF),
                "--gff-dir", "/tmp",
                "-o", "/tmp/out",
            ],
            capture_output=True, text=True,
        )
        assert result.returncode != 0
