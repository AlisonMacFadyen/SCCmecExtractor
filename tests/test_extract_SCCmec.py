#!/usr/bin/env python

"""
Tests for extract_SCCmec.py

This file contains tests that check if extract_SCCmec.py works correctly
"""

import pytest
import shutil
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch
from sccmecextractor.extract_SCCmec import (
    SCCmecExtractor, InputValidator, AttSite, AttSiteCollection,
    ExtractionReport, RlmHBlastAdapter,
)

HAS_BLAST = shutil.which("blastn") is not None

# Add the parent directory to Python's path for import/run scripts
sys.path.insert(0, str(Path(__file__).parent))

class TestSCCmecExtractor:
    """
    A class to group related tests together
    """

    def test_script_runs(self, test_genome, test_gff, test_tsv, temp_output_dir):
        """
        Test that extract_SCCmec.py runs successfully without errors using valid input
        """

        # Define output location
        output_dir = temp_output_dir

        # Specify the output file
        output_file = output_dir / ""

        # Run the script
        result = subprocess.run(
            [
                "python",
                "-m",
                "sccmecextractor.extract_SCCmec",
                "-f", str(test_genome),
                "-g", str(test_gff),
                "-a", str(test_tsv),
                "-s", str(output_dir)
            ],
            capture_output=True,
            text=True
        )

        # Check that the script didn't crash
        assert result.returncode == 0, f"Script failed with error:\n{result.stderr}"

        # Check if output file was created
        assert output_file.exists(), "Output file was not created"

    def test_name_matches_locate_att_output(self, test_genome, test_gff, temp_output_dir):
        """
        extract_sccmec should find entries produced by locate_att, even when the filename contains multiple dots.
        Regression test for name mismatch bug.
        """
        # Create symlinks with RefSeq-style dotted names
        dotted_fna = temp_output_dir / "GCF_000159575.1_ASM15957v1_genomic.fna"
        dotted_fna.symlink_to(test_genome.resolve())

        dotted_gff = temp_output_dir / "GCF_000159575.1_ASM15957v1_genomic.gff3"
        dotted_gff.symlink_to(test_gff.resolve())

        att_file = temp_output_dir / "att_sites.tsv"
        sccmec_dir = temp_output_dir / "sccmec_out"

        # Step 1: locate att sites
        locate_result = subprocess.run(
            ["python", "-m", "sccmecextractor.locate_att_sites",
             "-f", str(dotted_fna),
             "-g", str(dotted_gff),
             "-o", str(att_file)],
            capture_output=True, text=True
        )
        assert locate_result.returncode == 0, f"locate_att failed: {locate_result.stderr}"

        # Step 2: extract using that output
        extract_result = subprocess.run(
            ["python", "-m", "sccmecextractor.extract_SCCmec",
             "-f", str(dotted_fna),
             "-g", str(dotted_gff),
             "-a", str(att_file),
             "-s", str(sccmec_dir)],
            capture_output=True, text=True
        )
        assert extract_result.returncode == 0, (
            f"extract failed (likely name mismatch): {extract_result.stderr}"
        )

    def test_extract_handles_duplicate_headers(self, test_genome, test_gff, test_tsv, temp_output_dir):
        """
        extract_sccmec should tolerate duplicate header lines in the TSV (e.g. from concatenated batch output).
        Regression test for duplicate header parsing crash.
        """
        # Create a TSV with a duplicate header in the middle
        with open(test_tsv) as f:
            original = f.read()

        bad_tsv = temp_output_dir / "dup_header.tsv"
        lines = original.strip().split("\n")
        # Insert a duplicate header after the first data line
        with open(bad_tsv, 'w') as f:
            f.write(lines[0] + "\n")  # header
            f.write(lines[1] + "\n")  # data line 1
            f.write(lines[0] + "\n")  # DUPLICATE header
            if len(lines) > 2:
                f.write(lines[2] + "\n")  # data line 2

        result = subprocess.run(
            ["python", "-m", "sccmecextractor.extract_SCCmec",
             "-f", str(test_genome),
             "-g", str(test_gff),
             "-a", str(bad_tsv),
             "-s", str(temp_output_dir / "sccmec_out")],
            capture_output=True, text=True
        )
        assert result.returncode == 0, (
            f"extract crashed on duplicate header: {result.stderr}"
        )


class TestInputValidator:
    """Test InputValidator"""
    
    def test_validate_fasta_accepts_files(self, test_genome):
        """
        Test input FASTA file accepted if valid
        """
        validator = InputValidator()
        
        validator.validate_fasta_file(test_genome)
        
    def test_validate_fasta_rejects_files(self):
        """
        Test input FASTA file is rejected if invalid
        """
        validator = InputValidator()
        
        fake_path = Path("this_file_does_not_exist.fasta")
    
        # Test to ensure FileNotFoundError occurs
        with pytest.raises(FileNotFoundError):
            validator.validate_fasta_file(fake_path)
    
    def test_validate_gff_accepts_files(self, test_gff):
        """
        Test input GFF file is accepted if valid
        """
        validator = InputValidator()
            
        validator.validate_gff_file(test_gff)
        
    def test_validate_gff_rejects_files(self):
        """
        Test input GFF file is rejected if not valid
        """
        
        validator = InputValidator()
        
        fake_path = Path("this_file_does_not_exist.gff")
    
        # Test to ensure FileNotFoundError occurs
        with pytest.raises(FileNotFoundError):
            validator.validate_gff_file(fake_path)

    def test_validate_tsv_accepts_files(self, test_tsv):
        """
        Test input GFF file is accepted if valid
        """
        validator = InputValidator()
            
        validator.validate_tsv_file(test_tsv)
        
    def test_validate_gff_rejects_files(self):
        """
        Test input GFF file is rejected if not valid
        """

        validator = InputValidator()

        fake_path = Path("this_file_does_not_exist.tsv")

        # Test to ensure FileNotFoundError occurs
        with pytest.raises(FileNotFoundError):
            validator.validate_gff_file(fake_path)


class TestCompositeDetection:
    """Tests for composite element detection and extraction."""

    def test_composite_detection(self, test_genome, test_gff, composite_tsv, temp_output_dir):
        """Composite TSV with composite=False, no outer ccr → Is_Composite=False, downgraded."""
        report_file = temp_output_dir / "report.tsv"
        extractor = SCCmecExtractor(
            str(test_genome), str(test_gff), str(composite_tsv), composite=False
        )
        success = extractor.extract_sccmec(str(temp_output_dir / "out"), report_file=str(report_file))
        assert success

        lines = report_file.read_text().strip().split("\n")
        assert len(lines) == 2  # header + 1 data row
        cols = lines[1].split("\t")
        header = lines[0].split("\t")
        row = dict(zip(header, cols))

        assert row["Status"] == "extracted"
        assert row["Is_Composite"] == "False"
        assert row["Outer_AttL_Pattern"] == "cattL"
        assert row["Outer_AttL_Start"] == "2400000"

    def test_composite_extraction(self, test_genome, test_gff, composite_tsv, temp_output_dir):
        """Composite TSV with composite=True but no ccr in outer region → downgraded to extracted."""
        report_file = temp_output_dir / "report.tsv"
        out_dir = temp_output_dir / "out"

        extractor = SCCmecExtractor(
            str(test_genome), str(test_gff), str(composite_tsv), composite=True
        )
        success = extractor.extract_sccmec(str(out_dir), report_file=str(report_file))
        assert success

        lines = report_file.read_text().strip().split("\n")
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert row["Status"] == "extracted"
        assert "Composite downgraded" in row["Notes"]

        # The extracted FASTA should exist
        fasta_files = list(out_dir.glob("*.fasta"))
        assert len(fasta_files) == 1

    def test_simple_element_not_composite(self, test_genome, test_gff, test_tsv, temp_output_dir):
        """Standard test data (2 sites) → Is_Composite=False."""
        report_file = temp_output_dir / "report.tsv"
        extractor = SCCmecExtractor(
            str(test_genome), str(test_gff), str(test_tsv), composite=False
        )
        success = extractor.extract_sccmec(str(temp_output_dir / "out"), report_file=str(report_file))
        assert success

        lines = report_file.read_text().strip().split("\n")
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert row["Is_Composite"] == "False"
        assert row["Status"] == "extracted"
        assert row["Outer_AttL_Pattern"] == "-"

    def test_composite_cli_flag(self, test_genome, test_gff, composite_tsv, temp_output_dir):
        """CLI invocation with --composite --report: no ccr in outer → downgraded."""
        report_file = temp_output_dir / "cli_report.tsv"
        out_dir = temp_output_dir / "cli_out"

        result = subprocess.run(
            [
                "python", "-m", "sccmecextractor.extract_SCCmec",
                "-f", str(test_genome),
                "-g", str(test_gff),
                "-a", str(composite_tsv),
                "-s", str(out_dir),
                "--composite",
                "-r", str(report_file),
            ],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, f"CLI failed: {result.stderr}"
        assert report_file.exists()

        lines = report_file.read_text().strip().split("\n")
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert row["Status"] == "extracted"
        assert "Composite downgraded" in row["Notes"]


class TestCompositeCCRValidation:
    """Tests for composite CCR validation scenarios using mocked _detect_ccr_between."""

    @staticmethod
    def _ccr_side_effect(inner_result, outer_result):
        """Return a side_effect function that returns inner_result on first call,
        outer_result on second call."""
        call_count = [0]
        def side_effect(contig, pos_a, pos_b):
            call_count[0] += 1
            if call_count[0] == 1:
                return inner_result
            return outer_result
        return side_effect

    @patch.object(SCCmecExtractor, '_detect_ccr_between')
    def test_genuine_composite_extraction(self, mock_ccr, test_genome, test_gff,
                                          composite_tsv, temp_output_dir):
        """Both inner and outer regions have ccr → composite_extracted."""
        mock_ccr.side_effect = self._ccr_side_effect(True, True)
        report_file = temp_output_dir / "report.tsv"
        out_dir = temp_output_dir / "out"

        extractor = SCCmecExtractor(
            str(test_genome), str(test_gff), str(composite_tsv), composite=True
        )
        success = extractor.extract_sccmec(str(out_dir), report_file=str(report_file))
        assert success

        lines = report_file.read_text().strip().split("\n")
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert row["Status"] == "composite_extracted"
        assert row["Is_Composite"] == "True"
        assert row["Outer_AttL_Pattern"] == "cattL"

    @patch.object(SCCmecExtractor, '_detect_ccr_between')
    def test_composite_downgrade(self, mock_ccr, test_genome, test_gff,
                                 composite_tsv, temp_output_dir):
        """Inner ccr but no outer ccr → extracted with downgrade note."""
        mock_ccr.side_effect = self._ccr_side_effect(True, False)
        report_file = temp_output_dir / "report.tsv"
        out_dir = temp_output_dir / "out"

        extractor = SCCmecExtractor(
            str(test_genome), str(test_gff), str(composite_tsv), composite=True
        )
        success = extractor.extract_sccmec(str(out_dir), report_file=str(report_file))
        assert success

        lines = report_file.read_text().strip().split("\n")
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert row["Status"] == "extracted"
        assert "Composite downgraded" in row["Notes"]

    @patch.object(SCCmecExtractor, '_detect_ccr_between')
    def test_scar_bypass_extraction(self, mock_ccr, test_genome, test_gff,
                                    composite_tsv, temp_output_dir):
        """No inner ccr but outer ccr → extracted to outer attL (scar bypass)."""
        mock_ccr.side_effect = self._ccr_side_effect(False, True)
        report_file = temp_output_dir / "report.tsv"
        out_dir = temp_output_dir / "out"

        extractor = SCCmecExtractor(
            str(test_genome), str(test_gff), str(composite_tsv), composite=False
        )
        success = extractor.extract_sccmec(str(out_dir), report_file=str(report_file))
        assert success

        lines = report_file.read_text().strip().split("\n")
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert row["Status"] == "extracted"
        assert "Inner attL bypassed" in row["Notes"]
        assert row["Is_Composite"] == "False"
        # attL fields should point to outer attL
        assert row["AttL_Start"] == "2400000"
        assert row["Outer_AttL_Pattern"] == "-"

    @patch.object(SCCmecExtractor, '_detect_ccr_between')
    def test_no_ccr_with_outer_left(self, mock_ccr, test_genome, test_gff,
                                    composite_tsv, temp_output_dir):
        """No ccr in inner or outer region → no_ccr_element."""
        mock_ccr.side_effect = self._ccr_side_effect(False, False)
        report_file = temp_output_dir / "report.tsv"
        out_dir = temp_output_dir / "out"

        extractor = SCCmecExtractor(
            str(test_genome), str(test_gff), str(composite_tsv), composite=False
        )
        success = extractor.extract_sccmec(str(out_dir), report_file=str(report_file))
        assert not success

        lines = report_file.read_text().strip().split("\n")
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert row["Status"] == "no_ccr_element"
        assert "no ccr genes detected" in row["Notes"]


class TestPartialReporting:
    """Tests for partial element reporting when extraction fails."""

    @patch.object(SCCmecExtractor, '_has_ccr_in_genome', return_value=True)
    def test_diagnose_partial_right_only(self, mock_ccr, test_genome, test_gff, partial_right_tsv, temp_output_dir):
        """Right-only TSV with ccr present → partial_type='right_only'."""
        report_file = temp_output_dir / "report.tsv"
        extractor = SCCmecExtractor(
            str(test_genome), str(test_gff), str(partial_right_tsv)
        )
        success = extractor.extract_sccmec(str(temp_output_dir / "out"), report_file=str(report_file))
        assert not success

        lines = report_file.read_text().strip().split("\n")
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert row["Status"] == "failed"
        assert row["Partial_Type"] == "right_only"

    @patch.object(SCCmecExtractor, '_has_ccr_in_genome', return_value=True)
    def test_diagnose_partial_cross_contig(self, mock_ccr, test_genome, test_gff, cross_contig_tsv, temp_output_dir):
        """attR on contig_1, attL on contig_2 with ccr present → partial_type='cross_contig'."""
        report_file = temp_output_dir / "report.tsv"
        extractor = SCCmecExtractor(
            str(test_genome), str(test_gff), str(cross_contig_tsv)
        )
        success = extractor.extract_sccmec(str(temp_output_dir / "out"), report_file=str(report_file))
        assert not success

        lines = report_file.read_text().strip().split("\n")
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert row["Status"] == "failed"
        assert row["Partial_Type"] == "cross_contig"

    @patch.object(SCCmecExtractor, '_has_ccr_in_genome', return_value=False)
    def test_diagnose_partial_no_sites(self, mock_ccr, test_genome, test_gff, empty_sites_tsv, temp_output_dir):
        """Empty TSV, no ccr → partial_type='no_sites'."""
        report_file = temp_output_dir / "report.tsv"
        extractor = SCCmecExtractor(
            str(test_genome), str(test_gff), str(empty_sites_tsv)
        )
        success = extractor.extract_sccmec(str(temp_output_dir / "out"), report_file=str(report_file))
        assert not success

        lines = report_file.read_text().strip().split("\n")
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert row["Status"] == "failed"
        assert row["Partial_Type"] == "no_sites"


class TestNoCcrReporting:
    """Tests for no_ccr partial_type when ccr genes are absent."""

    @patch.object(SCCmecExtractor, '_has_ccr_in_genome', return_value=False)
    def test_no_ccr_right_only(self, mock_ccr, test_genome, test_gff, partial_right_tsv, temp_output_dir):
        """Right-only att sites but no ccr → partial_type='no_ccr', no ambiguous report."""
        report_file = temp_output_dir / "report.tsv"
        amb_file = temp_output_dir / "ambiguous.tsv"
        extractor = SCCmecExtractor(
            str(test_genome), str(test_gff), str(partial_right_tsv)
        )
        success = extractor.extract_sccmec(
            str(temp_output_dir / "out"), report_file=str(report_file),
            ambiguous_report_file=str(amb_file),
        )
        assert not success

        lines = report_file.read_text().strip().split("\n")
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert row["Status"] == "failed"
        assert row["Partial_Type"] == "no_ccr"
        assert "right_only" in row["Notes"]
        assert not amb_file.exists()

    @patch.object(SCCmecExtractor, '_has_ccr_in_genome', return_value=False)
    def test_no_ccr_cross_contig(self, mock_ccr, test_genome, test_gff, cross_contig_tsv, temp_output_dir):
        """Cross-contig att sites but no ccr → partial_type='no_ccr', no ambiguous report."""
        report_file = temp_output_dir / "report.tsv"
        amb_file = temp_output_dir / "ambiguous.tsv"
        extractor = SCCmecExtractor(
            str(test_genome), str(test_gff), str(cross_contig_tsv)
        )
        success = extractor.extract_sccmec(
            str(temp_output_dir / "out"), report_file=str(report_file),
            ambiguous_report_file=str(amb_file),
        )
        assert not success

        lines = report_file.read_text().strip().split("\n")
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert row["Status"] == "failed"
        assert row["Partial_Type"] == "no_ccr"
        assert "cross_contig" in row["Notes"]
        assert not amb_file.exists()

    @patch.object(SCCmecExtractor, '_has_ccr_in_genome', return_value=True)
    def test_ccr_present_right_only_writes_ambiguous(self, mock_ccr, test_genome, test_gff, partial_right_tsv, temp_output_dir):
        """Right-only with ccr present → partial_type='right_only', ambiguous report written."""
        report_file = temp_output_dir / "report.tsv"
        amb_file = temp_output_dir / "ambiguous.tsv"
        extractor = SCCmecExtractor(
            str(test_genome), str(test_gff), str(partial_right_tsv)
        )
        success = extractor.extract_sccmec(
            str(temp_output_dir / "out"), report_file=str(report_file),
            ambiguous_report_file=str(amb_file),
        )
        assert not success

        lines = report_file.read_text().strip().split("\n")
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert row["Status"] == "failed"
        assert row["Partial_Type"] == "right_only"
        assert amb_file.exists()

    @patch.object(SCCmecExtractor, '_has_ccr_in_genome', return_value=True)
    def test_no_sites_ccr_present_writes_ambiguous(self, mock_ccr, test_genome, test_gff, empty_sites_tsv, temp_output_dir):
        """No att sites but ccr present → partial_type='no_sites', ambiguous report written."""
        report_file = temp_output_dir / "report.tsv"
        amb_file = temp_output_dir / "ambiguous.tsv"
        extractor = SCCmecExtractor(
            str(test_genome), str(test_gff), str(empty_sites_tsv)
        )
        success = extractor.extract_sccmec(
            str(temp_output_dir / "out"), report_file=str(report_file),
            ambiguous_report_file=str(amb_file),
        )
        assert not success

        lines = report_file.read_text().strip().split("\n")
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert row["Status"] == "failed"
        assert row["Partial_Type"] == "no_sites"
        assert "ccr genes present" in row["Notes"]
        assert amb_file.exists()


class TestContigEdge:
    """Tests for contig-edge proximity flagging."""

    def test_contig_edge_near_start(self):
        """Site at position 100 on a 50000 bp contig → near start."""
        site = AttSite("attR", "c1", 100, 122)
        assert SCCmecExtractor._check_contig_edge(site, 50000) is True

    def test_contig_edge_near_end(self):
        """Site near the end of a 50000 bp contig → near end."""
        site = AttSite("attL", "c1", 49600, 49622)
        assert SCCmecExtractor._check_contig_edge(site, 50000) is True

    def test_contig_edge_not_near_boundary(self):
        """Site at 25000 on a 50000 bp contig → not near edge."""
        site = AttSite("attR", "c1", 25000, 25022)
        assert SCCmecExtractor._check_contig_edge(site, 50000) is False

    def test_edge_flags_in_report(self, test_genome, test_gff, edge_case_tsv, temp_output_dir):
        """Edge-case TSV with sites near contig boundaries → flags in report."""
        report_file = temp_output_dir / "report.tsv"
        extractor = SCCmecExtractor(
            str(test_genome), str(test_gff), str(edge_case_tsv)
        )
        # This will fail extraction (no rlmH on contig_2) but should still report edge flags
        success = extractor.extract_sccmec(str(temp_output_dir / "out"), report_file=str(report_file))
        assert not success

        lines = report_file.read_text().strip().split("\n")
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert "near_start" in row["Contig_Edge_Flags"]
        assert "near_end" in row["Contig_Edge_Flags"]


class TestReportGeneration:
    """Tests for report TSV generation."""

    def test_report_generated_on_success(self, test_genome, test_gff, test_tsv, temp_output_dir):
        """Successful extraction → report TSV with correct columns."""
        report_file = temp_output_dir / "report.tsv"
        extractor = SCCmecExtractor(str(test_genome), str(test_gff), str(test_tsv))
        success = extractor.extract_sccmec(str(temp_output_dir / "out"), report_file=str(report_file))
        assert success
        assert report_file.exists()

        lines = report_file.read_text().strip().split("\n")
        assert len(lines) == 2
        header_cols = lines[0].split("\t")
        assert header_cols[0] == "Input_File"
        assert header_cols[1] == "Status"
        assert len(header_cols) == 18

    @patch.object(SCCmecExtractor, '_has_ccr_in_genome', return_value=True)
    def test_report_generated_on_failure(self, mock_ccr, test_genome, test_gff, partial_right_tsv, temp_output_dir):
        """Failed extraction → report TSV with failed status."""
        report_file = temp_output_dir / "report.tsv"
        extractor = SCCmecExtractor(str(test_genome), str(test_gff), str(partial_right_tsv))
        success = extractor.extract_sccmec(str(temp_output_dir / "out"), report_file=str(report_file))
        assert not success
        assert report_file.exists()

        lines = report_file.read_text().strip().split("\n")
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert row["Status"] == "failed"

    @patch.object(SCCmecExtractor, '_has_ccr_in_genome', return_value=True)
    def test_report_appends_for_batch(self, mock_ccr, test_genome, test_gff, test_tsv, partial_right_tsv, temp_output_dir):
        """Two extractions appended to one report → 1 header + 2 data rows."""
        report_file = temp_output_dir / "batch_report.tsv"

        # First extraction — success
        ext1 = SCCmecExtractor(str(test_genome), str(test_gff), str(test_tsv))
        ext1.extract_sccmec(str(temp_output_dir / "out1"), report_file=str(report_file))

        # Second extraction — failure (right-only)
        ext2 = SCCmecExtractor(str(test_genome), str(test_gff), str(partial_right_tsv))
        ext2.extract_sccmec(str(temp_output_dir / "out2"), report_file=str(report_file))

        lines = report_file.read_text().strip().split("\n")
        assert len(lines) == 3  # 1 header + 2 data rows
        assert lines[0].startswith("Input_File")
        assert "extracted" in lines[1]
        assert "failed" in lines[2]

    def test_no_report_when_not_requested(self, test_genome, test_gff, test_tsv, temp_output_dir):
        """When report_file is None, no report is created."""
        extractor = SCCmecExtractor(str(test_genome), str(test_gff), str(test_tsv))
        success = extractor.extract_sccmec(str(temp_output_dir / "out"))
        assert success
        # No report files should exist in output dir
        assert not list(temp_output_dir.glob("*.tsv"))


class TestNoRlmHSource:
    """Test that missing GFF and BLAST raises an error."""

    def test_no_gff_no_blast_raises(self, test_genome, test_tsv):
        """ValueError raised when neither GFF nor blast_rlmh is provided."""
        with pytest.raises(ValueError, match="Either a GFF3 file"):
            SCCmecExtractor(str(test_genome), tsv_file=str(test_tsv))


@pytest.mark.skipif(not HAS_BLAST, reason="BLAST+ not installed")
class TestBlastRlmH:
    """Tests for BLAST-based rlmH detection in extraction."""

    def test_blast_adapter_finds_rlmh(self, test_genome):
        """RlmHBlastAdapter detects rlmH in the test genome."""
        adapter = RlmHBlastAdapter(str(test_genome))
        assert len(adapter.rlmH_positions) > 0, "BLAST should find rlmH"

    def test_blast_adapter_matches_gff(self, test_genome, test_gff):
        """BLAST and GFF find rlmH on the same contigs."""
        from sccmecextractor.extract_SCCmec import GeneAnnotations
        blast_adapter = RlmHBlastAdapter(str(test_genome))
        gff_genes = GeneAnnotations(str(test_gff))

        blast_contigs = set(blast_adapter.rlmH_positions.keys())
        gff_contigs = set(gff_genes.rlmH_positions.keys())
        assert gff_contigs == blast_contigs, (
            f"GFF found rlmH on {gff_contigs}, BLAST found on {blast_contigs}"
        )

    def test_extract_with_blast_rlmh(self, test_genome, test_tsv, temp_output_dir):
        """Extraction succeeds using BLAST rlmH (no GFF)."""
        report_file = temp_output_dir / "report.tsv"
        extractor = SCCmecExtractor(
            str(test_genome), tsv_file=str(test_tsv), blast_rlmh=True,
        )
        success = extractor.extract_sccmec(
            str(temp_output_dir / "out"), report_file=str(report_file),
        )
        assert success

        lines = report_file.read_text().strip().split("\n")
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert row["Status"] == "extracted"

    def test_blast_matches_gff_extraction(self, test_genome, test_gff, test_tsv, temp_output_dir):
        """BLAST-based and GFF-based extraction produce the same result."""
        # GFF-based
        gff_dir = temp_output_dir / "gff_out"
        gff_report = temp_output_dir / "gff_report.tsv"
        gff_ext = SCCmecExtractor(str(test_genome), str(test_gff), str(test_tsv))
        gff_ext.extract_sccmec(str(gff_dir), report_file=str(gff_report))

        # BLAST-based
        blast_dir = temp_output_dir / "blast_out"
        blast_report = temp_output_dir / "blast_report.tsv"
        blast_ext = SCCmecExtractor(
            str(test_genome), tsv_file=str(test_tsv), blast_rlmh=True,
        )
        blast_ext.extract_sccmec(str(blast_dir), report_file=str(blast_report))

        # Compare reports
        gff_lines = gff_report.read_text().strip().split("\n")
        blast_lines = blast_report.read_text().strip().split("\n")
        gff_row = dict(zip(gff_lines[0].split("\t"), gff_lines[1].split("\t")))
        blast_row = dict(zip(blast_lines[0].split("\t"), blast_lines[1].split("\t")))

        assert gff_row["Status"] == blast_row["Status"]
        assert gff_row["Contig"] == blast_row["Contig"]
        assert gff_row["AttR_Pattern"] == blast_row["AttR_Pattern"]
        assert gff_row["AttL_Pattern"] == blast_row["AttL_Pattern"]

        # Both should produce FASTA output
        gff_fastas = list(gff_dir.glob("*.fasta"))
        blast_fastas = list(blast_dir.glob("*.fasta"))
        assert len(gff_fastas) == 1
        assert len(blast_fastas) == 1

    def test_cli_blast_rlmh(self, test_genome, test_tsv, temp_output_dir):
        """CLI --blast-rlmh flag works for extraction."""
        out_dir = temp_output_dir / "cli_out"
        report_file = temp_output_dir / "cli_report.tsv"

        result = subprocess.run(
            [
                "python", "-m", "sccmecextractor.extract_SCCmec",
                "-f", str(test_genome),
                "-a", str(test_tsv),
                "-s", str(out_dir),
                "--blast-rlmh",
                "-r", str(report_file),
            ],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, f"CLI failed: {result.stderr}"
        assert report_file.exists()
        assert list(out_dir.glob("*.fasta"))

    def test_cli_auto_enables_blast(self, test_genome, test_tsv, temp_output_dir):
        """CLI auto-enables BLAST when no --gff provided."""
        out_dir = temp_output_dir / "auto_out"
        report_file = temp_output_dir / "auto_report.tsv"

        result = subprocess.run(
            [
                "python", "-m", "sccmecextractor.extract_SCCmec",
                "-f", str(test_genome),
                "-a", str(test_tsv),
                "-s", str(out_dir),
                "-r", str(report_file),
            ],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, f"CLI failed: {result.stderr}"
        assert "auto-enabling BLAST" in result.stdout
        assert report_file.exists()


class TestAttBOverlap:
    """Tests for att site overlap detection at chromosomal attB."""

    def test_overlap_only_fails(self, test_genome, test_gff, overlap_only_tsv, temp_output_dir):
        """Overlapping att sites with no alternative → fails with size_out_of_range."""
        report_file = temp_output_dir / "report.tsv"
        extractor = SCCmecExtractor(str(test_genome), str(test_gff), str(overlap_only_tsv))
        success = extractor.extract_sccmec(str(temp_output_dir / "out"), report_file=str(report_file))
        assert not success
        lines = report_file.read_text().strip().split("\n")
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert row["Status"] == "failed"
        assert row["Partial_Type"] == "size_out_of_range"
        assert "Overlapping" in row["Notes"]

    def test_overlap_skipped_for_valid_pair(self, test_genome, test_gff, overlap_with_valid_tsv, temp_output_dir):
        """Overlapping pair skipped, genuine attL selected → extraction succeeds."""
        report_file = temp_output_dir / "report.tsv"
        extractor = SCCmecExtractor(str(test_genome), str(test_gff), str(overlap_with_valid_tsv))
        success = extractor.extract_sccmec(str(temp_output_dir / "out"), report_file=str(report_file))
        assert success
        lines = report_file.read_text().strip().split("\n")
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
        assert row["Status"] == "extracted"
        assert row["AttL_Pattern"] == "cattL"  # Genuine pair, not the overlapping attL9

    def test_find_closest_pair_min_distance(self):
        """find_closest_pair with min_distance skips close pairs."""
        sites = AttSiteCollection.__new__(AttSiteCollection)
        sites.sites = [
            AttSite("attR", "c1", 1000, 1023),
            AttSite("cattL9", "c1", 1005, 1020),   # Overlapping (distance ~3)
            AttSite("cattL", "c1", 50000, 50022),   # Valid (distance ~49000)
        ]
        # Without min_distance: returns overlapping pair
        pair = sites.find_closest_pair()
        assert pair is not None
        assert pair[1].pattern == "cattL9"

        # With min_distance: skips overlapping, returns valid pair
        pair = sites.find_closest_pair(min_distance=1000)
        assert pair is not None
        assert pair[1].pattern == "cattL"