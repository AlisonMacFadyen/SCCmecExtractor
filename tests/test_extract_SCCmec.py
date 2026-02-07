#!/usr/bin/env python

"""
Tests for extract_SCCmec.py

This file contains tests that check if extract_SCCmec.py works correctly
"""

import pytest
import subprocess
import sys
from pathlib import Path
from sccmecextractor.extract_SCCmec import SCCmecExtractor
from sccmecextractor.extract_SCCmec import InputValidator

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