#!/usr/bin/env python

"""
Tests for locate_att_sites.py

This file contains tests that check if locate_att_sites.py works correctly
"""

import pytest
import subprocess
import sys
import pandas as pd
from pathlib import Path
from sccmecextractor.locate_att_sites import AttSiteFinder
from sccmecextractor.locate_att_sites import InputValidator

# Add the parent directory to Python's path for import/run scripts
sys.path.insert(0, str(Path(__file__).parent))

class TestLocateAttSites:
    """
    A class to group related tests together
    """

    def test_script_runs(self, test_genome, test_gff, temp_output_dir):
        """
        Test that locate_att_sites.py runs successfully without errors using valid input
        """

        # Define output location
        output_file = temp_output_dir / "att_sites.tsv"

        # Run the script
        result = subprocess.run(
            [
                "python",
                "-m", 
                "sccmecextractor.locate_att_sites",
                "-f", str(test_genome),
                "-g", str(test_gff),
                "-o", str(output_file)
            ],
            capture_output=True,
            text=True
        )

        # Check that the script didn't crash
        assert result.returncode == 0, f"Script failed with error:\n{result.stderr}"

        # Check if output file was created
        assert output_file.exists(), "Output file was not created"

    def test_output_format(self, test_genome, test_gff, temp_output_dir):
        """
        Test that the output TSV has the expected layout
        """

        output_file = temp_output_dir / "att_sites.tsv"

        # Run the script
        result = subprocess.run(
            ["python", 
             "-m", 
             "sccmecextractor.locate_att_sites",
             "-f", str(test_genome),
             "-g", str(test_gff),
             "-o", output_file
             ],
             capture_output=True
        )

        # Check the file was created
        assert result.returncode == 0, f"Script failed: {result.stderr}"
        assert output_file.exists(), "Output file not created"

        # Read the resulting TSV file
        with open(output_file) as out_file:
            header_line = out_file.readline().strip()

        # Check the header
        expected_columns = ['Input_File', 'Pattern', 'Contig', 'Start', 'End', 'Matching_Sequence']
        actual_columns = header_line.split('\t')

        assert actual_columns == expected_columns, \
            f"Expected columns {expected_columns}, got {actual_columns}"
        
    def test_att_sites_total(self, test_genome, test_gff, temp_output_dir):
        """
        Test that att sites are found and a realistic number
        """
        
        finder = AttSiteFinder(str(test_genome), str(test_gff))
        all_sites = finder.find_all_sites()
        filtered_sites = finder.filter_sites(all_sites)
        
        # Test we have obtained at least one site and no more than 10
        assert 1 <= len(filtered_sites) <= 9, \
            f"Expected 1-9 sites, found {len(filtered_sites)}"

    def test_output_content_matches_expected(self, test_genome, test_gff, temp_output_dir):
        """
        Test that output contains expected att sites
        """
        output_file = temp_output_dir / "att_sites.tsv"
        expected_file = Path("tests/test_data/expected_att_sites.tsv")
    
        # Run the script
        subprocess.run(
            ["python", "-m", "sccmecextractor.locate_att_sites",
             "-f", str(test_genome),
             "-g", str(test_gff),
             "-o", str(output_file)],
            capture_output=True
        )
    
        # Load both files
        actual_df = pd.read_csv(output_file, sep='\t')
        expected_df = pd.read_csv(expected_file, sep='\t')
    
        # Compare DataFrames
        pd.testing.assert_frame_equal(
            actual_df.sort_values(by=['Contig', 'Start']).reset_index(drop=True),
            expected_df.sort_values(by=['Contig', 'Start']).reset_index(drop=True),
            check_like=True  # Ignore column/row order
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