#!/usr/bin/env python

"""
Tests for locate_att_sites.py

This file contains tests that check if locate_att_sites.py works correctly
"""

import pytest
import subprocess
import sys
from pathlib import Path

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
                "./src/locate_att_sites.py",
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
        subprocess.run(
            ["python", 
             "./src/locate_att_sites.py",
             "-f", str(test_genome),
             "-g", str(test_gff),
             "-o", output_file
             ],
             capture_output=True
        )

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

        output_file = temp_output_dir / "att_sites.tsv"