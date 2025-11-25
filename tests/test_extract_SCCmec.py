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