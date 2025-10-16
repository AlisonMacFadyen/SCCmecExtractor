#!/usr/bin/env python

import pytest
import tempfile
import shutil
from pathlib import Path

"""
Pytest configuration and fixtures for SCCmecExtractor tests
"""

@pytest.fixture
def test_data_dir():
    """Return path to test data directory"""
    return Path(__file__).parent / "test_data"

@pytest.fixture
def test_genome(test_data_dir):
    """Return path to test genome file"""
    return test_data_dir / "test_genome.fna"

@pytest.fixture
def test_gff(test_data_dir):
    """Return path to test GFF file"""
    return test_data_dir / "test_genome.gff3"

@pytest.fixture
def expected_att_sites(test_data_dir):
    """Return path to expected att sites file"""
    return test_data_dir / "expected_att_sites.tsv"

@pytest.fixture
def temp_output_dir():
    """Create a temporary directory for test outputs"""
    temp_dir = tempfile.mkdtemp()
    yield Path(temp_dir)
    shutil.rmtree(temp_dir)
