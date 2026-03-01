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
def test_tsv(test_data_dir):
    """Return path to test TSV file"""
    return test_data_dir / "expected_att_sites.tsv"

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


@pytest.fixture
def composite_tsv(temp_output_dir):
    """Synthetic att sites TSV with 3 sites on contig_1 (composite element).

    Layout (reverse orientation — rlmH at 2584176, downstream of attR):
        attL_outer(2400000) ... attL_inner(2493654) ... attR(2584168) ... rlmH(2584176)
    """
    tsv_path = temp_output_dir / "composite_att_sites.tsv"
    tsv_path.write_text(
        "Input_File\tPattern\tContig\tStart\tEnd\tMatching_Sequence\n"
        "test_genome\tcattR\tcontig_1\t2584168\t2584190\tAAACCGCATCATTTATGATATGC\n"
        "test_genome\tcattL\tcontig_1\t2493654\t2493676\tGCTTATCGGTTAATGATGCGGTT\n"
        "test_genome\tcattL\tcontig_1\t2400000\t2400022\tGCTTATCAGTTGATGATGCGGTT\n"
    )
    return tsv_path


@pytest.fixture
def partial_right_tsv(temp_output_dir):
    """TSV with only right att sites (no left sites)."""
    tsv_path = temp_output_dir / "partial_right.tsv"
    tsv_path.write_text(
        "Input_File\tPattern\tContig\tStart\tEnd\tMatching_Sequence\n"
        "test_genome\tcattR\tcontig_1\t2584168\t2584190\tAAACCGCATCATTTATGATATGC\n"
    )
    return tsv_path


@pytest.fixture
def cross_contig_tsv(temp_output_dir):
    """TSV with attR on contig_1 and attL on contig_2."""
    tsv_path = temp_output_dir / "cross_contig.tsv"
    tsv_path.write_text(
        "Input_File\tPattern\tContig\tStart\tEnd\tMatching_Sequence\n"
        "test_genome\tcattR\tcontig_1\t2584168\t2584190\tAAACCGCATCATTTATGATATGC\n"
        "test_genome\tcattL\tcontig_2\t5000\t5022\tGCTTATCGGTTAATGATGCGGTT\n"
    )
    return tsv_path


@pytest.fixture
def empty_sites_tsv(temp_output_dir):
    """TSV with header only — no att sites for test_genome."""
    tsv_path = temp_output_dir / "empty_sites.tsv"
    tsv_path.write_text(
        "Input_File\tPattern\tContig\tStart\tEnd\tMatching_Sequence\n"
    )
    return tsv_path


@pytest.fixture
def edge_case_tsv(temp_output_dir):
    """TSV with att site near contig start (position 200) on contig_2 (27310 bp)."""
    tsv_path = temp_output_dir / "edge_case.tsv"
    tsv_path.write_text(
        "Input_File\tPattern\tContig\tStart\tEnd\tMatching_Sequence\n"
        "test_genome\tattR\tcontig_2\t200\t222\tGCATATCATAAATGATGCGGTTT\n"
        "test_genome\tattL\tcontig_2\t26900\t26922\tAACCGCATCATCAACTGATAAGC\n"
    )
    return tsv_path


@pytest.fixture
def overlap_only_tsv(temp_output_dir):
    """TSV with overlapping att sites at attB (attR contains cattL9)."""
    tsv_path = temp_output_dir / "overlap_only.tsv"
    tsv_path.write_text(
        "Input_File\tPattern\tContig\tStart\tEnd\tMatching_Sequence\n"
        "test_genome\tcattR\tcontig_1\t2584168\t2584190\tAAACCGCATCATTTATGATATGC\n"
        "test_genome\tattL9\tcontig_1\t2584174\t2584189\tCATCATTTATGATAAG\n"
    )
    return tsv_path


@pytest.fixture
def overlap_with_valid_tsv(temp_output_dir):
    """TSV with overlapping attB pair + genuine attL at valid distance."""
    tsv_path = temp_output_dir / "overlap_with_valid.tsv"
    tsv_path.write_text(
        "Input_File\tPattern\tContig\tStart\tEnd\tMatching_Sequence\n"
        "test_genome\tcattR\tcontig_1\t2584168\t2584190\tAAACCGCATCATTTATGATATGC\n"
        "test_genome\tattL9\tcontig_1\t2584174\t2584189\tCATCATTTATGATAAG\n"
        "test_genome\tcattL\tcontig_1\t2493654\t2493676\tGCTTATCGGTTAATGATGCGGTT\n"
    )
    return tsv_path
