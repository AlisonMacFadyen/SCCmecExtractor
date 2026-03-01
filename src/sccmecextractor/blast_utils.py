#!/usr/bin/env python

"""Shared BLAST utilities for SCCmecExtractor.

Provides BlastRunner for executing BLAST+ commands, BlastResult for parsed hits
and helper functions for filtering and resolving overlapping hits.
"""

import os
import csv
import shutil
import subprocess
import tempfile

from contextlib import contextmanager
from dataclasses import dataclass
from importlib.resources import files, as_file
from pathlib import Path
from typing import Dict, List, Optional, Tuple


class BlastNotFoundError(RuntimeError):
    """Raised when BLAST+ executables are not found on PATH."""

    pass


@dataclass
class BlastResult:
    """Represents a single BLAST outfmt 6 hit (all 12 standard columns)."""

    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float


def parse_blast_output(results_file: str) -> List[BlastResult]:
    """Parse a BLAST outfmt 6 results file into BlastResult objects."""
    hits = []
    path = Path(results_file)

    if not path.exists() or path.stat().st_size == 0:
        return hits

    with open(results_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 12:
                continue
            hits.append(
                BlastResult(
                    qseqid=row[0],
                    sseqid=row[1],
                    pident=float(row[2]),
                    length=int(row[3]),
                    mismatch=int(row[4]),
                    gapopen=int(row[5]),
                    qstart=int(row[6]),
                    qend=int(row[7]),
                    sstart=int(row[8]),
                    send=int(row[9]),
                    evalue=float(row[10]),
                    bitscore=float(row[11]),
                )
            )

    return hits


def filter_hits(
    hits: List[BlastResult],
    min_pident: float,
    min_coverage: float,
    ref_lengths: Dict[str, int],
) -> List[BlastResult]:
    """Filter BLAST hits by percent identity and query coverage.

    Args:
        hits: List of BlastResult objects.
        min_pident: Minimum percent identity threshold.
        min_coverage: Minimum query coverage fraction (0-1).
        ref_lengths: Dict mapping query sequence IDs to their lengths.

    Returns:
        Filtered list of BlastResult objects.
    """
    filtered = []

    for hit in hits:
        if hit.pident < min_pident:
            continue

        ref_len = ref_lengths.get(hit.qseqid)
        if ref_len is None:
            continue

        coverage = hit.length / ref_len
        if coverage < min_coverage:
            continue

        filtered.append(hit)

    return filtered


def get_best_non_overlapping_hits(
    hits: List[BlastResult], overlap_threshold: int = 5000
) -> List[BlastResult]:
    """Resolve overlapping BLAST hits, keeping the best (highest pident, then longest).

    Hits are sorted by bitscore descending, then greedily selected if they don't
    overlap with already-selected hits on the same subject sequence.

    Args:
        hits: List of BlastResult objects.
        overlap_threshold: Maximum allowed overlap in bp before hits are
            considered conflicting.

    Returns:
        Non-overlapping list of BlastResult objects.
    """
    if not hits:
        return []

    sorted_hits = sorted(hits, key=lambda h: (-h.bitscore, -h.length))

    selected = []

    for hit in sorted_hits:
        hit_start = min(hit.sstart, hit.send)
        hit_end = max(hit.sstart, hit.send)

        overlaps = False
        for sel in selected:
            if sel.sseqid != hit.sseqid:
                continue

            sel_start = min(sel.sstart, sel.send)
            sel_end = max(sel.sstart, sel.send)

            overlap = min(hit_end, sel_end) - max(hit_start, sel_start)
            if overlap > overlap_threshold:
                overlaps = True
                break

        if not overlaps:
            selected.append(hit)

    return selected


@contextmanager
def get_default_ref(filename: str):
    """Locate a bundled reference FASTA and yield a real filesystem path.

    Uses importlib.resources to locate data files bundled with the package.
    The context manager ensures temporary extractions are cleaned up.

    Args:
        filename: Name of the file in sccmecextractor/data/ (e.g. 'rlmH.fasta').

    Yields:
        Path to the reference file on the real filesystem.
    """
    ref = files("sccmecextractor").joinpath("data", filename)
    with as_file(ref) as path:
        yield path


class BlastRunner:
    """Wrapper for BLAST+ command-line tools."""

    def __init__(self):
        self._check_blast_installed()

    @staticmethod
    def _check_blast_installed():
        """Verify that blastn and makeblastdb are available on PATH."""
        for tool in ("blastn", "makeblastdb"):
            if shutil.which(tool) is None:
                raise BlastNotFoundError(
                    f"'{tool}' not found on PATH. "
                    "Please install BLAST+ (https://blast.ncbi.nlm.nih.gov/)"
                )

    def create_db(self, fasta_path: str, db_prefix: str) -> str:
        """Create a BLAST nucleotide database.

        Args:
            fasta_path: Path to the input FASTA file.
            db_prefix: Output database prefix path.

        Returns:
            The db_prefix path for use in subsequent blastn calls.
        """
        cmd = [
            "makeblastdb",
            "-in",
            str(fasta_path),
            "-dbtype",
            "nucl",
            "-out",
            str(db_prefix),
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return db_prefix

    def run_blastn(
        self,
        query: str,
        db: str,
        output: Optional[str] = None,
        evalue: float = 10,
        word_size: int = 28,
    ) -> str:
        """Run blastn with project-standard parameters.

        Args:
            query: Path to the query FASTA file.
            db: Path to the BLAST database prefix.
            output: Path for the output file. If None, a temp file is created.
            evalue: E-value threshold.
            word_size: Word size for blastn.

        Returns:
            Path to the output results file.
        """
        if output is None:
            fd, output = tempfile.mkstemp(suffix=".blast6", prefix="sccmec_")
            os.close(fd)

        cmd = [
            "blastn",
            "-query",
            str(query),
            "-db",
            str(db),
            "-out",
            str(output),
            "-evalue",
            str(evalue),
            "-word_size",
            str(word_size),
            "-gapopen",
            "0",
            "-gapextend",
            "2",
            "-penalty",
            "-2",
            "-reward",
            "1",
            "-task",
            "blastn",
            "-outfmt",
            "6",
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return output

    @staticmethod
    def cleanup_db(db_prefix: str):
        """Remove BLAST database files created by makeblastdb."""
        extensions = [".nin", ".nhr", ".nsq", ".ndb", ".nto", ".ntf", ".not"]
        for ext in extensions:
            path = Path(str(db_prefix) + ext)
            if path.exists():
                path.unlink()

    @staticmethod
    def cleanup_file(filepath: str):
        """Remove a single file if it exists."""
        path = Path(filepath)
        if path.exists():
            path.unlink()
