#!/usr/bin/env python

"""SCCmec Reference Comparison BLAST against reference database.

Aim is to determine if mosaic/hybrid SCCmec elements exist via BLAST against bundled reference databases.

"""

import argparse
import os
import sys
import logging
import tempfile

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
from pathlib import Path
from Bio import SeqIO
from collections import defaultdict

from sccmecextractor.blast_utils import (
    BlastRunner, get_default_ref, parse_blast_output, filter_hits
)

#You'll need a function that BLASTs each extracted element against the 32 references —
#  look at how SCCmecTyper._blast_ref() works in sccmec_type_classification.py for the pattern
#   (references as query, element as database)
#  - Think about what the function should return — a dict or dataclass per element with the
# - The 32 reference FASTAs need to be bundled into src/sccmecextractor/data/ as a single
##  hit details, so the pipeline can easily merge it into the report
# consolidated file (like how mec_class_reference.fasta is bundled)
#  - Consider the coverage direction question from our earlier plan — you want "what % of each
#   reference is found in this element", not the other way around

class MosaicTyper:
    """Orchestrate BLAST-based SCC element mosaic/hybrid detection.
    """

    def __init__(self):
        self.runner = BlastRunner()

        with get_default_ref("sccmec_type_references.fasta") as ref:
            self.ref_lengths = {}
            for record in SeqIO.parse(str(ref), "fasta"):
                self.ref_lengths[record.id] = len(record.seq) # obtain the lengths of the reference fasta to enable coverage calculations

    def type_file(self, input_fasta: str, db_prefix: str = None) -> dict:
        """Type a single FASTA file (extracted element or whole genome).

        Parameters
        ----------
        input_fasta : str
            Path to the input FASTA file.
        db_prefix : str, optional
            Pre-built BLAST database prefix.  When provided the existing
            database is reused, avoiding a redundant ``makeblastdb`` call.
            The caller is responsible for cleanup.

        Returns a dict with keys matching TYPING_HEADER.
        """
        input_name = Path(input_fasta).stem

        for record in SeqIO.parse(str(input_fasta), "fasta"):
                self.input_length = len(record.seq) # obtain the length of the extracted scc element

        if db_prefix is not None:
            # Reuse caller-provided BLAST DB
            tmp_dir = tempfile.mkdtemp(prefix=f"{input_name}_reference_")
            try:
                element_reference_hits = self._blast_ref(db_prefix)
            finally:
                try:
                    os.rmdir(tmp_dir)
                except OSError:
                    pass
        else:
            # Standalone mode: create a temporary BLAST DB
            tmp_dir = tempfile.mkdtemp(prefix=f"{input_name}_reference")
            db_prefix_local = os.path.join(tmp_dir, f"{input_name}_reference_db")
            try:
                self.runner.create_db(input_fasta, db_prefix_local)
                element_reference_hits = self._blast_ref( db_prefix_local)
            finally:
                self.runner.cleanup_db(db_prefix_local)
                try:
                    os.rmdir(tmp_dir)
                except OSError:
                    pass

        comparison_results = self._calculate_coverage(element_reference_hits)

        # Format output
        return self._format_result(input_name, comparison_results)

    def _blast_ref(self, db_prefix: str):
        """BLAST a reference set against the SCCmec database."""
    
        with get_default_ref("sccmec_type_references.fasta") as ref:
            results_file = self.runner.run_blastn(str(ref), db_prefix)
            hits = parse_blast_output(results_file)
            self.runner.cleanup_file(results_file)
            return hits
    
    def _merge_intervals(self, intervals):
        """Merge any overlapping BLAST hits for SCCmec References"""
        if not intervals:
            return 0
    
        sorted_intervals = sorted(intervals)
        merged_start, merged_end = sorted_intervals[0]
        total_covered = 0

        for start, end in sorted_intervals[1:]:
            if start <= merged_end:
                # Overlaps - extend the current interval
                merged_end = max(merged_end, end)
            else:
                # No overlap - save current, start new
                total_covered += merged_end - merged_start
                merged_start, merged_end = start, end

        total_covered += merged_end - merged_start
        
        return total_covered

    def _calculate_coverage(self, sccmec_reference_results, min_length=1000):
        """Calculate coverage for SCCmec Type reference against extracted SCC element"""

        hits_by_ref = defaultdict(list)

        for hit in sccmec_reference_results:
            if hit.length >= min_length:
                hits_by_ref[hit.qseqid].append(hit)

        results = {}

        for ref_name, ref_hits in hits_by_ref.items():
            intervals = [(hit.qstart, hit.qend) for hit in ref_hits]
            covered_bp = self._merge_intervals(intervals)
            ref_len = self.ref_lengths[ref_name]
            coverage = covered_bp / ref_len
            hsp_num = len(intervals)
            weighted_pid = sum(hit.pident * hit.length for hit in ref_hits) / sum(hit.length for hit in ref_hits)

            results[ref_name] = {
                "coverage": round(coverage * 100, 1),
                "ref_length": ref_len,
                "covered_bp": covered_bp,
                "num_hsps": hsp_num,
                "weighted_pident": weighted_pid,
            }

        return results

    def _format_result(
        self,
        input_name: str,
        comparison_results: dict,
    ) -> dict:
        """Format typing results into a dict for TSV output.

        For each SCCmec Type reference we need to group sub-types.  The reference types we have available are:

        I (Ia, Ib)
        II (IIa, IIb, IIc, IId, IIe)
        III
        IV (IVa, IVb, IVc, IVd, IVg, IVi, IVj, IVk, IVl, IVm, IVn)
        V (Va, Vb, Vc)
        VI, VII, VIII, IX, X, XI, XII, XIII, XIV, XV

        """

        _SCCMEC_REFERENCE_TYPE_MAP = {
            ("B", "1"):  "I",
            ("A", "2"):  "II",
        }

        return {
            "extracted_element_id": input_name,
            "reference_type": sccmec_reference_name,
            "reference_accession": sccmec_reference_accession,
            "nucleotide_identity": weighted_pid,
            "reference_coverage": reference_coverage_perc,
            "input_coverage": element_coverage,
            "total_alignmnet_length": reference_covered_bp,
            "hsps_num": reference_hsp_num,
            "reference_length": ref_length,
            "element_length": self.input_length,
        }


def main():
    parser = argparse.ArgumentParser(
        description="Examine if extracted SCC elements are mosaic or hybrids of reference SCCmec Types"
    )
    parser.add_argument(
        "-f",
        "--fasta",
        nargs="+",
        required=True,
        help="Input SCC element FASTA file(s) or directory",
    )
    
    args = parser.parse_args()

    # Collect input files
    input_files = collect_input_files(args.fasta)

    if not input_files:
        print("ERROR: No FASTA files found in the provided paths")
        return

    print(f"Found {len(input_files)} input file(s)")

    # Create typer
    typer = MosaicTyper(element_ref=input_files)

    # Type each file
    with open(args.outfile, "w") as f:
        f.write("\t".join(TYPING_HEADER) + "\n")

        for i, input_file in enumerate(input_files, 1):
            print(f"  [{i}/{len(input_files)}] Typing {Path(input_file).name}...")
            try:
                result = typer.type_file(input_file)
                line = "\t".join(str(result[col]) for col in TYPING_HEADER)
                f.write(line + "\n")
            except Exception as e:
                print(f"    ERROR: {e}")
                line = "\t".join(
                    [Path(input_file).stem] + ["ERROR"] * (len(TYPING_HEADER) - 1)
                )
                f.write(line + "\n")

    print(f"\nResults written to {args.outfile}")

if __name__ == "__main__":
    main()