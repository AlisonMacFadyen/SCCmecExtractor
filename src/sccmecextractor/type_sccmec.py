#!/usr/bin/env python

"""SCCmec typing by mec and ccr gene content using BLAST.

Classifies extracted SCCmec sequences by identifying mec and ccr genes
via BLAST against bundled reference databases.
"""

import argparse
import os
import tempfile

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from Bio import SeqIO

from sccmecextractor.blast_utils import (
    BlastRunner,
    filter_hits,
    get_best_non_overlapping_hits,
    get_default_ref,
    parse_blast_output,
)


@dataclass
class GeneHit:
    """A classified gene detection result.

    Classification values:
        - "full": confirmed identity (>=threshold) + high coverage (>=90%)
        - "partial": confirmed identity + moderate coverage (75-89.9%)
        - "novel_full": lower identity (below confirmed threshold) + high coverage
        - "novel_partial": lower identity + moderate coverage
    """

    gene_name: str
    pident: float
    coverage: float
    classification: str  # "full", "partial", "novel_full", "novel_partial"


class MecClassifier:
    """Classify mec gene hits from BLAST results.

    Thresholds:
        - Confirmed full: >=95% identity AND >=90% coverage
        - Confirmed partial: >=95% identity AND >=75% coverage
        - Novel full: 75-94.9% identity AND >=90% coverage
        - Novel partial: 75-94.9% identity AND >=75% coverage
    """

    CONFIRMED_PIDENT = 95.0
    NOVEL_PIDENT = 75.0
    FULL_COVERAGE = 0.90
    MIN_COVERAGE = 0.75

    def __init__(self, ref_fasta: str):
        self.ref_lengths = self._get_ref_lengths(ref_fasta)

    @staticmethod
    def _get_ref_lengths(fasta_path: str) -> Dict[str, int]:
        lengths = {}
        for record in SeqIO.parse(fasta_path, "fasta"):
            lengths[record.id] = len(record.seq)
        return lengths

    def classify(self, hits) -> List[GeneHit]:
        """Classify BLAST hits as confirmed/novel mec gene detections."""
        # Filter to minimum thresholds
        passing = filter_hits(
            hits,
            min_pident=self.NOVEL_PIDENT,
            min_coverage=self.MIN_COVERAGE,
            ref_lengths=self.ref_lengths,
        )

        # Resolve overlapping hits — keeps best hit per genomic location
        # This prevents cross-reactive allotype matches (e.g. mecA1/mecA2
        # appearing alongside mecA when only mecA is truly present)
        best_hits = get_best_non_overlapping_hits(passing, overlap_threshold=500)

        results = []
        for hit in best_hits:
            ref_len = self.ref_lengths.get(hit.qseqid, 1)
            coverage = hit.length / ref_len

            if hit.pident >= self.CONFIRMED_PIDENT:
                if coverage >= self.FULL_COVERAGE:
                    classification = "full"
                else:
                    classification = "partial"
            else:
                if coverage >= self.FULL_COVERAGE:
                    classification = "novel_full"
                else:
                    classification = "novel_partial"

            results.append(
                GeneHit(
                    gene_name=hit.qseqid,
                    pident=round(hit.pident, 1),
                    coverage=round(coverage * 100, 1),
                    classification=classification,
                )
            )

        return results


class CcrClassifier:
    """Classify ccr gene hits from BLAST results.

    Thresholds:
        - Confirmed full: >=85% identity AND >=90% coverage
        - Confirmed partial: >=85% identity AND >=75% coverage
        - Novel full: 70-84.4% identity AND >=90% coverage
        - Novel partial: 70-84.4% identity AND >=75% coverage
    """

    CONFIRMED_PIDENT = 84.5
    NOVEL_PIDENT = 70.0
    FULL_COVERAGE = 0.90
    MIN_COVERAGE = 0.75

    def __init__(self, ref_fasta: str):
        self.ref_lengths = self._get_ref_lengths(ref_fasta)

    @staticmethod
    def _get_ref_lengths(fasta_path: str) -> Dict[str, int]:
        lengths = {}
        for record in SeqIO.parse(fasta_path, "fasta"):
            lengths[record.id] = len(record.seq)
        return lengths

    def classify(self, hits) -> List[GeneHit]:
        """Classify BLAST hits as confirmed/novel ccr gene detections."""
        # Filter to minimum thresholds
        passing = filter_hits(
            hits,
            min_pident=self.NOVEL_PIDENT,
            min_coverage=self.MIN_COVERAGE,
            ref_lengths=self.ref_lengths,
        )

        # Resolve overlapping hits
        best_hits = get_best_non_overlapping_hits(passing, overlap_threshold=500)

        results = []
        for hit in best_hits:
            ref_len = self.ref_lengths.get(hit.qseqid, 1)
            coverage = hit.length / ref_len

            if hit.pident >= self.CONFIRMED_PIDENT:
                if coverage >= self.FULL_COVERAGE:
                    classification = "full"
                else:
                    classification = "partial"
            else:
                if coverage >= self.FULL_COVERAGE:
                    classification = "novel_full"
                else:
                    classification = "novel_partial"

            results.append(
                GeneHit(
                    gene_name=hit.qseqid,
                    pident=round(hit.pident, 1),
                    coverage=round(coverage * 100, 1),
                    classification=classification,
                )
            )

        return results


class CcrComplexLookup:
    """Map ccr allotype combinations to IWG-SCC ccr complex types (1-22).

    Types 1-9 follow original IWG-SCC designations.
    Type 10 from Xiao et al. 2023 (J Antimicrob Chemother 78:440-4).
    Types 11-22 from Huang et al. 2024 (J Infect Dis, doi:10.1093/infdis/jiae044).

    Standard pairings:
        Type 1:  ccrA1  + ccrB1
        Type 2:  ccrA2  + ccrB2
        Type 3:  ccrA3  + ccrB3
        Type 4:  ccrA4  + ccrB4
        Type 10: ccrA8  + ccrB9
        Type 13: ccrA10 + ccrB10

    Mixed pairings:
        Type 6:  ccrA5  + ccrB3
        Type 7:  ccrA1  + ccrB6
        Type 8:  ccrA1  + ccrB3
        Type 11: ccrA9  + ccrB3
        Type 12: ccrA10 + ccrB1
        Type 14: ccrA11 + ccrB7
        Type 15: ccrA11 + ccrB12
        Type 16: ccrA12 + ccrB1
        Type 17: ccrA12 + ccrB3
        Type 18: ccrA13 + ccrB3
        Type 19: ccrA14 + ccrB11

    Single genes:
        Type 5:  ccrC1
        Type 9:  ccrC2
        Type 20: ccrC3
        Type 21: ccrC4
        Type 22: ccrC5
    """

    # Ordered: pairs (size 2) first, standard pairings before mixed, then singles.
    # For shared ccrA allotypes, matched-number pairs precede mixed pairs
    # (e.g. type 13 ccrA10+ccrB10 before type 12 ccrA10+ccrB1).
    _COMPLEX_MAP = [
        # Established pairs (IWG-SCC types 1-4)
        (frozenset({"ccrA1", "ccrB1"}), "1"),
        (frozenset({"ccrA2", "ccrB2"}), "2"),
        (frozenset({"ccrA3", "ccrB3"}), "3"),
        (frozenset({"ccrA4", "ccrB4"}), "4"),
        # Established mixed pairs (IWG-SCC types 6-8)
        (frozenset({"ccrA5", "ccrB3"}), "6"),
        (frozenset({"ccrA1", "ccrB6"}), "7"),
        (frozenset({"ccrA1", "ccrB3"}), "8"),
        # Type 10 (Xiao et al. 2023)
        (frozenset({"ccrA8", "ccrB9"}), "10"),
        # Types 11-19 (Huang et al. 2024) — matched pairs before mixed
        (frozenset({"ccrA10", "ccrB10"}), "13"),
        (frozenset({"ccrA9", "ccrB3"}), "11"),
        (frozenset({"ccrA10", "ccrB1"}), "12"),
        (frozenset({"ccrA11", "ccrB7"}), "14"),
        (frozenset({"ccrA11", "ccrB12"}), "15"),
        (frozenset({"ccrA12", "ccrB1"}), "16"),
        (frozenset({"ccrA12", "ccrB3"}), "17"),
        (frozenset({"ccrA13", "ccrB3"}), "18"),
        (frozenset({"ccrA14", "ccrB11"}), "19"),
        # Single genes — established
        (frozenset({"ccrC1"}), "5"),
        (frozenset({"ccrC2"}), "9"),
        # Single genes — Huang et al. 2024
        (frozenset({"ccrC3"}), "20"),
        (frozenset({"ccrC4"}), "21"),
        (frozenset({"ccrC5"}), "22"),
    ]

    @classmethod
    def lookup(cls, ccr_results: List[GeneHit]) -> str:
        """Determine ccr complex type(s) from classified ccr gene hits.

        Uses greedy matching: pairs (size 2) are matched before singles,
        consuming matched genes from the remaining set. Composite elements
        with multiple ccr genes can yield multiple complex types (e.g. "2;5").

        Returns:
            Complex type string, e.g. "2", "2;5", "novel_combination", or "-".
            Appends " (novel_genes)" when any hit has a novel classification.
        """
        if not ccr_results:
            return "-"

        allotype_names = set(r.gene_name for r in ccr_results)
        has_novel = any(
            r.classification.startswith("novel") for r in ccr_results
        )

        remaining = set(allotype_names)
        matched_types = []

        for pattern, complex_type in cls._COMPLEX_MAP:
            if pattern.issubset(remaining):
                matched_types.append(complex_type)
                remaining -= pattern

        if not matched_types:
            label = "novel_combination"
        else:
            label = ";".join(matched_types)

        if has_novel:
            label += " (novel_genes)"

        return label


# Column header for typing output — shared with report_sccmec
TYPING_HEADER = [
    "Input_File",
    "mec_genes",
    "mec_identity",
    "mec_coverage",
    "ccr_genes",
    "ccr_allotypes",
    "ccr_identity",
    "ccr_complex_type",
]


class SCCmecTyper:
    """Orchestrate BLAST-based SCCmec typing for mec and ccr gene content."""

    def __init__(
        self,
        mec_ref: Optional[str] = None,
        ccr_ref: Optional[str] = None,
    ):
        self.mec_ref = mec_ref
        self.ccr_ref = ccr_ref
        self.runner = BlastRunner()

        # Cache classifiers — ref FASTAs parsed once, reused for all genomes
        self._mec_classifier = self._create_mec_classifier()
        self._ccr_classifier = self._create_ccr_classifier()

    def _create_mec_classifier(self) -> MecClassifier:
        """Create a MecClassifier, handling default ref resolution."""
        if self.mec_ref:
            return MecClassifier(self.mec_ref)
        with get_default_ref("mec_genes_allotypes.fasta") as ref:
            return MecClassifier(str(ref))

    def _create_ccr_classifier(self) -> CcrClassifier:
        """Create a CcrClassifier, handling default ref resolution."""
        if self.ccr_ref:
            return CcrClassifier(self.ccr_ref)
        with get_default_ref("ccr_genes.fasta") as ref:
            return CcrClassifier(str(ref))

    def type_file(self, input_fasta: str) -> dict:
        """Type a single SCCmec FASTA file.

        Returns a dict with keys:
            Input_File, mec_genes, mec_identity, mec_coverage,
            ccr_genes, ccr_allotypes, ccr_identity
        """
        input_name = Path(input_fasta).stem
        tmp_dir = tempfile.mkdtemp(prefix="sccmec_type_")
        db_prefix = os.path.join(tmp_dir, "sccmec_db")

        try:
            # Create BLAST DB from the input SCCmec sequence
            self.runner.create_db(input_fasta, db_prefix)

            # Type mec genes
            mec_hits = self._blast_ref("mec", db_prefix, tmp_dir)
            # Type ccr genes
            ccr_hits = self._blast_ref("ccr", db_prefix, tmp_dir)

        finally:
            self.runner.cleanup_db(db_prefix)
            try:
                os.rmdir(tmp_dir)
            except OSError:
                pass

        # Classify using cached classifiers (ref FASTAs parsed once in __init__)
        mec_results = self._mec_classifier.classify(mec_hits)
        ccr_results = self._ccr_classifier.classify(ccr_hits)

        # Format output
        return self._format_result(input_name, mec_results, ccr_results)

    def _blast_ref(self, ref_type: str, db_prefix: str, tmp_dir: str):
        """BLAST a reference set against the SCCmec database."""
        if ref_type == "mec":
            if self.mec_ref:
                ref_path = self.mec_ref
                results_file = self.runner.run_blastn(ref_path, db_prefix)
                hits = parse_blast_output(results_file)
                self.runner.cleanup_file(results_file)
                return hits
            else:
                with get_default_ref("mec_genes_allotypes.fasta") as ref:
                    results_file = self.runner.run_blastn(str(ref), db_prefix)
                    hits = parse_blast_output(results_file)
                    self.runner.cleanup_file(results_file)
                    return hits
        else:
            if self.ccr_ref:
                ref_path = self.ccr_ref
                results_file = self.runner.run_blastn(ref_path, db_prefix)
                hits = parse_blast_output(results_file)
                self.runner.cleanup_file(results_file)
                return hits
            else:
                with get_default_ref("ccr_genes.fasta") as ref:
                    results_file = self.runner.run_blastn(str(ref), db_prefix)
                    hits = parse_blast_output(results_file)
                    self.runner.cleanup_file(results_file)
                    return hits

    @staticmethod
    def _format_result(
        input_name: str,
        mec_results: List[GeneHit],
        ccr_results: List[GeneHit],
    ) -> dict:
        """Format typing results into a dict for TSV output."""
        if mec_results:
            mec_genes = ";".join(
                f"{r.gene_name}({r.classification})" for r in mec_results
            )
            mec_identity = ";".join(str(r.pident) for r in mec_results)
            mec_coverage = ";".join(str(r.coverage) for r in mec_results)
        else:
            mec_genes = "-"
            mec_identity = "-"
            mec_coverage = "-"

        if ccr_results:
            ccr_genes = ";".join(
                f"{r.gene_name}({r.classification})" for r in ccr_results
            )
            # Sort by gene name so allotypes and identity columns align
            sorted_ccr = sorted(ccr_results, key=lambda r: r.gene_name)
            ccr_allotypes = ";".join(r.gene_name for r in sorted_ccr)
            ccr_identity = ";".join(str(r.pident) for r in sorted_ccr)
        else:
            ccr_genes = "-"
            ccr_allotypes = "-"
            ccr_identity = "-"

        ccr_complex_type = CcrComplexLookup.lookup(ccr_results)

        return {
            "Input_File": input_name,
            "mec_genes": mec_genes,
            "mec_identity": mec_identity,
            "mec_coverage": mec_coverage,
            "ccr_genes": ccr_genes,
            "ccr_allotypes": ccr_allotypes,
            "ccr_identity": ccr_identity,
            "ccr_complex_type": ccr_complex_type,
        }


def collect_input_files(paths: List[str]) -> List[str]:
    """Collect FASTA files from file paths and/or directories.

    Args:
        paths: List of file paths or directories.

    Returns:
        List of resolved FASTA file paths.
    """
    fasta_extensions = {".fasta", ".fna", ".fa"}
    files = []

    for p in paths:
        path = Path(p)
        if path.is_dir():
            for ext in fasta_extensions:
                files.extend(str(f) for f in path.glob(f"*{ext}"))
        elif path.is_file():
            files.append(str(path))
        else:
            print(f"WARNING: Skipping {p} (not a file or directory)")

    return sorted(files)


def main():
    parser = argparse.ArgumentParser(
        description="Type extracted SCCmec sequences by mec and ccr gene content"
    )
    parser.add_argument(
        "-f",
        "--fasta",
        nargs="+",
        required=True,
        help="Input SCCmec FASTA file(s) or directory",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        required=True,
        help="Output TSV file for typing results",
    )
    parser.add_argument(
        "--mec-ref",
        help="Custom mec gene reference FASTA (default: bundled)",
    )
    parser.add_argument(
        "--ccr-ref",
        help="Custom ccr gene reference FASTA (default: bundled)",
    )
    args = parser.parse_args()

    # Collect input files
    input_files = collect_input_files(args.fasta)

    if not input_files:
        print("ERROR: No FASTA files found in the provided paths")
        return

    print(f"Found {len(input_files)} input file(s)")

    # Create typer
    typer = SCCmecTyper(mec_ref=args.mec_ref, ccr_ref=args.ccr_ref)

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
