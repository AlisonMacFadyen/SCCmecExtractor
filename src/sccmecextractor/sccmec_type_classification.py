#!/usr/bin/env python

"""SCCmec Type Classification by mec complex content and ccr gene presence using BLAST.

Assigns mec complex class (A-E), ccr complex type (1-22) and SCCmec type (I-XV)
via BLAST against bundled reference databases.

When custom --mec-ref or --ccr-ref references are provided, the tool operates
in gene content detection mode for that side, reporting detected genes without
class/complex/type assignment.  The other side still uses the bundled reference
for full typing.
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
    strand: str # "+" or "-"
    contig: str = "-"
    start: int = 0
    end: int = 0


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

            if hit.sstart < hit.send:
                orientation = "+"
            elif hit.sstart > hit.send:
                orientation = "-"
            else:
                orientation = "unknown"

            results.append(
                GeneHit(
                    gene_name=hit.qseqid,
                    pident=round(hit.pident, 1),
                    coverage=round(coverage * 100, 1),
                    classification=classification,
                    contig=hit.sseqid,
                    start=min(hit.sstart, hit.send),
                    end=max(hit.sstart, hit.send),
                    strand=orientation
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

            if hit.sstart < hit.send:
                orientation = "+"
            elif hit.sstart > hit.send:
                orientation = "-"
            else:
                orientation = "unknown"

            results.append(
                GeneHit(
                    gene_name=hit.qseqid,
                    pident=round(hit.pident, 1),
                    coverage=round(coverage * 100, 1),
                    classification=classification,
                    contig=hit.sseqid,
                    start=min(hit.sstart, hit.send),
                    end=max(hit.sstart, hit.send),
                    strand=orientation,
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
    "mec_locations",
    "mec_class_type",
    "ccr_genes",
    "ccr_allotypes",
    "ccr_identity",
    "ccr_locations",
    "ccr_complex_type",
    "SCCmec_Type",
    "SCCmec_Type_secondary"
]

class MecComplexLookup:
    """Map mec complex genes to IWG-SCC mec complex Classes
    
    Gene content for each class outlined below:

    Class A: IS431, mecA, mecR1, mecI
    Class B: IS431, mecA, mecR1, IS1272
    Class C: IS431, mecA, mecR1, IS431
    Class D: IS431, mecA, mecR1
    Class E: blaZ, mecC, mecR1, mecI
    
    Note Class C has two subtypes that are relevant to SCCmec Typing, C1 and C2.
    The distinguishing feature is the orientation of IS431 at each end, they are either
    the same or different:

    Class C1: <- IS431.........<- IS431
    Class C2: <- IS431.........-> IS431

    """

    _MEC_COMPLEX_RULES = [
    # Most specific first: Class E (mecC-based, distinct from all mecA classes)
    # Opted not to include "all" genes e.g. mecR1, as these can be truncated or absent
    {
        "class": "E",
        "requires": {"mecC", "mecI"},
        "absent": set(),
    },
    # Class A: mecA with intact regulation (mecI present)
    {
        "class": "A",
        "requires": {"mecA", "mecI"},
        "absent": set(),
    },
    # Class B: mecA with IS1272 (no mecI)
    {
        "class": "B",
        "requires": {"mecA", "IS1272"},
        "absent": {"mecI"},
    },
    # Class C: mecA flanked by IS431 (no mecI, no IS1272)
    # Class D: mecA with truncated mecR1, IS431 upstream, no mecI, no IS1272
    # NOTE: C1 vs C2 distinguished by IS431 orientation (as described above)
    # NOTE: C vs D can be distinguished by IS431 copy number (2 vs 1, respectively)
    {
        "class": "C_or_D",
        "requires": {"mecA", "IS431"},
        "absent": {"mecI", "IS1272"},
    },
    ]

    # Proximity threshold for IS elements relative to mecA (bp).
    # Genuine mec complex IS elements are within ~3.5 kb (canonical) or
    # ~7.5 kb (variant with extra IS431 insertion).  10 kb provides margin.
    IS_PROXIMITY_THRESHOLD = 10_000

    @staticmethod
    def _find_primary_meca(mec_results):
        """Find the primary mecA hit (highest coverage, then identity).

        Returns the best mecA GeneHit, or None if no mecA detected.
        """
        mec_hits = [r for r in mec_results if r.gene_name.startswith("mecA")
                     and not r.gene_name.startswith("mecA1")
                     and not r.gene_name.startswith("mecA2")]
        if not mec_hits:
            return None
        return max(mec_hits, key=lambda r: (r.coverage, r.pident))

    @classmethod
    def _is_near_meca(cls, mec_results, gene_name, threshold=None):
        """Check whether any hit for *gene_name* is within *threshold* bp of mecA.

        Only considers hits on the same contig as mecA.  The distance is
        measured as the gap between the two genes (not centre-to-centre).

        Returns True if at least one hit is within threshold.
        """
        if threshold is None:
            threshold = cls.IS_PROXIMITY_THRESHOLD

        meca = cls._find_primary_meca(mec_results)
        if meca is None:
            return False

        for r in mec_results:
            if r.gene_name != gene_name or r.contig != meca.contig:
                continue
            # Gap between the two genes
            if r.end < meca.start:
                dist = meca.start - r.end
            elif r.start > meca.end:
                dist = r.start - meca.end
            else:
                dist = 0  # overlapping
            if dist <= threshold:
                return True
        return False

    @classmethod
    def _get_flanking_is431(cls, mec_results, upstream_dist=5000, downstream_dist=10000):
        """Return IS431 hits that flank mecA (one upstream, one downstream).

        Only considers IS431 on the same contig as mecA and within the
        specified upstream/downstream distance thresholds.

        Returns a list of 0 or 2 GeneHit objects (closest on each side).
        """
        meca = cls._find_primary_meca(mec_results)
        if meca is None:
            return []

        upstream = []
        downstream = []
        for r in mec_results:
            if r.gene_name == "IS431" and r.contig == meca.contig:
                if (meca.start - r.end) > 0 and (meca.start - r.end) <= upstream_dist:
                    upstream.append(r)
                elif (r.start - meca.end) > 0 and (r.start - meca.end) <= downstream_dist:
                    downstream.append(r)

        if upstream and downstream:
            closest_upstream = max(upstream, key=lambda r: r.end)
            closest_downstream = min(downstream, key=lambda r: r.start)
            return [closest_upstream, closest_downstream]

        return []

    @staticmethod
    def subtype_class_c(flanking) -> str:
        """Distinguish mec complex C1 from C2 by IS431 orientation.

        C1: both IS431 copies in SAME orientation (← mecA ←)
        C2: IS431 copies in OPPOSITE orientation (← mecA →)

        Returns "C1" or "C2".
        """
        orientations = [hit.strand for hit in flanking]
        if len(set(orientations)) == 1:
            return "C1"
        else:
            return "C2"

    @classmethod
    def lookup(cls, mec_results: List[GeneHit]) -> str:
        """Classify mec complex from detected gene content near mecA.

        IS elements (IS1272, IS431) are only counted when they are within
        ``IS_PROXIMITY_THRESHOLD`` bp of mecA on the same contig.  This
        prevents distant transposase copies from causing false
        classification.

        Returns
        -------
        str
            Mec complex class (A, B, C1, C2, D, E) or ``"not_typeable"``
            or ``"-"`` (no mec gene detected).
        """
        mec_genes_list = [
            "mecA",
            "mecA1",
            "mecA2",
            "mecB",
            "mecC",
            "mecC1",
            "mecC2",
            "mecC3",
            "mecD",
        ]

        if not any(r.gene_name in mec_genes_list for r in mec_results):
            return "-"

        # Build detected-gene set with proximity validation for IS elements.
        # Resistance genes and regulatory genes are included genome-wide;
        # IS elements require proximity to mecA.
        detected = set()
        for r in mec_results:
            base_name = r.gene_name.split("_")[0]
            if base_name in ("IS1272", "IS431"):
                continue  # handled below via proximity check
            detected.add(base_name)

        if cls._is_near_meca(mec_results, "IS1272"):
            detected.add("IS1272")
        if cls._is_near_meca(mec_results, "IS431"):
            detected.add("IS431")

        for rule in cls._MEC_COMPLEX_RULES:
            has_required = rule["requires"].issubset(detected)
            lacks_excluded = rule["absent"].isdisjoint(detected)
            if has_required and lacks_excluded:
                if rule["class"] == "C_or_D":
                    flanking = cls._get_flanking_is431(mec_results)
                    if len(flanking) >= 2:
                        return cls.subtype_class_c(flanking)
                    else:
                        return "D"
                else:
                    return rule["class"]
        return "not_typeable"


class SCCmecTypeLookup:
    """Combine ccr type and mec complex information to assign SCCmec Type"""

    _SCCMEC_TYPE_MAP = {
        ("B", "1"):  "I",
        ("A", "2"):  "II",
        ("A", "3"):  "III",
        ("B", "2"):  "IV",
        ("C2", "5"): "V",
        ("B", "4"):  "VI",
        ("C1", "5"): "VII",
        ("A", "4"):  "VIII",
        ("C2", "1"): "IX",
        ("C1", "7"): "X",
        ("E", "8"):  "XI",
        ("C2", "9"): "XII",
        ("A", "9"):  "XIII",
        ("A", "5"):  "XIV",
        ("A", "7"):  "XV",
    }

    @classmethod
    def lookup(cls, mec_class_type, ccr_complex_type) -> str:
        """Assign SCCmec Type based on mec complex and ccr type"""

        if mec_class_type == "-" or ccr_complex_type == "-":
            return "not_typeable"

        # Direct match first before considering composites
        direct = cls._SCCMEC_TYPE_MAP.get((mec_class_type, ccr_complex_type))
        if direct:
            return direct
        
        # Composite check, try each ccr complex separately
        if ";" in ccr_complex_type:
            types = []
            for ccr_part in ccr_complex_type.split(";"):
                match = cls._SCCMEC_TYPE_MAP.get((mec_class_type, ccr_part))
                if match:
                    types.append(match)
            if types:
                return ";".join(types)
        
        return "novel_combination"

class SCCmecTyper:
    """Orchestrate BLAST-based SCCmec typing for mec and ccr gene content.

    Operating modes:

    **Typing mode** (default, no custom refs):
        Uses bundled mec_class_reference.fasta and ccr_genes.fasta to assign
        mec complex class (A-E), ccr complex type (1-22) and SCCmec type (I-XV).

    **Gene content mode** (custom --mec-ref and/or --ccr-ref):
        When a custom reference is provided for one side, that side reports
        gene content only (no class/complex assignment).  The other side
        still uses the bundled reference for full typing.  When both custom
        refs are provided, both sides report gene content only.
    """

    def __init__(
        self,
        mec_ref: Optional[str] = None,
        ccr_ref: Optional[str] = None,
    ):
        self.mec_ref = mec_ref
        self.ccr_ref = ccr_ref
        self.runner = BlastRunner()

        # Track which sides use custom refs (gene content only)
        self._custom_mec = mec_ref is not None
        self._custom_ccr = ccr_ref is not None

        # Cache classifiers — ref FASTAs parsed once, reused for all genomes
        self._mec_classifier = self._create_mec_classifier()
        self._ccr_classifier = self._create_ccr_classifier()

    def _create_mec_classifier(self) -> MecClassifier:
        """Create a MecClassifier, handling default ref resolution."""
        if self._custom_mec:
            return MecClassifier(self.mec_ref)
        with get_default_ref("mec_class_reference.fasta") as ref:
            return MecClassifier(str(ref))

    def _create_ccr_classifier(self) -> CcrClassifier:
        """Create a CcrClassifier, handling default ref resolution."""
        if self._custom_ccr:
            return CcrClassifier(self.ccr_ref)
        with get_default_ref("ccr_genes.fasta") as ref:
            return CcrClassifier(str(ref))

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

        if db_prefix is not None:
            # Reuse caller-provided BLAST DB
            tmp_dir = tempfile.mkdtemp(prefix="sccmec_type_")
            try:
                mec_hits = self._blast_ref("mec", db_prefix, tmp_dir)
                ccr_hits = self._blast_ref("ccr", db_prefix, tmp_dir)
            finally:
                try:
                    os.rmdir(tmp_dir)
                except OSError:
                    pass
        else:
            # Standalone mode: create a temporary BLAST DB
            tmp_dir = tempfile.mkdtemp(prefix="sccmec_type_")
            db_prefix_local = os.path.join(tmp_dir, "sccmec_db")
            try:
                self.runner.create_db(input_fasta, db_prefix_local)
                mec_hits = self._blast_ref("mec", db_prefix_local, tmp_dir)
                ccr_hits = self._blast_ref("ccr", db_prefix_local, tmp_dir)
            finally:
                self.runner.cleanup_db(db_prefix_local)
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
            ref_path = self.mec_ref
            default_name = "mec_class_reference.fasta"
        else:
            ref_path = self.ccr_ref
            default_name = "ccr_genes.fasta"

        if ref_path:
            results_file = self.runner.run_blastn(ref_path, db_prefix)
            hits = parse_blast_output(results_file)
            self.runner.cleanup_file(results_file)
            return hits
        else:
            with get_default_ref(default_name) as ref:
                results_file = self.runner.run_blastn(str(ref), db_prefix)
                hits = parse_blast_output(results_file)
                self.runner.cleanup_file(results_file)
                return hits
    
    @staticmethod
    def find_closest_ccr(mec_results, ccr_results):
        """Find the ccr gene(s) closest to mecA."""
        # Find mecA location
        mec_hits = [r for r in mec_results if r.gene_name.startswith("mecA")]
        if not mec_hits or not ccr_results:
            return ccr_results  # fallback: return all

        mec_hit = max(mec_hits, key=lambda r: (r.coverage, r.pident))

        # Only consider ccr on the same contig
        same_contig = [r for r in ccr_results if r.contig == mec_hit.contig]
        if not same_contig:
            return ccr_results  # fallback

        # Find the single closest ccr gene to mecA
        closest = min(same_contig, key=lambda r: min(
            abs(r.start - mec_hit.end),
            abs(mec_hit.start - r.end)
        ))
        
        # Include its pair partner if present (e.g. closest=ccrA2 → also grab ccrB2)
        # ccrC genes are unpaired so just return them alone
        if closest.gene_name.startswith("ccrC"):
            return [closest]
            
        # For ccrA/ccrB, grab the matching partner nearby
        primary = [closest]
        for r in same_contig:
            if r is not closest:
            # e.g. closest is ccrA2 → look for ccrB on same contig nearby
                if (closest.gene_name.startswith("ccrA") and r.gene_name.startswith("ccrB")) or \
                    (closest.gene_name.startswith("ccrB") and r.gene_name.startswith("ccrA")):
                    if abs(r.start - closest.end) < 5000 or abs(closest.start - r.end) < 5000:
                        primary.append(r)
                        break

        return primary

    def _format_result(
        self,
        input_name: str,
        mec_results: List[GeneHit],
        ccr_results: List[GeneHit],
    ) -> dict:
        """Format typing results into a dict for TSV output.

        When a custom reference is used for one side, that side reports
        gene content only (class/complex/type columns show '-').
        """
        if mec_results:
            mec_genes = ";".join(
                f"{r.gene_name}({r.classification})" for r in mec_results
            )
            mec_identity = ";".join(str(r.pident) for r in mec_results)
            mec_coverage = ";".join(str(r.coverage) for r in mec_results)
            mec_locations = ";".join(
                f"{r.contig}:{r.start}-{r.end}" for r in mec_results
            )
        else:
            mec_genes = "-"
            mec_identity = "-"
            mec_coverage = "-"
            mec_locations = "-"

        if ccr_results:
            ccr_genes = ";".join(
                f"{r.gene_name}({r.classification})" for r in ccr_results
            )
            # Sort by gene name so allotypes and identity columns align
            sorted_ccr = sorted(ccr_results, key=lambda r: r.gene_name)
            ccr_allotypes = ";".join(r.gene_name for r in sorted_ccr)
            ccr_identity = ";".join(str(r.pident) for r in sorted_ccr)
            ccr_locations = ";".join(
                f"{r.contig}:{r.start}-{r.end}" for r in sorted_ccr
            )
        else:
            ccr_genes = "-"
            ccr_allotypes = "-"
            ccr_identity = "-"
            ccr_locations = "-"

        # Typing lookups — only when using bundled references
        if self._custom_mec:
            mec_class_type = "-"
        else:
            mec_class_type = MecComplexLookup.lookup(mec_results)

        if self._custom_ccr:
            ccr_complex_type = "-"
        else:
            ccr_complex_type = CcrComplexLookup.lookup(ccr_results)
            closest_ccr = SCCmecTyper.find_closest_ccr(mec_results, ccr_results)
            primary_ccr_type = CcrComplexLookup.lookup(closest_ccr)

        # SCCmec type requires both mec class and ccr complex from bundled refs
        if self._custom_mec or self._custom_ccr:
            sccmec_type = "-"
            sccmec_type_secondary = "-"
        else:
            sccmec_type = SCCmecTypeLookup.lookup(mec_class_type, primary_ccr_type)

            # Secondary type(s) from non-primary ccr complexes
            if ";" in ccr_complex_type:
                primary_ids = {id(r) for r in closest_ccr}
                secondary_ccr = [r for r in ccr_results if id(r) not in primary_ids]
                secondary_ccr_type = CcrComplexLookup.lookup(secondary_ccr)
                secondary_types = []
                if ";" in secondary_ccr_type:
                    for part in secondary_ccr_type.split(";"):
                        secondary_types.append(
                            SCCmecTypeLookup.lookup(mec_class_type, part)
                        )
                else:
                    secondary_types.append(
                        SCCmecTypeLookup.lookup(mec_class_type, secondary_ccr_type)
                    )
                sccmec_type_secondary = ";".join(secondary_types)
            else:
                sccmec_type_secondary = "-"

        return {
            "Input_File": input_name,
            "mec_genes": mec_genes,
            "mec_identity": mec_identity,
            "mec_coverage": mec_coverage,
            "mec_locations": mec_locations,
            "mec_class_type": mec_class_type,
            "ccr_genes": ccr_genes,
            "ccr_allotypes": ccr_allotypes,
            "ccr_identity": ccr_identity,
            "ccr_locations": ccr_locations,
            "ccr_complex_type": ccr_complex_type,
            "SCCmec_Type": sccmec_type,
            "SCCmec_Type_secondary": sccmec_type_secondary,
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
        description="Type SCC elements by mec complex class (A-E), ccr complex type (1-22) and SCCmec type (I-XV)"
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
        help="Custom mec gene reference FASTA for gene content detection. "
             "When provided, mec genes are detected using this reference "
             "instead of the bundled database; mec complex class and SCCmec "
             "type will not be assigned. ccr typing still uses the bundled "
             "reference.",
    )
    parser.add_argument(
        "--ccr-ref",
        help="Custom ccr gene reference FASTA for gene content detection. "
             "When provided, ccr genes are detected using this reference "
             "instead of the bundled database; ccr complex type and SCCmec "
             "type will not be assigned. mec typing still uses the bundled "
             "reference.",
    )
    args = parser.parse_args()

    # Collect input files
    input_files = collect_input_files(args.fasta)

    if not input_files:
        print("ERROR: No FASTA files found in the provided paths")
        return

    print(f"Found {len(input_files)} input file(s)")

    # Report operating mode
    if args.mec_ref and args.ccr_ref:
        print("Mode: gene content detection (custom mec and ccr references)")
    elif args.mec_ref:
        print("Mode: custom mec gene detection + bundled ccr typing")
    elif args.ccr_ref:
        print("Mode: bundled mec typing + custom ccr gene detection")
    else:
        print("Mode: full SCCmec typing (mec class, ccr complex, SCCmec type)")

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