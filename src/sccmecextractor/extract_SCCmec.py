#!/usr/bin/env python

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

class InputValidator:
    """Check input files are valid"""
    
    def validate_fasta_file(self, input_fasta):
        """Validate that FASTA file exists and is readable."""
        fasta = Path(input_fasta)
            
        if not fasta.exists():
            raise FileNotFoundError(f"FASTA file not found: {input_fasta}")
    
        if not fasta.is_file():
            raise ValueError(f"Path to FASTA exists but is not a file: {input_fasta}")
            
        if fasta.stat().st_size == 0:
            raise ValueError(f"GFF file is empty: {input_fasta}")
        
    def validate_gff_file(self, input_gff):
        """Validate the GFF file exists and is readable"""
        gff = Path(input_gff)
        
        if not gff.exists():
            raise FileNotFoundError(f"GFF file not found: {input_gff}")
        
        if not gff.is_file():
            raise ValueError(f"Path to GFF file exists but is not a file: {input_gff}")
            
        if gff.stat().st_size == 0:
            raise ValueError(f"GFF file is empty: {input_gff}")
        
    def validate_tsv_file(self, input_tsv):
        """Validate the input TSV file that contains the att site locations"""
        tsv = Path(input_tsv)

        if not tsv.exists():
            raise FileNotFoundError(f"TSV file containing att sites not found: {input_tsv}")
        
        if not tsv.is_file():
            raise ValueError(f"Path exists for TSV but not a file: {input_tsv}")
        
        if tsv.stat().st_size == 0:
            raise ValueError(f"TSV file is empty: {input_tsv}")

class AttSite:
    """Represents a single att site with its properties."""
    
    def __init__(self, pattern: str, contig: str, start: int, end: int):
        self.pattern = pattern
        self.contig = contig
        self.start = start
        self.end = end
        self.is_right = pattern.lower().startswith(('attr', 'cattr'))
        self.is_left = pattern.lower().startswith(('attl', 'cattl'))
    
    def distance_to(self, other_site: 'AttSite') -> int:
        """Calculate distance to another att site on the same contig."""
        if self.contig != other_site.contig:
            return float('inf')
        return abs(self.start - other_site.end)
    
    def __str__(self):
        return f"{self.pattern}({self.start}-{self.end})"


@dataclass
class ExtractionReport:
    """Per-genome extraction metadata for reporting."""
    input_file: str = ""
    status: str = ""  # "extracted", "composite_extracted", "skipped", "failed"
    contig: str = "-"
    attr_pattern: str = "-"
    attr_start: str = "-"
    attr_end: str = "-"
    attl_pattern: str = "-"
    attl_start: str = "-"
    attl_end: str = "-"
    element_size: str = "-"
    is_composite: str = "False"
    outer_attl_pattern: str = "-"
    outer_attl_start: str = "-"
    outer_attl_end: str = "-"
    composite_size: str = "-"
    contig_edge_flags: str = "-"
    partial_type: str = "-"
    notes: str = "-"

    HEADER = (
        "Input_File\tStatus\tContig\tAttR_Pattern\tAttR_Start\tAttR_End\t"
        "AttL_Pattern\tAttL_Start\tAttL_End\tElement_Size_bp\tIs_Composite\t"
        "Outer_AttL_Pattern\tOuter_AttL_Start\tOuter_AttL_End\tComposite_Size_bp\t"
        "Contig_Edge_Flags\tPartial_Type\tNotes"
    )

    def to_tsv_row(self) -> str:
        return "\t".join([
            self.input_file, self.status, self.contig,
            self.attr_pattern, str(self.attr_start), str(self.attr_end),
            self.attl_pattern, str(self.attl_start), str(self.attl_end),
            str(self.element_size), str(self.is_composite),
            self.outer_attl_pattern, str(self.outer_attl_start), str(self.outer_attl_end),
            str(self.composite_size),
            self.contig_edge_flags, self.partial_type, self.notes,
        ])


@dataclass
class AmbiguousHitReport:
    """Per-genome report for actionable extraction failures."""
    input_file: str = ""
    reason: str = "-"
    right_att_sites: str = "-"
    left_att_sites: str = "-"
    rlmh_contigs: str = "-"
    element_size: str = "-"
    notes: str = "-"

    HEADER = (
        "Input_File\tReason\tRight_Att_Sites\tLeft_Att_Sites\t"
        "rlmH_Contigs\tElement_Size_bp\tNotes"
    )

    def to_tsv_row(self) -> str:
        return "\t".join([
            self.input_file, self.reason, self.right_att_sites,
            self.left_att_sites, self.rlmh_contigs,
            self.element_size, self.notes,
        ])


def _count_distinct_loci(coords):
    """Count distinct genomic loci from a collection of overlapping coordinates.

    Merges overlapping (start, end) intervals and returns the number of
    non-overlapping clusters.  This prevents multiple BLAST hits against
    the same gene (from different reference sequences) being counted as
    separate copies.
    """
    if not coords:
        return 0
    sorted_coords = sorted(coords)
    loci = 1
    _, current_end = sorted_coords[0]
    for start, end in sorted_coords[1:]:
        if start > current_end:
            loci += 1
            current_end = end
        else:
            current_end = max(current_end, end)
    return loci


class PrecomputedRlmH:
    """Wraps pre-computed rlmH positions to match GeneAnnotations interface.

    Used by the pipeline to pass Stage 1 rlmH results to Stage 2,
    eliminating a redundant BLAST search.
    """

    def __init__(self, rlmh_genes):
        self.rlmH_positions = {}
        self.multi_rlmH_contigs = {}
        for contig, coords in rlmh_genes.items():
            self.rlmH_positions[contig] = min(start for start, end in coords)
            n_loci = _count_distinct_loci(coords)
            if n_loci > 1:
                self.multi_rlmH_contigs[contig] = n_loci

    def has_rlmH(self, contig):
        return contig in self.rlmH_positions

    def get_rlmH_start(self, contig):
        return self.rlmH_positions.get(contig)

    def get_rlmH_copy_warning(self, contig):
        count = self.multi_rlmH_contigs.get(contig, 0)
        if count > 1:
            return f"Multiple rlmH copies detected ({count}) on {contig}"
        return ""


class RlmHBlastAdapter:
    """Adapts RlmHBlastDetector to the GeneAnnotations interface.

    Provides has_rlmH() and get_rlmH_start() using BLAST-based
    rlmH detection instead of GFF3 parsing.
    """

    def __init__(self, fasta_file: str, rlmh_ref: str = None):
        from sccmecextractor.locate_att_sites import RlmHBlastDetector
        self._detector = RlmHBlastDetector(fasta_file, rlmh_ref)
        self.rlmH_positions = self._build_positions()
        self.multi_rlmH_contigs = {
            contig: _count_distinct_loci(coords)
            for contig, coords in self._detector.rlmH_genes.items()
            if _count_distinct_loci(coords) > 1
        }

    def _build_positions(self) -> Dict[str, int]:
        """Convert BLAST detector format to simple contig -> start mapping."""
        positions = {}
        for contig, coords in self._detector.rlmH_genes.items():
            positions[contig] = min(start for start, end in coords)
        return positions

    def has_rlmH(self, contig: str) -> bool:
        """Check if a contig has an rlmH gene."""
        return contig in self.rlmH_positions

    def get_rlmH_start(self, contig: str) -> Optional[int]:
        """Get the start position of rlmH gene on a specific contig."""
        return self.rlmH_positions.get(contig)

    def get_rlmH_copy_warning(self, contig: str) -> str:
        count = self.multi_rlmH_contigs.get(contig, 0)
        if count > 1:
            return f"Multiple rlmH copies detected ({count}) on {contig}"
        return ""


class GeneAnnotations:
    """Handles parsing and storing gene annotation information."""

    RLMH_PRODUCT = 'ribosomal rna large subunit methyltransferase h'

    def __init__(self, gff3_file: str):
        self.gff3_file = gff3_file
        self.rlmH_positions = self._parse_rlmH_genes()

    @staticmethod
    def _is_rlmH_feature(attributes: dict) -> bool:
        """Check if GFF3 feature attributes identify rlmH.

        Handles both older Bakta (gene=rlmH) and newer Bakta where
        the gene name is dropped but the product name is retained.
        """
        if attributes.get('gene', '').lower() == 'rlmh':
            return True

        for attr in ('Name', 'product'):
            value = attributes.get(attr, '').lower()
            if 'rlmh' in value or GeneAnnotations.RLMH_PRODUCT in value:
                return True

        return False

    def _parse_rlmH_genes(self) -> Dict[str, int]:
        """Parse GFF3 file to extract rlmH gene start positions.

        Checks both gene and CDS feature types to handle annotation
        tools that may not assign a gene name to rlmH.
        """
        rlmH_info = {}
        rlmH_loci = {}  # contig -> set of unique start positions

        with open(self.gff3_file, 'r') as gff3:
            for line in gff3:
                if line.startswith("#"):
                    continue

                columns = line.strip().split("\t")
                if len(columns) != 9 or columns[2] not in ('gene', 'CDS'):
                    continue

                attributes = dict(
                    item.split('=', 1) for item in columns[8].split(';') if '=' in item
                )

                if self._is_rlmH_feature(attributes):
                    contig = columns[0]
                    start_position = int(columns[3])
                    rlmH_loci.setdefault(contig, set()).add(start_position)
                    # Prefer gene feature; don't overwrite with CDS
                    if contig not in rlmH_info:
                        rlmH_info[contig] = start_position

        if not rlmH_info:
            print(f"Warning: No rlmH genes found in {self.gff3_file}", file=sys.stderr)

        self.multi_rlmH_contigs = {
            contig: len(positions)
            for contig, positions in rlmH_loci.items()
            if len(positions) > 1
        }

        return rlmH_info
    
    def get_rlmH_start(self, contig: str) -> Optional[int]:
        """Get the start position of rlmH gene on a specific contig."""
        return self.rlmH_positions.get(contig)
    
    def has_rlmH(self, contig: str) -> bool:
        """Check if a contig has an rlmH gene."""
        return contig in self.rlmH_positions

    def get_rlmH_copy_warning(self, contig: str) -> str:
        count = self.multi_rlmH_contigs.get(contig, 0)
        if count > 1:
            return f"Multiple rlmH copies detected ({count}) on {contig}"
        return ""


class AttSiteCollection:
    """Manages a collection of att sites and provides analysis methods."""
    
    def __init__(self, tsv_file: str, target_file: str):
        self.tsv_file = tsv_file
        self.target_file = target_file
        self.sites = self._parse_att_sites()
    
    def _parse_att_sites(self) -> List[AttSite]:
        """Parse TSV file to extract att sites for the target file."""
        sites = []
        found_entries = False
        
        with open(self.tsv_file, 'r') as tsv:
            next(tsv)  # Skip header
            for line in tsv:
                columns = line.strip().split("\t")
                input_file = columns[0]

                # Skip duplicate header lines (from batch runs)
                if input_file == "Input_File":
                    continue
                
                # Only process entries for our target file
                if input_file == self.target_file:
                    found_entries = True
                    pattern = columns[1]
                    contig = columns[2]
                    start_pos = int(columns[3])
                    end_pos = int(columns[4])
                    
                    sites.append(AttSite(pattern, contig, start_pos, end_pos))
        
        if not found_entries:
            print(f"Warning: No entries found for {self.target_file} in {self.tsv_file}", file=sys.stderr)
        
        return sites
    
    def get_right_sites(self) -> List[AttSite]:
        """Get all right att sites (attR, cattR)."""
        return [site for site in self.sites if site.is_right]
    
    def get_left_sites(self) -> List[AttSite]:
        """Get all left att sites (attL, cattL)."""
        return [site for site in self.sites if site.is_left]
    
    def find_closest_pair(self, min_distance: int = 0) -> Optional[Tuple[AttSite, AttSite]]:
        """Find the closest attR-attL pair on the same contig.

        Parameters
        ----------
        min_distance : int
            Minimum element size (distance) between att sites. Pairs closer
            than this are skipped (e.g. overlapping patterns at attB).
        """
        right_sites = self.get_right_sites()
        left_sites = self.get_left_sites()

        if not right_sites or not left_sites:
            return None

        best_pair = None
        best_dist = float('inf')

        for right in right_sites:
            for left in left_sites:
                if right.contig == left.contig:
                    distance = right.distance_to(left)
                    if distance < min_distance:
                        continue
                    if distance < best_dist:
                        best_dist = distance
                        best_pair = (right, left)

        return best_pair
    
    def has_valid_sites(self) -> bool:
        """Check if we have both right and left sites."""
        return bool(self.get_right_sites() and self.get_left_sites())

    def get_sites_on_contig(self, contig: str) -> List[AttSite]:
        """Return all att sites on a specific contig."""
        return [site for site in self.sites if site.contig == contig]

    def find_outer_left_site(
        self, closest_pair: Tuple['AttSite', 'AttSite'], rlmH_start: int
    ) -> Optional['AttSite']:
        """Find the outermost left att site beyond the inner attL.

        For composite/nested elements there may be multiple attL sites.
        The inner attL is the one in the closest pair; the outer attL is
        the one furthest from rlmH on the same contig.

        Returns the outer AttSite, or None if the element is simple
        (only one left site on the contig, or no site beyond the inner).
        """
        att_right, att_left_inner = closest_pair
        contig = att_right.contig

        left_sites_on_contig = [
            s for s in self.get_sites_on_contig(contig) if s.is_left
        ]

        if len(left_sites_on_contig) <= 1:
            return None

        # Determine orientation: is rlmH upstream (forward) or downstream (reverse)?
        if rlmH_start < att_left_inner.start:
            # Forward orientation: rlmH ... attR ... attL_inner ... attL_outer
            # Outer = largest start coordinate, beyond inner
            candidates = [
                s for s in left_sites_on_contig
                if s.start > att_left_inner.start
            ]
            if not candidates:
                return None
            return max(candidates, key=lambda s: s.start)
        else:
            # Reverse orientation: attL_outer ... attL_inner ... attR ... rlmH
            # Outer = smallest start coordinate, beyond inner
            candidates = [
                s for s in left_sites_on_contig
                if s.start < att_left_inner.start
            ]
            if not candidates:
                return None
            return min(candidates, key=lambda s: s.start)

    def diagnose_partial(self) -> str:
        """Diagnose why extraction failed (no valid pair on same contig).

        Returns one of:
            'right_only'   — only right att sites found
            'left_only'    — only left att sites found
            'cross_contig' — both sides found but on different contigs
            'no_sites'     — no att sites at all
        """
        right_sites = self.get_right_sites()
        left_sites = self.get_left_sites()

        if not right_sites and not left_sites:
            return "no_sites"
        if right_sites and not left_sites:
            return "right_only"
        if left_sites and not right_sites:
            return "left_only"

        # Both exist — check if any share a contig
        right_contigs = {s.contig for s in right_sites}
        left_contigs = {s.contig for s in left_sites}
        if right_contigs.isdisjoint(left_contigs):
            return "cross_contig"

        # They share a contig but find_closest_pair still failed — shouldn't
        # normally reach here, but return a generic label.
        return "unknown"


class GenomeSequences:
    """Handles genome sequence data and extraction operations."""
    
    def __init__(self, fasta_file: str):
        self.fasta_file = fasta_file
        self.sequences = self._load_sequences()
    
    def _load_sequences(self) -> Dict:
        """Load all sequences from the FASTA file."""
        sequences = SeqIO.to_dict(SeqIO.parse(self.fasta_file, "fasta"))

        return sequences
    
    def get_contig_length(self, contig: str) -> Optional[int]:
        """Return the length of a contig, or None if not found."""
        if contig in self.sequences:
            return len(self.sequences[contig].seq)
        return None

    def get_sequence(self, contig: str):
        """Get sequence for a specific contig."""
        return self.sequences.get(contig)
    
    def extract_region(self, contig: str, start: int, end: int, reverse_complement: bool = False):
        """Extract a genomic region from a contig."""

        if contig not in self.sequences:
            raise ValueError(f"Contig {contig} not found in sequences")
        
        if reverse_complement:
            # Extract with end as start and reverse complement extracted sequence
            sequence = self.sequences[contig].seq[end:start]
            sequence = sequence.reverse_complement()
        
        else:
            sequence = self.sequences[contig].seq[start:end]

        return sequence


class SCCmecExtractor:
    """Main class that coordinates SCCmec extraction from genomic data."""

    def __init__(self, fasta_file: str, gff3_file: str = None, tsv_file: str = "",
                 composite: bool = False, blast_rlmh: bool = False,
                 rlmh_ref: str = None, rlmh_positions=None,
                 genome_sequences=None, genome_db_prefix: str = None):
        self.fasta_file = fasta_file
        self.target_file = self._get_input_filename(fasta_file)
        self.composite = composite
        self._genome_db_prefix = genome_db_prefix

        # Initialise component objects
        self.genome = genome_sequences if genome_sequences is not None else GenomeSequences(fasta_file)

        # rlmH detection: pre-computed > BLAST > GFF
        if rlmh_positions is not None:
            self.genes = PrecomputedRlmH(rlmh_positions)
        elif blast_rlmh:
            self.genes = RlmHBlastAdapter(fasta_file, rlmh_ref)
        elif gff3_file:
            self.genes = GeneAnnotations(gff3_file)
        else:
            raise ValueError(
                "Either a GFF3 file (--gff) or --blast-rlmh must be provided "
                "for rlmH detection"
            )

        self.att_sites = AttSiteCollection(tsv_file, self.target_file)

    def _get_input_filename(self, fna_path: str) -> str:
        """Extract the base filename without extension from the input path."""
        return Path(fna_path).stem

    @staticmethod
    def _check_contig_edge(site: AttSite, contig_length: int,
                           threshold: int = 500) -> bool:
        """Return True if the att site is within *threshold* bp of a contig boundary."""
        return site.start <= threshold or site.end >= (contig_length - threshold)

    def _build_edge_flags(self, sites: List[AttSite],
                          contig_length: int) -> str:
        """Build a semicolon-delimited string of edge-proximity flags."""
        flags = []
        for site in sites:
            if self._check_contig_edge(site, contig_length):
                if site.start <= 500:
                    flags.append(f"{site.pattern}_near_start")
                if site.end >= (contig_length - 500):
                    flags.append(f"{site.pattern}_near_end")
        return ";".join(flags) if flags else "-"

    def _annotate_edge_flags_all_sites(self) -> str:
        """Check ALL att sites across all contigs for edge proximity."""
        flags = []
        for site in self.att_sites.sites:
            contig_len = self.genome.get_contig_length(site.contig)
            if contig_len is None:
                continue
            if self._check_contig_edge(site, contig_len):
                if site.start <= 500:
                    flags.append(f"{site.pattern}@{site.contig}_near_start")
                if site.end >= (contig_len - 500):
                    flags.append(f"{site.pattern}@{site.contig}_near_end")
        return ";".join(flags) if flags else "-"

    @staticmethod
    def _format_sites_summary(sites):
        """Format a list of AttSite objects into a semicolon-delimited string."""
        if not sites:
            return "-"
        return "; ".join(f"{s.pattern}:{s.contig}:{s.start}-{s.end}" for s in sites)

    def _get_rlmh_contigs(self):
        """Return comma-delimited string of contigs with rlmH."""
        contigs = sorted(self.genes.rlmH_positions.keys())
        return ", ".join(contigs) if contigs else "-"

    @staticmethod
    def _write_ambiguous_report(report: 'AmbiguousHitReport', report_file: str):
        """Append one row to the ambiguous report TSV. Writes header if file is new/empty."""
        output_path = Path(report_file)
        write_header = not output_path.exists() or output_path.stat().st_size == 0
        with open(report_file, 'a') as f:
            if write_header:
                f.write(AmbiguousHitReport.HEADER + "\n")
            f.write(report.to_tsv_row() + "\n")

    def _build_ambiguous_report(self, reason: str, notes: str,
                                element_size: str = "-") -> 'AmbiguousHitReport':
        """Build an AmbiguousHitReport populated with all found att sites."""
        return AmbiguousHitReport(
            input_file=self.target_file,
            reason=reason,
            right_att_sites=self._format_sites_summary(self.att_sites.get_right_sites()),
            left_att_sites=self._format_sites_summary(self.att_sites.get_left_sites()),
            rlmh_contigs=self._get_rlmh_contigs(),
            element_size=element_size,
            notes=notes,
        )

    def _detect_ccr_between(self, contig: str, pos_a: int, pos_b: int) -> bool:
        """Check if ccr genes exist between two positions on a contig.

        Uses BLAST to detect ccr genes in the genome, then checks if any
        hit falls between pos_a and pos_b on the specified contig.

        Returns True if at least one valid ccr gene (ccrA+ccrB pair or ccrC)
        is found in the region.
        """
        region_start = min(pos_a, pos_b)
        region_end = max(pos_a, pos_b)

        runner = BlastRunner()

        # Reuse shared DB when provided, otherwise create a temporary one
        owns_db = self._genome_db_prefix is None
        if owns_db:
            tmp_dir = tempfile.mkdtemp(prefix="sccmec_ccr_")
            db_prefix = os.path.join(tmp_dir, "genome_db")
        else:
            tmp_dir = None
            db_prefix = self._genome_db_prefix

        try:
            if owns_db:
                runner.create_db(self.fasta_file, db_prefix)

            with get_default_ref("ccr_genes.fasta") as ref_path:
                results_file = runner.run_blastn(str(ref_path), db_prefix)
                hits = parse_blast_output(results_file)

                # Get reference lengths for coverage calculation
                ref_lengths = {}
                for record in SeqIO.parse(str(ref_path), "fasta"):
                    ref_lengths[record.id] = len(record.seq)

                runner.cleanup_file(results_file)

            # Filter hits: 70% identity (novel threshold), 75% coverage
            filtered = filter_hits(hits, min_pident=70.0, min_coverage=0.75,
                                   ref_lengths=ref_lengths)

            # Check which hits are on our contig and within the region
            ccr_in_region = set()
            for hit in filtered:
                if hit.sseqid != contig:
                    continue
                hit_start = min(hit.sstart, hit.send)
                hit_end = max(hit.sstart, hit.send)
                # Hit overlaps the region
                if hit_start <= region_end and hit_end >= region_start:
                    # Extract gene type: ccrA, ccrB, or ccrC
                    gene_type = hit.qseqid.rstrip('0123456789')
                    ccr_in_region.add(gene_type)

            # Valid ccr: (ccrA AND ccrB) OR ccrC
            has_AB = 'ccrA' in ccr_in_region and 'ccrB' in ccr_in_region
            has_C = 'ccrC' in ccr_in_region
            return has_AB or has_C

        finally:
            if owns_db:
                runner.cleanup_db(db_prefix)
                try:
                    os.rmdir(tmp_dir)
                except OSError:
                    pass

    MAX_ATTL_RLMH_DISTANCE = 120_000

    def _has_ccr_in_genome(self) -> bool:
        """Check if valid ccr genes exist anywhere in the genome.

        Uses BLAST to detect ccr genes without position or contig
        constraints.  This accommodates fragmented assemblies where ccr
        and rlmH may land on different contigs.

        Returns True if valid ccr (ccrA+ccrB or ccrC) found.
        """
        runner = BlastRunner()

        owns_db = self._genome_db_prefix is None
        if owns_db:
            tmp_dir = tempfile.mkdtemp(prefix="sccmec_ccr_")
            db_prefix = os.path.join(tmp_dir, "genome_db")
        else:
            tmp_dir = None
            db_prefix = self._genome_db_prefix

        try:
            if owns_db:
                runner.create_db(self.fasta_file, db_prefix)

            with get_default_ref("ccr_genes.fasta") as ref_path:
                results_file = runner.run_blastn(str(ref_path), db_prefix)
                hits = parse_blast_output(results_file)

                ref_lengths = {}
                for record in SeqIO.parse(str(ref_path), "fasta"):
                    ref_lengths[record.id] = len(record.seq)

                runner.cleanup_file(results_file)

            filtered = filter_hits(hits, min_pident=70.0, min_coverage=0.75,
                                   ref_lengths=ref_lengths)

            ccr_types = set()
            for hit in filtered:
                gene_type = hit.qseqid.rstrip('0123456789')
                ccr_types.add(gene_type)

            has_AB = 'ccrA' in ccr_types and 'ccrB' in ccr_types
            has_C = 'ccrC' in ccr_types
            return has_AB or has_C

        finally:
            if owns_db:
                runner.cleanup_db(db_prefix)
                try:
                    os.rmdir(tmp_dir)
                except OSError:
                    pass

    MIN_ELEMENT_SIZE = 1_000       # 1 kb — anything smaller is a pattern overlap artefact
    MAX_ELEMENT_SIZE = 200_000     # 200 kb — anything larger is a spurious cross-genome match
    MAX_COMPOSITE_SIZE = 250_000   # 250 kb — two tandem SCCs can reach ~170 kb; beyond this is mispaired outer attL

    def _attempt_left_only_recovery(self, report: ExtractionReport,
                                     output_dir: str, report_file: str = None,
                                     ambiguous_report_file: str = None) -> bool:
        """Attempt extraction when only attL is found (no attR).

        Assumes attR is within rlmH. Validates by checking:
        1. attL and rlmH are on the same contig
        2. Distance between attL and rlmH is <= 120kb
        3. ccr genes are present between attL and rlmH

        Returns True if extraction succeeded.
        """
        left_sites = self.att_sites.get_left_sites()

        for att_left in left_sites:
            contig = att_left.contig

            # Check 1: rlmH on same contig
            if not self.genes.has_rlmH(contig):
                continue

            rlmH_start = self.genes.get_rlmH_start(contig)
            rlmH_warning = self.genes.get_rlmH_copy_warning(contig)
            if rlmH_warning:
                print(f"Warning: {rlmH_warning} for {self.target_file}", file=sys.stderr)

            # Check 2: distance <= 120kb
            distance = abs(rlmH_start - att_left.end)
            if distance > self.MAX_ATTL_RLMH_DISTANCE:
                continue

            # Check 3: ccr genes between attL and rlmH
            if not self._detect_ccr_between(contig, rlmH_start, att_left.end):
                continue

            # All checks passed — proceed with extraction
            output_file = os.path.join(output_dir, f"{self.target_file}_SCCmec.fasta")
            if os.path.exists(output_file):
                report.status = "skipped"
                report.notes = "Output file already exists (left_only recovery)"
                if rlmH_warning:
                    report.notes = f"{rlmH_warning}; {report.notes}"
                if report_file:
                    self._write_report(report, report_file)
                return False

            # Create a synthetic attR representing the inferred rlmH position
            inferred_attr = AttSite(
                pattern="attR_inferred",
                contig=contig,
                start=rlmH_start,
                end=rlmH_start
            )

            # Populate report
            report.contig = contig
            report.attr_pattern = "attR_inferred"
            report.attr_start = str(rlmH_start)
            report.attr_end = str(rlmH_start)
            report.attl_pattern = att_left.pattern
            report.attl_start = str(att_left.start)
            report.attl_end = str(att_left.end)

            try:
                # Composite detection (reuse existing logic)
                best_pair = (inferred_attr, att_left)
                outer_left = self.att_sites.find_outer_left_site(best_pair, rlmH_start)
                if outer_left:
                    report.outer_attl_pattern = outer_left.pattern
                    report.outer_attl_start = str(outer_left.start)
                    report.outer_attl_end = str(outer_left.end)

                # Composite size sanity check — discard mispaired outer attL
                if outer_left:
                    composite_size = abs(rlmH_start - outer_left.end)
                    if composite_size > self.MAX_COMPOSITE_SIZE:
                        outer_left = None
                        report.notes = (f"Outer attL discarded: composite size {composite_size} bp "
                                        f"exceeds maximum ({self.MAX_COMPOSITE_SIZE} bp)")

                # Outer CCR check for composite validation
                outer_ccr = False
                if outer_left:
                    outer_ccr = self._detect_ccr_between(contig, att_left.end, outer_left.end)

                if outer_ccr and self.composite and outer_left:
                    extract_left = outer_left
                    report.status = "composite_fallback_extracted"
                    report.is_composite = "True"
                else:
                    extract_left = att_left
                    report.status = "fallback_extracted"
                    if self.composite and outer_left and not outer_ccr:
                        report.notes = "Composite downgraded: no ccr genes in outer region"

                # Determine coordinates (att_right param is unused internally)
                start_extract, end_extract, reverse_complement = \
                    self._determine_extraction_coordinates(
                        rlmH_start, inferred_attr, extract_left
                    )

                # Extract sequence
                extracted_seq = self.genome.extract_region(
                    contig, start_extract, end_extract, reverse_complement
                )

                # Element size
                element_size = abs(rlmH_start - att_left.end)
                report.element_size = str(element_size)
                if outer_left:
                    composite_size = abs(rlmH_start - outer_left.end)
                    report.composite_size = str(composite_size)

                # Size sanity check — always use element_size (inner element
                # must be plausible regardless of composite)
                if element_size < self.MIN_ELEMENT_SIZE or element_size > self.MAX_ELEMENT_SIZE:
                    report.status = "failed"
                    report.partial_type = "size_out_of_range"
                    report.notes = (f"Left-only element size {element_size} bp outside valid range "
                                    f"({self.MIN_ELEMENT_SIZE}-{self.MAX_ELEMENT_SIZE} bp)")
                    # Detect origin-spanning on circular chromosomes
                    if element_size > self.MAX_ELEMENT_SIZE:
                        contig_len = self.genome.get_contig_length(contig)
                        if contig_len and element_size > contig_len // 2:
                            wrap_around = contig_len - element_size
                            if self.MIN_ELEMENT_SIZE <= wrap_around <= self.MAX_ELEMENT_SIZE:
                                report.partial_type = "origin_spanning"
                                report.notes = (
                                    f"Probable origin-spanning element on circular chromosome "
                                    f"(contig {contig_len} bp). Linear distance {element_size} bp "
                                    f"but wrap-around distance ~{wrap_around} bp")
                    if rlmH_warning:
                        report.notes = f"{rlmH_warning}; {report.notes}" if report.notes != "-" else rlmH_warning
                    if report_file:
                        self._write_report(report, report_file)
                    if ambiguous_report_file:
                        amb = self._build_ambiguous_report(
                            report.partial_type, report.notes,
                            element_size=str(element_size))
                        self._write_ambiguous_report(amb, ambiguous_report_file)
                    return False

                # Edge flags
                contig_len = self.genome.get_contig_length(contig)
                if contig_len:
                    sites_on_contig = self.att_sites.get_sites_on_contig(contig)
                    report.contig_edge_flags = self._build_edge_flags(
                        sites_on_contig, contig_len
                    )

                # Create record and write
                record = self._create_sequence_record(
                    extracted_seq, inferred_attr, extract_left,
                    start_extract, end_extract
                )

                with open(output_file, "w") as fasta_output:
                    SeqIO.write([record], fasta_output, "fasta")

                left_only_note = "Left-only recovery: attR inferred from rlmH, validated by ccr"
                if report.notes != "-":
                    report.notes = f"{report.notes}; {left_only_note}"
                else:
                    report.notes = left_only_note
                if rlmH_warning:
                    report.notes = f"{rlmH_warning}; {report.notes}"
                print(f"Successfully processed {self.target_file} (left-only recovery): "
                      f"Extracted sequence of length {len(extracted_seq)} bp")

                if report_file:
                    self._write_report(report, report_file)
                return True

            except Exception as e:
                print(f"Error in left-only recovery for {self.target_file}: {str(e)}",
                      file=sys.stderr)
                report.status = "failed"
                report.notes = f"Left-only recovery failed: {str(e)}"
                if rlmH_warning:
                    report.notes = f"{rlmH_warning}; {report.notes}"
                if report_file:
                    self._write_report(report, report_file)
                return False

        # No valid left site passed all checks
        return False

    @staticmethod
    def _write_report(report: ExtractionReport, report_file: str):
        """Append one row to the report TSV. Writes header if file is new/empty."""
        output_path = Path(report_file)
        write_header = not output_path.exists() or output_path.stat().st_size == 0

        with open(report_file, 'a') as f:
            if write_header:
                f.write(ExtractionReport.HEADER + "\n")
            f.write(report.to_tsv_row() + "\n")
    
    def _determine_extraction_coordinates(self, rlmH_start: int, att_right: AttSite, att_left: AttSite) -> Tuple[int, int, bool]:
        """Determine the coordinates for SCCmec extraction."""
        # Default extraction with padding
        start_extract = rlmH_start - 30
        end_extract = att_left.end + 30
        reverse_complement = False
        
        # Check if we need reverse complement (SCCmec on reverse strand)
        if start_extract > end_extract:
            start_extract = rlmH_start + 600
            end_extract = att_left.end - 60
            reverse_complement = True
        
        return start_extract, end_extract, reverse_complement
    
    def _create_sequence_record(self, sequence, att_right: AttSite, att_left: AttSite, 
                              start: int, end: int):
        """Create a SeqRecord with appropriate ID and description."""
        record_id = f"{self.target_file}_{att_right.contig}_{start}_{end}"
        description = f"attR:{att_right}_attL:{att_left}"
        
        return SeqIO.SeqRecord(
            sequence,
            id=record_id,
            description=description
        )
    
    def extract_sccmec(self, output_dir: str, report_file: str = None,
                        ambiguous_report_file: str = None) -> bool:
        """Extract SCCmec sequence and save to file.

        Parameters
        ----------
        output_dir : str
            Directory to write extracted FASTA files.
        report_file : str, optional
            Path to a TSV report file. When provided, a row is appended for
            every genome processed (success, failure, or skip).
        ambiguous_report_file : str, optional
            Path to TSV file for actionable failure details. When provided,
            a row listing all found att sites is appended for every genome
            that fails extraction for an actionable reason.
        """
        os.makedirs(output_dir, exist_ok=True)
        report = ExtractionReport(input_file=self.target_file)

        # --- Failure: no valid att sites ---
        if not self.att_sites.has_valid_sites():
            partial_type = self.att_sites.diagnose_partial()

            # Try left-only recovery if we have attL but no attR
            if partial_type == "left_only":
                recovery_result = self._attempt_left_only_recovery(
                    report, output_dir, report_file,
                    ambiguous_report_file=ambiguous_report_file,
                )
                if recovery_result:
                    return True
                # If recovery returned False without writing a report, fall through
                if report.status:
                    return False

            print(f"Warning: Missing required att sites for {self.target_file}", file=sys.stderr)
            report.status = "failed"
            report.contig_edge_flags = self._annotate_edge_flags_all_sites()

            # Prioritise ccr absence over att site diagnosis
            write_ambiguous = False
            if partial_type == "no_sites" and self._has_ccr_in_genome():
                report.partial_type = partial_type
                report.notes = "Missing required att sites (right and/or left); ccr genes present"
                write_ambiguous = True
            elif partial_type == "no_sites":
                report.partial_type = partial_type
                report.notes = "Missing required att sites (right and/or left)"
            elif self._has_ccr_in_genome():
                report.partial_type = partial_type
                report.notes = "Missing required att sites (right and/or left)"
                write_ambiguous = True
            else:
                report.partial_type = "no_ccr"
                report.notes = (f"No ccr genes detected in genome (att site diagnosis: {partial_type})")

            if report_file:
                self._write_report(report, report_file)
            if ambiguous_report_file and write_ambiguous:
                amb = self._build_ambiguous_report(report.partial_type, report.notes)
                self._write_ambiguous_report(amb, ambiguous_report_file)
            return False

        # --- Try to find pair, skipping overlapping att sites at attB ---
        best_pair = self.att_sites.find_closest_pair(min_distance=self.MIN_ELEMENT_SIZE)
        if not best_pair:
            # Check if there WAS a pair that was filtered (attB overlap)
            unfiltered_pair = self.att_sites.find_closest_pair()
            if unfiltered_pair:
                # Pairs exist but all below min distance → attB overlap artefact
                att_r, att_l = unfiltered_pair
                overlap_size = abs(att_r.start - att_l.end)
                report.status = "failed"
                report.partial_type = "size_out_of_range"
                report.contig = att_r.contig
                report.attr_pattern = att_r.pattern
                report.attr_start = str(att_r.start)
                report.attr_end = str(att_r.end)
                report.attl_pattern = att_l.pattern
                report.attl_start = str(att_l.start)
                report.attl_end = str(att_l.end)
                report.element_size = str(overlap_size)
                report.notes = (f"Overlapping att site patterns at chromosomal attB "
                                f"({att_r.pattern}/{att_l.pattern}, {overlap_size} bp); "
                                f"no valid pair above {self.MIN_ELEMENT_SIZE} bp")
                if report_file:
                    self._write_report(report, report_file)
                if ambiguous_report_file:
                    amb = self._build_ambiguous_report(
                        "size_out_of_range", report.notes,
                        element_size=str(overlap_size))
                    self._write_ambiguous_report(amb, ambiguous_report_file)
                return False

            # Genuinely no pair — existing diagnosis logic
            partial_type = self.att_sites.diagnose_partial()

            # Try left-only recovery for cross-contig cases too
            if partial_type in ("left_only", "cross_contig"):
                recovery_result = self._attempt_left_only_recovery(
                    report, output_dir, report_file,
                    ambiguous_report_file=ambiguous_report_file,
                )
                if recovery_result:
                    return True
                if report.status:
                    return False

            print(f"Warning: No valid attR-attL pair found for {self.target_file}", file=sys.stderr)
            report.status = "failed"
            report.contig_edge_flags = self._annotate_edge_flags_all_sites()

            # Prioritise ccr absence over att site diagnosis
            if self._has_ccr_in_genome():
                report.partial_type = partial_type
                report.notes = "No attR-attL pair on same contig"
                write_ambiguous = True
            else:
                report.partial_type = "no_ccr"
                report.notes = (f"No ccr genes detected in genome (att site diagnosis: {partial_type})")
                write_ambiguous = False

            if report_file:
                self._write_report(report, report_file)
            if ambiguous_report_file and write_ambiguous:
                amb = self._build_ambiguous_report(report.partial_type, report.notes)
                self._write_ambiguous_report(amb, ambiguous_report_file)
            return False

        att_right, att_left = best_pair
        contig = att_right.contig
        report.contig = contig
        report.attr_pattern = att_right.pattern
        report.attr_start = str(att_right.start)
        report.attr_end = str(att_right.end)
        report.attl_pattern = att_left.pattern
        report.attl_start = str(att_left.start)
        report.attl_end = str(att_left.end)

        # --- Failure: no rlmH ---
        if not self.genes.has_rlmH(contig):
            print(f"Warning: No rlmH gene found for {self.target_file} on {contig}", file=sys.stderr)
            report.status = "failed"
            report.notes = f"No rlmH gene found on {contig}"
            contig_len = self.genome.get_contig_length(contig)
            if contig_len:
                sites_on_contig = self.att_sites.get_sites_on_contig(contig)
                report.contig_edge_flags = self._build_edge_flags(sites_on_contig, contig_len)
            if report_file:
                self._write_report(report, report_file)
            return False

        rlmH_start = self.genes.get_rlmH_start(contig)
        rlmH_warning = self.genes.get_rlmH_copy_warning(contig)
        if rlmH_warning:
            print(f"Warning: {rlmH_warning} for {self.target_file}", file=sys.stderr)
        output_file = os.path.join(output_dir, f"{self.target_file}_SCCmec.fasta")

        # --- Skip: output already exists ---
        if os.path.exists(output_file):
            print(f"Skipping {self.target_file}: Output file already exists", file=sys.stderr)
            report.status = "skipped"
            report.notes = "Output file already exists"
            if report_file:
                self._write_report(report, report_file)
            return False

        try:
            # Composite detection
            outer_left = self.att_sites.find_outer_left_site(best_pair, rlmH_start)
            if outer_left:
                report.outer_attl_pattern = outer_left.pattern
                report.outer_attl_start = str(outer_left.start)
                report.outer_attl_end = str(outer_left.end)

            # Composite size sanity check — discard mispaired outer attL
            if outer_left:
                composite_size = abs(att_right.start - outer_left.end)
                if composite_size > self.MAX_COMPOSITE_SIZE:
                    outer_left = None
                    report.notes = (f"Outer attL discarded: composite size {composite_size} bp "
                                    f"exceeds maximum ({self.MAX_COMPOSITE_SIZE} bp)")

            # CCR gene checks — validate inner and (if applicable) outer regions
            inner_ccr = self._detect_ccr_between(contig, att_right.start, att_left.end)
            outer_ccr = False
            if outer_left:
                outer_ccr = self._detect_ccr_between(contig, att_left.end, outer_left.end)

            # Determine extraction boundary based on CCR results
            if inner_ccr and outer_ccr and self.composite and outer_left:
                # Scenario 1: Genuine composite — both regions have ccr
                extract_left = outer_left
                report.status = "composite_extracted"
                report.is_composite = "True"
            elif inner_ccr:
                # Scenario 2: Standard extraction (or downgraded composite)
                extract_left = att_left
                report.status = "extracted"
                if outer_left and not outer_ccr:
                    report.notes = "Composite downgraded: no ccr genes in outer region"
            elif outer_left and outer_ccr:
                # Scenario 3: Inner attL is scar — extract to outer attL as single element
                extract_left = outer_left
                report.status = "extracted"
                report.notes = ("Inner attL bypassed (no ccr in inner region); "
                                "extracted to outer attL")
                # Update report attL fields to reflect actual boundary
                report.attl_pattern = outer_left.pattern
                report.attl_start = str(outer_left.start)
                report.attl_end = str(outer_left.end)
                report.outer_attl_pattern = "-"
                report.outer_attl_start = "-"
                report.outer_attl_end = "-"
            else:
                # Scenario 4: No CCR anywhere
                element_size = abs(att_right.start - att_left.end)
                report.element_size = str(element_size)
                report.status = "no_ccr_element"
                report.notes = "Valid att site pair but no ccr genes detected in element"
                if rlmH_warning:
                    report.notes = f"{rlmH_warning}; {report.notes}"
                if report_file:
                    self._write_report(report, report_file)
                if ambiguous_report_file:
                    amb = self._build_ambiguous_report(
                        "no_ccr_element", report.notes,
                        element_size=str(element_size))
                    self._write_ambiguous_report(amb, ambiguous_report_file)
                return False

            # Determine extraction coordinates
            start_extract, end_extract, reverse_complement = self._determine_extraction_coordinates(
                rlmH_start, att_right, extract_left
            )

            # Extract sequence
            extracted_seq = self.genome.extract_region(
                contig, start_extract, end_extract, reverse_complement
            )

            # Element sizes
            element_size = abs(att_right.start - extract_left.end)
            report.element_size = str(element_size)
            if report.status == "composite_extracted" and outer_left:
                composite_size = abs(att_right.start - outer_left.end)
                report.composite_size = str(composite_size)

            # Size sanity check
            if element_size < self.MIN_ELEMENT_SIZE or element_size > self.MAX_ELEMENT_SIZE:
                report.status = "failed"
                report.partial_type = "size_out_of_range"
                report.notes = (f"Element size {element_size} bp outside valid range "
                                f"({self.MIN_ELEMENT_SIZE}-{self.MAX_ELEMENT_SIZE} bp)")
                # Detect origin-spanning: element covers most of a circular
                # chromosome but the wrap-around distance is plausible
                if element_size > self.MAX_ELEMENT_SIZE:
                    contig_len = self.genome.get_contig_length(contig)
                    if contig_len and element_size > contig_len // 2:
                        wrap_around = contig_len - element_size
                        if self.MIN_ELEMENT_SIZE <= wrap_around <= self.MAX_ELEMENT_SIZE:
                            report.partial_type = "origin_spanning"
                            report.notes = (
                                f"Probable origin-spanning element on circular chromosome "
                                f"(contig {contig_len} bp). Linear distance {element_size} bp "
                                f"but wrap-around distance ~{wrap_around} bp")
                if rlmH_warning:
                    report.notes = f"{rlmH_warning}; {report.notes}" if report.notes != "-" else rlmH_warning
                if report_file:
                    self._write_report(report, report_file)
                if ambiguous_report_file:
                    amb = self._build_ambiguous_report(
                        report.partial_type, report.notes,
                        element_size=str(element_size))
                    self._write_ambiguous_report(amb, ambiguous_report_file)
                return False

            # Contig-edge flags
            contig_len = self.genome.get_contig_length(contig)
            if contig_len:
                sites_on_contig = self.att_sites.get_sites_on_contig(contig)
                report.contig_edge_flags = self._build_edge_flags(sites_on_contig, contig_len)

            # Create record and write
            record = self._create_sequence_record(
                extracted_seq, att_right, extract_left, start_extract, end_extract
            )

            with open(output_file, "w") as fasta_output:
                SeqIO.write([record], fasta_output, "fasta")

            print(f"Successfully processed {self.target_file}: Extracted sequence of length {len(extracted_seq)} bp")
            if rlmH_warning:
                report.notes = f"{rlmH_warning}; {report.notes}" if report.notes != "-" else rlmH_warning
            if report_file:
                self._write_report(report, report_file)
            return True

        except Exception as e:
            print(f"Error processing {self.target_file}: {str(e)}", file=sys.stderr)
            report.status = "failed"
            report.notes = str(e)
            if rlmH_warning:
                report.notes = f"{rlmH_warning}; {report.notes}"
            if report_file:
                self._write_report(report, report_file)
            return False


def main():
    parser = argparse.ArgumentParser(description="Extract SCCmec sequences based on att sites")
    parser.add_argument("-f", "--fna", required=True, help=".fasta or .fna file containing genome sequence")
    parser.add_argument("-g", "--gff", help=".gff3 file containing gene annotation information")
    parser.add_argument("-a", "--att", required=True, help=".tsv file containing att site location information")
    parser.add_argument("-s", "--sccmec", required=True, help="Output directory for SCCmec sequences")
    parser.add_argument("--composite", action="store_true",
                        help="Extract to outermost boundary for composite elements")
    parser.add_argument("-r", "--report", default=None,
                        help="Output TSV file for extraction report (appends for batch)")
    parser.add_argument(
        "--blast-rlmh",
        action="store_true",
        default=False,
        help="Use BLAST to detect rlmH gene (auto-enabled when no --gff provided)",
    )
    parser.add_argument(
        "--rlmh-ref",
        help="Custom rlmH reference FASTA for BLAST detection (optional)",
    )
    args = parser.parse_args()

    # Validate inputs
    validator = InputValidator()

    validator.validate_fasta_file(args.fna)

    if args.gff:
        validator.validate_gff_file(args.gff)

    if args.att:
        validator.validate_tsv_file(args.att)

    # Auto-enable BLAST rlmH when no GFF provided
    blast_rlmh = args.blast_rlmh
    if not args.gff and not blast_rlmh:
        try:
            from sccmecextractor.blast_utils import BlastRunner
            BlastRunner._check_blast_installed()
            blast_rlmh = True
            print("No GFF provided; auto-enabling BLAST-based rlmH detection")
        except Exception:
            print(
                "ERROR: No GFF provided and BLAST+ not available. "
                "Either provide --gff or install BLAST+ for --blast-rlmh.",
                file=sys.stderr,
            )
            sys.exit(1)

    # Derive ambiguous report path from report path
    ambiguous_report = None
    if args.report:
        report_dir = os.path.dirname(args.report) or "."
        ambiguous_report = os.path.join(report_dir, "ambiguous_att_sites.tsv")

    # Create extractor and process
    extractor = SCCmecExtractor(
        args.fna, gff3_file=args.gff, tsv_file=args.att,
        composite=args.composite, blast_rlmh=blast_rlmh,
        rlmh_ref=args.rlmh_ref,
    )
    success = extractor.extract_sccmec(
        args.sccmec, report_file=args.report,
        ambiguous_report_file=ambiguous_report,
    )

    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()