#!/usr/bin/env python

import argparse
import os
import sys

from typing import Dict, List, Optional, Tuple
from pathlib import Path
from Bio import SeqIO
from collections import defaultdict

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


class GeneAnnotations:
    """Handles parsing and storing gene annotation information."""
    
    def __init__(self, gff3_file: str):
        self.gff3_file = gff3_file
        self.rlmH_positions = self._parse_rlmH_genes()
    
    def _parse_rlmH_genes(self) -> Dict[str, int]:
        """Parse GFF3 file to extract rlmH gene start positions."""
        rlmH_info = {}
        
        with open(self.gff3_file, 'r') as gff3:
            for line in gff3:
                if line.startswith("#"):
                    continue
                    
                columns = line.strip().split("\t")
                if columns[2] == "gene" and "gene=rlmH" in columns[-1]:
                    contig = columns[0]
                    start_position = int(columns[3])
                    rlmH_info[contig] = start_position
        
        if not rlmH_info:
            print(f"Warning: No rlmH genes found in {self.gff3_file}", file=sys.stderr)
        
        return rlmH_info
    
    def get_rlmH_start(self, contig: str) -> Optional[int]:
        """Get the start position of rlmH gene on a specific contig."""
        return self.rlmH_positions.get(contig)
    
    def has_rlmH(self, contig: str) -> bool:
        """Check if a contig has an rlmH gene."""
        return contig in self.rlmH_positions


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
    
    def find_closest_pair(self) -> Optional[Tuple[AttSite, AttSite]]:
        """Find the closest attR-attL pair on the same contig."""
        right_sites = self.get_right_sites()
        left_sites = self.get_left_sites()
        
        if not right_sites or not left_sites:
            return None
        
        best_pair = None
        min_distance = float('inf')
        
        for right in right_sites:
            for left in left_sites:
                if right.contig == left.contig:
                    distance = right.distance_to(left)
                    if distance < min_distance:
                        min_distance = distance
                        best_pair = (right, left)
        
        return best_pair
    
    def has_valid_sites(self) -> bool:
        """Check if we have both right and left sites."""
        return bool(self.get_right_sites() and self.get_left_sites())


class GenomeSequences:
    """Handles genome sequence data and extraction operations."""
    
    def __init__(self, fasta_file: str):
        self.fasta_file = fasta_file
        self.sequences = self._load_sequences()
    
    def _load_sequences(self) -> Dict:
        """Load all sequences from the FASTA file."""
        sequences = SeqIO.to_dict(SeqIO.parse(self.fasta_file, "fasta"))

        return sequences
    
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
    
    def __init__(self, fasta_file: str, gff3_file: str, tsv_file: str):
        self.fasta_file = fasta_file
        self.target_file = self._get_input_filename(fasta_file)
        
        # Initialise component objects
        self.genome = GenomeSequences(fasta_file)
        self.genes = GeneAnnotations(gff3_file)
        self.att_sites = AttSiteCollection(tsv_file, self.target_file)
    
    def _get_input_filename(self, fna_path: str) -> str:
        """Extract the base filename without extension from the input path."""
        base = os.path.basename(fna_path)
        base_full = os.path.splitext(base)[0]
        return base_full.split(".")[0]
    
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
    
    def extract_sccmec(self, output_dir: str) -> bool:
        """Extract SCCmec sequence and save to file."""
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Check if we have valid att sites
        if not self.att_sites.has_valid_sites():
            print(f"Warning: Missing required att sites for {self.target_file}", file=sys.stderr)
            return False
        
        # Find the best att site pair
        best_pair = self.att_sites.find_closest_pair()
        if not best_pair:
            print(f"Warning: No valid attR-attL pair found for {self.target_file}", file=sys.stderr)
            return False
        
        att_right, att_left = best_pair
        contig = att_right.contig
        
        # Check for rlmH gene
        if not self.genes.has_rlmH(contig):
            print(f"Warning: No rlmH gene found for {self.target_file} on {contig}", file=sys.stderr)
            return False
        
        rlmH_start = self.genes.get_rlmH_start(contig)
        output_file = os.path.join(output_dir, f"{self.target_file}_SCCmec.fasta")
        
        # Skip if file already exists
        if os.path.exists(output_file):
            print(f"Skipping {self.target_file}: Output file already exists", file=sys.stderr)
            return False
        
        try:
            # Determine extraction coordinates
            start_extract, end_extract, reverse_complement = self._determine_extraction_coordinates(
                rlmH_start, att_right, att_left
            )
            
            # Extract sequence
            extracted_seq = self.genome.extract_region(
                contig, start_extract, end_extract, reverse_complement
            )
            
            # Create sequence record
            record = self._create_sequence_record(
                extracted_seq, att_right, att_left, start_extract, end_extract
            )
            
            # Write to file
            with open(output_file, "w") as fasta_output:
                SeqIO.write([record], fasta_output, "fasta")
            
            print(f"Successfully processed {self.target_file}: Extracted sequence of length {len(extracted_seq)} bp")
            return True
            
        except Exception as e:
            print(f"Error processing {self.target_file}: {str(e)}", file=sys.stderr)
            return False


def main():
    parser = argparse.ArgumentParser(description="Extract SCCmec sequences based on att sites")
    parser.add_argument("-f", "--fna", required=True, help=".fasta or .fna file containing genome sequence")
    parser.add_argument("-g", "--gff", required=True, help=".gff3 file containing gene annotation information")
    parser.add_argument("-a", "--att", required=True, help=".tsv file containing att site location information")
    parser.add_argument("-s", "--sccmec", required=True, help="Output directory for SCCmec sequences")
    args = parser.parse_args()
    
    # Validate inputs
    validator = InputValidator()
    
    validator.validate_fasta_file(args.fna)
    
    if args.gff:
        validator.validate_gff_file(args.gff)

    if args.att:
        validator.validate_tsv_file(args.att)

    # Create extractor and process
    extractor = SCCmecExtractor(args.fna, args.gff, args.att)
    success = extractor.extract_sccmec(args.sccmec)
    
    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()