#!/usr/bin/env python

import re
import argparse

from Bio import SeqIO
from typing import List, Dict, Set, Tuple
from pathlib import Path



class InputValidator:
    """Check input files are valid"""
    
    def validate_fasta_file(self, input_fasta):
        """Validate that FASTA file exists and is readable."""
        fasta = Path(input_fasta)
            
        if not fasta.exists():
            raise FileNotFoundError(f"FASTA file not found: {input_fasta}")
    
        if not fasta.is_file():
            raise ValueError(f"Path exists but is not a file: {input_fasta}")
            
        if fasta.stat().st_size == 0:
            raise ValueError(f"File is empty: {input_fasta}")
        
    def validate_gff_file(self, input_gff):
        """Validate the GFF file exists and is readable"""
        gff = Path(input_gff)
        
        if not gff.exists():
            raise FileNotFoundError(f"GFF file not found: {input_gff}")
        
        if not gff.is_file():
            raise ValueError(f"Path exists but is not a file: {input_gff}")
            
        if gff.stat().st_size == 0:
            raise ValueError(f"File is empty: {input_gff}")

class AttSite:
    """Represents a single att site match with its properties."""
    
    def __init__(self, pattern_name: str, contig: str, start: int, end: int, match_seq: str):
        self.pattern_name = pattern_name
        self.contig = contig
        self.start = start  # 1-based coordinates
        self.end = end
        self.match_seq = match_seq
        self.within_rlmH = False  # Will be set by gene checker
    
    def to_tsv_line(self, input_file_name: str) -> str:
        """Convert the AttSite to a TSV line format."""
        return f"{input_file_name}\t{self.pattern_name}\t{self.contig}\t{self.start}\t{self.end}\t{self.match_seq}"
    
    def __str__(self) -> str:
        return f"AttSite({self.pattern_name} at {self.contig}:{self.start}-{self.end})"


class GeneAnnotationParser:
    """Handles parsing and querying of gene annotations from GFF3 files."""
    
    def __init__(self, gff3_file: str):
        self.gff3_file = gff3_file
        self.rlmH_genes = self._parse_gff3()
    
    def _parse_gff3(self) -> Dict[str, Set[Tuple[int, int]]]:
        """Parse GFF3 file to extract rlmH gene locations."""
        genes = {}
        
        with open(self.gff3_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) == 9 and fields[2] == 'gene':
                    contig = fields[0]
                    start, end = int(fields[3]), int(fields[4])
                    attributes = dict(item.split('=', 1) for item in fields[8].split(';'))
                    
                    if attributes.get('gene') == 'rlmH':
                        genes.setdefault(contig, set()).add((start, end))
        
        return genes
    
    def is_within_rlmH(self, site: AttSite, tolerance: int = 10) -> bool:
        """Check if an AttSite falls within an rlmH gene (with tolerance)."""
        if site.contig not in self.rlmH_genes:
            return False
        
        for gene_start, gene_end in self.rlmH_genes[site.contig]:
            if gene_start <= site.start <= site.end <= (gene_end + tolerance):
                return True
        
        return False


class AttSiteFinder:
    """Main class that coordinates att site searching and analysis."""
    
    def __init__(self, fasta_file: str, gff3_file: str = None):
        self.fasta_file = fasta_file
        self.sequences = self._parse_fasta()
        self.gene_parser = GeneAnnotationParser(gff3_file) if gff3_file else None
        
        # Patterns for att sites
        self.patterns = {
            'attR': re.compile('GC[AG]TATCA[TC]AA[GA]TGATGCGGTTT'),
            'cattR': re.compile('AAACCGCATCA[CT]TT[GA]TGATA[CT]GC'),
            'attR2': re.compile('GC[GT]TA[TC]CA[TC]AAATAAAACTAAAA'),
            'cattR2': re.compile('TTTTAGTTTTATTT[GA]TG[AG]TA[AC]GC'),
            'attL': re.compile('AACC[TG]CATCA[TC][TC][AT][AC]C[TC]GATAAG[CT]'),
            'cattL': re.compile('[AG]CTTATC[GA]G[GT][AT][GA][GA]TGATG[CA]GGTT'),
            'attL2': re.compile('[TA][TA]TT[TA][AG][GC][TCA]AA[TA]AT[CA]ACT[GA][TGA][TC]A[AG]GG'),
            'cattL2': re.compile('CC[CT]T[GA][ACT][TC]AGT[TG]AT[TA]TT[TGA][CG][CT][TA]AA[TA][TA]'),
        }
    
    def _parse_fasta(self) -> Dict[str, str]:
        """Parse FASTA file to get sequences by contig."""
        sequences = {}
        
        for record in SeqIO.parse(self.fasta_file, 'fasta'):
            contig = record.id.split()[0]
            sequences[contig] = str(record.seq)
        
        return sequences
    
    def find_all_sites(self) -> List[AttSite]:
        """Search for all att sites in all sequences."""
        sites = []
        
        for contig, sequence in self.sequences.items():
            print(f"Processing contig: {contig}")
            
            for pattern_name, compiled_pattern in self.patterns.items():
                matches = list(compiled_pattern.finditer(sequence))
                
                if matches:
                    print(f"  Pattern {pattern_name}: {len(matches)} matches")
                
                for match in matches:
                    site = AttSite(
                        pattern_name=pattern_name,
                        contig=contig,
                        start=match.start() + 1,  # Convert to 1-based
                        end=match.end(),
                        match_seq=match.group()
                    )
                    
                    # Check if site is within rlmH gene
                    if self.gene_parser:
                        site.within_rlmH = self.gene_parser.is_within_rlmH(site)
                    
                    sites.append(site)
                    print(f"    Found at position {site.start}-{site.end}")
                    
                    if site.within_rlmH:
                        print(f"      Match falls within rlmH gene")
        
        return sites
    
    def filter_sites(self, sites: List[AttSite], require_rlmH: bool = False) -> List[AttSite]:
        """Filter sites based on criteria."""
        filtered_sites = []
        
        for site in sites:
            # For attR and attR2, only include if within rlmH (if gene parser available)
            if site.pattern_name in ['attR', 'attR2'] and self.gene_parser:
                if site.within_rlmH:
                    filtered_sites.append(site)
                else:
                    print(f"      Excluding {site.pattern_name} not in rlmH gene")
            else:
                # Include all other patterns
                filtered_sites.append(site)
        
        return filtered_sites
    
    def write_results(self, sites: List[AttSite], output_file: str):
        """Write results to TSV file."""
        input_file_name = Path(self.fasta_file).stem
        
        with open(output_file, 'a') as f:
            f.write("Input_File\tPattern\tContig\tStart\tEnd\tMatching_Sequence\n")
            
            for site in sites:
                f.write(site.to_tsv_line(input_file_name) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Find att sites in genomic sequences")
    parser.add_argument("-f", "--fna", required=True, help=".fna file for searching")
    parser.add_argument("-g", "--gff", help="gff3 file with the locations of genes")
    parser.add_argument("-o", "--outfile", required=True, help="Output file to save results")
    args = parser.parse_args()
    
    # Validate inputs
    validator = InputValidator()
    
    validator.validate_fasta_file(args.fna)
    
    if args.gff:
        validator.validate_gff_file(args.gff)
    
    # Create the finder
    finder = AttSiteFinder(args.fna, args.gff)
    
    # Find all sites
    all_sites = finder.find_all_sites()
    
    # Filter sites according to rules
    filtered_sites = finder.filter_sites(all_sites)
    
    # Write results
    finder.write_results(filtered_sites, args.outfile)
    
    print(f"\nFound {len(filtered_sites)} att sites total")
    print(f"Results written to {args.outfile}")


if __name__ == '__main__':
    main()