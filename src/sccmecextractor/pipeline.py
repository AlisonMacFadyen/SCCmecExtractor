#!/usr/bin/env python

"""Master pipeline orchestrating all four SCCmecExtractor stages.

Runs locate-att → extract → type → report in a single command.
"""

import argparse
import csv
import glob
import os
import sys
import tempfile
import threading

from pathlib import Path
from typing import Dict, List, Optional

from sccmecextractor.locate_att_sites import AttSiteFinder
from sccmecextractor.extract_SCCmec import SCCmecExtractor, ExtractionReport, AmbiguousHitReport, GenomeSequences
from sccmecextractor.type_sccmec import SCCmecTyper, TYPING_HEADER
from sccmecextractor.report_sccmec import (
    read_tsv,
    normalise_typing_keys,
    merge_reports,
    write_unified_report,
)


def resolve_gff(fasta_stem: str, gff_files: Optional[Dict[str, str]] = None,
                gff_dir: Optional[str] = None) -> Optional[str]:
    """Find a matching GFF3 file for a FASTA stem name.

    Parameters
    ----------
    fasta_stem : str
        Stem name of the FASTA file (e.g. ``Path("genome.fna").stem``).
    gff_files : dict, optional
        Pre-built mapping of ``{stem: path}`` from explicit ``--gff`` files.
    gff_dir : str, optional
        Directory to search for ``{fasta_stem}.gff3``.

    Returns
    -------
    str or None
        Path to the matching GFF3 file, or None if not found.
    """
    if gff_files is not None:
        return gff_files.get(fasta_stem)

    if gff_dir is not None:
        candidate = os.path.join(gff_dir, f"{fasta_stem}.gff3")
        if os.path.isfile(candidate):
            return candidate

    return None


def _write_typing_row(result: dict, output_file: str, write_header: bool):
    """Append one typing result dict as a TSV row."""
    with open(output_file, "a", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=TYPING_HEADER, delimiter="\t", extrasaction="ignore",
        )
        if write_header:
            writer.writeheader()
        writer.writerow(result)


def _process_genome(
    fasta_path: str,
    att_dir: str,
    sccmec_dir: str,
    extraction_report_file: str,
    ambiguous_report_file: str,
    typer: SCCmecTyper,
    gff_files: Optional[Dict[str, str]],
    gff_dir: Optional[str],
    blast_rlmh: bool,
    rlmh_ref: Optional[str],
    composite: bool,
    index: int,
    total: int,
    print_lock: Optional[threading.Lock] = None,
) -> dict:
    """Process a single genome through stages 1-3.

    Returns a result dict with keys:
        stem, status, typing_result, success
    where status is one of "extracted", "failed", "error_locate", "error_extract".
    """
    stem = Path(fasta_path).stem
    progress = f"[{index}/{total}] {stem}:"
    result = {
        "stem": stem,
        "status": "failed",
        "typing_result": None,
        "success": False,
        "extracted": False,
        "typed_sccmec": False,
        "typed_wgs": False,
    }

    def _print(*args, **kwargs):
        if print_lock is not None:
            with print_lock:
                print(*args, **kwargs)
        else:
            print(*args, **kwargs)

    # --- Resolve GFF for this genome ---
    gff_path = resolve_gff(stem, gff_files=gff_files, gff_dir=gff_dir)
    use_blast = blast_rlmh or (gff_path is None)

    # --- Parse genome once, share between stages ---
    genome = GenomeSequences(fasta_path)
    string_sequences = {cid: str(rec.seq) for cid, rec in genome.sequences.items()}

    # --- Create shared BLAST DB when using BLAST mode ---
    genome_db_prefix = None
    tmp_db_dir = None
    if use_blast:
        from sccmecextractor.blast_utils import BlastRunner
        tmp_db_dir = tempfile.mkdtemp(prefix="sccmec_pipeline_")
        genome_db_prefix = os.path.join(tmp_db_dir, "genome_db")
        BlastRunner().create_db(fasta_path, genome_db_prefix)

    try:
        # --- Stage 1: Locate att sites ---
        _print(f"{progress} locating att sites...", end="", file=sys.stderr, flush=True)
        att_output = os.path.join(att_dir, f"{stem}_att_sites.tsv")

        try:
            finder = AttSiteFinder(
                fasta_path, gff3_file=gff_path,
                blast_rlmh=use_blast, rlmh_ref=rlmh_ref,
                sequences=string_sequences,
                genome_db_prefix=genome_db_prefix,
            )
            all_sites = finder.find_all_sites()
            filtered_sites = finder.filter_sites(all_sites)
            finder.write_results(filtered_sites, att_output)
        except Exception as e:
            _print(f" ERROR (locate): {e}", file=sys.stderr)
            result["status"] = "error_locate"
            return result

        # --- Stage 2: Extract SCCmec ---
        _print(" extracting...", end="", file=sys.stderr, flush=True)

        # Pass pre-computed rlmH positions to avoid redundant BLAST
        rlmh_positions = getattr(finder.gene_parser, 'rlmH_genes', None)

        try:
            extractor = SCCmecExtractor(
                fasta_path, gff3_file=gff_path, tsv_file=att_output,
                composite=composite, blast_rlmh=use_blast, rlmh_ref=rlmh_ref,
                rlmh_positions=rlmh_positions,
                genome_sequences=genome,
                genome_db_prefix=genome_db_prefix,
            )
            success = extractor.extract_sccmec(
                sccmec_dir, report_file=extraction_report_file,
                ambiguous_report_file=ambiguous_report_file,
            )
        except Exception as e:
            _print(f" ERROR (extract): {e}", file=sys.stderr)
            result["status"] = "error_extract"
            return result
    finally:
        # Clean up shared BLAST DB
        if genome_db_prefix is not None:
            from sccmecextractor.blast_utils import BlastRunner as _BR
            _BR.cleanup_db(genome_db_prefix)
            try:
                os.rmdir(tmp_db_dir)
            except OSError:
                pass

    # --- Stage 3: Type ---
    sccmec_fasta = os.path.join(sccmec_dir, f"{stem}_SCCmec.fasta")
    if success and os.path.isfile(sccmec_fasta):
        _print(" typing (sccmec)...", end="", file=sys.stderr, flush=True)
        try:
            typing_result = typer.type_file(sccmec_fasta)
            result["typing_result"] = typing_result
            result["typed_sccmec"] = True
        except Exception as e:
            _print(f" typing ERROR: {e}", end="", file=sys.stderr)
        result["extracted"] = True
        result["success"] = True
    else:
        _print(" FAILED, typing (wgs)...", end="", file=sys.stderr, flush=True)
        try:
            typing_result = typer.type_file(fasta_path)
            result["typing_result"] = typing_result
            result["typed_wgs"] = True
        except Exception as e:
            _print(f" typing ERROR: {e}", end="", file=sys.stderr)

    _print(" done", file=sys.stderr)
    return result


def run_pipeline(
    fasta_files: List[str],
    outdir: str,
    gff_files: Optional[Dict[str, str]] = None,
    gff_dir: Optional[str] = None,
    blast_rlmh: bool = False,
    rlmh_ref: Optional[str] = None,
    composite: bool = False,
    threads: int = 1,
) -> dict:
    """Run the full SCCmecExtractor pipeline on one or more genomes.

    Parameters
    ----------
    fasta_files : list of str
        Paths to input FASTA/FNA files.
    outdir : str
        Root output directory.
    gff_files : dict, optional
        ``{stem: path}`` mapping for explicit GFF files.
    gff_dir : str, optional
        Directory of GFF3 files matched by stem.
    blast_rlmh : bool
        Force BLAST-based rlmH detection.
    rlmh_ref : str, optional
        Custom rlmH reference FASTA for BLAST.
    composite : bool
        Extract to outermost boundary for composite elements.
    threads : int
        Number of parallel threads (default 1 = sequential).

    Returns
    -------
    dict
        Summary with keys: total, extracted, failed, typed_sccmec, typed_wgs.
    """
    # Create output subdirectories
    att_dir = os.path.join(outdir, "att_sites")
    sccmec_dir = os.path.join(outdir, "sccmec")
    typing_dir = os.path.join(outdir, "typing")
    for d in (att_dir, sccmec_dir, typing_dir):
        os.makedirs(d, exist_ok=True)

    extraction_report_file = os.path.join(outdir, "extraction_report.tsv")
    ambiguous_report_file = os.path.join(outdir, "ambiguous_att_sites.tsv")
    typing_results_file = os.path.join(outdir, "typing_results.tsv")

    # Instantiate one typer (reuses BLAST runner across all genomes)
    typer = SCCmecTyper()

    total = len(fasta_files)

    # Pre-write report headers to avoid race condition in parallel mode
    if threads > 1 and not os.path.isfile(extraction_report_file):
        with open(extraction_report_file, 'w') as f:
            f.write(ExtractionReport.HEADER + "\n")
    if threads > 1 and not os.path.isfile(ambiguous_report_file):
        with open(ambiguous_report_file, 'w') as f:
            f.write(AmbiguousHitReport.HEADER + "\n")

    # Common kwargs for _process_genome
    common_kwargs = dict(
        att_dir=att_dir,
        sccmec_dir=sccmec_dir,
        extraction_report_file=extraction_report_file,
        ambiguous_report_file=ambiguous_report_file,
        typer=typer,
        gff_files=gff_files,
        gff_dir=gff_dir,
        blast_rlmh=blast_rlmh,
        rlmh_ref=rlmh_ref,
        composite=composite,
        total=total,
    )

    if threads <= 1:
        # Sequential processing
        results = []
        for i, fasta_path in enumerate(fasta_files, 1):
            result = _process_genome(
                fasta_path, index=i, **common_kwargs,
            )
            results.append(result)
    else:
        from concurrent.futures import ThreadPoolExecutor, as_completed

        print_lock = threading.Lock()
        results = [None] * total

        with ThreadPoolExecutor(max_workers=threads) as pool:
            futures = {}
            for i, fasta_path in enumerate(fasta_files):
                future = pool.submit(
                    _process_genome,
                    fasta_path,
                    index=i + 1,
                    print_lock=print_lock,
                    **common_kwargs,
                )
                futures[future] = i

            for future in as_completed(futures):
                idx = futures[future]
                try:
                    results[idx] = future.result()
                except Exception as e:
                    stem = Path(fasta_files[idx]).stem
                    print(f"ERROR: {stem}: {e}", file=sys.stderr)
                    results[idx] = {
                        "stem": stem,
                        "status": "error",
                        "typing_result": None,
                        "success": False,
                        "extracted": False,
                        "typed_sccmec": False,
                        "typed_wgs": False,
                    }

    # Write typing results in input order (deterministic output)
    typing_header_written = False
    for result in results:
        if result and result.get("typing_result"):
            _write_typing_row(
                result["typing_result"],
                typing_results_file,
                not typing_header_written,
            )
            typing_header_written = True

    # Tally results
    extracted_count = sum(1 for r in results if r and r.get("extracted"))
    failed_count = total - extracted_count
    typed_sccmec = sum(1 for r in results if r and r.get("typed_sccmec"))
    typed_wgs = sum(1 for r in results if r and r.get("typed_wgs"))

    # --- Stage 4: Unified report ---
    unified_report_file = os.path.join(outdir, "sccmec_unified_report.tsv")

    if os.path.isfile(extraction_report_file) and os.path.isfile(typing_results_file):
        extraction_rows = read_tsv(extraction_report_file)
        typing_rows = read_tsv(typing_results_file)
        typing_rows = normalise_typing_keys(typing_rows)
        merged = merge_reports(extraction_rows, typing_rows)
        write_unified_report(merged, unified_report_file)
    elif os.path.isfile(extraction_report_file):
        # Typing may have failed for all genomes; still produce a report
        extraction_rows = read_tsv(extraction_report_file)
        merged = merge_reports(extraction_rows, {})
        write_unified_report(merged, unified_report_file)

    # Summary
    summary = {
        "total": total,
        "extracted": extracted_count,
        "failed": failed_count,
        "typed_sccmec": typed_sccmec,
        "typed_wgs": typed_wgs,
    }

    print(
        f"\nPipeline complete: {total} genomes processed\n"
        f"  Extracted: {extracted_count} ({extracted_count/total*100:.1f}%)\n"
        f"  Failed: {failed_count} ({failed_count/total*100:.1f}%)\n"
        f"  Typed (SCCmec): {typed_sccmec}, Typed (WGS): {typed_wgs}\n"
        f"Report: {unified_report_file}",
        file=sys.stderr,
    )

    return summary


def main():
    parser = argparse.ArgumentParser(
        description=(
            "SCCmecExtractor pipeline: locate att sites, extract SCCmec, "
            "type by gene content, and generate a unified report"
        ),
    )
    fna_group = parser.add_mutually_exclusive_group(required=True)
    fna_group.add_argument(
        "-f", "--fna", nargs="+",
        help="One or more FASTA/FNA genome files",
    )
    fna_group.add_argument(
        "--fna-dir",
        help="Directory of FASTA/FNA genome files (all .fna files will be used)",
    )
    gff_group = parser.add_mutually_exclusive_group()
    gff_group.add_argument(
        "-g", "--gff", nargs="+",
        help="One or more GFF3 annotation files (matched to FASTA by stem name)",
    )
    gff_group.add_argument(
        "--gff-dir",
        help="Directory of GFF3 files (matched to FASTA by stem name)",
    )
    parser.add_argument(
        "--blast-rlmh", action="store_true", default=False,
        help="Use BLAST for rlmH detection (auto-enabled when no GFF provided)",
    )
    parser.add_argument(
        "--rlmh-ref",
        help="Custom rlmH reference FASTA for BLAST detection",
    )
    parser.add_argument(
        "--composite", action="store_true",
        help="Extract to outermost boundary for composite elements",
    )
    parser.add_argument(
        "-o", "--outdir", required=True,
        help="Output directory for all results",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1,
        help="Number of parallel threads (default: 1 = sequential)",
    )
    args = parser.parse_args()

    # Resolve FASTA files from --fna or --fna-dir
    if args.fna:
        fasta_files = args.fna
    else:
        if not os.path.isdir(args.fna_dir):
            print(f"ERROR: FNA directory not found: {args.fna_dir}", file=sys.stderr)
            sys.exit(1)
        fasta_files = sorted(glob.glob(os.path.join(args.fna_dir, "*.fna")))
        if not fasta_files:
            print(f"ERROR: No .fna files found in: {args.fna_dir}", file=sys.stderr)
            sys.exit(1)
        print(f"Found {len(fasta_files)} .fna files in {args.fna_dir}", file=sys.stderr)

    # Validate FASTA files exist
    for fna in fasta_files:
        if not os.path.isfile(fna):
            print(f"ERROR: FASTA file not found: {fna}", file=sys.stderr)
            sys.exit(1)

    # Build GFF lookup
    gff_files = None
    gff_dir = None
    if args.gff:
        gff_files = {Path(g).stem: g for g in args.gff}
    elif args.gff_dir:
        if not os.path.isdir(args.gff_dir):
            print(f"ERROR: GFF directory not found: {args.gff_dir}", file=sys.stderr)
            sys.exit(1)
        gff_dir = args.gff_dir

    # Auto-enable BLAST when no GFF source provided
    blast_rlmh = args.blast_rlmh
    if not args.gff and not args.gff_dir and not blast_rlmh:
        try:
            from sccmecextractor.blast_utils import BlastRunner
            BlastRunner._check_blast_installed()
            blast_rlmh = True
            print(
                "No GFF provided; auto-enabling BLAST-based rlmH detection",
                file=sys.stderr,
            )
        except Exception:
            print(
                "ERROR: No GFF provided and BLAST+ not available. "
                "Either provide --gff / --gff-dir or install BLAST+ for --blast-rlmh.",
                file=sys.stderr,
            )
            sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    run_pipeline(
        fasta_files=fasta_files,
        outdir=args.outdir,
        gff_files=gff_files,
        gff_dir=gff_dir,
        blast_rlmh=blast_rlmh,
        rlmh_ref=args.rlmh_ref,
        composite=args.composite,
        threads=args.threads,
    )


if __name__ == "__main__":
    main()
