#!/usr/bin/env python

"""Unified SCCmec report merging extraction metadata with typing results.

Produces two TSV reports from the pipeline:

**Summary report** (``sccmec_summary.tsv``)
    Concise one-line-per-genome overview with the most useful columns for
    screening large datasets (~15 columns).

**Full report** (``sccmec_unified_report.tsv``)
    Complete diagnostic detail including att site patterns, coordinates,
    composite boundaries, identity/coverage values, and WGS mec data.

The ``typing_source`` column indicates where typing data came from:
    - ``"sccmec"`` — extracted SCCmec element
    - ``"wgs"``    — whole-genome sequence (extraction failed)
    - ``"-"``      — no typing data available

The ``mec_context`` column indicates where mec resistance genes reside
relative to the extracted element:
    - ``"-"``                    — no mec detected in WGS, or not extracted
    - ``"in_element"``           — all mec genes within the element
    - ``"mec_adjacent"``         — same contig, <50 kb from element
    - ``"mec_chromosomal"``      — same contig, >=50 kb from element
    - ``"mec_different_contig"`` — different contig from the element
"""

import argparse
import csv
import re
import sys

from sccmecextractor.extract_SCCmec import ExtractionReport
from sccmecextractor.sccmec_type_classification import TYPING_HEADER

# Derive extraction columns from the canonical header string
EXTRACTION_HEADER = ExtractionReport.HEADER.split("\t")

# Typing columns excluding Input_File (already in extraction header)
TYPING_EXTRA_COLS = [c for c in TYPING_HEADER if c != "Input_File"]

# WGS mec columns added to the full report for extracted genomes
WGS_MEC_COLS = ["wgs_mec_genes", "wgs_mec_locations"]

# Full diagnostic report header
UNIFIED_HEADER = (
    EXTRACTION_HEADER
    + TYPING_EXTRA_COLS
    + ["typing_source", "element_type", "mec_context"]
    + WGS_MEC_COLS
)

# Concise summary report header
SUMMARY_HEADER = [
    "Input_File",
    "Status",
    "Element_Size_bp",
    "Is_Composite",
    "mec_genes",
    "mec_class_type",
    "ccr_allotypes",
    "ccr_complex_type",
    "SCCmec_Type",
    "typing_source",
    "element_type",
    "mec_context",
    "Failure_Reason",
    "Notes",
]

# Actual mec resistance genes (not IS elements or regulatory genes)
_MEC_GENE_PATTERN = re.compile(r"^(mecA|mecB|mecC|mecA1|mecA2|mecD)\(")

# Adjacency threshold in bp
_ADJACENT_THRESHOLD = 50_000


def read_tsv(filepath, key_column="Input_File"):
    """Read a TSV file into a dict keyed by *key_column*.

    Handles duplicate header lines (common from batch concatenation).
    Returns ``{key_value: {col: val, ...}, ...}``.

    Exits with an error message if *filepath* does not exist.
    """
    try:
        fh = open(filepath, "r", newline="")
    except FileNotFoundError:
        print(f"ERROR: File not found: {filepath}", file=sys.stderr)
        sys.exit(1)

    rows = {}
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        key = row.get(key_column, "")
        # Skip duplicate header lines
        if key == key_column:
            continue
        rows[key] = dict(row)
    fh.close()
    return rows


def normalise_typing_keys(typing_rows):
    """Strip trailing ``_SCCmec`` from typing keys so they match extraction keys.

    Returns a new dict with normalised keys.
    """
    normalised = {}
    for key, val in typing_rows.items():
        if key.endswith("_SCCmec"):
            new_key = key[: -len("_SCCmec")]
        else:
            new_key = key
        # Update the Input_File value inside the row too
        val = dict(val)
        val["Input_File"] = new_key
        normalised[new_key] = val
    return normalised


_EXTRACTED_STATUSES = {
    "extracted", "composite_extracted",
    "fallback_extracted", "composite_fallback_extracted",
}


def _classify_mec_context(ext_row, wgs_row):
    """Determine mec_context for an extracted element.

    Compares WGS mec gene locations against the extracted element boundaries
    to classify where mec genes reside relative to the element.

    Parameters
    ----------
    ext_row : dict
        Extraction report row with Contig, AttR_Start, AttL_End columns.
    wgs_row : dict or None
        WGS typing row with mec_genes and mec_locations columns.

    Returns
    -------
    str
        One of: "-", "in_element", "mec_adjacent", "mec_chromosomal",
        "mec_different_contig".
    """
    if not wgs_row:
        return "-"

    wgs_mec_genes = wgs_row.get("mec_genes", "-").strip()
    if wgs_mec_genes == "-" or wgs_mec_genes == "":
        return "-"

    # Check if any actual mec resistance gene is present (not just IS/regulatory)
    genes = wgs_mec_genes.split(";")
    if not any(_MEC_GENE_PATTERN.match(g.strip()) for g in genes):
        return "-"

    # Get element boundaries
    elem_contig = ext_row.get("Contig", "-").strip()
    attr_start = ext_row.get("AttR_Start", "-").strip()
    attl_end = ext_row.get("AttL_End", "-").strip()

    if elem_contig == "-" or attr_start == "-" or attl_end == "-":
        return "-"

    try:
        elem_start = int(attr_start)
        elem_end = int(attl_end)
    except ValueError:
        return "-"

    if elem_start > elem_end:
        elem_start, elem_end = elem_end, elem_start

    # Parse WGS mec locations and classify each actual mec gene
    wgs_locs = wgs_row.get("mec_locations", "-").strip()
    if wgs_locs == "-" or wgs_locs == "":
        return "-"

    loc_parts = wgs_locs.split(";")
    contexts = set()

    for gene, loc in zip(genes, loc_parts):
        gene = gene.strip()
        loc = loc.strip()

        # Only check actual mec resistance genes
        if not _MEC_GENE_PATTERN.match(gene):
            continue

        if ":" not in loc:
            continue

        mec_contig, coords = loc.rsplit(":", 1)
        try:
            mec_start, mec_end = coords.split("-")
            mec_start, mec_end = int(mec_start), int(mec_end)
        except ValueError:
            continue

        if mec_contig != elem_contig:
            contexts.add("mec_different_contig")
        elif mec_start >= elem_start and mec_end <= elem_end:
            contexts.add("in_element")
        else:
            dist = min(abs(mec_start - elem_end), abs(mec_end - elem_start))
            if dist < _ADJACENT_THRESHOLD:
                contexts.add("mec_adjacent")
            else:
                contexts.add("mec_chromosomal")

    if not contexts:
        return "-"

    # Priority: in_element > mec_adjacent > mec_chromosomal > mec_different_contig
    for label in ("in_element", "mec_adjacent", "mec_chromosomal", "mec_different_contig"):
        if label in contexts:
            return label

    return "-"


def merge_reports(extraction_rows, typing_rows, wgs_typing_rows=None):
    """Full outer join of extraction and typing dicts.

    The ``typing_source`` column is inferred from the extraction status:
        - ``"sccmec"`` — extraction succeeded, typing from extracted element
        - ``"wgs"``    — extraction failed/missing, typing from whole genome
        - ``"-"``      — no typing data available

    The ``mec_context`` column is derived by comparing WGS-level mec detection
    against the extracted element boundaries (only for extracted genomes).

    Parameters
    ----------
    extraction_rows : dict
        Extraction report keyed by Input_File.
    typing_rows : dict
        Typing results keyed by Input_File (element or WGS depending on status).
    wgs_typing_rows : dict, optional
        WGS typing results for extracted genomes (genome-wide mec screen).

    Missing columns are filled with ``"-"``. Returns a list of dicts
    sorted by Input_File.
    """
    if wgs_typing_rows is None:
        wgs_typing_rows = {}

    all_keys = sorted(set(extraction_rows) | set(typing_rows))
    merged = []

    for key in all_keys:
        row = {"Input_File": key}

        ext = extraction_rows.get(key, {})
        for col in EXTRACTION_HEADER:
            if col == "Input_File":
                continue
            row[col] = ext.get(col, "-")

        typ = typing_rows.get(key)
        if typ:
            for col in TYPING_EXTRA_COLS:
                row[col] = typ.get(col, "-")
            # Infer source from extraction status
            status = ext.get("Status", "-")
            row["typing_source"] = (
                "sccmec" if status in _EXTRACTED_STATUSES else "wgs"
            )
        else:
            for col in TYPING_EXTRA_COLS:
                row[col] = "-"
            row["typing_source"] = "-"

        # Classify element type based on extraction + typing
        status = ext.get("Status", "-")
        if status in _EXTRACTED_STATUSES:
            mec = row.get("mec_genes", "-")
            row["element_type"] = "SCCmec" if mec != "-" else "SCC"

            # Determine mec_context and WGS mec data
            wgs = wgs_typing_rows.get(key)
            row["mec_context"] = _classify_mec_context(ext, wgs)
            if wgs:
                row["wgs_mec_genes"] = wgs.get("mec_genes", "-")
                row["wgs_mec_locations"] = wgs.get("mec_locations", "-")
            else:
                row["wgs_mec_genes"] = "-"
                row["wgs_mec_locations"] = "-"
        else:
            row["element_type"] = "-"
            row["mec_context"] = "-"
            row["wgs_mec_genes"] = "-"
            row["wgs_mec_locations"] = "-"

        merged.append(row)

    return merged


def write_unified_report(merged, outfile):
    """Write the full diagnostic report as a TSV file."""
    with open(outfile, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=UNIFIED_HEADER, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(merged)


def write_summary_report(merged, outfile):
    """Write the concise summary report as a TSV file."""
    with open(outfile, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=SUMMARY_HEADER, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(merged)


def main():
    parser = argparse.ArgumentParser(
        description="Merge SCCmec extraction report with typing results"
    )
    parser.add_argument(
        "-e",
        "--extraction-report",
        required=True,
        help="TSV from sccmec-extract --report",
    )
    parser.add_argument(
        "-t",
        "--typing-results",
        required=True,
        help="TSV from sccmec-type (extracted SCCmec sequences, whole genomes, or both)",
    )
    parser.add_argument(
        "-w",
        "--wgs-typing-results",
        default=None,
        help="TSV from WGS mec screen (genome-wide typing for extracted genomes)",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        required=True,
        help="Output unified report TSV",
    )
    parser.add_argument(
        "-s",
        "--summary",
        default=None,
        help="Output summary report TSV (optional)",
    )
    args = parser.parse_args()

    extraction_rows = read_tsv(args.extraction_report)
    typing_rows = read_tsv(args.typing_results)
    typing_rows = normalise_typing_keys(typing_rows)

    wgs_typing_rows = {}
    if args.wgs_typing_results:
        wgs_typing_rows = read_tsv(args.wgs_typing_results)

    merged = merge_reports(extraction_rows, typing_rows, wgs_typing_rows)
    write_unified_report(merged, args.outfile)

    if args.summary:
        write_summary_report(merged, args.summary)

    # Summary counts
    sources = [r.get("typing_source", "-") for r in merged]
    n_sccmec = sources.count("sccmec")
    n_wgs = sources.count("wgs")
    n_none = sources.count("-")
    contexts = [r.get("mec_context", "-") for r in merged]
    n_adjacent = contexts.count("mec_adjacent")
    n_chromosomal = contexts.count("mec_chromosomal")
    n_diff_contig = contexts.count("mec_different_contig")
    print(
        f"Unified report: {len(merged)} entries written to {args.outfile}\n"
        f"  Typing source: {n_sccmec} sccmec, {n_wgs} wgs, {n_none} no typing\n"
        f"  Mec context: {n_adjacent} adjacent, {n_chromosomal} chromosomal, "
        f"{n_diff_contig} different contig"
    )


if __name__ == "__main__":
    main()
