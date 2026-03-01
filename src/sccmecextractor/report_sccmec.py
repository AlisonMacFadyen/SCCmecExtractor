#!/usr/bin/env python

"""Unified SCCmec report merging extraction metadata with typing results.

Produces a single TSV combining columns from sccmec-extract (ExtractionReport)
and sccmec-type (typing results) via a full outer join on Input_File.

The ``typing_source`` column is inferred automatically from the extraction
status: ``"sccmec"`` when the genome was successfully extracted (typing came
from an extracted SCCmec element), ``"wgs"`` when extraction failed (typing
came from a whole-genome sequence), or ``"-"`` when no typing data exists.
"""

import argparse
import csv
import sys

from sccmecextractor.extract_SCCmec import ExtractionReport
from sccmecextractor.type_sccmec import TYPING_HEADER

# Derive extraction columns from the canonical header string
EXTRACTION_HEADER = ExtractionReport.HEADER.split("\t")

# Typing columns excluding Input_File (already in extraction header)
TYPING_EXTRA_COLS = [c for c in TYPING_HEADER if c != "Input_File"]

# Unified header (typing_source tracks origin: "sccmec", "wgs", or "-")
UNIFIED_HEADER = EXTRACTION_HEADER + TYPING_EXTRA_COLS + ["typing_source"]


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


_EXTRACTED_STATUSES = {"extracted", "composite_extracted"}


def merge_reports(extraction_rows, typing_rows):
    """Full outer join of extraction and typing dicts.

    The ``typing_source`` column is inferred from the extraction status:
        - ``"sccmec"`` — extraction succeeded, typing from extracted element
        - ``"wgs"``    — extraction failed/missing, typing from whole genome
        - ``"-"``      — no typing data available

    Missing columns are filled with ``"-"``. Returns a list of dicts
    sorted by Input_File.
    """
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

        merged.append(row)

    return merged


def write_unified_report(merged, outfile):
    """Write merged rows as a TSV file."""
    with open(outfile, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=UNIFIED_HEADER, delimiter="\t", extrasaction="ignore"
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
        "-o",
        "--outfile",
        required=True,
        help="Output unified report TSV",
    )
    args = parser.parse_args()

    extraction_rows = read_tsv(args.extraction_report)
    typing_rows = read_tsv(args.typing_results)
    typing_rows = normalise_typing_keys(typing_rows)

    merged = merge_reports(extraction_rows, typing_rows)
    write_unified_report(merged, args.outfile)

    # Summary counts
    sources = [r.get("typing_source", "-") for r in merged]
    n_sccmec = sources.count("sccmec")
    n_wgs = sources.count("wgs")
    n_none = sources.count("-")
    print(
        f"Unified report: {len(merged)} entries written to {args.outfile}\n"
        f"  Typing source: {n_sccmec} sccmec, {n_wgs} wgs, {n_none} no typing"
    )


if __name__ == "__main__":
    main()
