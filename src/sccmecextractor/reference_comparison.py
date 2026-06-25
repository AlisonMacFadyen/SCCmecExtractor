#!/usr/bin/env python

"""SCCmec Reference Comparison — hybrid element detection.

BLASTs extracted SCC elements against 32 bundled SCCmec type references
(Types I–XV, including subtypes) to determine whether an element is a
canonical single-type match or a hybrid composed of regions from different
SCCmec types.

When a typing report is available (from sccmec-type or the pipeline), the
tool distinguishes true hybrids (single ccr complex, multi-type similarity)
from structural composites (multiple ccr complexes).

Outputs:
    - Always: hybrid_summary.tsv (one row per element)
    - With --detailed: per-element detail files in hybrid_detail/ subdirectory
"""

import argparse
import csv
import logging
import os
import sys
import tempfile

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO

from sccmecextractor.blast_utils import (
    BlastRunner,
    get_default_ref,
    parse_blast_output,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Reverse lookup: subtype name -> parent SCCmec type
_SUBTYPE_TO_TYPE = {
    "Ia": "I", "Ib": "I",
    "IIa": "II", "IIb": "II", "IIc": "II", "IId": "II", "IIe": "II",
    "III": "III",
    "IVa": "IV", "IVb": "IV", "IVc": "IV", "IVd": "IV",
    "IVg": "IV", "IVi": "IV", "IVj": "IV", "IVk": "IV",
    "IVl": "IV", "IVm": "IV", "IVn": "IV",
    "Va": "V", "Vb": "V", "Vc": "V",
    "VI": "VI", "VII": "VII", "VIII": "VIII", "IX": "IX",
    "X": "X", "XI": "XI", "XII": "XII", "XIII": "XIII",
    "XIV": "XIV", "XV": "XV",
}

# BLAST filtering thresholds
MIN_HIT_LENGTH = 1000   # bp — ignore short spurious alignments
MIN_HIT_PIDENT = 80.0   # % — ignore low-identity conserved-gene cross-matches

# Hybrid detection thresholds
NOVEL_BP_THRESHOLD = 1000      # bp — minimum novel territory for a secondary type
NOVEL_FRAC_THRESHOLD = 0.05    # fraction of element length — alternative minimum

# Summary report columns
HYBRID_SUMMARY_HEADER = [
    "element_id",
    "hybrid_best_match",
    "hybrid_best_subtype",
    "hybrid_best_coverage",
    "hybrid_best_identity",
    "hybrid_secondary_match",
    "hybrid_secondary_subtype",
    "hybrid_secondary_coverage",
    "hybrid_secondary_identity",
    "hybrid_call",
    "hybrid_components",
]

# Detailed evidence report columns
HYBRID_DETAIL_HEADER = [
    "element_id",
    "parent_type",
    "best_subtype",
    "ref_coverage",
    "element_coverage",
    "weighted_pident",
    "covered_bp",
    "element_covered_bp",
    "ref_length",
    "element_length",
    "num_hsps",
    "novel_bp",
]

# Reference FASTA bundled in data/
_REF_FILENAME = "sccmec_references.fasta"


# ---------------------------------------------------------------------------
# Typing report reader
# ---------------------------------------------------------------------------

def read_typing_report(filepath: str) -> Dict[str, dict]:
    """Read a sccmec_summary.tsv and return per-element typing context.

    Keyed by Input_File (accession).  For each element, extracts:
        - ccr_complex_type: e.g. "2", "2;5", "-"
        - Is_Composite: "True"/"False"
        - SCCmec_Type: e.g. "II", "IV (composite)", "not_typeable"
        - mec_genes: e.g. "mecA(full);IS431(full)", "-"

    Only includes rows where Status contains "extracted" (filters out
    failed genomes).
    """
    typing_context = {}
    with open(filepath, "r") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            status = row.get("Status", "")
            if "extracted" not in status:
                continue
            accession = row.get("Input_File", "")
            if not accession:
                continue

            ccr_raw = row.get("ccr_complex_type", "-")
            # Count distinct ccr complexes (e.g. "2;5" = 2 complexes)
            ccr_types = [c.strip() for c in ccr_raw.split(";")
                         if c.strip() and c.strip() != "-"]

            typing_context[accession] = {
                "ccr_complex_type": ccr_raw,
                "n_ccr_complexes": len(set(ccr_types)),
                "Is_Composite": row.get("Is_Composite", "False"),
                "SCCmec_Type": row.get("SCCmec_Type", "-"),
                "mec_genes": row.get("mec_genes", "-"),
            }

    return typing_context


def _extract_accession(element_id: str) -> str:
    """Extract accession from element filename stem.

    Handles patterns like 'GCF_000009585_SCCmec' or 'GCF_000009585_SCC'
    by stripping the _SCCmec/_SCC suffix.
    """
    for suffix in ("_SCCmec", "_SCC"):
        if element_id.endswith(suffix):
            return element_id[: -len(suffix)]
    return element_id


def classify_with_typing(
    hybrid_call: str,
    n_components: int,
    typing_info: Optional[dict],
) -> str:
    """Refine hybrid_call using typing context.

    Returns one of:
        - "canonical": single type, single ccr complex
        - "hybrid": multi-type match, single ccr complex
        - "multi_ccr": multi-type match, multiple ccr complexes
        - "multi_ccr_canonical": single type match, multiple ccr complexes
        - "multi_type": multi-type match, no typing context available
        - "no_match": no BLAST hits
    """
    if hybrid_call == "no_match":
        return "no_match"

    multi_type = n_components > 1
    has_typing = typing_info is not None

    if not has_typing:
        return "multi_type" if multi_type else "canonical"

    multi_ccr = typing_info["n_ccr_complexes"] > 1

    if multi_type and multi_ccr:
        return "multi_ccr"
    if multi_type and not multi_ccr:
        return "hybrid"
    if not multi_type and multi_ccr:
        return "multi_ccr_canonical"
    return "canonical"


# ---------------------------------------------------------------------------
# HybridTyper
# ---------------------------------------------------------------------------

class HybridTyper:
    """BLAST-based SCC element hybrid detection.

    BLASTs each extracted element against the 32 SCCmec type references,
    groups subtypes into parent types, and determines whether the element
    is canonical (single type) or a hybrid of multiple types.
    """

    def __init__(self):
        self.runner = BlastRunner()

        with get_default_ref(_REF_FILENAME) as ref:
            self.ref_lengths = {}
            for record in SeqIO.parse(str(ref), "fasta"):
                self.ref_lengths[record.id] = len(record.seq)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def type_file(
        self,
        input_fasta: str,
        db_prefix: str = None,
        typing_info: Optional[dict] = None,
    ) -> Tuple[dict, list]:
        """Analyse a single extracted SCC element for hybrid composition.

        Parameters
        ----------
        input_fasta : str
            Path to the extracted element FASTA.
        db_prefix : str, optional
            Pre-built BLAST database prefix.  When provided the existing
            database is reused, avoiding a redundant ``makeblastdb`` call.
        typing_info : dict, optional
            Per-element typing context from ``read_typing_report``.  When
            provided, ``hybrid_call`` distinguishes composite from hybrid.

        Returns
        -------
        summary : dict
            Keys matching ``HYBRID_SUMMARY_HEADER`` — for the unified report.
        detail_rows : list[dict]
            One dict per parent type with BLAST evidence — for the detailed
            hybrid evidence file.
        """
        input_name = Path(input_fasta).stem

        # Get element length
        element_length = 0
        for record in SeqIO.parse(str(input_fasta), "fasta"):
            element_length = len(record.seq)

        if db_prefix is not None:
            hits = self._blast_ref(db_prefix)
        else:
            tmp_dir = tempfile.mkdtemp(prefix=f"{input_name}_hybrid_")
            db_prefix_local = os.path.join(tmp_dir, f"{input_name}_db")
            try:
                self.runner.create_db(input_fasta, db_prefix_local)
                hits = self._blast_ref(db_prefix_local)
            finally:
                self.runner.cleanup_db(db_prefix_local)
                try:
                    os.rmdir(tmp_dir)
                except OSError:
                    pass

        per_subtype = self._calculate_coverage(hits, element_length)
        type_profiles = self._merge_subtypes(per_subtype)
        summary, detail_rows = self._classify_hybrid(
            input_name, element_length, type_profiles, typing_info
        )

        return summary, detail_rows

    def type_batch(
        self,
        input_fastas: List[str],
        typing_context: Optional[Dict[str, dict]] = None,
    ) -> List[Tuple[dict, list]]:
        """Batch-analyse multiple extracted elements with a single BLAST.

        Concatenates all element FASTAs into one BLAST database, runs one
        BLAST search, then splits results by subject (element) ID.

        Parameters
        ----------
        input_fastas : list of str
            Paths to extracted element FASTA files.
        typing_context : dict, optional
            Keyed by accession, from ``read_typing_report``.

        Returns
        -------
        list of (summary, detail_rows) tuples, one per input element,
        in the same order as input_fastas.
        """
        if not input_fastas:
            return []

        typing_context = typing_context or {}

        # Read element lengths and build combined FASTA
        element_lengths = {}
        tmp_dir = tempfile.mkdtemp(prefix="hybrid_batch_")
        combined_fasta = os.path.join(tmp_dir, "combined_elements.fasta")
        try:
            with open(combined_fasta, "w") as out_fh:
                for fasta_path in input_fastas:
                    for record in SeqIO.parse(str(fasta_path), "fasta"):
                        element_lengths[record.id] = len(record.seq)
                        out_fh.write(f">{record.id}\n{record.seq}\n")

            # Single BLAST: 32 references (query) vs all elements (DB)
            db_prefix = os.path.join(tmp_dir, "combined_db")
            self.runner.create_db(combined_fasta, db_prefix)
            all_hits = self._blast_ref(db_prefix)
            self.runner.cleanup_db(db_prefix)
        finally:
            try:
                os.remove(combined_fasta)
                os.rmdir(tmp_dir)
            except OSError:
                pass

        # Group hits by element (subject sequence ID)
        hits_by_element = defaultdict(list)
        for hit in all_hits:
            hits_by_element[hit.sseqid].append(hit)

        # Map FASTA path -> record ID (element name in the combined DB)
        path_to_record_id = {}
        for fasta_path in input_fastas:
            for record in SeqIO.parse(str(fasta_path), "fasta"):
                path_to_record_id[fasta_path] = record.id
                break  # first record only

        # Process each element
        results = []
        for fasta_path in input_fastas:
            input_name = Path(fasta_path).stem
            record_id = path_to_record_id.get(fasta_path, input_name)
            element_length = element_lengths.get(record_id, 0)
            element_hits = hits_by_element.get(record_id, [])

            per_subtype = self._calculate_coverage(
                element_hits, element_length
            )
            type_profiles = self._merge_subtypes(per_subtype)

            # Look up typing context
            accession = _extract_accession(input_name)
            typing_info = typing_context.get(accession)

            summary, detail_rows = self._classify_hybrid(
                input_name, element_length, type_profiles, typing_info
            )
            results.append((summary, detail_rows))

        return results

    # ------------------------------------------------------------------
    # BLAST
    # ------------------------------------------------------------------

    def _blast_ref(self, db_prefix: str) -> list:
        """BLAST the bundled references (query) against the element DB."""
        with get_default_ref(_REF_FILENAME) as ref:
            results_file = self.runner.run_blastn(str(ref), db_prefix)
            hits = parse_blast_output(results_file)
            self.runner.cleanup_file(results_file)
            return hits

    # ------------------------------------------------------------------
    # Coverage calculation (per subtype)
    # ------------------------------------------------------------------

    @staticmethod
    def _merge_intervals(intervals):
        """Merge overlapping intervals into a sorted list of [start, end].

        Returns an empty list when no intervals are provided.
        """
        if not intervals:
            return []

        sorted_intervals = sorted(intervals)
        merged = [list(sorted_intervals[0])]

        for start, end in sorted_intervals[1:]:
            if start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])

        return merged

    def _calculate_coverage(self, blast_hits, element_length):
        """Per-subtype coverage and interval data.

        Filters hits by MIN_HIT_LENGTH and MIN_HIT_PIDENT, then for each
        reference subtype calculates:
            - reference coverage (what % of the reference is found)
            - element intervals (where on the element the hits land)
            - weighted percent identity
        """
        hits_by_ref = defaultdict(list)
        for hit in blast_hits:
            if hit.length >= MIN_HIT_LENGTH and hit.pident >= MIN_HIT_PIDENT:
                hits_by_ref[hit.qseqid].append(hit)

        results = {}
        for ref_name, ref_hits in hits_by_ref.items():
            ref_intervals = [(hit.qstart, hit.qend) for hit in ref_hits]
            # Normalise subject coords for minus-strand hits
            ele_intervals = [
                (min(hit.sstart, hit.send), max(hit.sstart, hit.send))
                for hit in ref_hits
            ]

            ref_merged = self._merge_intervals(ref_intervals)
            ele_merged = self._merge_intervals(ele_intervals)

            ref_len = self.ref_lengths[ref_name]
            ref_covered_bp = sum(e - s for s, e in ref_merged)
            ele_covered_bp = sum(e - s for s, e in ele_merged)
            coverage = ref_covered_bp / ref_len if ref_len else 0
            weighted_pident = (
                sum(h.pident * h.length for h in ref_hits)
                / sum(h.length for h in ref_hits)
            )

            results[ref_name] = {
                "coverage": round(coverage * 100, 1),
                "ref_length": ref_len,
                "ref_covered_bp": ref_covered_bp,
                "ele_covered_bp": ele_covered_bp,
                "num_hsps": len(ref_hits),
                "weighted_pident": round(weighted_pident, 2),
                "ref_merged_intervals": ref_merged,
                "ele_merged_intervals": ele_merged,
                "element_intervals": ele_intervals,
                "raw_hits": ref_hits,
            }

        return results

    # ------------------------------------------------------------------
    # Merge subtypes into parent types
    # ------------------------------------------------------------------

    def _merge_subtypes(self, per_subtype: dict) -> dict:
        """Group subtype results by parent SCCmec type.

        For each parent type, pools element intervals from all subtypes,
        merges them into a single footprint, and computes aggregate stats.
        """
        grouped = defaultdict(list)
        for ref_name, stats in per_subtype.items():
            subtype = ref_name.split(" ")[0]
            parent_type = _SUBTYPE_TO_TYPE.get(subtype)
            if parent_type:
                grouped[parent_type].append((subtype, stats))

        type_profiles = {}
        for type_name, subtype_hits in grouped.items():
            # Pool all element intervals from every subtype
            all_ele_intervals = []
            all_raw_hits = []
            for _subtype, stats in subtype_hits:
                all_ele_intervals.extend(stats["element_intervals"])
                all_raw_hits.extend(stats["raw_hits"])

            element_footprint = self._merge_intervals(all_ele_intervals)
            element_covered_bp = sum(e - s for s, e in element_footprint)

            # Weighted pident across all hits for this parent type
            total_len = sum(h.length for h in all_raw_hits)
            weighted_pident = (
                sum(h.pident * h.length for h in all_raw_hits) / total_len
                if total_len > 0
                else 0.0
            )

            # Best subtype: highest coverage, then highest identity as tiebreak
            best_subtype, best_stats = max(
                subtype_hits,
                key=lambda x: (x[1]["coverage"], x[1]["weighted_pident"]),
            )

            type_profiles[type_name] = {
                "element_footprint": element_footprint,
                "element_covered_bp": element_covered_bp,
                "best_subtype": best_subtype,
                "best_coverage": best_stats["coverage"],
                "best_ref_length": best_stats["ref_length"],
                "weighted_pident": round(weighted_pident, 2),
                "num_hsps": sum(s["num_hsps"] for _, s in subtype_hits),
            }

        return type_profiles

    # ------------------------------------------------------------------
    # Novel coverage (interval subtraction)
    # ------------------------------------------------------------------

    @staticmethod
    def _calculate_novel_coverage(primary_intervals, secondary_intervals):
        """Calculate regions of secondary that fall outside primary.

        Parameters
        ----------
        primary_intervals : list of [start, end]
            Merged intervals already claimed (the primary type's footprint).
        secondary_intervals : list of [start, end]
            Merged intervals of the candidate secondary type.

        Returns
        -------
        novel : list of [start, end]
            Regions unique to the secondary type.
        novel_bp : int
            Total bp of novel territory.
        """
        if not secondary_intervals or not primary_intervals:
            return secondary_intervals or [], sum(
                e - s for s, e in (secondary_intervals or [])
            )

        novel = []
        for s_start, s_end in secondary_intervals:
            remaining = [[s_start, s_end]]
            for p_start, p_end in primary_intervals:
                new_remaining = []
                for r_start, r_end in remaining:
                    if r_end <= p_start or r_start >= p_end:
                        # No overlap — keep as is
                        new_remaining.append([r_start, r_end])
                    else:
                        # Overlap — keep parts outside primary
                        if r_start < p_start:
                            new_remaining.append([r_start, p_start])
                        if r_end > p_end:
                            new_remaining.append([p_end, r_end])
                remaining = new_remaining
            novel.extend(remaining)

        novel_bp = sum(e - s for s, e in novel)
        return novel, novel_bp

    # ------------------------------------------------------------------
    # Hybrid classification
    # ------------------------------------------------------------------

    def _classify_hybrid(
        self,
        input_name: str,
        element_length: int,
        type_profiles: dict,
        typing_info: Optional[dict] = None,
    ) -> Tuple[dict, list]:
        """Determine whether the element is canonical, hybrid, or composite.

        Ranks parent types by a combined score (coverage fraction x identity),
        then iteratively checks whether secondary types contribute novel
        territory on the element.  When typing_info is provided, uses ccr
        complex count to distinguish hybrid from composite.

        Returns (summary_dict, detail_rows).
        """
        # -- No hits at all --
        if not type_profiles:
            summary = {col: "-" for col in HYBRID_SUMMARY_HEADER}
            summary["element_id"] = input_name
            summary["hybrid_call"] = "no_match"
            return summary, []

        # Rank by combined score: element coverage fraction * identity
        def _score(profile):
            cov_frac = (
                profile["element_covered_bp"] / element_length
                if element_length > 0
                else 0
            )
            return cov_frac * profile["weighted_pident"] / 100

        ranked = sorted(
            type_profiles.items(), key=lambda x: _score(x[1]), reverse=True
        )

        # Primary type
        primary_name, primary_prof = ranked[0]
        claimed = list(primary_prof["element_footprint"])
        hybrid_components = [primary_name]
        novel_bp_per_type = {primary_name: 0}

        # Check remaining types for novel coverage
        for type_name, profile in ranked[1:]:
            novel, novel_bp = self._calculate_novel_coverage(
                claimed, profile["element_footprint"]
            )
            novel_bp_per_type[type_name] = novel_bp

            is_novel = (
                novel_bp >= NOVEL_BP_THRESHOLD
                and (novel_bp / element_length if element_length > 0 else 0)
                >= NOVEL_FRAC_THRESHOLD
            )
            if is_novel:
                hybrid_components.append(type_name)
                claimed = self._merge_intervals(claimed + novel)

        # -- Classify using typing context --
        raw_call = "multi_type" if len(hybrid_components) > 1 else "canonical"
        hybrid_call = classify_with_typing(
            raw_call, len(hybrid_components), typing_info
        )

        # -- Build summary dict --
        secondary_name = hybrid_components[1] if len(hybrid_components) > 1 else "-"
        secondary_prof = type_profiles.get(secondary_name, {})

        summary = {
            "element_id": input_name,
            "hybrid_best_match": primary_name,
            "hybrid_best_subtype": primary_prof["best_subtype"],
            "hybrid_best_coverage": primary_prof["best_coverage"],
            "hybrid_best_identity": primary_prof["weighted_pident"],
            "hybrid_secondary_match": secondary_name,
            "hybrid_secondary_subtype": secondary_prof.get(
                "best_subtype", "-"
            ),
            "hybrid_secondary_coverage": secondary_prof.get(
                "best_coverage", "-"
            ),
            "hybrid_secondary_identity": secondary_prof.get(
                "weighted_pident", "-"
            ),
            "hybrid_call": hybrid_call,
            "hybrid_components": ";".join(hybrid_components),
        }

        # -- Build detail rows (all types, for per-element detail files) --
        detail_rows = []
        for type_name, profile in ranked:
            ele_cov = (
                round(profile["element_covered_bp"] / element_length * 100, 1)
                if element_length > 0
                else 0.0
            )
            detail_rows.append(
                {
                    "element_id": input_name,
                    "parent_type": type_name,
                    "best_subtype": profile["best_subtype"],
                    "ref_coverage": profile["best_coverage"],
                    "element_coverage": ele_cov,
                    "weighted_pident": profile["weighted_pident"],
                    "covered_bp": profile["element_covered_bp"],
                    "element_covered_bp": profile["element_covered_bp"],
                    "ref_length": profile["best_ref_length"],
                    "element_length": element_length,
                    "num_hsps": profile["num_hsps"],
                    "novel_bp": novel_bp_per_type.get(type_name, 0),
                }
            )

        return summary, detail_rows


# ---------------------------------------------------------------------------
# Standalone CLI
# ---------------------------------------------------------------------------

def _collect_input_files(paths: list) -> list:
    """Resolve a mix of files and directories into FASTA file paths."""
    fasta_extensions = {".fasta", ".fa", ".fna", ".fsa"}
    collected = []
    for p in paths:
        path = Path(p)
        if path.is_dir():
            for child in sorted(path.iterdir()):
                if child.suffix.lower() in fasta_extensions:
                    collected.append(str(child))
        elif path.is_file():
            collected.append(str(path))
    return collected


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Examine extracted SCC elements for hybrid composition "
            "by BLAST comparison against SCCmec type references."
        )
    )
    parser.add_argument(
        "-f",
        "--fasta",
        nargs="+",
        required=True,
        help="Input SCC element FASTA file(s) or directory",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        default=".",
        help="Output directory (default: current directory)",
    )
    parser.add_argument(
        "-d",
        "--detailed",
        action="store_true",
        help="Write per-element detail files to hybrid_detail/ subdirectory",
    )
    parser.add_argument(
        "-t",
        "--typing-report",
        default=None,
        help=(
            "Path to sccmec_summary.tsv from sccmec-type or sccmec-pipeline. "
            "When provided, multi-type elements are classified as 'composite' "
            "(multiple ccr complexes) or 'hybrid' (single ccr complex). "
            "Without this, multi-type elements are reported as 'multi_type'."
        ),
    )

    args = parser.parse_args()

    input_files = _collect_input_files(args.fasta)
    if not input_files:
        print("ERROR: No FASTA files found in the provided paths",
              file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(input_files)} input file(s)", file=sys.stderr)

    # Load typing context if provided
    typing_context = {}
    if args.typing_report:
        if not os.path.isfile(args.typing_report):
            print(f"ERROR: Typing report not found: {args.typing_report}",
                  file=sys.stderr)
            sys.exit(1)
        typing_context = read_typing_report(args.typing_report)
        print(f"Loaded typing context for {len(typing_context)} elements",
              file=sys.stderr)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    summary_path = outdir / "hybrid_summary.tsv"

    detail_dir = None
    if args.detailed:
        detail_dir = outdir / "hybrid_detail"
        detail_dir.mkdir(parents=True, exist_ok=True)

    typer = HybridTyper()

    # Batch BLAST: single BLAST call for all elements
    print("Running batch BLAST comparison...", file=sys.stderr, flush=True)
    batch_results = typer.type_batch(
        input_files, typing_context=typing_context or None
    )

    with open(summary_path, "w") as sf:
        sf.write("\t".join(HYBRID_SUMMARY_HEADER) + "\n")

        for i, (input_file, (summary, detail_rows)) in enumerate(
            zip(input_files, batch_results), 1
        ):
            stem = Path(input_file).stem
            sf.write(
                "\t".join(str(summary[c]) for c in HYBRID_SUMMARY_HEADER)
                + "\n"
            )

            # Write per-element detail file if requested
            if detail_dir is not None and detail_rows:
                detail_path = detail_dir / f"{stem}_hybrid_detail.tsv"
                with open(detail_path, "w") as df:
                    df.write(
                        "\t".join(HYBRID_DETAIL_HEADER) + "\n"
                    )
                    for row in detail_rows:
                        df.write(
                            "\t".join(
                                str(row[c]) for c in HYBRID_DETAIL_HEADER
                            )
                            + "\n"
                        )

            print(
                f"  [{i}/{len(input_files)}] {stem}... "
                f"{summary['hybrid_call']}",
                file=sys.stderr,
            )

    print(f"\nSummary: {summary_path}", file=sys.stderr)
    if detail_dir is not None:
        print(f"Detail:  {detail_dir}/", file=sys.stderr)


if __name__ == "__main__":
    main()
