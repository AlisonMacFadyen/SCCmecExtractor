#!/usr/bin/env python

"""Tests for report_sccmec.py"""

import subprocess

import pytest

from sccmecextractor.report_sccmec import (
    EXTRACTION_HEADER,
    TYPING_EXTRA_COLS,
    UNIFIED_HEADER,
    merge_reports,
    normalise_typing_keys,
    read_tsv,
    write_unified_report,
    _extract_mec_resistance_gene,
    _format_ccr,
    _build_summary_row,
    _MEC_GENE_PATTERN,
    SUMMARY_HEADER,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_extraction_tsv(path, rows):
    """Write a minimal extraction report TSV."""
    header = "\t".join(EXTRACTION_HEADER)
    lines = [header]
    for row in rows:
        vals = [row.get(c, "-") for c in EXTRACTION_HEADER]
        lines.append("\t".join(vals))
    path.write_text("\n".join(lines) + "\n")


def _write_typing_tsv(path, rows):
    """Write a minimal typing results TSV."""
    from sccmecextractor.sccmec_type_classification import TYPING_HEADER

    header = "\t".join(TYPING_HEADER)
    lines = [header]
    for row in rows:
        vals = [row.get(c, "-") for c in TYPING_HEADER]
        lines.append("\t".join(vals))
    path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# TestReadTsv
# ---------------------------------------------------------------------------

class TestReadTsv:
    """Tests for read_tsv."""

    def test_reads_tsv(self, tmp_path):
        """Normal TSV is read and keyed by Input_File."""
        tsv = tmp_path / "data.tsv"
        tsv.write_text(
            "Input_File\tStatus\ngenome1\textracted\ngenome2\tfailed\n"
        )
        rows = read_tsv(str(tsv))
        assert "genome1" in rows
        assert rows["genome1"]["Status"] == "extracted"
        assert "genome2" in rows

    def test_duplicate_headers_skipped(self, tmp_path):
        """Duplicate header lines from batch concatenation are ignored."""
        tsv = tmp_path / "batch.tsv"
        tsv.write_text(
            "Input_File\tStatus\n"
            "genome1\textracted\n"
            "Input_File\tStatus\n"
            "genome2\tfailed\n"
        )
        rows = read_tsv(str(tsv))
        assert len(rows) == 2
        assert "genome1" in rows
        assert "genome2" in rows

    def test_missing_file_exits(self, tmp_path):
        """Missing file causes sys.exit."""
        with pytest.raises(SystemExit):
            read_tsv(str(tmp_path / "nonexistent.tsv"))


# ---------------------------------------------------------------------------
# TestNormaliseTypingKeys
# ---------------------------------------------------------------------------

class TestNormaliseTypingKeys:
    """Tests for normalise_typing_keys."""

    def test_strips_sccmec_suffix(self):
        """Keys ending with _SCCmec have suffix stripped."""
        rows = {
            "genome1_SCCmec": {"Input_File": "genome1_SCCmec", "mec_genes": "mecA(full)"},
        }
        result = normalise_typing_keys(rows)
        assert "genome1" in result
        assert "genome1_SCCmec" not in result
        assert result["genome1"]["Input_File"] == "genome1"

    def test_leaves_other_names_unchanged(self):
        """Keys without _SCCmec suffix are unchanged."""
        rows = {
            "genome1": {"Input_File": "genome1", "mec_genes": "-"},
        }
        result = normalise_typing_keys(rows)
        assert "genome1" in result


# ---------------------------------------------------------------------------
# TestMergeReports
# ---------------------------------------------------------------------------

class TestMergeReports:
    """Tests for merge_reports."""

    def test_matching_keys_all_filled(self):
        """Both extraction and typing present → all columns filled."""
        ext = {"g1": {"Input_File": "g1", "Status": "extracted"}}
        typ = {"g1": {"Input_File": "g1", "mec_genes": "mecA(full)"}}
        merged = merge_reports(ext, typ)
        assert len(merged) == 1
        assert merged[0]["Input_File"] == "g1"
        assert merged[0]["Status"] == "extracted"
        assert merged[0]["mec_genes"] == "mecA(full)"
        assert merged[0]["typing_source"] == "sccmec"
        assert merged[0]["element_type"] == "SCCmec"

    def test_extraction_only(self):
        """Extraction-only entry gets dash for typing columns."""
        ext = {"g1": {"Input_File": "g1", "Status": "extracted"}}
        typ = {}
        merged = merge_reports(ext, typ)
        assert len(merged) == 1
        assert merged[0]["mec_genes"] == "-"
        assert merged[0]["ccr_complex_type"] == "-"
        assert merged[0]["typing_source"] == "-"
        assert merged[0]["element_type"] == "SCC"

    def test_typing_only(self):
        """Typing-only entry gets dash for extraction columns, source is wgs."""
        ext = {}
        typ = {"g1": {"Input_File": "g1", "mec_genes": "mecA(full)"}}
        merged = merge_reports(ext, typ)
        assert len(merged) == 1
        assert merged[0]["Status"] == "-"
        assert merged[0]["mec_genes"] == "mecA(full)"
        # No extraction status → inferred as wgs
        assert merged[0]["typing_source"] == "wgs"
        assert merged[0]["element_type"] == "-"

    def test_full_outer_join(self):
        """Both sides present, some keys only in one side."""
        ext = {"g1": {"Input_File": "g1", "Status": "extracted"}}
        typ = {
            "g1": {"Input_File": "g1", "mec_genes": "mecA(full)"},
            "g2": {"Input_File": "g2", "mec_genes": "mecC(full)"},
        }
        merged = merge_reports(ext, typ)
        assert len(merged) == 2
        keys = [r["Input_File"] for r in merged]
        assert keys == ["g1", "g2"]

    def test_empty_inputs(self):
        """Both empty → empty result."""
        assert merge_reports({}, {}) == []

    def test_sorted_output(self):
        """Results are sorted by Input_File."""
        ext = {"z_genome": {"Input_File": "z_genome"}}
        typ = {"a_genome": {"Input_File": "a_genome"}}
        merged = merge_reports(ext, typ)
        keys = [r["Input_File"] for r in merged]
        assert keys == ["a_genome", "z_genome"]


# ---------------------------------------------------------------------------
# TestTypingSourceInference
# ---------------------------------------------------------------------------

class TestTypingSourceInference:
    """Tests for typing_source inference from extraction status."""

    def test_extracted_genome_is_sccmec(self):
        """Typing + Status=extracted → typing_source='sccmec'."""
        ext = {"g1": {"Input_File": "g1", "Status": "extracted"}}
        typ = {"g1": {"Input_File": "g1", "mec_genes": "mecA(full)"}}
        merged = merge_reports(ext, typ)
        assert merged[0]["typing_source"] == "sccmec"

    def test_composite_extracted_is_sccmec(self):
        """Typing + Status=composite_extracted → typing_source='sccmec'."""
        ext = {"g1": {"Input_File": "g1", "Status": "composite_extracted"}}
        typ = {"g1": {"Input_File": "g1", "mec_genes": "mecA(full)"}}
        merged = merge_reports(ext, typ)
        assert merged[0]["typing_source"] == "sccmec"

    def test_fallback_extracted_is_sccmec(self):
        """Typing + Status=fallback_extracted → typing_source='sccmec'."""
        ext = {"g1": {"Input_File": "g1", "Status": "fallback_extracted"}}
        typ = {"g1": {"Input_File": "g1", "mec_genes": "mecA(full)"}}
        merged = merge_reports(ext, typ)
        assert merged[0]["typing_source"] == "sccmec"

    def test_composite_fallback_extracted_is_sccmec(self):
        """Typing + Status=composite_fallback_extracted → typing_source='sccmec'."""
        ext = {"g1": {"Input_File": "g1", "Status": "composite_fallback_extracted"}}
        typ = {"g1": {"Input_File": "g1", "mec_genes": "mecA(full)"}}
        merged = merge_reports(ext, typ)
        assert merged[0]["typing_source"] == "sccmec"

    def test_failed_genome_is_wgs(self):
        """Typing + Status=failed → typing_source='wgs'."""
        ext = {"g1": {"Input_File": "g1", "Status": "failed"}}
        typ = {"g1": {"Input_File": "g1", "mec_genes": "mecA(full)"}}
        merged = merge_reports(ext, typ)
        assert merged[0]["typing_source"] == "wgs"

    def test_no_extraction_is_wgs(self):
        """Typing but no extraction row → typing_source='wgs'."""
        ext = {}
        typ = {"g1": {"Input_File": "g1", "mec_genes": "mecA(full)"}}
        merged = merge_reports(ext, typ)
        assert merged[0]["typing_source"] == "wgs"

    def test_mixed_sources(self):
        """Different genomes get different typing sources based on status."""
        ext = {
            "g1": {"Input_File": "g1", "Status": "extracted"},
            "g2": {"Input_File": "g2", "Status": "failed"},
            "g3": {"Input_File": "g3", "Status": "failed"},
        }
        typ = {
            "g1": {"Input_File": "g1", "mec_genes": "mecA(full)"},
            "g2": {"Input_File": "g2", "mec_genes": "mecC(full)"},
        }
        merged = merge_reports(ext, typ)

        by_key = {r["Input_File"]: r for r in merged}
        assert by_key["g1"]["typing_source"] == "sccmec"
        assert by_key["g2"]["typing_source"] == "wgs"
        assert by_key["g3"]["typing_source"] == "-"

    def test_no_typing_data(self):
        """Failed extraction + no typing → typing_source='-'."""
        ext = {"g1": {"Input_File": "g1", "Status": "failed"}}
        typ = {}
        merged = merge_reports(ext, typ)
        assert merged[0]["typing_source"] == "-"
        assert merged[0]["mec_genes"] == "-"


# ---------------------------------------------------------------------------
# TestElementType
# ---------------------------------------------------------------------------

class TestElementType:
    """Tests for element_type classification."""

    def test_extracted_with_mec_is_sccmec(self):
        """Extracted element with mec genes → element_type='SCCmec'."""
        ext = {"g1": {"Input_File": "g1", "Status": "extracted"}}
        typ = {"g1": {"Input_File": "g1", "mec_genes": "mecA(full)"}}
        merged = merge_reports(ext, typ)
        assert merged[0]["element_type"] == "SCCmec"

    def test_extracted_without_mec_is_scc(self):
        """Extracted element without mec genes → element_type='SCC'."""
        ext = {"g1": {"Input_File": "g1", "Status": "extracted"}}
        typ = {"g1": {"Input_File": "g1", "mec_genes": "-"}}
        merged = merge_reports(ext, typ)
        assert merged[0]["element_type"] == "SCC"

    def test_extracted_no_typing_is_scc(self):
        """Extracted but no typing data → mec_genes defaults to '-' → 'SCC'."""
        ext = {"g1": {"Input_File": "g1", "Status": "extracted"}}
        merged = merge_reports(ext, {})
        assert merged[0]["element_type"] == "SCC"

    def test_fallback_extracted_with_mec_is_sccmec(self):
        """Fallback-extracted element with mec → element_type='SCCmec'."""
        ext = {"g1": {"Input_File": "g1", "Status": "fallback_extracted"}}
        typ = {"g1": {"Input_File": "g1", "mec_genes": "mecA(full)"}}
        merged = merge_reports(ext, typ)
        assert merged[0]["element_type"] == "SCCmec"

    def test_fallback_extracted_without_mec_is_scc(self):
        """Fallback-extracted element without mec → element_type='SCC'."""
        ext = {"g1": {"Input_File": "g1", "Status": "fallback_extracted"}}
        typ = {"g1": {"Input_File": "g1", "mec_genes": "-"}}
        merged = merge_reports(ext, typ)
        assert merged[0]["element_type"] == "SCC"

    def test_composite_extracted_with_mec_is_sccmec(self):
        """Composite-extracted element with mec → element_type='SCCmec'."""
        ext = {"g1": {"Input_File": "g1", "Status": "composite_extracted"}}
        typ = {"g1": {"Input_File": "g1", "mec_genes": "mecA(full)"}}
        merged = merge_reports(ext, typ)
        assert merged[0]["element_type"] == "SCCmec"

    def test_composite_fallback_extracted_without_mec_is_scc(self):
        """Composite-fallback-extracted without mec → element_type='SCC'."""
        ext = {"g1": {"Input_File": "g1", "Status": "composite_fallback_extracted"}}
        typ = {"g1": {"Input_File": "g1", "mec_genes": "-"}}
        merged = merge_reports(ext, typ)
        assert merged[0]["element_type"] == "SCC"

    def test_failed_is_dash(self):
        """Failed extraction → element_type='-'."""
        ext = {"g1": {"Input_File": "g1", "Status": "failed"}}
        typ = {"g1": {"Input_File": "g1", "mec_genes": "mecA(full)"}}
        merged = merge_reports(ext, typ)
        assert merged[0]["element_type"] == "-"

    def test_no_ccr_element_is_dash(self):
        """no_ccr_element status → element_type='-'."""
        ext = {"g1": {"Input_File": "g1", "Status": "no_ccr_element"}}
        typ = {"g1": {"Input_File": "g1", "mec_genes": "-"}}
        merged = merge_reports(ext, typ)
        assert merged[0]["element_type"] == "-"

    def test_no_extraction_row_is_dash(self):
        """No extraction row at all → element_type='-'."""
        typ = {"g1": {"Input_File": "g1", "mec_genes": "mecA(full)"}}
        merged = merge_reports({}, typ)
        assert merged[0]["element_type"] == "-"

    def test_mixed_element_types(self):
        """Multiple genomes with different element types."""
        ext = {
            "g1": {"Input_File": "g1", "Status": "extracted"},
            "g2": {"Input_File": "g2", "Status": "extracted"},
            "g3": {"Input_File": "g3", "Status": "failed"},
        }
        typ = {
            "g1": {"Input_File": "g1", "mec_genes": "mecA(full)"},
            "g2": {"Input_File": "g2", "mec_genes": "-"},
            "g3": {"Input_File": "g3", "mec_genes": "mecA(full)"},
        }
        merged = merge_reports(ext, typ)
        by_key = {r["Input_File"]: r for r in merged}
        assert by_key["g1"]["element_type"] == "SCCmec"
        assert by_key["g2"]["element_type"] == "SCC"
        assert by_key["g3"]["element_type"] == "-"


# ---------------------------------------------------------------------------
# TestWriteUnifiedReport
# ---------------------------------------------------------------------------

class TestWriteUnifiedReport:
    """Tests for write_unified_report."""

    def test_correct_header_and_data(self, tmp_path):
        """Output file has correct header and data rows."""
        merged = [
            {col: "val" for col in UNIFIED_HEADER},
        ]
        merged[0]["Input_File"] = "genome1"
        outfile = tmp_path / "report.tsv"
        write_unified_report(merged, str(outfile))

        lines = outfile.read_text().strip().split("\n")
        assert len(lines) == 2  # header + 1 data row
        assert lines[0].split("\t") == UNIFIED_HEADER
        assert lines[1].split("\t")[0] == "genome1"


# ---------------------------------------------------------------------------
# TestCLI
# ---------------------------------------------------------------------------

class TestCLI:
    """Tests for CLI entry point."""

    def test_help(self):
        """--help runs without error."""
        result = subprocess.run(
            ["python", "-m", "sccmecextractor.report_sccmec", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "extraction-report" in result.stdout
        assert "typing-results" in result.stdout

    def test_full_merge_run(self, tmp_path):
        """End-to-end: merge extraction + typing → unified report."""
        ext_file = tmp_path / "extraction.tsv"
        _write_extraction_tsv(ext_file, [
            {"Input_File": "genome1", "Status": "extracted", "Contig": "c1"},
        ])

        typ_file = tmp_path / "typing.tsv"
        _write_typing_tsv(typ_file, [
            {
                "Input_File": "genome1_SCCmec",
                "mec_genes": "mecA(full)",
                "mec_identity": "99.0",
                "mec_coverage": "100.0",
                "ccr_genes": "ccrA2(full);ccrB2(full)",
                "ccr_allotypes": "ccrA2;ccrB2",
                "ccr_identity": "95.0;93.0",
                "ccr_complex_type": "2",
            },
        ])

        out_file = tmp_path / "unified.tsv"
        result = subprocess.run(
            [
                "python", "-m", "sccmecextractor.report_sccmec",
                "-e", str(ext_file),
                "-t", str(typ_file),
                "-o", str(out_file),
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"Failed:\n{result.stderr}"
        assert out_file.exists()

        lines = out_file.read_text().strip().split("\n")
        header = lines[0].split("\t")
        assert header == UNIFIED_HEADER
        # Data row should have genome1 (suffix stripped from typing)
        assert lines[1].startswith("genome1\t")
        data = lines[1].split("\t")
        # Check extraction col
        assert data[header.index("Status")] == "extracted"
        # Check typing col
        assert data[header.index("mec_genes")] == "mecA(full)"
        # Check typing source
        assert data[header.index("typing_source")] == "sccmec"
        # Check element type
        assert data[header.index("element_type")] == "SCCmec"

    def test_cli_mixed_sources(self, tmp_path):
        """End-to-end: typing source inferred from extraction status."""
        ext_file = tmp_path / "extraction.tsv"
        _write_extraction_tsv(ext_file, [
            {"Input_File": "genome1", "Status": "extracted", "Contig": "c1"},
            {"Input_File": "genome2", "Status": "failed", "Failure_Reason": "right_only"},
        ])

        # Single typing file with results from both sources
        typ_file = tmp_path / "typing.tsv"
        _write_typing_tsv(typ_file, [
            {
                "Input_File": "genome1_SCCmec",
                "mec_genes": "mecA(full)",
                "mec_identity": "99.0",
                "mec_coverage": "100.0",
                "ccr_genes": "ccrA2(full);ccrB2(full)",
                "ccr_allotypes": "ccrA2;ccrB2",
                "ccr_identity": "95.0;93.0",
                "ccr_complex_type": "2",
            },
            {
                "Input_File": "genome2",
                "mec_genes": "mecA(full)",
                "mec_identity": "98.5",
                "mec_coverage": "99.0",
                "ccr_genes": "ccrC1(full)",
                "ccr_allotypes": "ccrC1",
                "ccr_identity": "91.0",
                "ccr_complex_type": "5",
            },
        ])

        out_file = tmp_path / "unified.tsv"
        result = subprocess.run(
            [
                "python", "-m", "sccmecextractor.report_sccmec",
                "-e", str(ext_file),
                "-t", str(typ_file),
                "-o", str(out_file),
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"Failed:\n{result.stderr}"

        lines = out_file.read_text().strip().split("\n")
        header = lines[0].split("\t")
        assert header == UNIFIED_HEADER

        # genome1: extracted → typing_source='sccmec'
        d1 = dict(zip(header, lines[1].split("\t")))
        assert d1["Input_File"] == "genome1"
        assert d1["typing_source"] == "sccmec"
        assert d1["Status"] == "extracted"
        assert d1["element_type"] == "SCCmec"

        # genome2: failed → typing_source='wgs'
        d2 = dict(zip(header, lines[2].split("\t")))
        assert d2["Input_File"] == "genome2"
        assert d2["typing_source"] == "wgs"
        assert d2["Status"] == "failed"
        assert d2["mec_genes"] == "mecA(full)"
        assert d2["ccr_complex_type"] == "5"
        assert d2["element_type"] == "-"

    def test_missing_input_error(self, tmp_path):
        """Missing input file causes non-zero exit."""
        result = subprocess.run(
            [
                "python", "-m", "sccmecextractor.report_sccmec",
                "-e", str(tmp_path / "missing.tsv"),
                "-t", str(tmp_path / "missing2.tsv"),
                "-o", str(tmp_path / "out.tsv"),
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0


# ---------------------------------------------------------------------------
# element_type classification bug fix tests
# ---------------------------------------------------------------------------

class TestElementTypeClassification:
    """Ensure element_type uses only actual mec resistance genes, not IS/regulatory."""

    def test_is431_alone_is_scc_not_sccmec(self):
        """IS431 without a mec resistance gene should be SCC, not SCCmec."""
        mec_genes = "IS431(full);IS431(novel_full)"
        genes_list = mec_genes.split(";")
        has_mec = any(_MEC_GENE_PATTERN.match(g.strip()) for g in genes_list)
        assert has_mec is False

    def test_mecA_is_sccmec(self):
        mec_genes = "mecA(full);IS431(full)"
        genes_list = mec_genes.split(";")
        has_mec = any(_MEC_GENE_PATTERN.match(g.strip()) for g in genes_list)
        assert has_mec is True

    def test_mecC_is_sccmec(self):
        mec_genes = "mecC(full);mecI(full)"
        genes_list = mec_genes.split(";")
        has_mec = any(_MEC_GENE_PATTERN.match(g.strip()) for g in genes_list)
        assert has_mec is True

    def test_mecR1_alone_is_not_sccmec(self):
        mec_genes = "mecR1(full);IS1272(full)"
        genes_list = mec_genes.split(";")
        has_mec = any(_MEC_GENE_PATTERN.match(g.strip()) for g in genes_list)
        assert has_mec is False

    def test_empty_is_not_sccmec(self):
        mec_genes = "-"
        genes_list = mec_genes.split(";") if mec_genes != "-" else []
        has_mec = any(_MEC_GENE_PATTERN.match(g.strip()) for g in genes_list)
        assert has_mec is False

    def test_mecA1_is_sccmec(self):
        mec_genes = "mecA1(full)"
        genes_list = mec_genes.split(";")
        has_mec = any(_MEC_GENE_PATTERN.match(g.strip()) for g in genes_list)
        assert has_mec is True

    def test_mecD_is_sccmec(self):
        mec_genes = "mecD(novel_full)"
        genes_list = mec_genes.split(";")
        has_mec = any(_MEC_GENE_PATTERN.match(g.strip()) for g in genes_list)
        assert has_mec is True


# ---------------------------------------------------------------------------
# Summary helper tests
# ---------------------------------------------------------------------------

class TestExtractMecResistanceGene:
    def test_mecA(self):
        assert _extract_mec_resistance_gene("mecA(full);IS431(full)") == "mecA"

    def test_mecC(self):
        assert _extract_mec_resistance_gene("mecC(full);mecI(full)") == "mecC"

    def test_is431_only(self):
        assert _extract_mec_resistance_gene("IS431(full);IS431(novel_full)") == "-"

    def test_no_genes(self):
        assert _extract_mec_resistance_gene("-") == "-"

    def test_empty(self):
        assert _extract_mec_resistance_gene("") == "-"

    def test_mecA1(self):
        assert _extract_mec_resistance_gene("mecA1(full)") == "mecA1"

    def test_mecA_with_regulatory(self):
        assert _extract_mec_resistance_gene(
            "mecA(full);mecR1(full);IS431(full);mecI(full)"
        ) == "mecA"


class TestFormatCcr:
    def test_known_complex(self):
        assert _format_ccr("2", "ccrA2;ccrB2") == "2"

    def test_multiple_known(self):
        assert _format_ccr("2;5", "ccrA2;ccrB2;ccrC1") == "2;5"

    def test_novel_combination(self):
        assert _format_ccr("novel_combination", "ccrA5;ccrB3") == "ccrA5;ccrB3"

    def test_no_ccr(self):
        assert _format_ccr("-", "-") == "-"

    def test_empty(self):
        assert _format_ccr("", "") == "-"


class TestBuildSummaryRow:
    def test_extracted_sccmec(self):
        row = {
            "Input_File": "GCF_001",
            "Status": "extracted",
            "mec_genes": "mecA(full);IS431(full)",
            "mec_context": "in_element",
            "ccr_complex_type": "2",
            "ccr_allotypes": "ccrA2;ccrB2",
            "SCCmec_Type": "IV",
            "element_type": "SCCmec",
            "typing_source": "sccmec",
            "Failure_Reason": "-",
        }
        summary = _build_summary_row(row)
        assert summary["mec_gene"] == "mecA"
        assert summary["ccr"] == "2"
        assert summary["element_type"] == "SCCmec"
        assert summary["closest_ref"] == "-"
        assert summary["hybrid_call"] == "-"

    def test_extracted_scc_with_adjacent_mec(self):
        row = {
            "Input_File": "GCF_002",
            "Status": "extracted",
            "mec_genes": "IS431(full);IS431(novel_full)",
            "mec_context": "mec_adjacent",
            "ccr_complex_type": "5",
            "ccr_allotypes": "ccrC1",
            "SCCmec_Type": "not_typeable",
            "element_type": "SCC",
            "typing_source": "sccmec",
            "Failure_Reason": "-",
            "wgs_mec_genes": "mecA(full);IS1272(full)",
        }
        summary = _build_summary_row(row)
        assert summary["mec_gene"] == "mecA"  # from WGS fallback
        assert summary["mec_context"] == "mec_adjacent"
        assert summary["element_type"] == "SCC"  # element itself is SCC

    def test_wgs_typed(self):
        row = {
            "Input_File": "GCF_003",
            "Status": "failed",
            "mec_genes": "mecA(full)",
            "mec_context": "-",
            "ccr_complex_type": "2",
            "ccr_allotypes": "ccrA2;ccrB2",
            "SCCmec_Type": "IV",
            "element_type": "-",
            "typing_source": "wgs",
            "Failure_Reason": "cross_contig",
        }
        summary = _build_summary_row(row)
        assert summary["element_type"] == "WGS"
        assert summary["mec_gene"] == "mecA"

    def test_novel_ccr(self):
        row = {
            "Input_File": "GCF_004",
            "Status": "extracted",
            "mec_genes": "-",
            "mec_context": "-",
            "ccr_complex_type": "novel_combination",
            "ccr_allotypes": "ccrA5;ccrB3",
            "SCCmec_Type": "-",
            "element_type": "SCC",
            "typing_source": "sccmec",
            "Failure_Reason": "-",
        }
        summary = _build_summary_row(row)
        assert summary["ccr"] == "ccrA5;ccrB3"

    def test_with_hybrid_info(self):
        row = {
            "Input_File": "GCF_005",
            "Status": "extracted",
            "mec_genes": "mecA(full)",
            "mec_context": "in_element",
            "ccr_complex_type": "2",
            "ccr_allotypes": "ccrA2;ccrB2",
            "SCCmec_Type": "IV",
            "element_type": "SCCmec",
            "typing_source": "sccmec",
            "Failure_Reason": "-",
        }
        hybrid_info = {
            "hybrid_best_match": "IV",
            "hybrid_best_subtype": "IVa",
            "hybrid_best_coverage": 94.9,
            "hybrid_best_identity": 98.5,
            "hybrid_call": "canonical",
            "hybrid_components": "IV",
        }
        summary = _build_summary_row(row, hybrid_info)
        assert summary["closest_ref"] == "IVa (94.9%)"
        assert summary["hybrid_call"] == "canonical"

    def test_with_hybrid_no_match(self):
        row = {
            "Input_File": "GCF_006",
            "Status": "extracted",
            "mec_genes": "-",
            "mec_context": "-",
            "ccr_complex_type": "5",
            "ccr_allotypes": "ccrC1",
            "SCCmec_Type": "-",
            "element_type": "SCC",
            "typing_source": "sccmec",
            "Failure_Reason": "-",
        }
        hybrid_info = {
            "hybrid_best_match": "-",
            "hybrid_best_subtype": "-",
            "hybrid_best_coverage": "-",
            "hybrid_best_identity": "-",
            "hybrid_call": "no_match",
            "hybrid_components": "-",
        }
        summary = _build_summary_row(row, hybrid_info)
        assert summary["closest_ref"] == "-"
        assert summary["hybrid_call"] == "no_match"

    def test_summary_has_all_headers(self):
        row = {
            "Input_File": "GCF_007",
            "Status": "extracted",
            "mec_genes": "mecA(full)",
            "mec_context": "in_element",
            "ccr_complex_type": "2",
            "ccr_allotypes": "ccrA2;ccrB2",
            "SCCmec_Type": "IV",
            "element_type": "SCCmec",
            "typing_source": "sccmec",
            "Failure_Reason": "-",
        }
        summary = _build_summary_row(row)
        for col in SUMMARY_HEADER:
            assert col in summary
