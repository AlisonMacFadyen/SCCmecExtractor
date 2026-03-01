# SCCmecExtractor

A Python toolkit for extracting and typing SCC*mec* and non-*mec* SCC elements from *Staphylococcus* and *Mammaliicoccus* whole genome sequences.  The tool identifies attachment (*att*) sites, extracts complete SCC elements and performs gene-level typing of *mec* and *ccr* gene complexes.

**Note: extraction requires both *att* sites to be on the same contig as *rlmH*.  This is a same-contig requirement with assembly fragmentation being the main source of extraction failure, particularly in non-*aureus* species.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![PyPI version](https://img.shields.io/pypi/v/sccmecextractor)](https://pypi.org/project/sccmecextractor/)
[![Docker Image Version](https://img.shields.io/docker/v/alisonmacfadyen/sccmecextractor?sort=semver)](https://hub.docker.com/r/alisonmacfadyen/sccmecextractor)

[![CD - PyPI](https://github.com/AlisonMacFadyen/SCCmecExtractor/actions/workflows/cd-pypi.yaml/badge.svg)](https://github.com/AlisonMacFadyen/SCCmecExtractor/actions/workflows/cd-pypi.yaml)
[![CD - Docker](https://github.com/AlisonMacFadyen/SCCmecExtractor/actions/workflows/cd_docker.yaml/badge.svg)](https://github.com/AlisonMacFadyen/SCCmecExtractor/actions/workflows/cd_docker.yaml)


## Overview

SCC*mec*Extractor provides five CLI commands that work together to identify, extract and type SCC*mec* and non-*mec* SCC elements:

| Command | Description |
|---|---|
| `sccmec-pipeline` | Master pipeline orchestrating all steps (recommended) |
| `sccmec-locate-att` | Locate attachment (*att*) sites in genomic sequences |
| `sccmec-extract` | Extract SCC elements bounded by *att* site pairs |
| `sccmec-type` | Type extracted elements or WGS by *mec* and *ccr* gene content |
| `sccmec-report` | Merge extraction and typing results into a unified report |

### Key Capabilities

- **FASTA-only mode** — no GFF annotation required; *rlmH* detected via BLAST against a 70-species reference database
- **Non-*mec* SCC detection** — extracts SCC elements that carry *ccr* genes but lack *mec* genes
- **Composite element detection** — identifies tandem/nested SCC elements with multiple *att* site pairs
- **Gene-level typing** — classifies *mec* complex (*mecA*, *mecB*, *mecC*, *mecD* allotypes) and *ccr* complex (*ccrA/B*, *ccrC* allotypes) via BLAST
- **Cross-genus support** — validated on *Staphylococcus* (64 species) and *Mammaliicoccus* (6 species)

## Table of Contents

- [Installation](#installation)
  - [Using Conda/Mamba](#using-condamamba-recommended)
  - [Using pip](#using-pip)
  - [Using Docker](#using-docker)
  - [Using Singularity](#using-singularity)
- [Requirements](#requirements)
- [Usage](#usage)
  - [FASTA-only Mode (Recommended)](#fasta-only-mode-recommended)
  - [GFF Mode](#gff-mode)
  - [Command Reference](#command-reference)
- [Complete Workflow Examples](#complete-workflow-examples)
- [How It Works](#how-it-works)
- [Output Format](#output-format)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)
- [Contributing](#contributing)
- [Contact](#contact)

## Installation

### Using Conda/Mamba (Recommended)

```bash
# Create a new environment
conda create -n sccmecextractor python=3.11
conda activate sccmecextractor

# Install dependencies
conda install -c conda-forge -c bioconda biopython blast

# Install Bakta only if you need GFF-based workflow
conda install -c conda-forge -c bioconda bakta

# Install SCCmecExtractor
pip install sccmecextractor

# Test that commands are available
sccmec-pipeline --help
sccmec-locate-att --help
sccmec-extract --help
sccmec-type --help
sccmec-report --help
```

If commands do not run, make sure the environment's `bin/` directory is in your PATH:

```bash
export PATH="$CONDA_PREFIX/bin:$PATH"
```

### Using pip

Note: installation with `pip` does not provide BLAST+ or Bakta. BLAST+ is required for *ccr* gene checks, typing and FASTA-only mode.

```bash
# Install SCCmecExtractor
pip install sccmecextractor

# Test that commands are available
sccmec-pipeline --help
```

### Using Docker

Docker provides a containerised environment with all dependencies pre-installed, including BLAST+ and Bakta.

```bash
# Pull the pre-built image
docker pull alisonmacfadyen/sccmecextractor:latest

# Or build from source
git clone https://github.com/AlisonMacFadyen/SCCmecExtractor.git
cd SCCmecExtractor
docker build -t sccmecextractor:latest -f containers/Dockerfile .
```

**Quick Start with Docker:**

```bash
# FASTA-only mode (no Bakta database needed)
docker run --rm \
  -v $PWD:/work \
  sccmecextractor:latest \
  sccmec-pipeline --fna-dir genomes/ -o results/

# With GFF annotation
docker run --rm \
  -v $PWD:/work \
  sccmecextractor:latest \
  sccmec-pipeline -f genome.fna -g genome.gff3 -o results/
```

See [CONTAINER_GUIDE.md](CONTAINER_GUIDE.md) for detailed Docker usage instructions.

### Using Singularity

Singularity/Apptainer is ideal for HPC environments where Docker is not available.

```bash
# Pull from Docker Hub
singularity pull docker://alisonmacfadyen/sccmecextractor:latest

# Or build from definition file
singularity build sccmecextractor.sif containers/sccmecextractor.def
```

**Quick Start with Singularity:**

```bash
# FASTA-only mode
singularity exec \
  --bind $PWD:/work \
  sccmecextractor.sif \
  sccmec-pipeline --fna-dir genomes/ -o results/

# With GFF annotation
singularity exec \
  --bind $PWD:/work \
  sccmecextractor.sif \
  sccmec-pipeline -f genome.fna -g genome.gff3 -o results/
```

See [CONTAINER_GUIDE.md](CONTAINER_GUIDE.md) for detailed Singularity usage instructions.

## Requirements

### Dependencies

* Python 3.9+
* Biopython
* BLAST+ (bundled in containers)
* Bakta (optional, only needed for GFF-based workflow) and automatically included in containers

### Input Files

* **Genome sequence**: `.fasta` or `.fna` file containing the assembled genome

When using GFF mode:
* **Gene annotations**: `.gff3` file with gene annotations (we recommend [Bakta](https://github.com/oschwengers/bakta) for annotation)

When using FASTA-only mode (`--blast-rlmh`): No annotation file is needed *rlmH* is detected via BLAST against a bundled 70-species reference database

### Bakta Database

Only required if using Bakta for GFF annotation:

```bash
# Light database (faster, smaller)
bakta_db download --output bakta_db --type light

# Full database
bakta_db download --output bakta_db
```

## Usage

### FASTA-only Mode (Recommended)

The simplest way to run SCC*mec*Extractor.  No GFF annotation or Bakta database is required — *rlmH* is detected via BLAST.  If using a novel species, recommended to supply *rlmH* gene reference.

```bash
# Single genome
sccmec-pipeline -f genome.fna -o results/

# Multiple genomes
sccmec-pipeline -f genome1.fna genome2.fna genome3.fna -o results/

# Directory of genomes
sccmec-pipeline --fna-dir genomes/ -o results/

# With composite element detection and multithreading
sccmec-pipeline --fna-dir genomes/ --composite -o results/ -t 4
```

### GFF Mode

If you have pre-computed GFF3 annotations (e.g. from Bakta), you can provide them.  GFF files are matched to FASTA files by stem name.

```bash
# Single genome with GFF
sccmec-pipeline -f genome.fna -g genome.gff3 -o results/

# Directories of genomes and GFFs
sccmec-pipeline --fna-dir genomes/ --gff-dir annotations/ -o results/
```

### Command Reference

#### `sccmec-pipeline`

Master pipeline that runs all steps: *att* site location, extraction, typing and report generation.

```
sccmec-pipeline [-h] (-f FNA [FNA ...] | --fna-dir FNA_DIR)
                [-g GFF [GFF ...] | --gff-dir GFF_DIR] [--blast-rlmh]
                [--rlmh-ref RLMH_REF] [--composite] -o OUTDIR [-t THREADS]
```

| Argument | Description |
|---|---|
| `-f`, `--fna` | One or more FASTA/FNA genome files |
| `--fna-dir` | Directory of FASTA/FNA genome files |
| `-g`, `--gff` | One or more GFF3 annotation files (matched by stem name) |
| `--gff-dir` | Directory of GFF3 files (matched by stem name) |
| `--blast-rlmh` | Use BLAST for *rlmH* detection (auto-enabled when no GFF provided) |
| `--rlmh-ref` | Custom *rlmH* reference FASTA for BLAST detection |
| `--composite` | Extract to outermost boundary for composite elements |
| `-o`, `--outdir` | Output directory for all results |
| `-t`, `--threads` | Number of parallel threads (default: 1) |

#### `sccmec-locate-att`

Identifies attachment sites in genomic sequences.

```
sccmec-locate-att [-h] -f FNA [-g GFF] -o OUTFILE [--blast-rlmh] [--rlmh-ref RLMH_REF]
```

| Argument | Description |
|---|---|
| `-f`, `--fna` | Input genome file (.fasta or .fna) |
| `-g`, `--gff` | Gene annotation file (.gff3 format, optional) |
| `-o`, `--outfile` | Output TSV file containing *att* site locations |
| `--blast-rlmh` | Use BLAST for *rlmH* detection (auto-enabled when no GFF provided) |
| `--rlmh-ref` | Custom *rlmH* reference FASTA |

#### `sccmec-extract`

Extracts SCC elements bounded by *att* site pairs.

```
sccmec-extract [-h] -f FNA [-g GFF] -a ATT -s SCCMEC [--composite] [-r REPORT]
               [--blast-rlmh] [--rlmh-ref RLMH_REF]
```

| Argument | Description |
|---|---|
| `-f`, `--fna` | Input genome file (.fasta or .fna) |
| `-g`, `--gff` | Gene annotation file (.gff3 format, optional) |
| `-a`, `--att` | TSV file from `sccmec-locate-att` with *att* site locations |
| `-s`, `--sccmec` | Output directory for extracted SCC sequences |
| `--composite` | Extract to outermost boundary for composite elements |
| `-r`, `--report` | Output TSV file for extraction report (appends for batch) |
| `--blast-rlmh` | Use BLAST for *rlmH* detection (auto-enabled when no GFF provided) |
| `--rlmh-ref` | Custom *rlmH* reference FASTA |

#### `sccmec-type`

Types extracted SCC elements (or whole genomes) by *mec* and *ccr* gene content using BLAST.

```
sccmec-type [-h] -f FASTA [FASTA ...] -o OUTFILE [--mec-ref MEC_REF] [--ccr-ref CCR_REF]
```

| Argument | Description |
|---|---|
| `-f`, `--fasta` | Input FASTA file(s) or directory of extracted SCC elements |
| `-o`, `--outfile` | Output TSV file for typing results |
| `--mec-ref` | Custom *mec* gene reference FASTA (default: bundled) |
| `--ccr-ref` | Custom *ccr* gene reference FASTA (default: bundled) |

#### `sccmec-report`

Merges extraction and typing results into a unified report.

```
sccmec-report [-h] -e EXTRACTION_REPORT -t TYPING_RESULTS -o OUTFILE
```

| Argument | Description |
|---|---|
| `-e`, `--extraction-report` | TSV from `sccmec-extract --report` |
| `-t`, `--typing-results` | TSV from `sccmec-type` |
| `-o`, `--outfile` | Output unified report TSV |

## Complete Workflow Examples

### FASTA-only Mode (Recommended)

```bash
# One-liner for a directory of genomes
sccmec-pipeline --fna-dir genomes/ -o results/ -t 4
```

### GFF Mode with Bakta Annotation

```bash
# 1. Annotate your genome with Bakta
bakta --db bakta_db genome.fna --output bakta_output

# 2. Run the pipeline with GFF
sccmec-pipeline -f genome.fna -g bakta_output/genome.gff3 -o results/
```

### Step-by-Step (Individual Commands)

```bash
# 1. Locate att sites
sccmec-locate-att -f genome.fna -o att_sites/att_sites.tsv

# 2. Extract SCC elements
sccmec-extract -f genome.fna -a att_sites/att_sites.tsv -s sccmec/ -r extraction_report.tsv

# 3. Type extracted elements and whole genome
sccmec-type -f sccmec/ genome.fna -o typing_results.tsv

# 4. Generate unified report
sccmec-report -e extraction_report.tsv -t typing_results.tsv -o unified_report.tsv
```

### Docker

```bash
# FASTA-only mode (simplest)
docker run --rm \
  -v $PWD:/work \
  sccmecextractor:latest \
  sccmec-pipeline --fna-dir genomes/ -o results/

# With Bakta annotation
docker run --rm \
  -v $PWD:/work \
  -v ~/bakta_db:/data/bakta_db \
  sccmecextractor:latest \
  bash -c "
    bakta --db /data/bakta_db genome.fna --output bakta_output --prefix genome && \
    sccmec-pipeline -f genome.fna -g bakta_output/genome.gff3 -o results/
  "
```

### Singularity

```bash
# FASTA-only mode (simplest)
singularity exec \
  --bind $PWD:/work \
  sccmecextractor.sif \
  sccmec-pipeline --fna-dir genomes/ -o results/

# With Bakta annotation
singularity exec \
  --bind $PWD:/work \
  --bind ~/bakta_db:/data/bakta_db \
  sccmecextractor.sif \
  bash -c "
    bakta --db /data/bakta_db genome.fna --output bakta_output --prefix genome && \
    sccmec-pipeline -f genome.fna -g bakta_output/genome.gff3 -o results/
  "
```

## How It Works

### Attachment Site Detection

The tool searches for 24 DNA motif patterns (8 attR/cattR + 16 attL/cattL) representing the attachment sites that flank SCC elements.  These patterns use regex with degeneracy to account for sequence variation across species.  *attR* sites are anchored within the *rlmH* gene, which serves as the chromosomal integration site.

- **attR/cattR**: Right attachment sites (within *rlmH*), including CcrC-associated attR2/cattR2 and attR3/cattR3 variants
- **attL/cattL**: Left attachment sites, including multiple bridge-length variants (attL2-attL8) covering diverse *ccr*-mediated recombination products

### Extraction Logic

1. **rlmH Detection**: Locates the *rlmH* gene via GFF annotation or BLAST against a 70-species reference database
2. **Site Validation**: *attR* sites must fall within or near *rlmH* (within 100 bp)
3. **Pair Selection**: Identifies the closest *attR*-*attL* pair on the same contig
4. **ccr Validation**: Verifies *ccr* genes are present between the *att* sites.  Elements without *ccr* are classified as `no_ccr_element` (not extracted)
5. **Size Filtering**: Rejects artefacts < 1,000 bp (pattern overlaps) and spurious matches > 200,000 bp
6. **Composite Detection**: Identifies tandem/nested elements with multiple *att* site pairs (with `--composite`)
7. **Origin-Spanning Detection**: Flags elements where chromosome linearisation splits the SCC across contig ends
8. **Fallback Extraction**: When standard extraction fails, but *attL* is identified, fallback extraction is utilised using the location of *rlmH* as a proxy for *attR*
9. **Strand Awareness**: Automatically handles reverse complement extraction when necessary

### Gene-Level Typing

`sccmec-type` carries out gene-typing by BLAST-based detection of:

- **_mec_ gene**: *mecA*, *mecB*, *mecC* *mecD* allotypes
- **_ccr_ complex**: *ccrA/ccrB* pairs and *ccrC* allotypes, incorporating all 22 *ccr* complex types.

## Output Format

### Pipeline Output Directory Structure

```
results/
├── ambiguous_att_sites.tsv # collates extraction limited genomes that may be of interest
├── att_sites/              # att site locations (one TSV per genome)
│   └── *.tsv
├── sccmec/                 # Extracted SCC element FASTAs
│   └── *_SCCmec.fasta
├── typing/                 # Typing results for extracted + whole genomes
│   └── *.tsv
├── extraction_report.tsv   # Per-genome extraction status and coordinates
├── typing_results.tsv      # Gene-level typing for all inputs
└── sccmec_unified_report.tsv  # Merged extraction + typing report
```

### Unified Report

The unified report (`sccmec_unified_report.tsv`) includes per-genome columns for extraction status, *att* site coordinates, element size, *mec*/*ccr* gene content, allotype classifications and typing method e.g. `sccmec` or `wgs`.

## Troubleshooting

### Common Issues

**No *att* sites found:**

May be genuine due to non-matching patterns, however:

- Verify that the input FASTA file is properly formatted
- If using GFF mode, ensure the GFF3 file corresponds to the same genome assembly

**No *rlmH* gene found:**

- In FASTA-only mode, *rlmH* is detected via BLAST.  If your species has a highly divergent *rlmH*, provide a custom reference with `--rlmh-ref`
- In GFF mode, verify that gene annotation was performed correctly and *rlmH* is annotated as such
 	- The tool expects either the gene annotated as *rlmH* or the product is named "Ribosomal RNA large subunit methyltransferase H" (case-insenstive)
- *rlmH* is a conserved housekeeping gene — its absence may indicate an assembly issue

**Missing *attR*-*attL* pairs (cross_contig):**

- Caused by assembly fragmentation, the *att* sites are on different contigs
- Check the *att* site TSV output to see which sites were detected and on which contigs
- Long-read sequencing or hybrid assembly can improve extraction rates

**Missing *attR*-*attL* pairs (right_only):**

- Missing *att* pairs will most commonly occur due to a missing *attL*, therefore SCC elements cannot be extracted.
- This is a current limitation of the tool

**BLAST+ not found:**

- BLAST+ is required for FASTA-only mode (`--blast-rlmh`), typing and *ccr* location checks
- Install via conda: `conda install -c bioconda blast`
- BLAST+ is pre-installed in Docker and Singularity containers

**Container-specific issues:**

* See [CONTAINER_GUIDE.md](CONTAINER_GUIDE.md) for troubleshooting Docker and Singularity problems

### Warning Messages

The tools provide informative warning messages to help diagnose issues:

- Missing gene annotations
- Incomplete *att* site pairs
- File processing errors

## Citation

If you use SCC*mec*Extractor in your research, please cite this repository:

```
MacFadyen, A.C. SCCmecExtractor: A toolkit for extracting and typing SCCmec elements
from Staphylococcus and Mammaliicoccus genomes.
GitHub repository: https://github.com/AlisonMacFadyen/SCCmecExtractor
```

## License

[MIT License](https://github.com/AlisonMacFadyen/SCCmecExtractor/tree/main?tab=MIT-1-ov-file#readme)

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Contact

Email: [alison.macfadyen86@gmail.com](mailto:alison.macfadyen86@gmail.com)

## Acknowledgments

- [Bakta](https://github.com/oschwengers/bakta) for bacterial genome annotation
- The Biopython project for sequence manipulation tools
- [BLAST+](https://doi.org/10.1016/s0022-2836(05)80360-2) for sequence comparison/alignment
