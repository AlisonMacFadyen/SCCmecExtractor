# SCCmecExtractor

A Python toolkit for extracting SCC*mec* (Staphylococcal Cassette Chromosome *mec*) sequences from *Staphylococcus* whole genome sequences.  This tool identifies attachment (*att*) sites and extracts the complete SCC*mec* element based on genomic context.

**Note the tool is quite stringent and requires the *att* sites to be located on the same contig as each other and as the gene *rlmH* in order to extract the DNA sequence of the SCC*mec***

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![PyPI](https://img.shields.io/pypi/v/:sccmecextractor)](https://pypi.org/project/sccmecextractor/1.1.0/)
[![Docker Image Version](https://img.shields.io/docker/alisonmacfadyen/sccmecextractor?sort=semver&label=docker))](https://hub.docker.com/repository/docker/alisonmacfadyen/sccmecextractor/general)

## Overview

SCCmecExtractor consists of two main scripts that work together to identify and extract SCC*mec* sequences:

1. **`locate_att_sites.py`** - Identifies attachment sites in genomic sequences
    - Canonical *attR* sites: *attR* and the complement, *cattR*
    - Divergent CcrC associated *attR2* and the complement, *cattR2*
    - Canonical *attL* and the complement, *cattL*
    - Divergent CcrC associated *attL2* and the complement, *cattL2*
2. **`extract_SCCmec.py`** - Extracts the SCC*mec* sequence based on identified *att* sites and gene annotations

## Table of Contents

- [Installation](#installation)
  - [Using Conda/Mamba](#using-condamamba-recommended)
  - [Using pip](#using-pip)
  - [Using Docker](#using-docker)
  - [Using Singularity](#using-singularity)
- [Requirements](#requirements)
- [Usage](#usage)
- [Complete Workflow](#complete-workflow-example)
- [How It Works](#how-it-works)
- [Output Format](#output-format)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)
- [Contributing](#contributing)
- [Contact](#contact)

## Installation

### Using Conda/Mamba

```bash
# Create a new environment
conda create -n sccmecextractor python=3.11
conda activate sccmecextractor

# Install dependencies
conda install -c conda-forge -c bioconda biopython bakta

# Install SCCmecExtractor
pip install sccmecextractor

# Test that scripts are available
sccmec-locate-att --help
sccmec-extract --help
```

If scripts do not run, make sure the environment’s `bin/` directory is in your PATH:

```bash
export PATH="$CONDA_PREFIX/bin:$PATH"
```

### Using pip

Note, installation with `pip` does not provide Bakta.

```bash
# Install SCCmecExtractor
pip install sccmecextractor

# Test that scripts are available
sccmec-locate-att --help
sccmec-extract --help
```

### Using Docker

Docker provides a containerised environment with all dependencies pre-installed, including Bakta.

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
# Download Bakta Database (light in this example)

# Create a directory for the Bakta database
mkdir -p ~/bakta_db

# Download using Docker
docker run --rm -v ~/bakta_db/:/data/bakta_db \
  sccmecextractor:latest \
  bakta_db download --output /data/bakta_db --type light

# Run the complete pipeline
docker run --rm \
  -v $PWD:/work \
  -v ~/bakta_db:/data/bakta_db \
  sccmecextractor:latest \
  bash -c "bakta --db /data/bakta_db genome.fna.gz --output bakta_out && \
           sccmec-locate-att -f genome.fna -g bakta_out/genome.gff3 -o att_sites.tsv && \
           sccmec-extract -f genome.fna -g bakta_out/genome.gff3 -a att_sites.tsv -s output"
```

See [CONTAINER_GUIDE.md](CONTAINER_GUIDE.md) for detailed Docker usage instructions.

### Using Singularity

Singularity/Apptainer is ideal for HPC environments where Docker is not available.

```bash
# Build from definition file
singularity build sccmecextractor.sif containers/sccmecextractor.def

# Or pull from Docker Hub
singularity pull docker://alisonmacfadyen/sccmecextractor:latest
```

**Quick Start with Singularity:**

```bash
# Download Bakta Database (light in this example)
singularity exec \
  --bind $PWD:/work \
  sccmecextractor.sif \
  bakta_db download --output ~/bakta_db --type light

# Run the complete pipeline
singularity exec \
  --bind $PWD:/work \
  --bind ~/bakta_db:/data/bakta_db \
  sccmecextractor.sif \
  bash -c "bakta --db /data/bakta_db genome.fna --output bakta_out && \
           sccmec-locate-att -f genome.fna -g bakta_out/genome.gff3 -o att_sites.tsv && \
           sccmec-extract -f genome.fna -g bakta_out/genome.gff3 -a att_sites.tsv -s output"
```

See [CONTAINER_GUIDE.md](CONTAINER_GUIDE.md) for detailed Singularity usage instructions.

## Requirements

### Dependencies

* Python 3.9+
* Biopython (`pip install biopython`)
* Bakta (for genome annotation) - automatically included in containers

### Input Files

* **Genome sequence**: `.fasta` or `.fna` file containing the assembled genome.  Note a compressed version is required to run Bakta.
* **Gene annotations**: `.gff3` file with gene annotations (we recommend using [bakta](https://github.com/oschwengers/bakta) for annotation)

### Bakta Database

If using Bakta for annotation, you'll need to download the Bakta database:

```bash
# Light database (faster, smaller)
bakta_db download --output bakta_db --type light

# Full database
bakta_db download --output bakta_db
```

## Usage

### Step 1: Locate Attachment Sites

First, identify *att* sites in your genome:

```bash
sccmec-locate-att -f genome.fna -g genome.gff3 -o att_sites.tsv
```

Or using the Python script directly:

```bash
python src/sccmecextractor/locate_att_sites.py -f genome.fna -g genome.gff3 -o att_sites.tsv
```

**Parameters:**
- `-f, --fna`: Input genome file (.fasta or .fna)
- `-g, --gff`: Gene annotation file (.gff3 format)
- `-o, --outfile`: Output TSV file containing *att* site locations

**Output:**
The script generates a TSV file with the following columns:
- Input\_File
- Pattern (*attR*, *attL*, *cattR*, *cattL*, *attR2*, *cattR2*)
- Contig
- Start position
- End position
- Matching\_Sequence

### Step 2: Extract SCC*mec* Sequences

Extract the SCC*mec* sequence using the identified *att* sites:

```bash
python extract_SCCmec.py -f genome.fna -g genome.gff3 -a att_sites.tsv -s output_directory
```

Or using the Python script directly:

```bash
python src/sccmecextractor/extract_SCCmec.py -f genome.fna -g genome.gff3 -a att_sites.tsv -s output_directory
```


**Parameters:**
- `-f, --fna`: Input genome file (.fasta or .fna)
- `-g, --gff`: Gene annotation file (.gff3 format)
- `-a, --att`: TSV file from step 1 containing att site locations
- `-s, --sccmec`: Output directory for extracted SCCmec sequences

**Output:**
The script creates a FASTA file named `{genome}_SCCmec.fasta` in the specified output directory containing the extracted SCC*mec* sequence.

## Complete Workflow Example

### Local Installation

```bash
# 1. Annotate your genome with bakta (recommended)
bakta --db bakta_db genome.fna --output bakta_output

# 2. Locate att sites
sccmec-locate-att -f genome.fna -g bakta_output/genome.gff3 -o att_sites.tsv

# 3. Extract SCCmec sequence
sccmec-extract -f genome.fna -g bakta_output/genome.gff3 -a att_sites.tsv -s sccmec_output
```

### Docker

```bash
# Complete pipeline in one command
docker run --rm \
  -v $PWD:/work \
  -v ~/bakta_db:/data/bakta_db \
  sccmecextractor:latest \
  bash -c "
    bakta --db /data/bakta_db genome.fna --output bakta_output --prefix genome && \
    sccmec-locate-att -f genome.fna -g bakta_output/genome.gff3 -o att_sites.tsv && \
    sccmec-extract -f genome.fna -g bakta_output/genome.gff3 -a att_sites.tsv -s sccmec_output
  "
```

### Singularity

```bash
# Complete pipeline in one command
singularity exec \
  --bind $PWD:/work \
  --bind ~/bakta_db:/data/bakta_db \
  sccmecextractor.sif \
  bash -c "
    bakta --db /data/bakta_db genome.fna --output bakta_output --prefix genome && \
    sccmec-locate-att -f genome.fna -g bakta_output/genome.gff3 -o att_sites.tsv && \
    sccmec-extract -f genome.fna -g bakta_output/genome.gff3 -a att_sites.tsv -s sccmec_output
  "
```

## How It Works

I hope to publish this tool someday but until then here is an overview of how the tool performs its functions.

### Attachment Site Detection
The tool searches for specific DNA motifs that represent attachment sites:

- **attR**: Right attachment site patterns
- **attL**: Left attachment site patterns  
- **cattR**: Complementary right attachment sites
- **cattL**: Complementary left attachment sites
- **attR2/cattR2**: Alternative right attachment site patterns
- **attL2/cattL2**: Alternative left attachment site patterns

The script uses regex patterns with degeneracy to account for sequence variation in these sites.

### SCCmec Extraction Logic
1. **Site Validation**: Identifies the closest *attR*-*attL* pair on the same contig
2. **Gene Context**: Locates the *rlmH* gene, which is used as a reference point
3. **Coordinate Determination**: Calculates extraction coordinates based on *rlmH* position and *att* sites
4. **Sequence Extraction**: Extracts the region between *att* sites with appropriate padding
5. **Orientation Handling**: Automatically handles reverse complement extraction when necessary

### Key Features
- **Intelligent Filtering**: *attR* and *attR2* sites are only considered if they fall within *rlmH* genes
- **Distance Optimisation**: Selects the closest *attR*-*attL* pair to minimise extraction of non-SCC*mec* sequences
- **Strand Awareness**: Automatically detects and handles SCC*mec* elements on reverse strands
- **Quality Control**: Validates presence of required genes and *att* sites before extraction

## Output Format

The extracted SCC*mec* sequence is saved as a FASTA file with:
- **ID**: `{input_file}_{contig}_{start}_{end}`
- **Description**: `attR:{right_att_info}_attL:{left_att_info}`

## Troubleshooting

### Common Issues

**No *att* sites found:**
- Check that your genome contains SCC*mec* elements
- Verify that the input FASTA file is properly formatted
- Ensure the GFF3 file corresponds to the same genome assembly

**No *rlmH* gene found:**
- This may indicate there is an issue with your input genome as *rlmH* is a conserved gene for *Staphylococcus*
- Verify that gene annotation was performed correctly
- Check that the GFF3 file contains gene features with proper naming - *rlmH* must be annotated as such

**Missing *attR*-*attL* pairs:**
- Some genomes may have incomplete or atypical SCC*mec* elements
- Check the att_sites.tsv output to see which sites were detected

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
MacFadyen, A.C. SCCmecExtractor: A toolkit for extracting SCCmec sequences from Staphylococcus genomes. 
GitHub repository: https://github.com/AlisonMacFadyen/SCCmecExtractor
```

## Work in Progress

I aim to add in `bakta` annotation as part of the pipeline, as well as to include information on SCC*mec* gene carriage and Typing information.  In the meantime, for typing, I recommend checking out this tool:  [sccmec](https://github.com/rpetit3/sccmec)

If you have any additional ideas, please let me know.

## License

[MIT License](https://github.com/AlisonMacFadyen/SCCmecExtractor/tree/main?tab=MIT-1-ov-file#readme)

## Contributing

Contributions are welcome!  Please feel free to submit issues or pull requests.

## Contact

Email: [alison.macfadyen86@gmail.com](mailto:alison.macfadyen86@gmail.com)

## Acknowledgments

- [Bakta](https://github.com/oschwengers/bakta) for bacterial genome annotation
- The Biopython project for sequence manipulation tools