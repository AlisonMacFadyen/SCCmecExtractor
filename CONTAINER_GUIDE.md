# SCCmecExtractor Container Guide

This guide explains how to build and use Docker and Singularity containers for SCCmecExtractor.

## Table of Contents

- [Repository Structure](#repository-structure)
- [Docker Setup](#docker-setup)
- [Singularity Setup](#singularity-setup)
- [Usage Examples](#usage-examples)
- [Troubleshooting](#troubleshooting)

## Repository Structure

```
SCCmecExtractor/
├── src/
│   └── sccmecextractor/
│       ├── __init__.py
│       ├── pipeline.py
│       ├── locate_att_sites.py
│       ├── extract_SCCmec.py
│       ├── type_sccmec.py
│       ├── report_sccmec.py
│       ├── blast_utils.py
│       └── data/
│           ├── mec_genes_allotypes.fasta
│           ├── ccr_genes.fasta
│           └── rlmH.fasta
├── tests/
│   └── test_*.py
├── containers/
│   ├── Dockerfile
│   ├── sccmecextractor.def
│   └── environment.yml
├── docker-entrypoint.sh
├── .dockerignore
├── .gitignore
├── pyproject.toml
├── requirements.txt
├── VERSION
├── LICENSE
├── README.md
└── CONTAINER_GUIDE.md
```

## Docker Setup

### Prerequisites

- Docker installed on your system
- Sufficient disk space (~2-3 GB for the image)

### Building the Docker Image

```bash
# Clone the repository
git clone https://github.com/AlisonMacFadyen/SCCmecExtractor.git
cd SCCmecExtractor

# Build the Docker image
docker build -t sccmecextractor:latest -f containers/Dockerfile .
```

### Downloading the Bakta Database

Only required if using Bakta for GFF annotation.  Not needed for FASTA-only mode.

```bash
# Create a directory for the Bakta database
mkdir -p ~/bakta_db

# Download using Docker
docker run --rm -v ~/bakta_db/:/data/bakta_db \
  sccmecextractor:latest \
  bakta_db download --output /data/bakta_db --type light
```

For the full database:

```bash
docker run --rm -v ~/bakta_db/:/data/bakta_db \
  sccmecextractor:latest \
  bakta_db download --output /data/bakta_db
```

### Running the Docker Container

#### Running the Complete Pipeline (FASTA-only Mode)

The simplest workflow, no GFF annotation or Bakta database required:

```bash
# Single genome
docker run --rm \
  -v $PWD:/work \
  sccmecextractor:latest \
  sccmec-pipeline -f genome.fna -o results/

# Directory of genomes with multithreading
docker run --rm \
  -v $PWD:/work \
  sccmecextractor:latest \
  sccmec-pipeline --fna-dir genomes/ -o results/ -t 4

# With composite element detection
docker run --rm \
  -v $PWD:/work \
  sccmecextractor:latest \
  sccmec-pipeline --fna-dir genomes/ --composite -o results/ -t 4
```

#### Running the Complete Pipeline (GFF Mode)

```bash
docker run --rm \
  -v $PWD:/work \
  -v ~/bakta_db:/data/bakta_db \
  sccmecextractor:latest \
  bash -c "
    bakta --db /data/bakta_db genome.fna --output bakta_output --prefix genome && \
    sccmec-pipeline -f genome.fna -g bakta_output/genome.gff3 -o results/
  "
```

#### Running Individual Commands

```bash
# Locate attachment sites
docker run --rm \
  -v $PWD:/work \
  sccmecextractor:latest \
  sccmec-locate-att -f genome.fna -o att_sites.tsv

# Extract SCC elements
docker run --rm \
  -v $PWD:/work \
  sccmecextractor:latest \
  sccmec-extract -f genome.fna -a att_sites.tsv -s sccmec_output -r extraction_report.tsv

# Gene level typing extracted elements
docker run --rm \
  -v $PWD:/work \
  sccmecextractor:latest \
  sccmec-type -f sccmec_output/ -o typing_results.tsv

# Generate unified report
docker run --rm \
  -v $PWD:/work \
  sccmecextractor:latest \
  sccmec-report -e extraction_report.tsv -t typing_results.tsv -o unified_report.tsv

# Run Bakta annotation
docker run --rm \
  -v $PWD:/work \
  -v ~/bakta_db:/data/bakta_db \
  sccmecextractor:latest \
  bakta --db /data/bakta_db/db-light genome.fna.gz --output bakta_output --prefix sample
```

## Singularity Setup

### Prerequisites

- Singularity/Apptainer installed on your system
- Root/sudo access for building (or use `--fakeroot`)

### Building the Singularity Image

#### Option 1: Pull from Docker Hub (Recommended)

```bash
singularity pull docker://alisonmacfadyen/sccmecextractor:latest
```

#### Option 2: Build from Definition File

The definition file pulls the published Docker image, so the result is identical:

```bash
# Build the image (requires sudo or --fakeroot)
sudo singularity build sccmecextractor.sif containers/sccmecextractor.def

# Or with fakeroot (no sudo required)
singularity build --fakeroot sccmecextractor.sif containers/sccmecextractor.def
```

#### Option 3: Build from Local Docker Image

```bash
singularity build sccmecextractor.sif docker-daemon://sccmecextractor:latest
```

### Running the Singularity Container

#### Running the Complete Pipeline (FASTA-only Mode)

```bash
# Single genome
singularity exec \
  --bind $PWD:/work \
  sccmecextractor.sif \
  sccmec-pipeline -f genome.fna -o results/

# Directory of genomes with multithreading
singularity exec \
  --bind $PWD:/work \
  sccmecextractor.sif \
  sccmec-pipeline --fna-dir genomes/ -o results/ -t 8
```

#### Running the Complete Pipeline (GFF Mode)

```bash
singularity exec \
  --bind $PWD:/work \
  --bind ~/bakta_db:/data/bakta_db \
  sccmecextractor.sif \
  bash -c "
    bakta --db /data/bakta_db genome.fna --output bakta_output --prefix genome && \
    sccmec-pipeline -f genome.fna -g bakta_output/genome.gff3 -o results/
  "
```

#### Running Individual Commands

```bash
# Locate attachment sites
singularity exec \
  --bind $PWD:/work \
  sccmecextractor.sif \
  sccmec-locate-att -f genome.fna -o att_sites.tsv

# Extract SCC elements
singularity exec \
  --bind $PWD:/work \
  sccmecextractor.sif \
  sccmec-extract -f genome.fna -a att_sites.tsv -s sccmec_output -r extraction_report.tsv

# Gene level typing of extracted elements
singularity exec \
  --bind $PWD:/work \
  sccmecextractor.sif \
  sccmec-type -f sccmec_output/ -o typing_results.tsv

# Generate unified report
singularity exec \
  --bind $PWD:/work \
  sccmecextractor.sif \
  sccmec-report -e extraction_report.tsv -t typing_results.tsv -o unified_report.tsv

# Download Bakta Database
singularity exec \
  --bind ~/bakta_db:/data/bakta_db \
  sccmecextractor.sif \
  bakta_db download --output /data/bakta_db --type light

# Run Bakta annotation
singularity exec \
  --bind $PWD:/work \
  --bind ~/bakta_db:/data/bakta_db \
  sccmecextractor.sif \
  bakta --db /data/bakta_db genome.fna --output bakta_output --prefix sample
```

## Troubleshooting

### Common Issues

#### Issue: "Bakta database not found"

**Solution:** This only applies to the GFF workflow. Ensure you've downloaded the Bakta database and mounted it correctly:

```bash
# Docker
-v ~/bakta_db:/data/bakta_db

# Singularity
--bind ~/bakta_db:/data/bakta_db
```

For FASTA-only mode, no Bakta database is needed.

#### Issue: "Permission denied" errors

**Docker Solution:**

- Ensure the mounted directories have appropriate permissions
- The container runs as a non-root user by default

**Singularity Solution:**

- Check file permissions in bound directories
- Singularity runs as your user by default, so permissions should match

#### Issue: "No such file or directory" when accessing files

**Solution:** Make sure you're mounting the correct directory and using paths relative to `/work` inside the container:

```bash
# If your genome is in $PWD/data/genome.fna
docker run --rm -v $PWD:/work sccmecextractor:latest \
  sccmec-pipeline -f data/genome.fna -o output/results/
```

#### Issue: Container build fails

**Docker Solution:**

- Check Docker daemon is running
- Ensure sufficient disk space
- Try pruning old images: `docker system prune -a`

**Singularity Solution:**

- Ensure SINGULARITY_TMPDIR has sufficient space
- Try: `export SINGULARITY_TMPDIR=/path/to/large/tmp`
- Use `--fakeroot` if you don't have sudo access

### Getting Help

If you encounter issues not covered here:

1. Check the GitHub issues: https://github.com/AlisonMacFadyen/SCCmecExtractor/issues
2. Review the Bakta documentation: https://github.com/oschwengers/bakta
3. Open a new issue with:
   - Container type (Docker/Singularity)
   - Error message
   - Command used
   - System information

## Additional Resources

- [Docker Documentation](https://docs.docker.com/)
- [Singularity Documentation](https://docs.sylabs.io/guides/latest/user-guide/)
- [Bakta GitHub](https://github.com/oschwengers/bakta)
- [SCCmecExtractor GitHub](https://github.com/AlisonMacFadyen/SCCmecExtractor)
