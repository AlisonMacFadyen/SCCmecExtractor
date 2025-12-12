# SCCmecExtractor Container Guide

This guide explains how to build and use Docker and Singularity containers for SCCmecExtractor.

## Table of Contents

- [Repository Structure](#repository-structure)
- [Docker Setup](#docker-setup)
- [Singularity Setup](#singularity-setup)
- [Usage Examples](#usage-examples)
- [Troubleshooting](#troubleshooting)

## Repository Structure

Your repository should be organised as follows:

```
SCCmecExtractor/
├── src/
│   └── sccmecextractor/
│       ├── __init__.py
│       ├── locate_att_sites.py
│       └── extract_SCCmec.py
├── tests/
│   └── test_*.py
├── containers/
│   ├── Dockerfile
│   ├── sccmecextractor.def    # Singularity definition
│   ├── environment.yml
│   └── docker-entrypoint.sh
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

# Or build with a specific tag
docker build -t sccmecextractor:1.0.0 -f containers/Dockerfile .
```

### Downloading the Bakta Database

Before using Bakta within the container, you need to download the database:

```bash
# Create a directory for the Bakta database
mkdir -p ~/bakta_db

# Download using Docker
docker run --rm -v ~/bakta_db/:/data/bakta_db \
  sccmecextractor:latest \
  bakta_db download --output /data/bakta_db --type light
```

For the full database (recommended for production):

```bash
docker run --rm -v ~/bakta_db/:/data/bakta_db \
  sccmecextractor:latest \
  bakta_db download --output /data/bakta_db
```

### Running the Docker Container

#### Running Bakta Annotation

```bash
docker run --rm \
  -v $PWD:/work \
  -v ~/bakta_db:/data/bakta_db \
  sccmecextractor:latest \
  bakta --db /data/bakta_db genome.fna --output bakta_output --prefix sample
```

#### Locating Attachment Sites

```bash
docker run --rm \
  -v $PWD:/work \
  sccmecextractor:latest \
  sccmec-locate-att -f genome.fna -g genome.gff3 -o att_sites.tsv
```

#### Extracting SCCmec Sequences

```bash
docker run --rm \
  -v $PWD:/work \
  sccmecextractor:latest \
  sccmec-extract -f genome.fna -g genome.gff3 -a att_sites.tsv -s sccmec_output
```

#### Complete Pipeline

```bash
# 1. Annotate genome with Bakta
docker run --rm \
  -v $PWD:/work \
  -v ~/bakta_db:/data/bakta_db \
  sccmecextractor:latest \
  bakta --db /data/bakta_db genome.fna --output bakta_output --prefix genome

# 2. Locate att sites
docker run --rm \
  -v $PWD:/work \
  sccmecextractor:latest \
  sccmec-locate-att -f genome.fna -g bakta_output/genome.gff3 -o att_sites.tsv

# 3. Extract SCCmec
docker run --rm \
  -v $PWD:/work \
  sccmecextractor:latest \
  sccmec-extract -f genome.fna -g bakta_output/genome.gff3 -a att_sites.tsv -s sccmec_sequences
```

## Singularity Setup

### Prerequisites

- Singularity/Apptainer installed on your system
- Root/sudo access for building (or use `--fakeroot`)

### Building the Singularity Image

#### Option 1: Build from Definition File

```bash
# Build the image (requires sudo or --fakeroot)
sudo singularity build sccmecextractor.sif containers/sccmecextractor.def

# Or with fakeroot (no sudo required)
singularity build --fakeroot sccmecextractor.sif containers/sccmecextractor.def
```

#### Option 2: Build from Docker Image

```bash
# Pull from Docker Hub (if published)
singularity pull docker://yourusername/sccmecextractor:latest

# Or build from local Docker image
singularity build sccmecextractor.sif docker-daemon://sccmecextractor:latest
```

### Running the Singularity Container

#### Running Commands with exec

```bash
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

# Locate attachment sites
singularity exec \
  --bind $PWD:/work \
  sccmecextractor.sif \
  sccmec-locate-att -f genome.fna -g genome.gff3 -o att_sites.tsv

# Extract SCCmec sequences
singularity exec \
  --bind $PWD:/work \
  sccmecextractor.sif \
  sccmec-extract -f genome.fna -g genome.gff3 -a att_sites.tsv -s sccmec_output
```

## Troubleshooting

### Common Issues

#### Issue: "Bakta database not found"

**Solution:** Ensure you've downloaded the Bakta database and mounted it correctly:

```bash
# Docker
-v ~/bakta_db:/data/bakta_db

# Singularity
--bind ~/bakta_db:/data/bakta_db
```

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
  sccmec-locate-att -f data/genome.fna -g data/genome.gff3 -o output/att_sites.tsv
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
