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
docker build -t sccmecextractor:latest .

# Or build with a specific tag
docker build -t sccmecextractor:1.0.0 .
```

### Downloading the Bakta Database

Before using Bakta within the container, you need to download the database:

```bash
# Create a directory for the Bakta database
mkdir -p ~/bakta_db

# Download using Docker
docker run --rm -v ~/bakta_db:/data/bakta_db \
  oschwengers/bakta:latest \
  bakta_db download --output /data/bakta_db --type light
```

For the full database (recommended for production):

```bash
docker run --rm -v ~/bakta_db:/data/bakta_db \
  oschwengers/bakta:latest \
  bakta_db download --output /data/bakta_db
```

### Running the Docker Container

#### Interactive Mode

```bash
docker run -it --rm \
  -v $PWD:/work \
  -v ~/bakta_db:/data/bakta_db \
  sccmecextractor:latest bash
```

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

### Publishing to Docker Hub

```bash
# Tag the image
docker tag sccmecextractor:latest yourusername/sccmecextractor:latest
docker tag sccmecextractor:latest yourusername/sccmecextractor:1.0.0

# Login to Docker Hub
docker login

# Push the image
docker push yourusername/sccmecextractor:latest
docker push yourusername/sccmecextractor:1.0.0
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

#### Interactive Shell

```bash
singularity shell \
  --bind $PWD:/work \
  --bind ~/bakta_db:/data/bakta_db \
  sccmecextractor.sif
```

#### Running Commands with exec

```bash
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

#### Complete Pipeline with Singularity

```bash
# Set up bind paths
export WORK_DIR=$PWD
export BAKTA_DB=~/bakta_db

# 1. Annotate with Bakta
singularity exec \
  --bind $WORK_DIR:/work \
  --bind $BAKTA_DB:/data/bakta_db \
  sccmecextractor.sif \
  bakta --db /data/bakta_db genome.fna --output bakta_output --prefix genome

# 2. Locate att sites
singularity exec \
  --bind $WORK_DIR:/work \
  sccmecextractor.sif \
  sccmec-locate-att -f genome.fna -g bakta_output/genome.gff3 -o att_sites.tsv

# 3. Extract SCCmec
singularity exec \
  --bind $WORK_DIR:/work \
  sccmecextractor.sif \
  sccmec-extract -f genome.fna -g bakta_output/genome.gff3 -a att_sites.tsv -s sccmec_sequences
```

### HPC/Cluster Usage

For running on HPC systems with job schedulers:

#### SLURM Example

```bash
#!/bin/bash
#SBATCH --job-name=sccmec_extraction
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00

module load singularity

# Define paths
WORK_DIR=/path/to/working/directory
BAKTA_DB=/path/to/bakta_db
GENOME=genome.fna

cd $WORK_DIR

# Run Bakta
singularity exec \
  --bind $WORK_DIR:/work \
  --bind $BAKTA_DB:/data/bakta_db \
  sccmecextractor.sif \
  bakta --db /data/bakta_db --threads $SLURM_CPUS_PER_TASK \
  $GENOME --output bakta_output --prefix $GENOME

# Run SCCmecExtractor
singularity exec \
  --bind $WORK_DIR:/work \
  sccmecextractor.sif \
  sccmec-locate-att -f $GENOME -g bakta_output/${GENOME}.gff3 -o att_sites.tsv

singularity exec \
  --bind $WORK_DIR:/work \
  sccmecextractor.sif \
  sccmec-extract -f $GENOME -g bakta_output/${GENOME}.gff3 \
  -a att_sites.tsv -s sccmec_output
```

## Usage Examples

### Processing Multiple Genomes

#### Docker

```bash
# Create a script to process multiple genomes
for genome in genomes/*.fna; do
  basename=$(basename $genome .fna)
  
  docker run --rm \
    -v $PWD:/work \
    -v ~/bakta_db:/data/bakta_db \
    sccmecextractor:latest \
    bakta --db /data/bakta_db $genome --output bakta_${basename} --prefix ${basename}
  
  docker run --rm \
    -v $PWD:/work \
    sccmecextractor:latest \
    sccmec-locate-att -f $genome -g bakta_${basename}/${basename}.gff3 \
    -o att_sites_${basename}.tsv
  
  docker run --rm \
    -v $PWD:/work \
    sccmecextractor:latest \
    sccmec-extract -f $genome -g bakta_${basename}/${basename}.gff3 \
    -a att_sites_${basename}.tsv -s sccmec_${basename}
done
```

#### Singularity

```bash
for genome in genomes/*.fna; do
  basename=$(basename $genome .fna)
  
  singularity exec \
    --bind $PWD:/work \
    --bind ~/bakta_db:/data/bakta_db \
    sccmecextractor.sif \
    bakta --db /data/bakta_db $genome --output bakta_${basename} --prefix ${basename}
  
  singularity exec --bind $PWD:/work sccmecextractor.sif \
    sccmec-locate-att -f $genome -g bakta_${basename}/${basename}.gff3 \
    -o att_sites_${basename}.tsv
  
  singularity exec --bind $PWD:/work sccmecextractor.sif \
    sccmec-extract -f $genome -g bakta_${basename}/${basename}.gff3 \
    -a att_sites_${basename}.tsv -s sccmec_${basename}
done
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
