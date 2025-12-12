#!/bin/bash
set -e

# -----------------------------
# Activate micromamba environment
# -----------------------------
# Ensure micromamba environment is loaded
eval "$(micromamba shell hook -s bash)"
micromamba activate base

# -----------------------------
# Check for Bakta database if needed
# -----------------------------
# Only warn if user is trying to run Bakta
if [ "$1" = "bakta" ]; then
    if [ -z "${BAKTA_DB}" ]; then
        export BAKTA_DB=/data/bakta_db
    fi

    if [ ! -d "$BAKTA_DB" ] || [ -z "$(ls -A "$BAKTA_DB")" ]; then
        echo "WARNING: Bakta database not found at $BAKTA_DB"
        echo "Please mount the Bakta database using: -v /path/to/bakta_db:/data/bakta_db"
        echo ""
        echo "To download the Bakta database, run:"
        echo "  micromamba run -n base bakta_db download --output /path/to/bakta_db"
        echo ""
    fi
fi

# -----------------------------
# Handle no arguments / interactive shell
# -----------------------------
if [ $# -eq 0 ]; then
    echo "SCCmecExtractor Container"
    echo "========================="
    echo ""
    echo "Available commands:"
    echo "  sccmec-locate-att    - Locate attachment sites"
    echo "  sccmec-extract       - Extract SCCmec sequences"
    echo "  bakta                - Run Bakta annotation"
    echo "  bash                 - Interactive shell"
    echo ""
    echo "Example usage:"
    echo "  docker run -v \$PWD:/work -v /path/to/bakta_db:/data/bakta_db sccmecextractor \\"
    echo "    bakta --db /data/bakta_db genome.fna --output bakta_output"
    echo ""
    echo "  docker run -v \$PWD:/work sccmecextractor \\"
    echo "    sccmec-locate-att -f genome.fna -g genome.gff3 -o att_sites.tsv"
    echo ""
    exec /bin/bash
fi

# -----------------------------
# Run user-provided command
# -----------------------------
# Pass all arguments to the shell in the activated environment
exec "$@"
