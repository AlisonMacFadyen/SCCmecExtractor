#!/bin/bash
set -e

# Check if Bakta database is mounted/available
if [ ! -d "$BAKTA_DB" ] || [ -z "$(ls -A $BAKTA_DB)" ]; then
    echo "WARNING: Bakta database not found at $BAKTA_DB"
    echo "Please mount the Bakta database using: -v /path/to/bakta_db:/data/bakta_db"
    echo ""
    echo "To download the Bakta database, run:"
    echo "  bakta_db download --output /path/to/bakta_db"
    echo ""
fi

# If no arguments provided, print usage
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
    exec /bin/bash
else
    exec "$@"
fi
