#!/bin/bash

# Function to print usage
usage() {
    echo "Usage: docker run [docker-options] image-name COMMAND [options]"
    echo
    echo "Commands:"
    echo " cellranger Run cellranger ATAC analysis"
    echo " clustering Run clustering analysis"
    echo " python Run Python command"
    echo " R Run R command"
    echo " version Show version information"
    echo " shell Start an interactive shell"
    echo
    echo "Options for cellranger:"
    echo " -b SAMPLE Sample name"
    echo " -f DIR FASTQ directory"
    echo " -c CORES Number of cores (default: 24)"
    echo " -m MEM Memory in GB (default: 200)"
}

# Function to show version information
show_versions() {
    echo "Installed software versions:"
    echo "Python: $(python3 --version)"
    echo "R: $(R --version | head -n 1)"
    echo "Cellranger ATAC: $(cellranger-atac --version)"
}

# Check if command is provided
if [ $# -eq 0 ]; then
    usage
    exit 1
fi

# Get the command
CMD=$1
shift

# Execute appropriate script based on command
case "$CMD" in
    "cellranger")
        exec /scripts/runCellranger.sh "$@"
        ;;
    "clustering")
        exec /scripts/ClusteringScript.R "$@"
        ;;
    "python")
        exec python3 "$@"
        ;;
    "R")
        exec R "$@"
        ;;
    "version")
        show_versions
        ;;
    "shell")
        exec /bin/bash "$@"
        ;;
    *)
        echo "Error: Unknown command '$CMD'"
        usage
        exit 1
        ;;
esac