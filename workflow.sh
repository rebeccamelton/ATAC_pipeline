#!/bin/bash

# Get script directory for relative paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
SCRIPTS_SOURCE_DIR="${SCRIPT_DIR}/scripts"

# Docker image names
CELLRANGER_IMAGE="cellranger:latest"
CLUSTERING_IMAGE="snatac-clustering:latest"

# Default values
CORES=24
TOTAL_MEMORY=200    # Total memory for cellranger (GB)
WORKING_DIR=""
START_PHASE="cellranger"  # Default start phase

# Function to show usage
show_usage() {
    echo "Usage: $0 COMMAND -w WORKING_DIR [options]"
    echo
    echo "Commands:"
    echo "  setup     Create directory structure and copy scripts"
    echo "  run       Run complete analysis pipeline"
    echo
    echo "Required for all commands:"
    echo "  -w DIR    Working directory where analysis will be performed"
    echo
    echo "Options for run:"
    echo "  -s NAME   Sample name"
    echo "  -c CORES  Number of cores (default: 24)"
    echo "  -t MEM    Total memory in GB for cellranger (default: 200)"
    echo "  -p PHASE  Start from specific phase (cellranger|clustering)"
    echo "               Default: cellranger"
    echo "  -T TMPDIR TMP directory where analysis will be performed"
    echo
    echo "Example:"
    echo "  $0 setup -w /path/to/010825_test"
    echo "  $0 run -w /path/to/010825_test -s sample_name -c 24 -t 200"
    echo "  $0 run -w /path/to/010825_test -s sample_name -p clustering  # Start from clustering phase"
}

# Function to log messages
log_message() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] $1" | tee -a "${LOG_FILE}"
}

# Function to process common arguments
process_args() {
    local OPTIND
    while getopts "w:s:c:t:p:x:h:" flag; do
        case ${flag} in
            w) WORKING_DIR=${OPTARG};;
            s) SAMPLE_NAME=${OPTARG};;
            c) CORES=${OPTARG};;
            t) TOTAL_MEMORY=${OPTARG};;
            p) START_PHASE=${OPTARG};;
            x) TMP_DIR=${OPTARG};;
            h) show_usage; exit 0;;
            ?) show_usage; exit 1;;
        esac
    done

    if [ -z "$WORKING_DIR" ]; then
        echo "Error: Working directory (-w) must be specified"
        show_usage
        exit 1
    fi

    # Convert working directory to absolute path
    WORKING_DIR=$(realpath "$WORKING_DIR")
}

# Function to set up directory structure
setup_directories() {
    echo "==================================================================="
    echo "          Setting up snATAC-seq Analysis Directory Structure         "
    echo "==================================================================="
    echo "Working directory: ${WORKING_DIR}"

    # Create directories and set permissions
    local dirs=(
        "${WORKING_DIR}/01_code"
        "${WORKING_DIR}/02_fastq"
        "${WORKING_DIR}/03_cellranger"
        "${WORKING_DIR}/04_clustering"
        "${WORKING_DIR}/logs"
        "${WORKING_DIR}/tmp"
    )

    for dir in "${dirs[@]}"; do
        if [ ! -d "$dir" ]; then
            echo "Creating directory: $dir"
            mkdir -p "$dir"
            chmod 777 "$dir"
        else
            echo "Directory already exists: $dir"
            chmod 777 "$dir"
        fi
    done

    # Copy scripts from scripts directory to 01_code
    echo
    echo "Copying scripts to ${WORKING_DIR}/01_code..."
    if [ -d "$SCRIPTS_SOURCE_DIR" ]; then
        cp -v "${SCRIPTS_SOURCE_DIR}"/* "${WORKING_DIR}/01_code/"
        # Make all shell scripts executable
        find "${WORKING_DIR}/01_code" -type f -name "*.sh" -exec chmod +x {} \;
        echo "Scripts copied successfully!"
    else
        echo "Warning: Scripts directory not found at: ${SCRIPTS_SOURCE_DIR}"
        echo "Please ensure the 'scripts' directory exists alongside this script."
    fi

    echo
    echo "Directory structure created successfully!"
    echo
    echo "Next steps:"
    echo "1. Place your FASTQ files in: ${WORKING_DIR}/02_fastq"
    echo "2. Run analysis using: $0 run -w ${WORKING_DIR} -b <sample_name>"
    echo
    echo "For more information, refer to the README file."
    echo "==================================================================="
}

# Function to run the complete pipeline
run_complete_pipeline() {
    # Validate required parameters
    if [ -z "$SAMPLE_NAME" ]; then
        echo "Error: Sample name (-b) must be specified"
        show_usage
        exit 1
    fi

    # Validate start phase
    case "$START_PHASE" in
        cellranger|clustering) ;;
        *)
            echo "Error: Invalid phase specified. Must be one of: cellranger, clustering"
            exit 1
            ;;
    esac

    # Set up paths relative to working directory
    local FASTQ_DIR="${WORKING_DIR}/02_fastq"
    local CELLRANGER_OUT="${WORKING_DIR}/03_cellranger"
    local CLUSTERING_OUT="${WORKING_DIR}/04_clustering"
    local LOG_DIR="${WORKING_DIR}/logs"
    local TMP_DIR="${WORKING_DIR}/tmp"

    # Initialize logging
    local timestamp=$(date '+%Y%m%d_%H%M%S')
    local LOG_FILE="${LOG_DIR}/pipeline_${timestamp}.log"
    mkdir -p "$LOG_DIR"
    chmod 777 "$LOG_DIR"

    log_message "Starting pipeline analysis from phase: ${START_PHASE}"
    log_message "Working directory: ${WORKING_DIR}"
    log_message "Sample name: ${SAMPLE_NAME}"
    log_message "Cores: ${CORES}"
    log_message "Total memory for cellranger: ${TOTAL_MEMORY}GB"

    # 1. Run Cellranger analysis
    if [ "$START_PHASE" = "cellranger" ]; then
        log_message "Starting Cellranger analysis..."
        docker run --rm \
            -v "${FASTQ_DIR}:/fastq:ro" \
            -v "${CELLRANGER_OUT}:/output" \
            ${CELLRANGER_IMAGE} cellranger \
            -s ${SAMPLE_NAME} \
            -f /fastq \
            -c ${CORES} \
            -m ${TOTAL_MEMORY} \
            -o /output

        if [ $? -ne 0 ]; then
            log_message "Error: Cellranger analysis failed"
            exit 1
        fi
        log_message "Cellranger analysis completed successfully"
    fi

     # 2. Run Clustering analysis
     if [[ "$START_PHASE" = "cellranger" || "$START_PHASE" = "clustering" ]]; then
         log_message "Starting Clustering analysis..."
         docker run --rm \
             -v "${CELLRANGER_OUT}:/03_cellranger" \
             -v "${CLUSTERING_OUT}:/04_clustering" \
             ${CLUSTERING_IMAGE} clustering \
             -s ${SAMPLE_NAME} \
             -i /03_cellranger/ \
             -o /04_clustering/  \
             --amulet

         if [ $? -ne 0 ]; then
             log_message "Error: Clustering analysis failed"
             exit 1
         fi
         log_message "Clustering analysis completed successfully"
     fi

    log_message "Complete pipeline finished successfully"
}

# Main script logic
case "$1" in
    "setup")
        shift  # Remove 'setup' from arguments
        process_args "$@"
        setup_directories
        ;;
    "run")
        shift  # Remove 'run' from arguments
        process_args "$@"
        run_complete_pipeline
        ;;
    *)
        show_usage
        exit 1
        ;;
esac