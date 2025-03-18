#!/bin/bash
### bash script for snATAC cellranger docker

# Default values
cores=24
mem=200
ref="/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/"

# Parse command line arguments
while getopts "b:f:r:c:m:o:" flag
do
    case "${flag}" in
        b) sample=${OPTARG};;      # Sample name
        f) fastqDir=${OPTARG};;    # FASTQ directory
        r) ref=${OPTARG};;         # Reference directory
        c) cores=${OPTARG};;       # Number of cores
        m) mem=${OPTARG};;         # Memory in GB
        o) output=${OPTARG}
    esac
done

# Check if required parameters are provided
if [ -z "$sample" ] || [ -z "$fastqDir" ]; then
    echo "Error: Sample name (-b) and FASTQ directory (-f) are required"
    echo "Usage: $0 -b <sample_name> -f <fastq_dir> [-r <reference_dir>] [-c <cores>] [-m <memory_GB>]"
    exit 1
fi

cd ${output}

# Run cellranger-atac command
cellranger-atac count \
    --id ${sample} \
    --fastqs ${fastqDir} \
    --sample ${sample} \
    --reference ${ref} \
    --localcores ${cores} \
    --disable-ui \
    --localmem ${mem}