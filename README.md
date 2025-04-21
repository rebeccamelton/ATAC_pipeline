# Single-cell ATAC-seq Analysis Pipeline

## Overview
This pipeline provides an automated workflow for processing single-cell ATAC-seq data through three main phases:
1. Cellranger ATAC processing
2. Window-based clustering using Seurat
---
# TSCC
Directions for running this pipeline on TSCC's HPC system using slurm and singularity


---
# Gaulton Lab servers
## Prerequisites
- Clone this repo 
- Docker installed and running
- Pull Docker images from Dockerhub:
  - `docker pull rlmelton1112/cellranger:latest`
  - `docker pull rlmelton1112/snatac-clustering:latest`
- Sufficient disk space for analysis
- Minimum 16GB RAM recommended
- Multi-core processor recommended

## Directory Structure
```
.
├── workflow.sh         # Main workflow script
├── scripts/            # Analysis scripts directory
├── [analysis_dir]/     # Your analysis directory containing:
    ├── 01_code/        # Analysis scripts (copied automatically)
    ├── 02_fastq/       # Raw FASTQ input files
    ├── 03_cellranger/  # Cellranger analysis output
    ├── 04_clustering/  # Clustering analysis results
    └── logs/           # Log files from all steps
```

## Quick Start

### 1. Setup Directory Structure
```bash
./workflow.sh setup -w /path/to/analysis_dir
```

### 2. Place FASTQ Files
Put (either copy or move) your FASTQ files into the `02_fastq` directory:
```bash
cp /path/to/fastq/files/* /path/to/analysis_dir/02_fastq/
```

### 3. Run Analysis
```bash
./workflow.sh run \
    -w /path/to/analysis_dir \
    -b sample_name \
    -c 24 \           # Number of cores
    -t 200 \          # Total memory for cellranger (GB)
    -m 4              # Memory per core for pipeline (GB)
```

## Detailed Usage

### Setup Command
```bash
./workflow.sh setup -w PATH
```
Options:
- `-w PATH`: Working directory where analysis will be performed (required)

### Run Command
```bash
./workflow.sh run [options]
```
Required options:
- `-w PATH`: Working directory where analysis will be performed
- `-b NAME`: Sample name

Optional parameters:
- `-c CORES`: Number of cores to use (default: 24)
- `-t MEM`: Total memory in GB for cellranger (default: 200)
- `-m MEM`: Memory per core in GB for pipeline (default: 4)
- `-p PHASE`: Start from specific phase (cellranger|pipeline|clustering)

## Resuming Failed Runs
If the pipeline fails at any point, you can resume from a specific phase using the `-p` option:

```bash
# Resume from pipeline phase
./workflow.sh run -w /path/to/analysis_dir -b sample_name -p pipeline

# Resume from clustering phase
./workflow.sh run -w /path/to/analysis_dir -b sample_name -p clustering
```

## Output Structure

### Cellranger Output (03_cellranger/)
- Contains the standard 10x Genomics Cellranger ATAC output
- Key files include aligned BAM files and initial analyses
  
### Clustering Output (04_clustering/)
- Final clustering results
- Visualization files
- Analysis summaries

### Logs (logs/)
- Detailed logs from each step of the pipeline
- Timestamped for tracking multiple runs
- Includes error messages and progress information

## Troubleshooting

### Common Issues

1. "No such file or directory" for FASTQ files
   - Verify FASTQ files are in the 02_fastq directory
   - Check file permissions

2. Memory Issues
   - Adjust -t parameter for cellranger memory
   - Adjust -m parameter for pipeline memory per core

3. Docker Issues
   - Ensure Docker is running
   - Verify both Docker images are built and available

### Log Files
Check the logs directory for detailed error messages:
```bash
ls -l /path/to/analysis_dir/logs/
```

## Resource Requirements

### Minimum Requirements
- RAM: 16GB
- CPU: 8 cores
- Storage: 100GB free space

### Recommended Requirements
- RAM: 32GB+
- CPU: 24 cores
- Storage: 500GB+ free space

## Support
For issues and questions, please contact Rebecca Melton (rlmelton@health.ucsd.edu)

## Features to add
[] Additional QC plots to the clustering rscript  
[] HTML report with all QC plots and quick metrics  
[] AMULET docker, script, and doublet removal for clustering  
[] Store final in s3 bucket   
[] Provide a detailed BP cells notebook for merging multiple atac samples together   
