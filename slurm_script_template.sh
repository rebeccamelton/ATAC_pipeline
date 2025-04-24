#!/bin/bash
#SBATCH --job-name=snATAC-seq         # Job name
#SBATCH --output=snATAC_%j.out        # Standard output log (%j expands to jobID)
#SBATCH --error=snATAC_%j.err         # Standard error log
#SBATCH --time=48:00:00               # Request runtime of 48 hours (adjust as needed)
#SBATCH --nodes=1                     # Request 1 node
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --account=csd854              # Account name (adjust as needed)
#SBATCH --partition=platinum             # Partition name (adjust as needed)
#SBATCH --qos=hcp-csd854                    
#SBATCH --cpus-per-task=30            # Request 30 CPUs per task (matches your -c parameter)
#SBATCH --mem=200G                    # Request 200GB memory (adjust based on your -t parameter)
#SBATCH --mail-type=BEGIN,END,FAIL    # Email notifications (optional)
#SBATCH --mail-user=rlmelton@health.ucsd.com  # User email (optional)

# Load any necessary modules
module load singularitypro
singularity cache list
singularity cache clean

# Set working directory (crucial to ensure the job runs in the correct directory)
WORKING_DIR="/tscc/projects/ps-gaultonlab/rlmelton/github/ATAC_pipeline/042125"
cd "${WORKING_DIR}"

# Define parameters
SAMPLE_NAME="MM_12" 
CORES=30
MEMORY=150
PHASE="cellranger"

# Path to workflow script
WORKFLOW_SCRIPT="/tscc/projects/ps-gaultonlab/rlmelton/github/ATAC_pipeline/TSCC_workflow.sh"

# Run the workflow
echo "Starting snATAC-seq workflow at $(date)"
${WORKFLOW_SCRIPT} run -w ${WORKING_DIR} -s ${SAMPLE_NAME} -c ${CORES} -t ${MEMORY} -p ${PHASE}
echo "Completed snATAC-seq workflow at $(date)"

# Optional: generate summary or cleanup
echo "Job completed successfully"