# TSCC: Tutorial for processing a single sample of 10X single nuclear ATAC data

If needed, download practice data from 10X website ([fastq](https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_5k_nextgem/atac_pbmc_5k_nextgem_fastqs.tar)). 

### (optional) Create your singularity images from docker hub images
 **Not advised since they take up quite a bit of space (~15gb for cellranger image).** <br> 

This process requests more memory then the default TSCC setting, start an interactive job.
Examples code ```srun --pty --nodes=1 --ntasks-per-node=1 --cpus-per-task=4 --mem=2G --account=csd854 -t 2:00:00 -p platinum -q hcp-csd854 --wait 0 /bin/bash```

Activate the correct module and clear any singularity cache:
```
module load singularitypro
singularity cache list
singularity cache clean
```
Create sif files: <br>
```
singularity pull --tmpdir {working_dir}/tmp/ cellranger.sif docker://rlmelton1112/cellranger:latest 

singularity pull --tmpdir {working_dir}/tmp/ snatac-clustering.sif docker://rlmelton1112/snatac-clustering:latest
```
If you do not do this, then you will call the singularity containers from ```/tscc/projects/ps-gaultonlab/rlmelton/github/ATAC_pipeline```
## Running the pipeline

### Step 1: Get necessary files
There's 2 ways to retrieve all the necessary files for running the snATAC pipeline, either git clone the repo or navigate to the lab's copy of the repo `/tscc/projects/ps-gaultonlab/rlmelton/github/ATAC_pipeline`

### Step 2: Set up working directory
In order to make sure everything runs smoothly, we will want our directories to be set up the same. To do this, run:

```bash
sh TSCC_workflow.sh setup -w /path/to/analysis_dir
```
#### example
```bash
sh TSCC_workflow.sh setup -w 042125
```

I would advise making a single directory per project, and run all samples out of that.

### Step 3: Move your files to appropiate directories
If you're planning to run cellranger then all fastqs will need to located in `/path/to/analysis_dir/02_fastq/`

If cellranger has already been run for your samples, then move the full cellranger sample directory in `/path/to/analysis_dir/03_cellranger/`

### Step 4: Run pipeline
#### <u> Basic usage </u>

```bash
sh TSCC_workflow.sh run \
    -w /path/to/analysis_dir \
    -s sample_name \
    -c 24 \           # Number of cores
    -t 200 \          # Total memory for cellranger (GB)
    -m 4              # Memory per core for pipeline (GB)
```
#### <u> Using SLURM on TSCC </u>
In the github repo there is a file `slurm_script_template.sh`, use this to submit a SLURM job to run a single sample

Copy that file and **edit the sample name and working directory!**

Non-gaulton lab users: edit the account and qos flags in `slurm_script_template.sh `

### Understanding the outputs
```
.
├── TSCC_workflow.sh        
├── scripts/            
├── [analysis_dir]/                   
    ├── 01_code/        
    ├── 02_fastq/       
    ├── 03_cellranger/  
    ├── 04_clustering/  
        ├──AMULET_{SAMPLE ID}/        
            ├──{SAMPLE ID}_amulet_doublet_plots.pdf
        ├──log/
        ├──{SAMPLE ID}__Initial_QC_Plots.pdf
        ├──{SAMPLE ID}_initial.RDS
        ├──{SAMPLE ID}_postAmulet.RDS
    └── logs/           
```
* `AMULET_{SAMPLE ID}/` contains all the amulet outputs and plots
* `{SAMPLE ID}_amulet_doublet_plots.pdf` is a pdf containing multiple plots re. AMULET including number of doublets detected and post doublet removal UMAP
* `{SAMPLE ID}__Initial_QC_Plots.pdf` is a pdf containing volin plots for QC metrics and initial clustering UMAP
* `{SAMPLE ID}_initial.RDS` an RDS object post initial QC (only removing barcodes that do not meet the min fragment threshold)
* `{SAMPLE ID}_postAmulet.RDS` an RDS object post Amulet doublet removal (removes any called as a doublet)

### Optional changes to analysis 
* You can skip the cellranger step by using the start phase `-p clustering `
* If you do NOT want AMULET ran, then duplicate `TSCC_workflow.sh` and remove `--amulet` from the clustering command
* A minimum number of ATAC fragments for QC can be set by adding the flag `-f` to the clustering command