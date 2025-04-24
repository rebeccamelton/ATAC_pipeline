# TSCC: Tutorial for processing a single sample of 10X single nuclear ATAC data

If needed, download practice data from 10X website ([fastq](https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_5k_nextgem/atac_pbmc_5k_nextgem_fastqs.tar)). 

### (optional) Create your singularity images from docker hub images
Not advised since they take up quite a bit of space (~15gb for cellranger image). <br>

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