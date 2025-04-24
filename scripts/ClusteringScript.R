#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(Matrix)
library(stringr)
library(gridExtra)
library(logr)
library(BSgenome.Hsapiens.UCSC.hg38)})

# Source helper functions
source("/scripts/atac_clustering_helper_functions.R")

# Set global options
options("logr.on" = TRUE, "logr.notes" = FALSE)

# Trap any unhandled errors
options(error = function() {
    cat("An unexpected error occurred:\n")
    traceback()
    quit(status = 1)
})

# Argument parsing
parse_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    
    # Default values
    opts <- list(
        sample = NULL,
        cellranger_dir = '/03_cellranger/',
        output_dir = '/04_clustering/',
        pdf_dir = '/04_clustering/',
        resolution = 0.5,
        min_fragments = 1000,
        run_amulet = FALSE  # Add AMULET flag

    )
    
    # Parse arguments
    i <- 1
    while (i <= length(args)) {
        if (args[i] == "-s" || args[i] == "--sample") {
            opts$sample <- args[i+1]
            i <- i + 2
        } else if (args[i] == "-i" || args[i] == "--cellranger") {
            opts$cellranger_dir <- args[i+1]
            i <- i + 2
        } else if (args[i] == "-o" || args[i] == "--output") {
            opts$output_dir <- args[i+1]
            i <- i + 2
        } else if (args[i] == "-r" || args[i] == "--resolution") {
            opts$resolution <- as.numeric(args[i+1])
            i <- i + 2
        } else if (args[i] == "-f" || args[i] == "--min-fragments") {
            opts$min_fragments <- as.numeric(args[i+1])
            i <- i + 2
        } else if (args[i] == "--amulet") {
            opts$run_amulet <- TRUE
            i <- i + 1
        } else {
            cat("Unknown argument:", args[i], "\n")
            cat("Usage: Rscript cluster.r -s <sample_name> [options]\n")
            cat("Options:\n")
            cat("  -s, --sample        Sample name (required)\n")
            cat("  -c, --cellranger    Cellranger directory (default: /03_cellranger/)\n")
            cat("  -o, --output        Output directory (default: /04_clustering/)\n")
            cat("  -r, --resolution    Clustering resolution (default: 0.5)\n")
            cat("  -f, --min-fragments Minimum fragments per cell (default: 1000)\n")
            cat("  --amulet            Run AMULET doublet detection (default: disabled)\n")
            quit(status = 1)
        }
    }
    
    # Validate required arguments
    if (is.null(opts$sample)) {
        cat("Error: Sample name is required\n")
        cat("Usage: Rscript cluster.r -s <sample_name> [options]\n")
        quit(status = 1)
    }
    
    return(opts)
}

# Main clustering function
run_clustering <- function(opts) {
    # Open the log
    log_file <- paste(format(Sys.time(), "%Y_%m_%d"), "log", sep = ".")
    log_open(paste0(opts$output_dir, log_file))
    
    # Set directory paths
    wd <- paste0(opts$cellranger_dir, opts$sample, "/outs/")
    
    # Log start of processing
    log_print(paste("Starting analysis for sample:", opts$sample))
    log_print(paste("analysis wd:", wd))
    log_print(paste("analysis outdir:", opts$output_dir))

    # Load 10x data
    log_print("Loading h5")
    inputdata.10x <- Read10X_h5(file.path(wd, 'raw_peak_bc_matrix.h5'))
    
    atac.counts <- inputdata.10x
    
    # Load ATAC object
    chrom_assay <- load_atac_object(atac.counts, fragments=paste0(wd, 'fragments.tsv.gz'))
    
    adata <- CreateSeuratObject(
        counts = chrom_assay,
        assay = "ATAC",
        meta.data=read.table(paste0(wd, 'singlecell.csv'), sep=',', header=T, row.names=1)
    )
    
    # Get all barcodes
    amulet_barcodes <- data.frame(barcode=colnames(adata))
    
    # Print details and add sample metadata
    log_print(paste0(opts$sample, ": ", length(colnames(adata[["ATAC"]]))," total BCs"))
    log_print(paste0(opts$sample, " metadata dimensions: ",  stringr::str_flatten(dim(adata@meta.data), collapse='x')))
    log_print(paste0(opts$sample, " Seurat summary:"))
    print(adata)
    gc()
    
    # Filter cells by fragment count
    adata <- subset(
      x = adata,
      subset = passed_filters >= opts$min_fragments
    )
    log_print(paste(length(colnames(adata[["ATAC"]])),"total BCs after ATAC fragment threshold"))
    log_print(paste("ATAC median peaks per cell:",median(adata[[]][,'nFeature_ATAC'])))
    log_print(paste("ATAC median fragments per cell:",median(adata[[]][,'nCount_ATAC'])))
   
    # Make windows 
    genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
    fragments <- Fragments(adata)
    log_print('Starting windows...')

    windows <- GenomeBinMatrix(
         fragments = fragments,
         genome = genome,
         binsize = 5000,
         process_n = 2000
    )
    adata[["windows"]]  <- CreateChromatinAssay(
         counts = windows,
         sep = c("-", "-"),
         fragments = fragments,
         min.cells = 0,
         min.features = 0
    )
    log_print('Starting clustering...')

    # Cluster with specified resolution
    adata <- atac_cluster_single_samp(adata, ATAC_assay = 'windows', res = opts$resolution)
    
    # Save initial RDS
    saveRDS(adata, paste0(opts$output_dir, opts$sample,"_initial", ".RDS"))
    
    # Generate plots
    plots <- list()
    plots[[opts$sample]] <- atac_single_samp_plots(adata, opts$sample, reduction = 'windows')
    print_single_sample_pdf_atac(
        opts$sample, 
        opts$pdf_dir, 
        plots, 
        '_Initial_QC_Plots.pdf', 
        reduction = 'windows'
    ) 
    if (opts$run_amulet) {
        # Prepare AMULET barcode table
        log_print("Running AMULET doublet detection...")
        # Prepare AMULET barcode table
        amulet_table <- dplyr::left_join(
        amulet_barcodes, 
        data.frame(
            barcode=colnames(adata), 
            cell_id=colnames(adata), 
            is__cell_barcode=TRUE
        )
            )
    amulet_table$is__cell_barcode[is.na(amulet_table$is__cell_barcode)] <- FALSE
    amulet_table$is__cell_barcode <- as.numeric(amulet_table$is__cell_barcode)
    amulet_table$cell_id[is.na(amulet_table$cell_id)] <- 'None'
    
    # Write AMULET barcode table
    write.table(
        amulet_table, 
        paste0(opts$output_dir, opts$sample,  "_AMULET_barcode_table.csv"), 
        sep=',', col.names=T, row.names=T, quote=F
    )
    # Open the log
    Amulet_script_dir <- '/scripts/AMULET_v1.1/AMULET.sh'
    Amulet_out_dir_base <- paste0(opts$output_dir, 'AMULET_')
    fragment_base_path <- paste0(wd, 'fragments.tsv.gz')
    Amulet_table_path <- paste0(opts$output_dir, opts$sample,  "_AMULET_barcode_table.csv")
    autosome_file <- '/scripts/AMULET_v1.1/human_autosomes.txt'

    restriction_file <- '/scripts/AMULET_v1.1/restrictionlist_repeats_segdups_rmsk_hg38.bed'
    Amulet_package_path <- '/scripts/AMULET_v1.1'

    amulet_command1 <- c(paste0('mkdir ', Amulet_out_dir_base, opts$sample, ' -p'))
    amulet_command2 <- c(paste0('bash ',Amulet_script_dir,' --forcesorted --bambc CB --bcidx 1 --cellidx 2 --iscellidx 3 \ ',fragment_base_path,' ',Amulet_table_path,' ' ,autosome_file,' \ ',restriction_file,' \ ',Amulet_out_dir_base, opts$sample,' ', Amulet_package_path))

    log_print('Starting AMULET...')
    system(amulet_command1)   
    
    log_print('Run AMULET...')
    system(amulet_command2, intern = TRUE)

    
    log_print('Done with AMULET...')
    amulet_dir <- paste0(opts$output_dir, 'AMULET_')

    reduction <- 'windows'
    pdf(paste0(amulet_dir, opts$sample, '/',opts$sample, "_amulet_doublet_plots.pdf"))
    adata <- readRDS(paste0(opts$output_dir, opts$sample,"_initial",".RDS"))
    gc()
    
    multiplets <- read.table(paste0(amulet_dir, opts$sample, '/MultipletBarcodes_01.txt'))
    
    adata[[]]$amulet_doublet <- FALSE
    adata[[]][multiplets$V1,]$amulet_doublet <- TRUE
    
    #plot UMAPs of ATAC
    p1 <- DimPlot(adata, reduction=sprintf('umap.%s', reduction), group.by='seurat_clusters', label=FALSE, label.size=6, repel=TRUE)
    p1 <- p1 + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle(paste0('ATAC - ', opts$sample))
    p1 <- p1 + NoLegend() & theme(plot.title=element_text(hjust=0.5))
    print(p1)
    
    #plot UMAPs of ATAC
    p1 <- DimPlot(adata, reduction=sprintf('umap.%s', reduction), group.by='amulet_doublet', label=FALSE, label.size=6, repel=TRUE)
    p1 <- p1 + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle(paste0('ATAC - ', opts$sample))
    p1 <- p1 + NoLegend() & theme(plot.title=element_text(hjust=0.5)) 
    p1 <- p1 + annotate('text',x=Inf, y=Inf, hjust=1, vjust=1, label=paste0("# Doublets: ",sum(adata[[]]$amulet_doublet)))  + 
         annotate('text',x=Inf, y=Inf, hjust=1, vjust=3, 
             label=paste0("% Doublets: ", format(sum(adata[[]]$amulet_doublet) * 100 / nrow(adata[[]]), digits=2, nsmall=2), "%"))
    print(p1)
    
    p1 <- VlnPlot(adata, features=sprintf('nFeature_%s',reduction), group.by='amulet_doublet', pt.size=0, log=TRUE) + 
        geom_boxplot(width=.6, fill='white', alpha=.6, pt.size=0) + 
        geom_hline(yintercept=median(adata$nFeature_ATAC), linetype='dashed', lw=2) + 
        ggtitle(paste0('ATAC Features -  ', opts$sample))
    
    p2 <- VlnPlot(adata, features=sprintf('nCount_%s',reduction), group.by='amulet_doublet', pt.size=0, log=TRUE) + 
        geom_boxplot(width=.6, fill='white', alpha=.6, pt.size=0) + 
        geom_hline(yintercept=median(adata$nCount_ATAC), linetype='dashed', lw=2) + 
        ggtitle(paste0('ATAC Counts -  ', opts$sample))
    print(p1 / p2)
    
    adata <- subset(
      x = adata,
      subset = amulet_doublet==FALSE
    )
    
    adata <- atac_cluster_single_samp(adata)
    
    #plot UMAPs of ATAC
    p1 <- DimPlot(adata, reduction=sprintf('umap.%s', reduction), group.by='seurat_clusters', label=TRUE, label.size=6, repel=TRUE)
    p1 <- p1 + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle(paste0('ATAC - ', opts$sample))
    p1 <- p1 + NoLegend() & theme(plot.title=element_text(hjust=0.5))
    print(p1)
    print_single_sample_pdf_atac(
        opts$sample, 
        opts$pdf_dir, 
        plots, 
        '_PostAmulet_Plots.pdf', 
        reduction = 'windows'
    )
    saveRDS(adata, paste0(opts$output_dir, opts$sample,"_postAmulet",".RDS"))
    
    log_print("Clustering analysis including AMULET completed successfully.")
    } else {
        log_print("Clustering analysis completed successfully.")

    }
}

# Main execution
main <- function() {
    # Parse command-line arguments
    opts <- parse_args()
    
    # Run clustering
    run_clustering(opts)
}

# Call main function
main()