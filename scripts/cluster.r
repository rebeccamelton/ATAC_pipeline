
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(Matrix)
library(stringr)
library(gridExtra)
library(logr)

#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library('BSgenome.Hsapiens.UCSC.hg38')
options("logr.on" = TRUE, "logr.notes" = FALSE)

### testing the new R clustering code
outdir <- paste0('/05_clustering/')
samp <- 'atac_pbmc_5k_nextgem'
samp_dir <- '/03_cellranger/'
pdf_dir <- '/05_clustering/'

atac_cluster_single_samp <- function(adata, ATAC_assay='ATAC', min.cutoff='q0', res=0.5) {
    # ATAC analysis
    # We exclude the first dimension as this is typically correlated with sequencing depth
    DefaultAssay(adata) <- ATAC_assay
    adata <- RunTFIDF(adata)
    adata <- FindTopFeatures(adata, min.cutoff=min.cutoff)
    adata <- RunSVD(adata)
    adata <- RunUMAP(adata, reduction='lsi', dims=2:50,
                        reduction.name=paste0('umap.', str_to_lower(ATAC_assay)), 
                        reduction.key=paste0(str_to_lower(ATAC_assay),'UMAP_'))
    
    adata <- FindNeighbors(adata, assay=ATAC_assay, reduction='lsi', dims=2:50, verbose = FALSE)
    adata <- FindClusters(adata, algorithm=4, verbose=FALSE, resolution=res)
    
    return(adata)
}

load_atac_object <- function(atac.counts, meta=NULL, fragments=NULL) {
    if (!is.null(meta)) {
        meta <- read.table(meta, sep=',', header=T, row.names=1)
    }
    
    print("Creating ATAC Object")
    annotation <- suppressWarnings(GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86))
    seqlevels(annotation) <- suppressWarnings(paste0('chr', seqlevels(annotation)))
    
    chrom_assay <- CreateChromatinAssay(
        counts = atac.counts,
        sep = c(":", "-"),
        fragments = fragments,
        annotation = annotation,
        max.lines = NULL
    )
    
    return(chrom_assay)
}

atac_single_samp_plots <- function(adata, samp, marker.list,reduction) {
    samp_plots <- list()
    
    #plot UMAPs of ATAC
    p1 <- DimPlot(adata, reduction=sprintf('umap.%s', reduction), group.by='seurat_clusters', label=TRUE, label.size=6, repel=TRUE)
    p1 <- p1 + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle(paste0('ATAC - ', samp))
    p1 <- p1 + NoLegend() & theme(plot.title=element_text(hjust=0.5))
    samp_plots[[sprintf("UMAP_%s",reduction)]] <- p1
    
    #plot ATAC, RNA, and SCT(normalized RNA) metrics by cluster
    p3 <- VlnPlot(adata, features=sprintf('nCount_%s',reduction), group.by='seurat_clusters', pt.size=0, log=TRUE) + 
        geom_boxplot(width=.6, fill='white', alpha=.6, pt.size=0) + 
        geom_hline(yintercept=median(adata$nCount_ATAC), linetype='dashed', lw=2) + 
        ggtitle(paste0('ATAC Counts -  ', samp))
    p4 <- VlnPlot(adata, features=sprintf('nFeature_%s',reduction), group.by='seurat_clusters', pt.size=0, log=TRUE) + 
        geom_boxplot(width=.6, fill='white', alpha=.6, pt.size=0) + 
        geom_hline(yintercept=median(adata$nFeature_ATAC), linetype='dashed', lw=2) + 
        ggtitle(paste0('ATAC Features -  ', samp))
    samp_plots[[sprintf("Vln_%s_Count",reduction)]] <- p3
    samp_plots[[sprintf("Vln_%s_Feature",reduction)]] <- p4
    
    #LSI elbow plot
    samp_plots[["Elbow_LSI"]] <- ElbowPlot(adata, reduction='lsi') + ggtitle(paste0('LSI Elbow Plot -  ', samp))
    
    samp_plots[['DepthCor']] <- DepthCor(adata, assay = reduction, reduction = "lsi", n = 10) + 
        ggtitle(paste0(samp, ' - Correlation between depth and reduced dimension components'))
    
    return(samp_plots)
}

print_single_sample_pdf_atac <- function(samp, pdf_dir, plots, suffix,reduction) {
    pdf(paste0(pdf_dir, samp,"_",suffix), width=9, height=8) 
    
    options(repr.plot.width=6, repr.plot.height=6)
    print(plots[[samp]][[sprintf("UMAP_%s",reduction)]])
    
    #plot ATAC, RNA, and SCT(normalized RNA) metrics by cluster
    options(repr.plot.width=10, repr.plot.height=10)
    print(plots[[samp]][[sprintf("Vln_%s_Count",reduction)]] / plots[[samp]][[sprintf("Vln_%s_Feature",reduction)]])

    #Plot elbow plot
    print(plots[[samp]][["Elbow_LSI"]])

    options(repr.plot.width=12, repr.plot.height=10)
    print(plots[[samp]][['DepthCor']])
    
    dev.off()
}


# Load in 10x and calculate QC metrics

    # Open the log
   log_file <- paste(format(Sys.time(), "%Y_%m_%d"), "log", sep = ".")

    log_open(paste0(outdir,log_file))
    #Filter sample QC metrics, make sure it is one sample
    
    log_print("Loading h5")
    # Load in 10x data
    wd <- paste0(samp_dir, samp, "/outs/")
    inputdata.10x <- Read10X_h5(file.path(wd, 'raw_peak_bc_matrix.h5'))
    
    atac.counts <- inputdata.10x
    
    chrom_assay <- load_atac_object(atac.counts, fragments=paste0(wd, 'fragments.tsv.gz'))
    
    adata <- CreateSeuratObject(
        counts = chrom_assay,
        assay = "ATAC",
        meta.data=read.table(paste0(wd, 'singlecell.csv'), sep=',', header=T, row.names=1)
    )
    ## get all barcodes
    amulet_barcodes <- data.frame(barcode=colnames(adata))
    
    #Print some details and add sample metadata
    log_print(paste0(samp, ": ", length(colnames(adata[["ATAC"]]))," total BCs"))
    log_print(paste0(samp, " metadata dimensions: ",  stringr::str_flatten(dim(adata@meta.data), collapse='x')))
    log_print(paste0(samp, " Seurat summary:"))
    print(adata)
    gc()
    
    adata <- subset(
      x = adata,
      subset = passed_filters >= 1000
    )
    log_print(paste(length(colnames(adata[["ATAC"]])),"total BCs after ATAC fragment threshold"))
    log_print(paste("ATAC median peaks per cell:",median(adata[[]][,'nFeature_ATAC'])))
    log_print(paste("ATAC median fragments per cell:",median(adata[[]][,'nCount_ATAC'])))
   #make windows 
    genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
    fragments <- Fragments(adata)
    log_print('Starting windows...')

    windows <- GenomeBinMatrix(
         fragments = fragments,
         genome = genome,
         binsize = 5000,
        process_n = 2000,
       )
    adata[["windows"]]  <- CreateChromatinAssay(
         counts = windows,
         sep = c("-", "-"),
         fragments = fragments,
         min.cells = 0,
         min.features = 0
       )
    log_print('Starting clustering...')

    adata <- atac_cluster_single_samp(adata, ATAC_assay = 'windows' )
    
    saveRDS(adata, paste0(outdir, samp,"_initial", ".RDS"))
    
    
   amulet_table <- dplyr::left_join(amulet_barcodes, data.frame(barcode=colnames(adata), 
                                             cell_id=colnames(adata), 
                                             is__cell_barcode=TRUE))
    amulet_table$is__cell_barcode[is.na(amulet_table$is__cell_barcode)] <- FALSE
    amulet_table$is__cell_barcode <- as.numeric(amulet_table$is__cell_barcode)
    amulet_table$cell_id[is.na(amulet_table$cell_id)] <- 'None'
    
    write.table(amulet_table, paste0(outdir, samp,  "_AMULET_barcode_table.csv"), 
                sep=',', col.names=T, row.names=T, quote=F)

    log_print('Ready for amulet...')
    # Open the log
    Amulet_script_dir <- '/scripts/AMULET_v1.1/AMULET.sh'
    Amulet_out_dir_base <- paste0(outdir, 'AMULET_')
    fragment_base_path <- paste0(wd, 'fragments.tsv.gz')
    Amulet_table_path <- paste0(outdir, samp,  "_AMULET_barcode_table.csv")
    autosome_file <- '/scripts/AMULET_v1.1/human_autosomes.txt'

    restriction_file <- '/scripts/AMULET_v1.1/restrictionlist_repeats_segdups_rmsk_hg38.bed'
    Amulet_package_path <- '/scripts/AMULET_v1.1'

    amulet_command1 <- c(paste0('mkdir ', Amulet_out_dir_base, samp, ' -p'))
    amulet_command2 <- c(paste0('bash ',Amulet_script_dir,' --forcesorted --bambc CB --bcidx 1 --cellidx 2 --iscellidx 3 \ ',fragment_base_path,' ',Amulet_table_path,' ' ,autosome_file,' \ ',restriction_file,' \ ',Amulet_out_dir_base, samp,' ', Amulet_package_path))

    log_print('Starting AMULET...')
   # log_print('AMULET output directory: ',amulet_command1)
    #log_print('AMULET command: ',amulet_command2)

    system(amulet_command1)   
    
    log_print('Run AMULET...')
    system(amulet_command2, intern = TRUE)

    
    log_print('Done with AMULET...')
    amulet_dir <- paste0(outdir, 'AMULET_')
    log_print('Plot initial QC plots')
    plots <- list()
    plots[[samp]] <- atac_single_samp_plots(adata, samp, marker.genes, reduction = 'windows')
    print_single_sample_pdf_atac(samp, pdf_dir, plots, '_Initial_QC_Plots.pdf', reduction = 'windows')

    reduction <- 'windows'
    pdf(paste0(amulet_dir, samp, '/',samp, "_amulet_doublet_plots.pdf"))
    adata <- readRDS(paste0(outdir, samp,"_initial",".RDS"))
    gc()
    
    multiplets <- read.table(paste0(amulet_dir, samp, '/MultipletBarcodes_01.txt'))
    
    adata[[]]$amulet_doublet <- FALSE
    adata[[]][multiplets$V1,]$amulet_doublet <- TRUE
    
    #plot UMAPs of ATAC
    p1 <- DimPlot(adata, reduction=sprintf('umap.%s', reduction), group.by='seurat_clusters', label=FALSE, label.size=6, repel=TRUE)
    p1 <- p1 + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle(paste0('ATAC - ', samp))
    p1 <- p1 + NoLegend() & theme(plot.title=element_text(hjust=0.5))
    print(p1)
    
    #plot UMAPs of ATAC
    p1 <- DimPlot(adata, reduction=sprintf('umap.%s', reduction), group.by='amulet_doublet', label=FALSE, label.size=6, repel=TRUE)
    p1 <- p1 + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle(paste0('ATAC - ', samp))
    p1 <- p1 + NoLegend() & theme(plot.title=element_text(hjust=0.5)) 
    p1 <- p1 + annotate('text',x=Inf, y=Inf, hjust=1, vjust=1, label=paste0("# Doublets: ",sum(adata[[]]$amulet_doublet)))  + 
         annotate('text',x=Inf, y=Inf, hjust=1, vjust=3, 
             label=paste0("% Doublets: ", format(sum(adata[[]]$amulet_doublet) * 100 / nrow(adata[[]]), digits=2, nsmall=2), "%"))
    print(p1)
    
    p1 <- VlnPlot(adata, features=sprintf('nFeature_%s',reduction), group.by='amulet_doublet', pt.size=0, log=TRUE) + 
        geom_boxplot(width=.6, fill='white', alpha=.6, pt.size=0) + 
        geom_hline(yintercept=median(adata$nFeature_ATAC), linetype='dashed', lw=2) + 
        ggtitle(paste0('ATAC Features -  ', samp))
    
    p2 <- VlnPlot(adata, features=sprintf('nCount_%s',reduction), group.by='amulet_doublet', pt.size=0, log=TRUE) + 
        geom_boxplot(width=.6, fill='white', alpha=.6, pt.size=0) + 
        geom_hline(yintercept=median(adata$nCount_ATAC), linetype='dashed', lw=2) + 
        ggtitle(paste0('ATAC Counts -  ', samp))
    print(p1 / p2)
    
    adata <- subset(
      x = adata,
      subset = amulet_doublet==FALSE
    )
    
    adata <- atac_cluster_single_samp(adata)
    
    #plot UMAPs of ATAC
    p1 <- DimPlot(adata, reduction=sprintf('umap.%s', reduction), group.by='seurat_clusters', label=TRUE, label.size=6, repel=TRUE)
    p1 <- p1 + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle(paste0('ATAC - ', samp))
    p1 <- p1 + NoLegend() & theme(plot.title=element_text(hjust=0.5))
    print(p1)
    
    saveRDS(adata, paste0(outdir, samp,"_postAmulet",".RDS"))