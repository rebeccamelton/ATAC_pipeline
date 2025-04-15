# Helper functions for ATAC clustering analysis

# Function to find files with more robust searching
find_file <- function(search_dir, file_patterns, recursive = TRUE) {
    # Ensure search directory exists
    if (!dir.exists(search_dir)) {
        cat("Error: Search directory does not exist:", search_dir, "\n")
        return(NULL)
    }
    
    # Print directory contents for debugging
    cat("Searching for files in directory:", search_dir, "\n")
    cat("Directory contents:\n")
    print(list.files(search_dir, recursive = recursive, full.names = TRUE))
    
    # Try each pattern
    for (pattern in file_patterns) {
        files <- list.files(
            path = search_dir, 
            pattern = pattern, 
            full.names = TRUE, 
            recursive = recursive
        )
        
        if (length(files) > 0) {
            cat("Found file(s):", files, "\n")
            return(files[1])
        }
    }
    
    # If no file found
    cat("No files found matching patterns:", paste(file_patterns, collapse = ", "), "\n")
    return(NULL)
}

# Find H5 file
find_h5_file <- function(search_dir) {
    h5_patterns <- c(
        "raw_peak_bc_matrix.h5", 
        "filtered_peak_bc_matrix.h5", 
        "*_peak_bc_matrix.h5",
        "*.h5"
    )
    
    file <- find_file(search_dir, h5_patterns)
    
    if (is.null(file)) {
        stop("Could not find H5 file for peak matrix in directory: ", search_dir)
    }
    
    return(file)
}

# Find fragments file
find_fragments_file <- function(search_dir) {
    frag_patterns <- c(
        "fragments.tsv.gz", 
        "*fragments.tsv.gz"
    )
    
    file <- find_file(search_dir, frag_patterns)
    
    if (is.null(file)) {
        stop("Could not find fragments file in directory: ", search_dir)
    }
    
    return(file)
}

# Find single cell metadata file
find_singlecell_file <- function(search_dir) {
    metadata_patterns <- c(
        "singlecell.csv", 
        "*singlecell.csv"
    )
    
    file <- find_file(search_dir, metadata_patterns)
    
    if (is.null(file)) {
        stop("Could not find singlecell metadata file in directory: ", search_dir)
    }
    
    return(file)
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

atac_single_samp_plots <- function(adata, samp, reduction) {
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

print_single_sample_pdf_atac <- function(samp, pdf_dir, plots, suffix, reduction) {
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