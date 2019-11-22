#' @title Generates a doublet score for each cell via doubletCells
#' @description Uses doubletCells from the scran package to determine
#'  cells within the dataset suspected to be doublets.
#' @param sce SingleCellExperiment object. Must contain a counts matrix
#' @param assayType Character. The name of the counts matrix within your
#'  singleCellExperiment object that you wish to run doubletCells.
#' @return SingleCellExperiment object containing the 'scran_doublet_score'.
#' @examples
#' sce <- runScranDoubletScore(sce, assayType = "rsemcounts")
#' @export
#' @import scran
runScranDoubletScore <- function(sce, assayType){
    if(!assayType %in% names(SummarizedExperiment::assays(sce))){
        stop("'assayType' must be one of the assays in the SingleCellExperiment
            object.")
    }

    doubletScore <- rep(NA, ncol(sce))
    allSampleNumbers <- sort(unique(sce@colData$sample))
    for(num in allSampleNumbers){
        sceSubIx <- which(sce@colData$sample == num)
        sceSub <- sce[,sceSubIx]
        dblDens <- scran::doubletCells(sceSub, assay.type = assayType)
        doubletScore[sceSubIx] <- dblDens
    }
    colData(sce)$scran_doublet_score <- doubletScore
    return(sce)
}

#' @title Generates a doublet score for each cell via doubletFinder
#' @description Uses doubletFinder to determine cells within the dataset
#'  suspected to be doublets.
#' @param sce SingleCellExperiment object. Must contain a counts matrix
#' @return SingleCellExperiment object containing the
#'  'doublet_finder_doublet_score'.
#' @examples
#' sce <- runDoubletFinder(sce, assayType = "rsemcounts")
#' @export
#' @import SummarizedExperiment
#' @import Seurat
#' @import DoubletFinder
runDoubletFinder <- function(sce, assayType){
    if(!assayType %in% names(SummarizedExperiment::assays(sce))){
        stop("'assayType' must be one of the assays in the SingleCellExperiment
            object.")
    }

    doubletScore <- rep(NA, ncol(sce))
    doubletLabel <- rep(NA, ncol(sce))
    allSampleNumbers <- sort(unique(sce@colData$sample))
    for(num in allSampleNumbers){
        sceSubIx <- which(sce@colData$sample == num)
        sceCounts <- SummarizedExperiment::assays(sce)[[assayType]]
        sceCountsSub <- sceCounts[,sceSubIx]

        seurat <- Seurat::CreateSeuratObject(counts = sceCountsSub,
            project = "seurat", min.features = 200)
        seurat <- Seurat::NormalizeData(object = seurat,
            normalization.method="LogNormalize", scale.factor = 10000,
            verbose = F)

        seurat <- Seurat::FindVariableFeatures(seurat, selection.method = "vst",
            nfeatures = 2000, verbose = F)

        allGenes <- rownames(seurat)
        seurat <- Seurat::ScaleData(seurat, features = allGenes, verbose = F)

        numPc <- min(ncol(seurat@assays$RNA@scale.data) - 1, 50)
        seurat <- Seurat::RunPCA(seurat, features =
                Seurat::VariableFeatures(object = seurat),
                npcs = numPc, verbose = F)

        seurat <- Seurat::FindNeighbors(seurat, dims = 1:15, verbose = F)
        seurat <- Seurat::FindClusters(seurat, resolution = 1.2, verbose = F)

        #seurat <- RunUMAP(seurat, dims = 1:10)
        #seurat <- Seurat::RunTSNE(seurat, dims = 1:15)

        sweepResListSeurat <- DoubletFinder::paramSweep_v3(seurat,
            PCs = 1:15, sct = FALSE, verbose = F)
        sweepStatsSeurat <- DoubletFinder::summarizeSweep(sweepResListSeurat,
            GT = FALSE, verbose = F)
        bcmvnSeurat <- DoubletFinder::find.pK(sweepStatsSeurat)
        pkOptimal <- as.numeric(as.matrix(bcmvnSeurat$pK[
            which.max(bcmvnSeurat$MeanBC)]))

        annotations <- seurat@meta.data$seurat_clusters
        homotypicProp <- DoubletFinder::modelHomotypic(annotations)
        nExpPoi <- round(0.075*ncol(seurat@assays$RNA))
        seurat <- DoubletFinder::doubletFinder_v3(seurat,
            PCs = 1:15, pN = 0.25, pK = pkOptimal, nExp = nExpPoi,
            reuse.pANN = FALSE, sct = FALSE)

        names(seurat@meta.data)[6] <- "doubletFinderAnnScore"
        names(seurat@meta.data)[7] <- "doubletFinderLabel"

        doubletScore[sceSubIx] <- seurat@meta.data$doubletFinderAnnScore
        doubletLabel[sceSubIx] <- seurat@meta.data$doubletFinderLabel
    }
    colData(sce)$doublet_finder_doublet_score <- doubletScore
    colData(sce)$doublet_finder_doublet_label <- doubletLabel
    return(sce)
}


#' @title Generates an empty drops score for each cell via EmptyDrops
#' @description Uses EmptyDrops to determine cells within the dataset
#'  suspected to be doublets.
#' @param sce SingleCellExperiment object. Must contain a counts matrix
#' @param lower Number. Specifies the lower bound on the total UMI count.
#' All cells containing UMIs less than this will be considered empty droplets.
#' @param assayType Character. The name of the counts matrix within your
#'  singleCellExperiment object that you wish to run EmptyDrops.
#' @return SingleCellExperiment object containing the
#'  'droplet_utils_empty_drops_score'.
#' @examples
#' sce <- runEmptyDrops(sce, lower = 100, assayType = "rsemcounts")
#' @export
#' @import SummarizedExperiment
#' @import DropletUtils
runEmptyDrops <- function(sce, lower, assayType){
 if(!assayType %in% names(SummarizedExperiment::assays(sce))){
   stop("'assayType' must be one of the assays in the SingleCellExperiment
        object.")
 }
    emptyDropsScore <- rep(NA, ncol(sce))
    allSampleNumbers <- sort(unique(sce@colData$sample))
 for(num in allSampleNumbers){
    sceSubIx <- which(sce@colData$sample == num)
    sceCounts <- SummarizedExperiment::assays(sce)[[assayType]]
    sceCountsSub <- sceCounts[,sceSubIx]
    emptyDropsScore[sceSubIx] <- DropletUtils::emptyDrops(m = sceCountsSub, lower = lower)
    }
    colData(sce)$droplet_utils_empty_drops_score <- emptyDropsScore
    return(sce)
}

