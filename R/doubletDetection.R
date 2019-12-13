#' @title Generates a doublet score for each cell via doubletCells
#' @description Uses doubletCells from the scran package to determine
#'  cells within the dataset suspected to be doublets.
#' @param sce SingleCellExperiment object. Must contain a counts matrix.
#' @param sampleLabels Numeric vector. Each cell will be assigned a sample number.
#' @return SingleCellExperiment object containing the 'scran_doublet_score'.
#' @examples
#' sce <- runScranDoubletScore(sce)
#' @export
#' @import scran
runScranDoubletScore <- function(sce, sampleLabels){

    doubletScore <- rep(NA, ncol(sce))
    allSampleNumbers <- sort(unique(sampleLabels))
    for(num in allSampleNumbers){
        sceSubIx <- which(sampleLabels == num)
        sceSub <- sce[,sceSubIx]
        dblDens <- scran::doubletCells(sceSub)
        doubletScore[sceSubIx] <- dblDens
    }
    colData(sce)$scran_doublet_score <- doubletScore
    return(sce)
}

#' @title Generates a doublet score for each cell via doubletFinder
#' @description Uses doubletFinder to determine cells within the dataset
#'  suspected to be doublets.
#' @param sce SingleCellExperiment object. Must contain a counts matrix
#' @param sampleLabels Numeric vector. Each cell will be assigned a sample number.
#' @return SingleCellExperiment object containing the
#'  'doublet_finder_doublet_score'.
#' @examples
#' sce <- runDoubletFinder(sce)
#' @export
#' @import SummarizedExperiment
#' @import Seurat
#' @import DoubletFinder
runDoubletFinder <- function(sce, sampleLabels){

    doubletScore <- rep(NA, ncol(sce))
    doubletLabel <- rep(NA, ncol(sce))
    allSampleNumbers <- sort(unique(sampleLabels))
    for(num in allSampleNumbers){
        sceSubIx <- which(sampleLabels == num)
        sceCounts <- SingleCellExperiment::counts(sce)
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

        sweepResListSeurat <- DoubletFinder::paramSweep_v3(seurat,
            PCs = 1:15, sct = FALSE)
        sweepStatsSeurat <- DoubletFinder::summarizeSweep(sweepResListSeurat,
            GT = FALSE)
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

#' @title Generates a SingleCellExperiment object containing the output of
#'  specified QC functions
#' @description A wrapper function for the individual QC algorithms.
#' @param sce SingleCellExperiment object. Must contain a counts matrix
#' @param sampleLabels Numeric vector. Label for each sample within the
#'  SingleCellExperiment. If NULL, all cells within SingleCellExperiment
#'  are assumed to be of the same sample. Default: NULL.
#' @param algorithm Character vector. Specifies which vector to run within
#'  function. Possible options are: "doubletscore", "doubletfinder",
#'  "emptydrops".
#' @return SingleCellExperiment object containing the outputs of the
#'  specified algorithms.
#' @examples
#' sce <- runQC(sce, algorithm = c("doubletscore", "doubletfinder"))
#' @export
#' @import SummarizedExperiment
runQC <- function(sce,
    sampleLabels = NULL,
    algorithm = c("doubletscore", "doubletfinder")){
    #Label samples
    if(is.null(sampleLabels)){
        sampleLabels <- rep(1,ncol(sce))
    }

    #Run algorithms
    if("doubletscore" %in% algorithm){
        sce = runScranDoubletScore(sce = sce, sampleLabels = sampleLabels)
    }
    if("doubletfinder" %in% algorithm){
        sce = runDoubletFinder(sce = sce, sampleLabels = sampleLabels)
    }
        return(sce)
}

