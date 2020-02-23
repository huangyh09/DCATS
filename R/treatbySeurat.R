# Functions used for Farbehi et al 2019 data, maybe adapted for other real datasets 

#' Main function to filter and normalise the large data matrix 
#' 
#' @param tip.data A large data matrix that may be unfiltered and not normalised. Preferably this should be imported with
#' Read10X() function available from Seurat 
#' 
#' @return A list containing the filtered and normalised count matrix and the processed seurat object 
treatbySeurat <- function(tip.data) {
  
  # all preprocessing steps modelled from 
  tip.aggr <- CreateSeuratObject(tip.data, min.cells = 10, min.features = 200, project = "TIP aggregate")
  mito.genes <- grep("^mt-", rownames(tip.aggr@data), value = TRUE)
  percent.mito <- Matrix::colSums(tip.aggr@raw.data[mito.genes, ])/Matrix::colSums(tip.aggr@raw.data)
  tip.aggr <- AddMetaData(tip.aggr, percent.mito, "percent.mito")
  batch.id <- sub(".*-(.*)","\\1", tip.aggr@cell.names)
  batch.id <- replace(batch.id, batch.id=="1", "Sham")
  batch.id <- replace(batch.id, batch.id=="2", "MI-day 3")
  batch.id <- replace(batch.id, batch.id=="3", "MI-day 7")
  tip.aggr <- AddMetaData(tip.aggr, batch.id, "batch")
  tip.aggr <- FilterCells(object = tip.aggr, subset.names = c("nGene", "nUMI", "percent.mito"),
                          low.thresholds = c(200, 500, -Inf), 
                          high.thresholds = c(4000, 20000, 0.05))
  tip.aggr <- NormalizeData(object = tip.aggr, normalization.method = "LogNormalize", 
                            scale.factor = 10000, display.progress = F)
  tip.aggr <- FindVariableGenes(object = tip.aggr, mean.function = ExpMean, 
                                dispersion.function = LogVMR, display.progress = F, 
                                x.low.cutoff = 0.01, x.high.cutoff = 5, 
                                y.cutoff = 0.5, y.high.cutoff = 20)
  tip.aggr <- ScaleData(tip.aggr, vars.to.regress = c("nUMI"))
  tip.aggr <- RunPCA(object = tip.aggr, pc.genes = tip.aggr@var.genes, do.print = FALSE, pcs.compute=60)
  tip.aggr <- FindClusters(object = tip.aggr, reduction.type = "pca", dims.use = 1:54, 
                           resolution =   seq(0.2, 1.8, 0.2), print.output = 0, save.SNN = TRUE)
  
  list(tip.aggr[["RNA"]]@data, tip.aggr)
}


