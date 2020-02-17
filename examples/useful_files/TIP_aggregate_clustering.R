## We're going to be using the Seurat R package for processing and clustering the scRNA-seq data.
## Seurat 2.1.0 was used for our analysis.
## First load the relevant packages we'll be using.

library(ggplot2)
library(Seurat)
library(plyr)
library(dplyr)
library(Matrix)
library(cowplot)

##Read in an aggregate of the TIP data containing Sham, MI-day 3 and MI-day 7
tip.data <- Read10X("data/TIP_ShamVsMI_Days3_and_7")

tip.aggr <- CreateSeuratObject(tip.data, min.cells = 10, min.genes = 200, project = "TIP aggregate")
remove(tip.data) #  remove the raw data

## Add some meta-data to the Seurat object. Will add percent of RNA mapped to mitochondrial genes 
## and the batch ID so we can compare conditions later.

mito.genes <- grep("^mt-", rownames(tip.aggr@data), value = TRUE)
percent.mito <- Matrix::colSums(tip.aggr@raw.data[mito.genes, ])/Matrix::colSums(tip.aggr@raw.data)
tip.aggr <- AddMetaData(tip.aggr, percent.mito, "percent.mito")

## Add batch ID to meta data
## batches are in the following order:
## 1 - Sham
## 2 - MI-day 3
## 3 - MI-day 7
batch.id <- sub(".*-(.*)","\\1", tip.aggr@cell.names)
batch.id <- replace(batch.id, batch.id=="1", "Sham")
batch.id <- replace(batch.id, batch.id=="2", "MI-day 3")
batch.id <- replace(batch.id, batch.id=="3", "MI-day 7")
table(batch.id)
names(batch.id) = tip.aggr@cell.names
tip.aggr <- AddMetaData(tip.aggr, batch.id, "batch")

## Generate plots showing differences in UMI/gene number distributions between conditions before filtering
VlnPlot(object = tip.aggr, features.plot = c("nGene", "nUMI", "percent.mito"), 
        nCol = 3, point.size.use = 0.5)

par(mfrow = c(1, 2))
GenePlot(tip.aggr, "nUMI", "percent.mito", cex.use=1)
GenePlot(tip.aggr, "nUMI", "nGene",  cex.use=1)
par(mfrow = c(1, 1))

## Filtering out genes with outlier numbers of genes, or total UMI
## Also will filter out cells with > 5% RNA mapped to mitochondrial genes
tip.aggr <- FilterCells(object = tip.aggr, subset.names = c("nGene", "nUMI", "percent.mito"),
                        low.thresholds = c(200, 500, -Inf), 
                        high.thresholds = c(4000, 20000, 0.05))

## Normalise data
tip.aggr <- NormalizeData(object = tip.aggr, normalization.method = "LogNormalize", 
                          scale.factor = 10000, display.progress = F)

## Find higly variable genes to be used for PCA
tip.aggr <- FindVariableGenes(object = tip.aggr, mean.function = ExpMean, 
                              dispersion.function = LogVMR, display.progress = F, 
                              x.low.cutoff = 0.01, x.high.cutoff = 5, 
                              y.cutoff = 0.5, y.high.cutoff = 20)

print(paste0("Number of highly variable genes: ", length(x = tip.aggr@var.genes)))

## Scale the data, regressing out number of UMIs 
tip.aggr <- ScaleData(tip.aggr, vars.to.regress = c("nUMI"))

## Run PCA up to the first 60 components
tip.aggr <- RunPCA(object = tip.aggr, pc.genes = tip.aggr@var.genes, do.print = FALSE, pcs.compute=60)

## Check what genes have highest correlation with the top PCs
PCHeatmap(object = tip.aggr, pc.use = c(1:6), cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

## We use the JackStraw function with 1000 replicates to identify PCs for clustering. 
## This takes a long time to run, but the output is saved as a pdf
tip.aggr <- JackStraw(object = tip.aggr, num.pc = 60, num.replicate = 1000, do.print = FALSE)

## The JackStraw plot indiciates that up to the first 53 PCs are significant with P < 0.001
JackStrawPlot(object = tip.aggr, PCs = 1:60)

## For clustering, we use the PCs identified by the JackStraw test, but have found that modification of
## the PCs only causes minor changes in the clustering solutions. Will run clustering with a range of resolutions
tip.aggr <- FindClusters(object = tip.aggr, reduction.type = "pca", dims.use = 1:54, 
                         resolution =   seq(0.2, 1.8, 0.2), print.output = 0, save.SNN = TRUE)

## Found that a resolution of 1.2 returns a sensible clustering of the data.
## But can increase this to get more clusters
tip.aggr <- SetAllIdent(tip.aggr, id="res.1.2")

# t-SNE analysis.
tip.aggr <- RunTSNE(tip.aggr, dims.use = 1:54, do.fast = T, seed.use=1)

new.labels <- c("M1M\u03A6", "EC1", "F-SL", "F-Act", "F-SH", "BC", "M2M\u03A6", "M1Mo", "MYO",
              "EC3", "DC", "EC2", "TC1-Cd8", "TC2-Cd4", "MAC-TR", "Mural", "M2M\u03A6-EC", "Cyc",
              "MAC6", "MAC-IFNIC", "MAC7", "MAC8", "F-EC", "F-WntX", "EC-L1", "NKC", "EC-L2", "BC-TC", "Glial")
tip.aggr@ident = plyr::mapvalues(x = tip.aggr@ident, from = names(table(tip.aggr@ident)), to = new.labels)

col.set <- c("#c10023", "#008e17", "#fb8500", "#f60000", "#fde800", "#bc9000","#4ffc00", "#00bcac", "#0099cc",
             "#D35400", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#ffb600", "#0053c8",
             "#f2a287","#ffb3ff", "#800000", "#77a7b7", "#630099", "#00896e", "#ffba4f", "#00cc99", "#a81c0d", "#00ffae", "#FE0092")
TSNEPlot(tip.aggr, do.label = TRUE, pt.size = 0.75, colors.use = col.set)

## Our analysis showed that the M2MAC-EC (M2M\u03A6-EC), F-EC, EC-L1, EC-L2 and BC-TC population 
## represent hybrid populations that cannot be discounted as doublets.
## We therefore remove these populations

hybrid.cl = c("M2M\u03A6-EC", "F-EC", "EC-L1", "EC-L2", "BC-TC")
hybrid.cells = names(tip.aggr@ident)[as.character(tip.aggr@ident) %in% hybrid.cl]
length(hybrid.cells)
cl.keep = setdiff(new.labels, hybrid.cl)

cells.use = setdiff(tip.aggr@cell.names, hybrid.cells)

tip.aggr = FilterCells(tip.aggr, subset.names = NULL, cells.use = cells.use)

col.set.update <- c("#c10023", "#008e17", "#fb8500", "#f60000", "#fde800", "#bc9000","#4ffc00", "#00bcac", "#0099cc",
                    "#D35400", "#00eefd", "#5f777f", "#cf6bd6", "#99cc00", "#aa00ff", "#ff00ff", "#0053c8",
                    "#f2a287","#ffb3ff", "#800000", "#77a7b7", "#00896e", "#00cc99", "#FE0092")
tip.aggr@ident = factor(tip.aggr@ident, levels = cl.keep)

TSNEPlot(tip.aggr, do.label = TRUE, colors.use = col.set.update, pt.size = 0.75)


## Finally can save the object as an R data file
save(tip.aggr, file = "data/TIP_aggregate_Seurat.RData")

