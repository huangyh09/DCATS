#' A script to reproduce the experimental results from testing on splatter simulated dataset (approach 3)
#' 
#' @description 
#' 
#' 
#' @data 
#' 1. count matrix generated from splatter (no input is necessary)
#' @return 
#' 
#' 
#' 
#### Required packages #### 

library(tidyr)
library(stringr)
library(nabor)
library(aod)
library(DCATS)
library(dplyr)
library(clue)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(plotly)
library(svMisc)
library(gtools)
library(DCATS)
library(lme4)
library(scDC)
library(splatter)
library(scater)
library(Seurat)

#### Function to generate simulated data using splatter ####
#' @param batchCells, a positive integer, Number of cells to simulate
#' @param nGenes, a positive integer, Number of genes measured for each cell 
#' @param group.prob A numeric vector of n dimension, whose i-th entry indicates the proportion of cells belonging to i-th cluster that exist in the simulated data 
#' @param de.prob A numeric vector of the same dimension as group.prob, whose i-th entry indicates the proportion of nGenes that are differentially expressed
#' for the i-th cluster compared with the rest of the clusters 
#' 
#' @return A list containing
#' (1) sim1count_mat, A matrix, the gene-cell count matrix of the simulated data
#' (2) celllabels_orig, A numeric vector of batchCells dimension whose i-th entry indicates the cluster identity of the i-th cell 
#' 
#' @example (uncomment to run)
# param.groups <- newSplatParams(batchCells = 10000, nGenes =3000)
# generating data with two closely associated groups
# sim1 <- splatSimulateGroups(param.groups, group.prob = c(1/3,1/3,1/3),
#                                     de.prob = c(0.01,0.01,0.3),
#                                     verbose = FALSE)
# sim1 <- logNormCounts(sim1)
# sim1 <- runPCA(sim1)
# plotPCA(sim1, colour_by = "Group") + ggtitle("Few DE genes")
# celllabels_orig <- sim1@colData@listData$Group
# sim1count_mat <- counts(sim1)
generateSimData <- function(batchCells, nGenes, group.prob, de.prob){
  param.groups <- newSplatParams(batchCells = batchCells, nGenes = nGenes)
  sim1 <- splatSimulateGroups(param.groups, group.prob = group.prob,
                                      de.prob = de.prob,
                                      verbose = FALSE)
  celllabels_orig <- sim1@colData@listData$Group
  sim1count_mat <- counts(sim1)
  return(list(sim1count_mat,celllabels_orig))
}
#### Function to cluster by Seurat and to generate PCA loadings and KNN graph #### 
#' @example 
#' treated <- SeuratTreat(sim1count_mat)
#' ElbowPlot(treated)
#' treated <- FindNeighbors(treated, dims = 1:5)
#' treated <- FindClusters(treated, resolution = 0.1)
#' PCAloadings <- treated@reductions$pca@cell.embeddings
#' KNNgraph <- treated@graphs$RNA_snn
SeuratTreat <- function(countdata) {
  seuratobject <- CreateSeuratObject(counts = countdata, project="FunctionEnv")
  seuratobject <- NormalizeData(seuratobject)
  seuratobject <- FindVariableFeatures(seuratobject, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(seuratobject)
  seuratobject <- ScaleData(seuratobject,features = all.genes)
  seuratobject <- RunPCA(seuratobject,features = VariableFeatures(object = seuratobject))
  return(seuratobject)
}
#### Function to process seuratobject created from SeuratTreat ####
#' @param seuratobject, the seurat object output from SeuratTreat
#' @param dimensions, a vector of seq_len(X) where X is the number of dimensions to choose for finding the KNN graph
#' @param resolution, a parameter to control for the clustering algorithms used in Seurat 
#' @param
#' 
#' @return processed_seuratobject, which is a list containing (1) PCAloadings, the PCA loadings of all the cells input
#' and (2) KNN graph found, and (3) the seuratobject processed 

processSeuratObject <- function(seuratobject, dimensions, resolution,K) {
  seuratobject_inst <- seuratobject
  seuratobject_inst <- FindNeighbors(seuratobject_inst, k.params=K, dims = seq_len(dimensions))
  seuratobject_inst <- FindClusters(seuratobject_inst, resolution = resolution)
  PCAloadings <- seuratobject_inst@reductions$pca@cell.embeddings
  KNNgraph <- seuratobject_inst@graphs$RNA_snn 
  return(list(PCAloadings,KNNgraph,seuratobject_inst))
}
#### Function to confuse the original cell type labels ####
#' @example
#' two_to_one_prop <- runif(1,0,0.5)
#' one_to_two_prop <- runif(1,0,0.5)
#' celltypeINDEX <- lapply(c("Group1","Group2","Group3"),function(t){which(celllabels_orig == t)})
#' confusion_labels <- celllabels_orig
#' two_to_one_confused <- sample(celltypeINDEX[[2]],round(two_to_one_prop*length(celltypeINDEX[[2]])),replace = FALSE)
#' one_to_two_confused <- sample(celltypeINDEX[[1]],round(two_to_one_prop*length(celltypeINDEX[[1]])),replace = FALSE)
#' confusion_labels[two_to_one_confused] <- "Group1"
#' confusion_labels[one_to_two_confused] <- "Group2"
#' DCATS::KNN_transition(treated@graphs$RNA_snn,confusion_labels)
confuse_original_labels <- function(celllabels_orig) {
  two_to_one_prop <- runif(1,0,0.5)
  one_to_two_prop <- runif(1,0,0.5)
  celltypeINDEX <- lapply(c("Group1","Group2","Group3"),function(t){which(celllabels_orig == t)})
  confusion_labels <- celllabels_orig
  two_to_one_confused <- sample(celltypeINDEX[[2]],round(two_to_one_prop*length(celltypeINDEX[[2]])),replace = FALSE)
  one_to_two_confused <- sample(celltypeINDEX[[1]],round(one_to_two_prop*length(celltypeINDEX[[1]])),replace = FALSE)
  confusion_labels[two_to_one_confused] <- "Group1"
  confusion_labels[one_to_two_confused] <- "Group2"
  answer <- list(confusion_labels,two_to_one_prop,one_to_two_prop,celltypeINDEX)
  return(answer)
}
#### Function to subsample cells from simulated data and annotate with subject labels ####
#' @param newinputdata, which should be the transpose of sim1count_mat produced by generateSimData with an added column indicating the 
#' original cluster identity of each cell
#' @param random, a Boolean indicating whether the cell type proportion vectors should be randomly drawn from a dirichlet distribution
#' @param prop1, a proportion vector specified for condition1 if random is FALSE
#' @param prop2, a proportion vector specified for condition2 if random is FALSE
#' @param n_subject, number of subjects to simulate
#' @param n_cells, number of cells to be contributed by each subject
#' @param first_prop, the proportion of Group3 cells which is the same across condition. It must be specified if random is TRUE
#' 
#' @return a list containing (1) the cell-(gene,cluster,subject) matrix for condition1, (2) the cell-(gene,cluster,subject) matrix for
#' condition2, (3) proportion vector simulated for condition 1 and (4) proportion vector simulated for condition 2
samplingfromInputData_approach3 <- function(newinputdata, random=TRUE,prop1,prop2,n_subject,n_cells,conc1,conc2,first_prop) {
  
  clusterorder <- c("Group1","Group2","Group3")
  n_columns <- dim(newinputdata)[2]
  INDICES <- lapply(clusterorder,function(t){which(newinputdata[,n_columns]==t)})
  
  if (random== TRUE) {
    
    prop1 <- rdirichlet(1,rep(conc1,2))*(1-first_prop)
    prop2 <- rdirichlet(1,rep(conc2,2))*(1-first_prop)
    prop1 <- c(prop1,first_prop)
    prop2 <- c(prop2,first_prop)
    
    subjectcountVectorscond1 <- lapply(1:n_subject,function(t){rmultinom(1,n_cells,prop1)})
    subjectcountVectorscond2 <- lapply(1:n_subject,function(t){rmultinom(1,n_cells,prop2)})
    subjectcountVectorS <- append(subjectcountVectorscond1,subjectcountVectorscond2)
    
    clusterorder <- c(1:3)
    RESULT <- list()
    
    for (i in 1:(n_subject*2)) { #i is for each subject 
      stINDICES <- lapply(clusterorder,function(t) {sample(INDICES[[t]],subjectcountVectorS[[i]][t],replace = FALSE)} )
      INDICES <- lapply(clusterorder,function(t){INDICES[[t]][-which(INDICES[[t]] %in% stINDICES[[t]])]})
      RESULT[[i]] <- stINDICES #RESULT is a list of list; 
    }
    
    RESULT <- lapply(1:(n_subject*2), function(t){unlist(RESULT[[t]])})
    RESULT <- lapply(1:(n_subject*2), function(t){newinputdata[RESULT[[t]],]})
    subjectnames <- paste("subject",1:(n_subject*2),sep="")
    
    for (i in 1:(n_subject*2)) {
      RESULT[[i]] <- cbind(RESULT[[i]],rep(subjectnames[i],dim(RESULT[[i]])[1]))
      colnames(RESULT[[i]])[n_columns+1] <- "subjects"
    }
    
    RESULT1 <- do.call(rbind,RESULT[1:n_subject])
    RESULT2 <- do.call(rbind,RESULT[(n_subject+1):(n_subject*2)])
    
    return(list(RESULT1,RESULT2,prop1,prop2))
    
  }
  else {
    
    subjectcountVectorscond1 <- lapply(1:n_subject,function(t){rmultinom(1,n_cells,prop1)})
    subjectcountVectorscond2 <- lapply(1:n_subject,function(t){rmultinom(1,n_cells,prop2)})
    subjectcountVectorS <- append(subjectcountVectorscond1,subjectcountVectorscond2)
    
    clusterorder <- c(1:3)
    RESULT <- list()
    
    for (i in 1:(n_subject*2)) {
      stINDICES <- lapply(clusterorder,function(t){sample(INDICES[[t]],subjectcountVectorS[[i]][t],replace=FALSE)})
      INDICES <- lapply(clusterorder,function(t){INDICES[[t]][-which(INDICES[[t]] %in% stINDICES[[t]])]})
      RESULT[[i]] <- stINDICES 
    }
    RESULT <- lapply(1:(n_subject*2), function(t){unlist(RESULT[[t]])})
    RESULT <- lapply(1:(n_subject*2), function(t){bigM[RESULT[[t]],]})
    subjectnames <- paste("subject",1:(n_subject*2),sep="")
    
    for (i in 1:(n_subject*2)) {
      RESULT[[i]] <- cbind(RESULT[[i]],rep(subjectnames[i],dim(RESULT[[i]])[1]))
      colnames(RESULT[[i]])[n_columns+1] <- "subjects"
    }
    
    RESULT1 <- do.call(rbind,RESULT[1:n_subject])
    RESULT2 <- do.call(rbind,RESULT[n_subject:(n_subject*2)])
    
    return(list(RESULT1,RESULT2,prop1,prop2))
  }
}

#### Function to generate subject-cell-type count matrix and KNN-based similarity matrix ####
#' @param PCAloadings, which is the output of the @function samplingfromInputData_approach3
#' @param cofused_labels_RESULT, which is the output of the @function confuse_original_labels
#' 
#' @example 
#'PCAloadings <- cbind(as.data.frame(PCAloadings),as.data.frame(as.character(celllabels_orig)))
# newPCAloadings <- samplingfromInputData_approach3(PCAloadings,random = TRUE,n_subject=3,n_cells=1000,conc1=5,conc2=10,first_prop=0.3)
# 
# PCAload_sampled <- rbind(newPCAloadings[[1]],newPCAloadings[[2]])
# colnames(PCAload_sampled)[51] <- "celllabels_orig"
# subjectinfomat <- table(PCAload_sampled[,52],PCAload_sampled[,51])
# 
# confused_labels_RESULT <- confuse_original_labels(PCAload_sampled[,51])
# confused_labels <- confused_labels_RESULT[[1]]
# confused_subjectinfomat <- table(PCAload_sampled[,52],confused_labels)
# 
# subKNNgraph <- treated@graphs$RNA_snn[rownames(PCAload_sampled),rownames(PCAload_sampled)]
# misclassM <- KNN_transition(subKNNgraph,confused_labels)
testDataGenerate <- function(PCAloadings, confused_labels_RESULT) {
  newPCAloadings <- PCAloadings
  PCAload_sampled <- rbind(newPCAloadings[[1]],newPCAloadings[[2]])
  colnames(PCAload_sampled)[51] <- "celllabels_orig"
  subjectinfomat <- table(PCAload_sampled[,52],PCAload_sampled[,51])

  confused_labels <- confused_labels_RESULT[[1]]
  confused_subjectinfomat <- table(PCAload_sampled[,52],confused_labels)
  
  subKNNgraph <- treated@graphs$RNA_snn[rownames(PCAload_sampled),rownames(PCAload_sampled)]
  misclassM <- KNN_transition(subKNNgraph,confused_labels)
  
  return(list(confused_subjectinfomat,subjectinfomat,misclassM))
}

#### EXAMPLE - Analysis of the relationship between clustering confusion rates and off-diagonal entries of KNN-based similarity matrix ####
# RESULT <- list()
# for (i in seq_len(1000)){
#   confused_labels_RESULT <- confuse_original_labels(PCAload_sampled[,51])
#   confused_labels <- confused_labels_RESULT[[1]]
#   misclassM <- KNN_transition(subKNNgraph,confused_labels)
#   answer <- c(confused_labels_RESULT[[2]],confused_labels_RESULT[[3]],misclassM[1,1],misclassM[2,2],misclassM[1,2],misclassM[2,1])
#   names(answer) <- c("two_to_one_prop","one_to_two_prop","one_one_entry","two_two_entry","one_two_entry","two_one_entry")
#   RESULT[[i]] <- answer
# }
# 
# RESULT <- do.call(rbind,RESULT)
# RESULT <- as.data.frame(RESULT)
# plot(RESULT$one_to_two_prop,RESULT$one_one_entry)


#### Function for testing DCATS algorithms ####
#' A function to test DCATS algorithms
#' 
#' @param algorithm DCATS algorithms which accepts algorithms that admit similarity_matrix
#' @param result1 the tsne profiles from the sampling function above
DCATtest <- function(algorithm,subtinfomat,similarity_mat,p_threshold,diff_size_threshold,n_subject) {
  subject_number <- dim(subtinfomat)[1]
  decisionmat <- algorithm(subtinfomat[1:n_subject,],subtinfomat[(n_subject+1):subject_number,],similarity_mat = similarity_mat,n_samples=100)
  decision <- decisionmat[9] < p_threshold
  truth <- abs(result1[[3]]-result1[[4]]) >= diff_size_threshold
  confusion_tab <- table(truth,decision)
  return(list(subtinfomat,decisionmat,confusion_tab))
}


#### Function for applying Fisher's exact test to a subject-cell-type count matrix ####
test_Fisher_exact <- function(subtinfomat,p,q) {
  subtinfomat <- confused_subjectinfomat
  subtinfomat_cond1 <- subtinfomat[1:p,]
  subtinfomat_cond1 <- colSums(subtinfomat_cond1)
  subtinfomat_cond2 <- subtinfomat[(p+1):(p+q),]
  subtinfomat_cond2 <- colSums(subtinfomat_cond2)
  test_mat <- rbind(subtinfomat_cond1,subtinfomat_cond2)
  RESULT <- c(NA,NA,NA)
  for (i in 1:dim(subtinfomat)[2]) {
    vec1 <- test_mat[,i]
    vec2 <- rowSums(test_mat)-vec1
    tested_mat <- cbind(vec1,vec2)
    test <- fisher.test(tested_mat)
    p_value <- test$p.value
    RESULT[i] <- p_value 
  }
  return(RESULT)
}