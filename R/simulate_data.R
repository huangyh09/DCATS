# Functions used to simulate data using splatter

#' Main simulate data function 
#' 
#' @param x A vector of proportions signifying the desired proportions of cell types for
#' subjects drawn from a hypothetical condition A. The number of clusters created by the function 
#' is equal to the length of x
#' @param y A vector of proportions signifying the desired proportions of cell types for 
#' subjects drawn from a hypothetical condition B.
#' @param replicates_number An integer signifying the number of subjects drawn from each condition
#' @param subject_cellnumber An integer signifying the number of cells drawn from a subject 
#' 
#' @return a list containing \code{primary_count_matrix}, a large gene expression count matrix simulated from
#' splatter, containing equal proportions of cells from different clusters; \code{subjects}, a list of
#' gene expression profiles of different subjects; \code{cluster_lists}, a list of gene expression profiles of cells from different clusters.
simulate_data <- function(x, y, replicates_number, subject_cellnumber){
  
  # extracting the number of clusters 
  cluster_number <- length(x)
  
  # ensure that x and y are of the same dimensions 
  if (length(x) != length(y)) {
    warning("input cell type proportions are not of the same dimension")
  }
  else {
    params.groups <- newSplatParams(batchCells = 22500, nGenes=1000) # the data will contain 22500 cells and 1000 genes 
    simdc <- splatSimulateGroups(params.groups, # de.prob controls the probability a gene selected is differetially expressed; for the meaning of other parameters run "browseVignettes(splatter)".
                                 group.prob = rep(1/cluster_number,cluster_number),
                                 de.prob = rep(0.01,cluster_number),
                                 de.facLoc = rep(0.01,cluster_number),
                                 de.facScale = rep(0.8,cluster_number),
                                 verbose = FALSE)
  }
  
  # extract the count matrix from simulated data 
  countdc <- counts(simdc)
  primary_count_matix <- countdc
  
  # k-means clustering applied to the primary count matrix
  primary_sce <- SingleCellExperiment(assays=list(counts=primary_count_matix))
  primary_sce <- logNormCounts(primary_sce)
  primary_kout <- kmeans((t(logcounts(primary_sce))),centers = cluster_number)
  cluster_lists <- lapply(seq_len(cluster_number),function(t){countdc[,primary_kout$cluster==t]})
  
  # create gene expression profiles for condition A and condition B
  conditionA <- lapply(seq_len(cluster_number),function(i){cluster_lists[[i]][,sample(ncol(cluster_lists[[i]]),size=round(5500*x[i])),drop=FALSE]})
  conditionB <- lapply(seq_len(cluster_number),function(i){cluster_lists[[i]][,sample(ncol(cluster_lists[[i]]),size=round(5500*y[i])),drop=FALSE]})
  conditionA <- do.call(cbind,conditionA) 
  conditionB <- do.call(cbind,conditionB)
  conditionA <- conditionA[,sample(ncol(conditionA))] #permute the columns of the matrix 
  conditionB <- conditionB[,sample(ncol(conditionB))]
  
  # create gene expression profiles for subjects in different conditions 
  subjectsA <- lapply(seq_len(replicates_number),function(i){conditionA[,sample(ncol(conditionA),size=subject_cellnumber),drop=FALSE]})
  subjectsB <- lapply(seq_len(replicates_number),function(i){conditionB[,sample(ncol(conditionB),size=subject_cellnumber),drop=FALSE]})
  subjects <- append(subjectsA,subjectsB)
  subjects <- renamecountmatrix(subjects)
  result <- list(countdc,subjects,cluster_lists)
  names(result) <- c("primary_count_matrix","subjects","cluster_lists")
  result 
}


#' Accessory function to simulate_data that renames the columns of the gene expression count 
#' matrices of subjects by the order of cell and the index of the subject 
renamecountmatrix <- function(x) {
  
  g <- x #Admit a count matrix of a subject, possible to admit multiple subjects
  
  for (i in 1:length(g)) {
    newnames <- paste("subject",i,"Cell",1:dim(g[[i]])[2],sep="")
    colnames(g[[i]]) <- newnames
  }
  
  g
}

