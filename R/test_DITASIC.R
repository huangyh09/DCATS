# Main function that tests DITASIC with similarity matrix created by create_misclass_matrix_01

#' @param N An integer specifying number of iterations the function is to be called
#' @param x A vector of proportions signifying the desired proportions of cell types for
#' subjects drawn from a hypothetical condition A. The number of clusters created by the function 
#' is equal to the length of x
#' @param y A vector of proportions signifying the desired proportions of cell types for 
#' subjects drawn from a hypothetical condition B.
#' @param replicates_number An integer signifying the number of subjects drawn from each condition
#' @param subject_cellnumber An integer signifying the number of cells drawn from a subject 
#' @param misclassfunction A function parameter specifying the function used to create misclassification,
#' with default being the identity matrix 
#' 
#' @return A list whose first element is a vector containing the number of cell types which exhibit differential
#' proportions between conditions, and whose second element is a vector containing the number of cell types which
#' do not exhibit such property.
test_DITASIC <- function(N, x, y, replicates_number, subject_cellnumber, misclassfunction) {
  
  # create a vector that stores the result 
  positive_result <- c()
  negative_result <- c()
  
  # create simulated data and extract required objects 
  result <- simulateData(x,y,replicatesNumber,subjectCellNumber)
  countdc <- result[[1]]
  cluster_lists <- result[[3]]
  cluster_number <- length(x)

  # iteration steps 
  for (i in 1:N) {
    # create gene expression profiles for condition A and condition B
    conditionA <- lapply(seq_len(cluster_number),function(j){cluster_lists[[j]][,sample(ncol(cluster_lists[[j]]),size=round(5500*x[j])),drop=FALSE]})
    conditionB <- lapply(seq_len(cluster_number),function(j){cluster_lists[[j]][,sample(ncol(cluster_lists[[j]]),size=round(5500*y[j])),drop=FALSE]})
    conditionA <- do.call(cbind,conditionA) 
    conditionB <- do.call(cbind,conditionB)
    conditionA <- conditionA[,sample(ncol(conditionA))] #permute the columns of the matrix 
    conditionB <- conditionB[,sample(ncol(conditionB))]
    
    # create gene expression profiles for subjects in different conditions 
    subjectsA <- lapply(seq_len(replicates_number),function(j){conditionA[,sample(ncol(conditionA),size=subject_cellnumber),drop=FALSE]})
    subjectsB <- lapply(seq_len(replicates_number),function(j){conditionB[,sample(ncol(conditionB),size=subject_cellnumber),drop=FALSE]})
    subjects <- append(subjectsA,subjectsB)
    subjects <- renamecountmatrix(subjects)
    
    # extract condition-specific gene expression count matrix 
    condAindex <- seq(replicates_number)
    condBindex <- sapply(seq_len(replicates_number),function(j) {j+replicates_number})
    condAmat <- lapply(condAindex,function(j) {subjects[[j]]})
    condBmat <- lapply(condBindex,function(j) {subjects[[j]]})
    condAmat <- do.call(cbind,condAmat)
    condBmat <- do.call(cbind,condBmat)
    
    # collect all cells from all subjects into a gene expression count matrix and perform k-means clustering
    SubjectsMatrix <- do.call(cbind,subjects)
    sce <- SingleCellExperiment(assays=list(counts=SubjectsMatrix))
    sce <- logNormCounts(sce)
    kout <- kmeans((t(logcounts(sce))), centers = cluster_number)
    
    # extract cluster labels output by kmeans clustering 
    clust <- as.factor(kout$cluster)
    levels(clust) <- sapply(seq_len(cluster_number),function(j) {j-1})
    misclassM <- misclassfunction(SubjectsMatrix,clust)
    
    # extract the names of the cells in different conditions 
    condAnames <- colnames(condAmat)
    condBnames <- colnames(condBmat)
    
    # extract the names of the cells in different clusters 
    subjectClusterNames <- lapply(seq_len(cluster_number),function(t){names(kout$cluster[kout$cluster==t])})
    
    # estimate the cell type counts for each condition base on name matches 
    condAcounts <- unlist(lapply(subjectClusterNames,function(i){detectcondA(i)}))
    condBcounts <- unlist(lapply(subjectClusterNames,function(i){detectcondB(i)}))
    
    # imputation of counts for absenct cell types 
    if (0 %in% condAcounts) {
      condAcounts[condAcounts==0] <- 1
    } else if (0 %in% condBcounts) {
      condBcounts[condBcounts==0] <- 1
    }
    
    # all the steps in DITASIC 
    abundA <- EstimateCount(misclassM,condAcounts)
    coeffA <- abundA[[1]]
    abundB <- EstimateCount(misclassM,condBcounts)
    coeffB <- abundB[[1]]
    empdistA <- create.empdist(abundA)
    empdistB <- create.empdist(abundB)
    AR <- DiffExpr(empdistA,empdistB,coeffA,coeffB)
    
    # recording results 
    positive_result[i] <- sum(AR[,4] < 0.05)
    negative_result[i] <- sum(AR[,4] >= 0.05)
  }
  
  list(positive_result,negative_result)
}


#' Accessory function for test_DITASIC that find the number of matches of names of cells 
#' in a cluster that belong to conditionA 
detectcondA <- function(x) {
  
  condA_subjectnames <- paste("subject",condAindex)
  matches <- sapply(condA_subjectnames,function(j) {sum(str_detect(x,j))})
  
  matches <- sum(matches)
}


#' Accessory function for test_DITASIC that find the number of matches of names of cells 
#' in a cluster that belong to conditionA 
detectcondB <- function(x) {
  
  condB_subjectnames <- paste("subject",condBindex)
  matches <- sapply(condB_subjectnames,function(j) {sum(str_detect(x,j))})
  
  matches <- sum(matches)
}