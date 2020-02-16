# A function that creates the identity matrix as the similarity matrix 

#' @param bigM the gene expression count matrix from which the similarity matrix is to be constructed
#' @param clust a vector containing the cluster labels of each cell in bigM 
#' 
#' @return \code{misclassM}, an identity matrix of dimension equal to the number of unique clusters 
create_misclass_matrix_00 <- function(bigM, clust) {
  
  cluster_number <- length(unique(clust))
  
  misclassM <- matrix(0,nrow = cluster_number,ncol = cluster_number)
  diag(misclassM) <- 1
  
  misclassM
}