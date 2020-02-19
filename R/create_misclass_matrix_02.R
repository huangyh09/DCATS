# Creating similarity matrix by bootstrapping (Hennig, 2007)

#' Main function to create similarity matrix 
#' 
#' 
create_misclass_matrix_02 <- function(bigM, clust) {
  
  # define the number of iterations 
  iterations <- 10
  
  # define the number of clusters 
  cluster_number <- length(unique(clust))
  
  # cluster_lists contains the gene-cell count matrices of each cluster
  cluster_lists <- lapply(seq_len(cluster_number), function(i) {bigM[,clust==i]} )
  
  # extract the names of the cells in each cluster 
  cluster_namelists <- lapply(cluster_lists, function(i) {colnames(i)} )
  
  # define a matrix called misclassM that will store the result 
  misclassM <- matrix(0, cluster_number, cluster_number)
  
  # generate copies of misclassM to store results during iteration 
  misclassMlists <- lapply(seq_len(iterations), function(i) {misclassM} )
  
  # begin iteration 
  for (i in 1:iterations) {
    
    # bootstrap the starting matrix and apply kmeans clustering 
    bootstrappedM <- bigM[,sample(1:ncol(bigM),1000)]
    scee <- SingleCellExperiment(assays=list(counts=bootstrappedM))
    scee <- logNormCounts(scee)
    kout <- kmeans((t(logcounts(scee))), centers = cluster_number)
    
    # extract the gene-cell count matrices of each cluster of the bootstrapped matrix
    boot_cluster_lists <- lapply(seq_len(cluster_number), function(j) {bootstrappedM[, kout$cluster==j]} )
    
    # extract the names of the cells in each cluster of the bootstrapped matrix
    boot_cluster_namelists <- lapply(boot_cluster_lists, function(j) {colnames(j)} )
    
    # for each parent cluster of the parent matrix, return the number of cells of each daughter cluster that exists in the parent cluster
    permuteset <- lapply(seq_len(cluster_number), function(k) {sapply(boot_cluster_namelists, function(j) {sum(j %in% cluster_namelists[[k]])} )} )
    
    # store the results of reindexing 
    memoriseind <- c()
    answerind <- c()
    
    # for each vector in permuteset, identify the index that returns the maximising element, no repeating is allowed 
    for (j in 1:cluster_number) {
      if (j == 1){
        answerind[j] <- which(permuteset[[j]]==max(permuteset[[j]]))[1]
        memoriseind[j] <- answerind[j]
      } else {
        answerind[j] <- which(permuteset[[j]]==max(permuteset[[j]][-memoriseind]))[1]
        memoriseind[j] <- answerind[j]
      }
    }
    
    # reindexing 
    newclusterlabels <- as.factor(kout$cluster)
    levels(newclusterlabels) <- answerind
    
    # new gene-cell count matrices of each cluster that is reindexed (the results is, say, cluster1 in newboot_cluster_lists is most similar to parent cluster 1)
    newboot_cluster_lists <- lapply(seq_len(cluster_number), function(j) {bootstrappedM[, newclusterlabels==j]} )
    
    # generate misclassification matrix 
    for (m in 1:cluster) {
      for (n in 1:cluster){
        misclassMlists[[i]][m,n] <- sum(colnames(newboot_cluster_lists[[m]]) %in% cluster_namelists[[n]])/(length(colnames(newboot_cluster_lists[[m]])) + length(cluster_namelists[[n]]))
      }
    }
  }
  
  # taking the average of the matrices in the list 
  misclassM <- Reduce("+", my.list) / length(misclassMlists)
  
  misclassM
  
}


