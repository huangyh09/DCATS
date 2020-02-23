# Functions used to create misclassification matrix 

#' Main create misclassification function 
#' 
#' @param bigM A gene expression count matrix
#' @param clust A vector labelling each cell to a cluster 
#' 
#' @return A sqaure matrix of a dimension equal to the number of clusters present in the data
create_misclass_matrix_01 <- function(bigM,clust) {
  
  #filtering bigM for the genes that are larger than the first quantiles 
  getVar <- apply(bigM, 1, var)
  bigM <- bigM[getVar >= quantile(getVar)[4], ]
  bigM <- apply(bigM, 1, function(i){log(i/sum(i)*10000)})
  bigM <- t(bigM)
  
  #normalisation of bigM 
  summarymatrix <- apply(bigM,1,function(i){summary(i)})
  summarymatrix <- t(summarymatrix)
  
  #trinarising bigM 
  trinarisedM <- apply(bigM, 2, function(i){interrogate(i,summarymatrix)})
  
  #create misclassification matrix 
  mat <- misclassification_01(trinarisedM,clust)
  
  mat
}


#' Accessory function to create_misclass_matrix_01 that produces the misclassification matrix
misclassification_01 <- function(trinarisedM, clust) {
  
  #determining the dimension of the misclassifcation matrix
  howlong <- length(unique(clust))
  misclassM <- matrix(0,howlong,howlong)
  
  #create representatives of cluster 
  RES <- createRepresentatives(trinarisedM,clust)
  clusters <- RES[[1]]
  representatives <- RES[[2]]
  
  for (i in 1:howlong) {
    for (j in 1:howlong) {
      J <- representatives[[i]] != representatives[[j]]
      if (sum(J) ==0) {
        misclassM[i,j] <- dim(clusters[[i]])[2]
      }
      else {
        shortrepresentative_j <- representatives[[j]][J]
        correctScore <- admission(clusters[[j]][J,],shortrepresentative_j)
        testScore <- admission(clusters[[i]][J,],shortrepresentative_j)
        probability <- sapply(testScore,function(i){ecdf(correctScore)(i)})
        misclassM[i,j] <- sum(probability>0.5)
      }
    }
  }
  
  #normalisation of row of misclassification matrix by diagonal elements
  misclassM <- lapply(seq_len(nrow(misclassM)),function(i){misclassM[i,]/misclassM[i,][i]})
  misclassM <- do.call(rbind,misclassM)
  
  misclassM
}


#' Accessory function to misclassification_01 that creates representatives of a cluster 
createRepresentatives <- function(trinarisedM, clust) {
  
  listObj <- lapply(seq_len(length(unique(clust))), function(i){trinarisedM[,clust==i-1]})
  representatives <- lapply(listObj, function(x){apply(x,1,function(i){getmode(i)})})
  
  result <- list(listObj,representatives)
}


#' Accessory function to create_misclass_matrix_01 that helps with trinarisation 
interrogate <- function(x,summarymatrix) {
  
  A1 <- x > summarymatrix[,2]
  A2 <- x > summarymatrix[,5]
  A3 <- 1*A1 + 1*A2
  
  A3
}


#' Accessory function to misclassification_01 that compares the cells of a trinarised matrix to a cluster's
#' representative 
admission <- function(trinarisedM2, representative1) {
  
  examresult <- apply(trinarisedM2, 2, function(i) {i==representative1} )
  examresult <- 1*examresult 
  scores <- colSums(examresult)
  
  scores
}


#' Accessory function to createRepresentatives 
getmode <- function(v) {
  
  uniqv <- unique(v)
  
  uniqv[which.max(tabulate(match(v, uniqv)))]
}





