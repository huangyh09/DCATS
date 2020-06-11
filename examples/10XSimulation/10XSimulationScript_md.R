#' A script to reproduce the experimental results of the simulation with 10X neuron data 
#' 
#' @data
#' 1. tsne10X.txt
#' 2. meta10X.txt 
#' which should have been stored in the working directory, and which should not have their names altered
#' 
#' @return
#' 1. Similarity matrix based on KNN
#' 2. Similarity matrix based on K-means clustering 
#' 3. Confusion matrix representing the confusion of k-means cluster with the original cluster 
#' which will all be exported to the working directory once the codes below have been run 

#### Required packages ####

library(dplyr)
library(clue)

#### Data import ####

wdname <- getwd()
data1 <- paste(wdname,"/tsne10X.txt",sep="")
data2 <- paste(wdname,"/meta10X.txt",sep="")
tsne10X <- read.csv(data1)
meta10X <- read.csv(data2)
meta10X <- meta10X[-1,]
tsne10X <- tsne10X[-1,]
rownames(meta10X) <- paste("Cell",1:1306127,sep="")
rownames(tsne10X) <- paste("Cell",1:1306127,sep="")
meta10X[,2] <- as.numeric(meta10X[,2])
meta10X[,2] <- as.factor(meta10X[,2])
levels(meta10X[,2]) <- paste("Clust",1:20,sep="")

#### Sampling ####

set.seed(7)
indices <- sample(1:1306127,100000)
smalltsne <- tsne10X[indices,]
smallmeta <- meta10X[indices,]

#### Similarity matrix based on K-means clustering #### 

subindices <- sample(1:100000,50000)
smalltsne_toclust <- smalltsne[subindices,-1]
smallmeta_toclust <- smallmeta[subindices,]
nclust <- length(unique(smallmeta_toclust[,2]))
km <- kmeans(smalltsne_toclust,nclust,10000)
kmclust <- km$cluster
kmclust <- as.factor(kmclust)
levels(kmclust) <- paste("expClust",1:nclust,sep="")
X <- table(kmclust,as.character(smallmeta_toclust[,2]))
permuteorder <- as.vector(solve_LSAP(X,maximum = TRUE))
sapply(permuteorder, function(s){X[,s]} )

