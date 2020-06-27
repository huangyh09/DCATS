#' A script to reproduce the experimental results from testing on Kang 2018 SLE datasets
#' 
#' @data
#' 1. GSE96583_batch2.total.tsne.df.tsv
#' which should have been stored in the working directory, and which should not have their names altered 
#' 
#' @return 
#' 1. All the objects listed in "Output data" in https://www.notion.so/hephaes/Kang-2018-Results-Reproduction-f8cd983ad6c548569bb808d0f7a0033f
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

#### Data import and processing ####
datapath <- paste(getwd(),"GSE96583_batch2.total.tsne.df.tsv",sep="/")
GSE96583_batch2.total.tsne.df <- read.delim(dapath)

# data preprocessing 
namesbatch <- GSE96583_batch2.total.tsne.df[,3:4]
pastetwo_names <- function(a,b){
  c <- paste(a,b,sep="_")
  c
}
namesbatch <- apply(namesbatch,1,function(t){pastetwo_names(t[1],t[2])})

batchlabels <- GSE96583_batch2.total.tsne.df[,6][!is.na(GSE96583_batch2.total.tsne.df[,6])]
namesbatch <- namesbatch[!is.na(GSE96583_batch2.total.tsne.df[,6])]
batch2info_subt <- table(namesbatch,batchlabels)
batch2info_subt <- rbind(batch2info_subt[c(1,3,5,7,9,11,13,15),],batch2info_subt[c(2,4,6,8,10,12,14,16),])
batch2infoheatmat <- heat_matrix(batch2info_subt,show_value = TRUE)

# export data 

write.csv(batch2info_subt,"batch2subject_celltypecount_mat.csv")
ggsave("batch2infoheatmat.png", plot=batch2infoheatmat,width=30,height=7,units="cm")

#### Similarity matrix creation #### 

createKNN <- function(vec,trynn){
  ans <- rep(0, dim(trynn)[1])
  
  for (i in vec){
    ans[i] <- 1
  }
  ans
}


tsnemat <- GSE96583_batch2.total.tsne.df[,1:2][!is.na(GSE96583_batch2.total.tsne.df[,6]),]
myK <- 101
nn <-knn(data=tsnemat, query=tsnemat, k=myK+1);
trynn <- nn$nn.idx
KNN <- apply(trynn, 1, function(vec){createKNN(vec,trynn)})
rownames(KNN) <- colnames(KNN) <- rownames(tsnemat)
names(batchlabels) <- rownames(tsnemat)
misclassM <- KNN_transition(KNN,batchlabels)
names(batchlabels) <- rownames(tsnemat)
misclassM <- KNN_transition(KNN,batchlabels)

clusters_uniq <- unique(batchlabels)
n_cluster <- length(clusters_uniq)

properorder <- sapply(colnames(misclassM),function(t){which(colnames(batch2info)==t)})
robatch2info_subt <- sapply(properorder, function(t){batch2info_subt[,t]})

#### Experiment A ####

uM <- get_similarity_mat(8,0.05)
counts1 <- robatch2info_subt[1:8,]
counts2 <- robatch2info_subt[9:16,]
result_binom_wM <- dcats_fit(counts1,counts2,similarity_mat = misclassM, n_samples = 100)
result_betabin_wM <- dcats_betabin(counts1,counts2,similarity_mat = misclassM,n_samples = 100)
result_binom_woM <- dcats_fit(counts1,counts2,similarity_mat = diag(8), n_samples = 100)
result_betabin_woM <- dcats_betabin(counts1,counts2,similarity_mat = diag(8), n_samples = 100)
result_binom_wuM <- dcats_fit(counts1,counts2,similarity_mat = uM, n_samples = 100)
result_betabin_wuM <- dcats_betabin(counts1,counts2,similarity_mat = uM, n_samples = 100)
result_betaLRT <- betabinLRT(counts1,counts2)
RESULT <- list(result_binom_wM,result_binom_woM,result_binom_wuM,result_betabin_wM,result_betabin_woM,result_betabin_wuM,result_betaLRT)

# export data 

write.csv(result_binom_wM,"result_binom_wM.csv")
write.csv(result_betabin_wM,"result_betabin_wM.csv")
write.csv(result_binom_woM,"result_binom_woM.csv")
write.csv(result_betabin_woM,"result_betabin_woM.csv")
write.csv(result_binom_wuM,"result_binom_wuM.csv")
write.csv(result_betabin_wuM,"result_betabin_wuM.csv")
write.csv(result_betaLRT,"result_betaLRT.csv")

# organisation of p-values vectors 

PVAL <- lapply(RESULT,function(t){t[,9]})
pvalmat <- do.call(rbind,PVAL)
methodnames <- c("result_binom_wM","result_binom_woM","result_binom_wuM","result_betabin_wM","result_betabin_woM","result_betabin_wuM","result_betaLRT")
rownames(pvalmat) <- methodnames
colnames(pvalmat) <- colnames(misclassM)
pvalmat <- t(pvalmat)

# plot as bar graphs 

pvalmattoplot <- as.data.frame(-log(pvalmat))
pvalmattoplot$celltypes <- rownames(pvalmattoplot)
pvalmattoplot <- gather(pvalmattoplot,"methods","negativeLogP",1:7)

pvalwholeexpAbar <- ggplot(pvalmattoplot,aes(x=celltypes, y=negativeLogP, fill=methods)) + geom_bar(stat="identity", position=position_dodge()) +
  theme_minimal() + geom_hline(yintercept=3, linetype="dashed",  color = "red", size=0.5)
ggsave("pvalwholeexpAbar.png", plot=pvalwholeexpAbar,width=30,height=10,units="cm")

# new p-value plots, with woM methods all removed 

logindex <- sapply(pvalmattoplot[,2], function(t){str_detect(t,"woM")})
alterpvalmattoplot <- pvalmattoplot[!logindex,]
pvalaltexpAbar <- ggplot(alterpvalmattoplot,aes(x=celltypes, y=negativeLogP, fill=methods)) + geom_bar(stat="identity", position=position_dodge()) +
  theme_minimal() + geom_hline(yintercept=3, linetype="dashed",  color = "red", size=0.5)
ggsave("pvalaltexpAbar.png", plot=pvalaltexpAbar,width=30,height=10,units="cm")


#### Experiment B ####
#' Experiment description:
#' Modify the confusion rate of the uniform similarity matrix on a scale from 0.1 to 0.9 with 1000 intervals between them, 
#' the modified uniform similarity matrix is then fed into DCATS binom and DCATS betabin methods 
#' the p-values that result are plotted 

create_progression <- function(leftend, rightend,nintervals){
  progression <- sapply(1:nintervals+1,function(t){leftend+(rightend-leftend)*(t-1)/nintervals})
  progression
}


progression <- create_progression(0.01,0.9,1000)
confusionmat_progression <- lapply(progression,function(t){get_similarity_mat(8,t)})
methodbinomp <- lapply(confusionmat_progression,function(t){dcats_fit(counts1,counts2,similarity_mat = t, n_samples = 100)[,9]})
methodbinomp <- do.call(rbind,methodbinomp)
colnames(methodbinomp) <- colnames(misclassM)
write.csv(methodbinomp,"methodbinompExpB.csv") # export data 
methodbinomp <- as.data.frame(methodbinomp)
methodbinomp$confusion_rate <- progression
methodbinomptp <- gather(methodbinomp,key="CellTypes",value="p-value",1:8)
Plotmethodbinomptp <- ggplot(methodbinomptp,aes(x=confusion_rate, y=methodbinomptp$`p-value`, color=CellTypes)) + geom_point() +
  theme_minimal()
ggsave("Plotmethodbinomptp.png", plot=Plotmethodbinomptp,width=30,height=15,units="cm") # save plot 

methodbetabinomp <- lapply(confusionmat_progression,function(t){dcats_betabin(counts1,counts2,similarity_mat = t, n_samples = 100)[,9]})
methodBETAbinomp <- methodbetabinomp
methodBETAbinomp <- do.call(rbind,methodBETAbinomp)
colnames(methodBETAbinomp) <- colnames(misclassM)
write.csv(methodBETAbinomp,"methodBETAbinompExpB.csv") # export data 
methodBETAbinomp <- as.data.frame(methodBETAbinomp)
methodBETAbinomp$confusion_rate <- progression
methodBETAbinomptp <- gather(methodBETAbinomp,key="CellTypes",value="p-value",1:8)
PlotmethodBETAbinomptp <- ggplot(methodBETAbinomptp,aes(x=confusion_rate, y=methodBETAbinomptp$`p-value`, color=CellTypes)) + geom_point() +
  theme_minimal()
ggsave("PlotmethodBETAbinomptp.png", plot=PlotmethodBETAbinomptp,width=30,height=15,units="cm") # save plot 


binombetacompare <- lapply(seq_len(8),function(t){cbind(methodbinomp[,t],methodBETAbinomp[,t])})
binombetacompare <- lapply(binombetacompare,function(t){cbind(t,methodbinomp$confusion_rate)})
cmbinomp <- colnames(methodbinomp)
cmbetap <- colnames(methodBETAbinomp)
cmbinomp <- sapply(cmbinomp, function(t){paste(t,"Binom",sep="_")})
cmbetap <-  sapply(cmbetap, function(t){paste(t,"BetaBinom",sep="_")})

binombetacomparecopy <- binombetacompare

for (i in 1:8){
  names <- append(cmbinomp[i],cmbetap[i])
  names <- append(names,'confusion_rate')
  colnames(binombetacomparecopy[[i]]) <- names
}

for (i in 1:8){
  binombetacomparecopy[[i]] <- gather(as.data.frame(binombetacomparecopy[[i]]),"methods","p_values",1:2)
}

outputnames <- paste("binombetacompare_cluster",1:8,sep="_") 

for (i in 1:8){
  write.csv(binombetacomparecopy[[i]],outputnames[i]) # export data 
}

plotlists <- lapply(binombetacomparecopy,function(t){ggplot(t,aes(x=confusion_rate, y=p_values, color=methods)) + geom_point() +
    theme_minimal()})
ml <- marrangeGrob(plotlists, nrow=2,ncol=4)
ggsave("ml.png", plot=ml,width=100,height=30,units="cm") # save plot 


