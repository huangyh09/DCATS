library(DCATS)
library(ggplot2)
count_mat <- as.matrix(read.table("examples/ALM_data/subjectinfomat_AML.txt"))
simMM <- as.matrix(read.table("examples/ALM_data/numemisclassM_AML.txt"))
simMM_uni <- get_similarity_mat(K=21, confuse_rate = 0.05)

heat_matrix(simMM)
heat_matrix(simMM_uni, digits = 3)

dcats_fit(count_mat[1, , drop=FALSE], count_mat[2:5, ])
dcats_fit(count_mat[1, , drop=FALSE], count_mat[2:5, ], simMM_uni)
dcats_fit(count_mat[1, , drop=FALSE], count_mat[2:5, ], simMM)
