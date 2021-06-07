library(testthat)
library(DCATS)

test_that("create_simMat generate a matrix we want", {
  stdRes = matrix(c(0.9, rep(0.1/3, 4), 0.9, rep(0.1/3, 4), 0.9, rep(0.1/3, 4), 0.9), ncol = 4)
  
  expect_equal(create_simMat(4, 0.1), stdRes)
  
})

test_that("knn_simMat generate a matrix we want", {
  data("simulation")
  stdResV = c(0.9502, 0.0364, 0.0128, 0.0010, 0.0174, 0.0032, 0.0012, 0.0008, 0.0209, 0.9490, 0.0054, 0.0000, 0.0070, 0.0001, 0.0000, 0.0000, 0.0090, 0.0066, 0.9767, 0.0001, 0.0046, 0.0004, 0.0000, 0.0000, 0.0017, 0.0000, 0.0002, 0.9984, 0.0000, 0.0003, 0.0001, 0.0001, 0.0108, 0.0076, 0.0041, 0.0000, 0.9695, 0.0006, 0.0000, 0.0000, 0.0049, 0.0004, 0.0010, 0.0003, 0.0014, 0.9950, 0.0002, 0.0002, 0.0013, 0.0000, 0.0000, 0.0001, 0.0000, 0.0002, 0.9982, 0.0001, 0.0013, 0.0000, 0.0000, 0.0001, 0.0000, 0.0002, 0.0002, 0.9987)
  stdRes = matrix(stdResV, ncol = 8)
  colnames(stdRes) = c("D", "H", "G", "A", "F", "B", "E", "C")
  rownames(stdRes) = c("D", "H", "G", "A", "F", "B", "E", "C")
  
  knn_mat = knn_simMat(simulation$knnGraphs, simulation$labels)
  
  expect_equal(round(knn_mat, 4), stdRes)
  
})

#> Test passed 🌈