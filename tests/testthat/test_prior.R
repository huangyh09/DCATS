library(testthat)
library(DCATS)

test_that("whether simulator works", {
  K <- 2
  totals1 = c(100, 800, 1300, 600)
  totals2 = c(250, 700, 1100)
  diri_s1 = diri_s2 = rep(1, K) * 20
  simil_mat = create_simMat(K, 0.2)
  sim_dat <- simulator_base(totals1, totals2, diri_s1, diri_s2, simil_mat)
  
  expect_equal(dim(sim_dat$numb_cond1)[2], K)
  expect_equal(dim(sim_dat$numb_cond2)[2], K)
  expect_equal(rowSums(sim_dat$numb_cond1), totals1)
  expect_equal(rowSums(sim_dat$numb_cond2), totals2)
})

test_that("whether EM process works", {
  X = c(100, 300, 1500, 500, 1000)
  simMM = create_simMat(5, confuse_rate=0.2)
  res = multinom_EM(X, simMM)
  
  expect_equal(sum(res$mu), 1)
})
