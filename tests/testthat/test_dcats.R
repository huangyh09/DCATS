library(testthat)
library(DCATS)

test_that("whether dcats_betabin works as we expected", {
  set.seed(3)
  K <- 2
  totals1 = rep(300, 4)
  totals2 = rep(300, 3)
  diri_s1 = rep(1, K) * 20
  diri_s2 = rep(1, K) * 20
  simil_mat = create_simMat(K, confuse_rate=0.2)
  sim_dat <- simulator_base(totals1, totals2, diri_s1, diri_s2, simil_mat)
  res = dcats_betabin(sim_dat[[1]], sim_dat[[2]], simil_mat, n_samples = 100)
  
  expect_equal(res$pvals > 0.05, rep(TRUE, K))
  
})

test_that("whether dcats_betabin works as we expected", {
  set.seed(3)
  K <- 2
  totals1 = rep(300, 4)
  totals2 = rep(300, 3)
  diri_s1 = c(1/6, 5/6) * 20
  diri_s2 = rep(1/K, K) * 20
  simil_mat = create_simMat(K, confuse_rate=0.2)
  sim_dat <- simulator_base(totals1, totals2, diri_s1, diri_s2, simil_mat)
  res = dcats_betabin(sim_dat[[1]], sim_dat[[2]], simil_mat, n_samples = 100)
  
  expect_equal(res$pvals > 0.05, rep(FALSE, K))
  
})

test_that("whether dcats_betabin works without similarity matrix", {
  set.seed(3)
  K <- 2
  totals1 = rep(300, 4)
  totals2 = rep(300, 3)
  diri_s1 = c(1/6, 5/6) * 20
  diri_s2 = rep(1/K, K) * 20
  simil_mat = create_simMat(K, confuse_rate=0.2)
  sim_dat <- simulator_base(totals1, totals2, diri_s1, diri_s2, simil_mat)
  res = dcats_betabin(sim_dat[[1]], sim_dat[[2]], n_samples = 100)
  
  expect_equal(res$pvals > 0.05, rep(FALSE, K))
  
})

# test_that("whether dcats_betabin works as we expected", {
#   set.seed(3)
#   K <- 3
#   totals1 = rep(300, 4)
#   totals2 = rep(300, 3)
#   diri_s1 = rep(1/K, K) * 200
#   diri_s2 = c(1/3, 1/2, 1/6) * 200
#   simil_mat = create_simMat(K, confuse_rate=0.1)
#   sim_dat <- DCATS::simulator_base(totals1, totals2, diri_s1, diri_s2, simil_mat)
#   sim_count = rbind(sim_dat$numb_cond1, sim_dat$numb_cond2)
#   sim_design = matrix(c("g1", "g1", "g1", "g2", "g2", "g2", "g2"), ncol = 1)
#   res = dcats_GLM(sim_count, sim_design, similarity_mat = NULL)
#   
#   expect_equal(as.vector(res$pvals) > 0.05, c(TRUE, FALSE, FALSE))
#   
# })

