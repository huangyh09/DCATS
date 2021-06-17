#' A Generalised linear model based likelihood ratio testing
#'
#' GLM supports both beta-binomial and negative binomial from aod package.
#'
#' @param count_mat A matrix of composition sizes (n_sample, n_cluster) for each
#'   cluster in each sample
#' @param design_mat A matrix of testing candidate factors (n_sample, n_factor)
#'   with same sample order as count_mat
#' @param similarity_mat A matrix of floats (n_cluster, n_cluster) for the
#'   similarity matrix between cluster group pair. The order of cluster should
#'   be consistent with those in `count_mat`.
#' @param pseudo_count A pseudo count to add for counts in all cell types.
#'   Default NULL means 0 except if a cell type is empty in one condition,
#'   otherwise pseudo_count will be: 0.01 * rowMeans for each condition
#' @param n_samples An integer for number samples in sampling for estimating the
#'   variance of the weights
#' @param base_model A string value: `NULL` for 1 factor vs NULL factor testing;
#'   `FULL` for FULL factors vs n-1 factors testing. Only `NULL` is supported
#'   now.
#'
#' @return a list of significance p values for each cluster
#'
#' @export
#' @import matrixStats
#'
#' @examples
#' K <- 3
#' totals1 = c(100, 800, 1300, 600)
#' totals2 = c(250, 700, 1100)
#' diri_s1 = rep(1, K) * 20
#' diri_s2 = rep(1, K) * 20
#' simil_mat = create_simMat(K, confuse_rate=0.2)
#' sim_dat <- DCATS::simulator_base(totals1, totals2, diri_s1, diri_s2, simil_mat)
#' sim_count = rbind(sim_dat$numb_cond1, sim_dat$numb_cond2)
#' sim_design = matrix(c("g1", "g1", "g1", "g1", "g2", "g2", "g2"), ncol = 1)
#' dcats_GLM(sim_count, sim_design, similarity_mat = simil_mat)
#' 
dcats_GLM <- function(count_mat, design_mat, similarity_mat=NULL, n_samples=50,
                      pseudo_count=NULL,  base_model='NULL') {
  # Output matrices
  coeffs     <- matrix(NA, ncol(count_mat), ncol(design_mat))
  coeffs_err <- matrix(NA, ncol(count_mat), ncol(design_mat))
  LR_vals    <- matrix(NA, ncol(count_mat), ncol(design_mat))
  LRT_pvals  <- matrix(NA, ncol(count_mat), ncol(design_mat))
  LRT_fdr    <- matrix(NA, ncol(count_mat), ncol(design_mat))

  # Check colnames
  if (is.null(colnames(count_mat)))
    colnames(count_mat) <- paste0('cell_type_', seq(ncol(count_mat)))
  if (is.null(colnames(design_mat)))
    colnames(design_mat) <- paste0('factor_', seq(ncol(design_mat)))

  # Add rownames and colnames
  rownames(LR_vals) <- rownames(LRT_pvals) <- rownames(LRT_fdr) <-
    rownames(coeffs) <- rownames(coeffs_err) <- colnames(count_mat)
  colnames(LR_vals) <- colnames(LRT_pvals) <- colnames(LRT_fdr) <-
    colnames(coeffs) <- colnames(coeffs_err) <- colnames(design_mat)
  
  
  ## using estimated the latent cell counts
  count_latent = count_mat
  if(!is.null(similarity_mat)) {
    for (i in seq_len(nrow(count_mat))) {
      count_latent[i, ] <- sum(count_mat[i, ]) *
        multinom_EM(count_mat[i, ], similarity_mat, verbose = FALSE)$mu
    }
  }
  
  K <- ncol(count_mat) ## number of cell types
  if (is.null(similarity_mat)) {
    n_samples <- 1
  }
  
  if (!is.null(n_samples) && !is.null(similarity_mat)) {
    count_use <- matrix(0, nrow(count_mat) * n_samples, K)
    for (i in seq_len(nrow(count_mat))) {
      idx <- seq((i - 1) * n_samples + 1, i * n_samples)
      for (j in seq_len(K)) {
        count_use[idx, ] <- (
          count_use[idx, ] + t(rmultinom(n_samples, count_latent[i, j], similarity_mat[j, ])))
      }
    }
  } else{
    count_use <- count_mat
  }
  
  # adding pseudo counts
  if (is.null(pseudo_count)) {
    if (any(colMeans(count_mat) == 0)) {
      print(paste("Empty cell type exists in at least one conidtion;",
                  "adding replicate & condition specific pseudo count:"))
      count_use <- count_use + 1
    }
  } else {
    count_use = count_use + pseudo_count
  }
  
  count_use = round(count_use)
  
  # Test each factor
  for (m in seq_len(ncol(count_use))) {
    for (k in seq_len(ncol(design_mat))) {
      df_use <- data.frame(n1 = count_use[, m], total=rowSums(count_use))
      df_use <- cbind(df_use, as.data.frame(design_mat)[, k, drop=FALSE])

      df_tmp <- df_use[!is.na(design_mat[, k]), ]
      
      ## model fitting using betabin
      fm0 <- aod::betabin(cbind(n1, total-n1) ~ 1, ~ 1, data = df_tmp)
      formula_fix <- as.formula(paste0('cbind(n1, total-n1)', '~ 1+',
                                         colnames(design_mat)[k], sep=''))
      fm1 <- aod::betabin(formula_fix, ~ 1, data = df_tmp)

      ## ignore the fitting if the hessian matrix is singular
      if (length(fm1@varparam) < 4 || is.na(fm1@varparam[2, 2])) {next}

      LR_vals[m, k] <- fm0@dev - fm1@dev
      LRT_pvals[m, k] <- pchisq(LR_vals[m, k], df=1, lower.tail = FALSE, log.p = FALSE)

      coeffs[m, k] <- fm1@param[2]
      coeffs_err[m, k] <- fm1@varparam[2, 2]
    }
  }

  # Return list
  LRT_fdr[,] <- p.adjust(LRT_pvals, method = 'fdr')
  res <- list('ceoffs'=coeffs, 'coeffs_err'=coeffs_err,
              'LR_vals'=LR_vals, 'pvals'=LRT_pvals, 'fdr'=LRT_fdr)
  res
}
