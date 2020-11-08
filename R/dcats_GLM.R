#' A Generalised linear model based likelihood ratio testing
#'
#' GLM supports both beta-binomial and negative binomial from aod package.
#'
#' @param count_mat A matrix of compsition sizes (n_sample, n_cluster) for each
#' cluster in each sample
#' @param design_mat A matrix of testing candidate factors (n_sample, n_factor)
#' with same sample order as count_mat
#' @param model A string value: `betabin` for beta-binomial and `negbin` for
#' negative binomial distribution
#' @param base_model A string value: `NULL` for 1 factor vs NULL factor testing;
#' `FULL` for FULL facotrs vs n-1 factors testing. Only `NULL` is supported now.
#'
#' @return a list of significance p values for each cluster
#'
#' @export
#' @import matrixStats
#'
dcats_GLM <- function(count_mat, design_mat, model='betabin', base_model='NULL') {
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

  # Test each factor
  for (m in seq_len(ncol(count_mat))) {
    for (k in seq_len(ncol(design_mat))) {
      df_use <- data.frame(n1 = count_mat[, m], total=rowSums(count_mat))
      df_use <- cbind(df_use, as.data.frame(design_mat)[, k, drop=FALSE])

      df_tmp <- df_use[!is.na(design_mat[, k]), ]

      if (model == 'betabin') {
        fm0 <- aod::betabin(cbind(n1, total-n1) ~ 1, ~ 1, data = df_tmp)

        formula_fix <- as.formula(paste0('cbind(n1, total-n1)', '~ 1+',
                                         colnames(design_mat)[k], sep=''))
        fm1 <- aod::betabin(formula_fix, ~ 1, data = df_tmp)
      } else {
        fm0 <- aod::negbin(n1 ~ total + 1, ~ 1, data = df_tmp)

        formula_fix <- as.formula(paste0('n1 ~  ', colnames(design_mat)[k],
                                         '+ total + 1', sep=''))
        fm1 <- aod::negbin(formula_fix, ~ 1, data = df_tmp)
      }

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
