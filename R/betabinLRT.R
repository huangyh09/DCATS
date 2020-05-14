#' A Likelihood ratio test based on beta-binomial regression for cell counts
#'
#' A GLM test with binomial distribution. In order to estimate the variance of
#' the weight, a boostrapping based on the composition similarity is performed.
#'
#' @param counts1 A matrix of compsition sizes (n_rep1, n_cluster) for each
#' replicate in each cluster for codition 1 as case
#' @param counts2 A matrix of compsition sizes (n_rep2, n_cluster) for each
#' replicate in each cluster for codition 2 as control
#' @param pseudo_count A pseudo count to add for counts in all cell types.
#' Default NULL means 0 except if a cell type is emplty in one condition,
#' otherwise pseudo_count will be: 0.01 * rowMeans for each condition
#' @param binom_only A bool value to turn the model into binoimal only and no
#' beta over dispersion
#'
#' @return a vector of significance p values for each cluster
#'
#' @export
#' @import matrixStats
#'
#' @examples
#' K <- 2
#' totals1 = c(100, 800, 1300, 600)
#' totals2 = c(250, 700, 1100)
#' diri_s1 = rep(1, K) * 20
#' diri_s2 = rep(1, K) * 20
#' simil_mat = get_similarity_mat(K, confuse_rate=0.2)
#' sim_dat <- DCATS::simulator_base(totals1, totals2, diri_s1, diri_s2, simil_mat)
#' DCATS::betabinLRT(sim_dat[[1]], sim_dat[[2]])
#'
betabinLRT <- function(counts1, counts2, pseudo_count=NULL, binom_only=FALSE) {
    ## Check counts1 and counts2 shape
    if (length(counts1) == 1 || is.null(dim(counts1)) ||
        length(dim(counts1)) < 2) {
        counts1 = matrix(counts1, nrow=1)
    }
    if (length(counts2) == 1 || is.null(dim(counts2)) ||
        length(dim(counts2)) < 2) {
        counts2 = matrix(counts2, nrow=1)
    }

    ## add base pseudo count if zero cells for all replicate
    if (is.null(pseudo_count)) {
        if (any(colMeans(counts1) == 0) || any(colMeans(counts2) == 0) ) {
            print(paste("Empty cell type exists in at least one conidtion;",
                        "adding replicate & condition specific pseudo count:"))
            print(0.01 * rowMeans(counts1))
            print(0.01 * rowMeans(counts2))

            counts1 = counts1 + 0.01 * rowMeans(counts1)
            counts2 = counts2 + 0.01 * rowMeans(counts2)
        }
    } else {
        counts1 = counts1 + pseudo_count
        counts2 = counts2 + pseudo_count
    }
    prop1 <- counts1 / rowSums(counts1)
    prop2 <- counts2 / rowSums(counts2)

    ## number of cell types
    K <- ncol(counts1)

    ## Binomial regression for each sampling
    intercept <- rep(NA, K)
    intercept_err <- rep(NA, K)

    coeffs_val <- rep(NA, K)
    coeffs_err <- rep(NA, K)
    LR_vals <- rep(NA, K)
    LRT_pvals <- rep(NA, K)

    totals <- c(rowSums(counts1), rowSums(counts2))
    labels <- c(rep(1, nrow(counts1)), rep(0, nrow(counts2)))

    if (binom_only) {
        for (i in seq_len(K)) {
            n1 <- c(counts1[, i], counts2[, i])
            df_tmp <- data.frame(n1 = n1, n2 = totals - n1, label=labels)
            model0 <- glm(cbind(n1, n2) ~ 1, family = binomial(),
                          data = df_tmp)
            model1 <- glm(cbind(n1, n2) ~ labels + 1, family = binomial(),
                          data = df_tmp)

            intercept[i] <- summary(model1)$coefficients[1, 1]
            intercept_err[i] <- summary(model1)$coefficients[1, 2]

            coeffs_val[i] <- summary(model1)$coefficients[2, 1]
            coeffs_err[i] <- summary(model1)$coefficients[2, 2]

            LR_vals[i] <- model0$deviance - model1$deviance
            LRT_pvals[i] <-  pchisq(LR_vals[i], df=1, lower.tail = FALSE,
                                    log.p = FALSE)
        }
    } else{
        for (i in seq_len(K)) {
            n1 <- c(counts1[, i], counts2[, i])
            df_tmp <- data.frame(n1 = n1, n2 = totals - n1, label=labels)

            fm0 <- aod::betabin(cbind(n1, n2) ~ 1, ~ 1, data = df_tmp)
            fm1 <- aod::betabin(cbind(n1, n2) ~ labels, ~ 1, data = df_tmp)

            intercept[i] <- fm1@param[1]      # summary(fm1)@Coef[1, 1]
            intercept_err[i] <- fm1@varparam[1, 1]

            coeffs_val[i] <- fm1@param[2]
            coeffs_err[i] <-fm1@varparam[2, 2]

            LR_vals[i] <- fm0@dev - fm1@dev
            LRT_pvals[i] <-  pchisq(LR_vals[i], df=1, lower.tail = FALSE,
                                    log.p = FALSE)
        }
    }

    data.frame(
        "prop1_mean" = colMeans(prop1),
        "prop1_std"  = matrixStats::colSds(prop1),
        "prop2_mean" = colMeans(prop2),
        "prop2_std"  = matrixStats::colSds(prop2),
        "coeff_mean" = coeffs_val,
        "coeff_std"  = sqrt(coeffs_err),
        "intecept_mean" = intercept,
        "intecept_std"  = sqrt(intercept_err),
        "pvals" = LRT_pvals,
        "LR" = LR_vals,
        row.names = colnames(counts1)
    )
}
