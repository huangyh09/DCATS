
#' Generate similarity matrix with uniform confusion rate to none-self clusters
#'
#' @param K A integer for number of cluster
#' @param confuse_rate A float for confusion rate, uniformly to none-self
#' clusters
#'
#' @return a similarity matrix with uniform confusion with other cluster
#'
#' @export
#'
#' @examples
#' get_similarity_mat(4, 0.1)
get_similarity_mat <- function(K, confuse_rate) {
    diag(K) * (1 - confuse_rate) + confuse_rate * (1 - diag(K)) / (K - 1)
}


#' A binomial regression test with similarity based bootstrapping
#'
#' A GLM test with binomial distribution. In order to estimate the variance of
#' the weight, a boostrapping based on the composition similarity is performed.
#'
#' @param counts1 A matrix of compsition sizes (n_rep1, n_cluster) for each
#' replicate in each cluster for codition 1 as case
#' @param counts2 A matrix of compsition sizes (n_rep2, n_cluster) for each
#' replicate in each cluster for codition 2 as control
#' @param similarity_mat A matrix of floats (n_cluster, n_cluster) for the
#' similarity matrix between cluster group pair
#' @param pseudo_count A pseudo count to add for counts in all cell types.
#' Default NULL means 0 except if a cell type is emplty in one condition,
#' otherwise pseudo_count will be: 0.01 * rowMeans for each condition
#' @param n_samples An integer for number samples in sampling for estimating the
#' variance of the weights
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
#' dcats_fit(sim_dat[[1]], sim_dat[[2]], confuse_mat, n_samples = 100)
#'
dcats_fit <- function(counts1, counts2, similarity_mat=NULL, n_samples=50,
                      pseudo_count=NULL) {
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
    if (is.null(similarity_mat)) {
        n_samples <- 1
    }

    if (!is.null(n_samples) && !is.null(similarity_mat)) {
        counts1_use <- matrix(0, nrow(counts1) * n_samples, K)
        counts2_use <- matrix(0, nrow(counts2) * n_samples, K)
        for (i in seq_len(nrow(counts1))) {
            idx <- seq((i - 1) * n_samples + 1, i * n_samples)
            for (j in seq_len(K)) {
                counts1_use[idx, ] <- (counts1_use[idx, ] +
                                           t(rmultinom(n_samples, counts1[i, j],
                                                       similarity_mat[j, ])))
            }
        }
        for (i in seq_len(nrow(counts2))) {
            idx <- seq((i - 1) * n_samples + 1, i * n_samples)
            for (j in seq_len(K)) {
                counts2_use[idx, ] <- (counts2_use[idx, ] +
                                           t(rmultinom(n_samples, counts2[i, j],
                                                       similarity_mat[j, ])))
            }
        }
    } else{
        counts1_use <- counts1
        counts2_use <- counts2
    }

    ## Binomial regression for each sampling
    coeffs_val <- matrix(NA, n_samples, K)
    coeffs_err <- matrix(NA, n_samples, K)
    intercept_val <- matrix(NA, n_samples, K)
    intercept_err <- matrix(NA, n_samples, K)
    total_all <- c(rowSums(counts1_use), rowSums(counts2_use))
    label_all <- c(rep(1, nrow(counts1_use)), rep(0, nrow(counts2_use)))
    for (ir in seq_len(n_samples)) {
        idx <- seq(1, length(total_all), n_samples) + ir - 1
        for (i in seq_len(K)) {
            n1 <- c(counts1_use[, i], counts2_use[, i])[idx]
            df <- data.frame(n1 = n1, n2 = total_all[idx] - n1,
                             label = label_all[idx])

            model1 <- glm(cbind(n1, n2) ~ label + 1,
                          family = binomial(), data = df)
            coeffs_val[ir, i] <- summary(model1)$coefficients[2, 1]
            coeffs_err[ir, i] <- summary(model1)$coefficients[2, 2]
            intercept_val[ir, i] <- summary(model1)$coefficients[1, 1]
            intercept_err[ir, i] <- summary(model1)$coefficients[1, 2]
        }
    }

    ## Averaging the coeffcients errors
    if (is.null(n_samples) || is.null(similarity_mat) || n_samples == 1) {
        coeff_val_mean <- colMeans(coeffs_val)
        coeff_err_pool <- colMeans(coeffs_err**2)
        intercept_val_mean <- colMeans(intercept_val)
        intercept_err_pool <- colMeans(intercept_err**2)
    } else{
        coeff_val_mean <- colMeans(coeffs_val)
        coeff_err_pool <- colMeans(coeffs_err**2) +
            matrixStats::colSds(coeffs_val) +
            matrixStats::colSds(coeffs_val) / n_samples

        intercept_val_mean <- colMeans(intercept_val)
        intercept_err_pool <- colMeans(intercept_err**2) +
            matrixStats::colSds(intercept_val) +
            matrixStats::colSds(intercept_val) / n_samples
    }

    pvals <- pnorm(-abs(coeff_val_mean) / sqrt(coeff_err_pool))  * 2

    # variance of two independent random variables:
    # Taylor expansion: http://www.stat.cmu.edu/~hseltman/files/ratio.pdf
    # Note, the variance is devided by n_replicates so it's the variance on mean
    # fold_var <- (colMeans(prop1)^2 / colMeans(prop2)^2 *
    #         (matrixStats::colSds(prop1)^2 / colMeans(prop1)^2 / nrow(counts1) +
    #          matrixStats::colSds(prop2)^2 / colMeans(prop2)^2) / nrow(counts2))
    # fold_mean = colMeans(prop1) / colMeans(prop2)

    data.frame(
        "prop1_mean" = colMeans(prop1),
        "prop1_std"  = matrixStats::colSds(prop1),
        "prop2_mean" = colMeans(prop2),
        "prop2_std"  = matrixStats::colSds(prop2),
        "coeff_mean" = coeff_val_mean,
        "coeff_std"  = sqrt(coeff_err_pool),
        "intecept_mean" = intercept_val_mean,
        "intecept_std"  = sqrt(intercept_err_pool),
        "pvals" = pvals,
        row.names = colnames(counts1)
    )
}


