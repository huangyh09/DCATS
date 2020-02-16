# Functions borrowed / inspired from DiTaSIC

#' Estimate count function adapted from DiTaSIC
#'
#' @param mat A matrix of floats (n_group, n_group) for similarity matrix
#' between clusters
#' @param observedCounts A vector of integers (n_group) for the observed sample
#' count in each group
#' @param filtering A bool value for ...
#' @param pvalThreshold A float for ...
#'
#' @return a list containing \code{glm.coeffs}, a vector of estimated
#' coefficients; \code{error.abs}, a vector of standard divation of estimated
#' coefficients; \code{mod.exclude}, an indictor verctor of a cluster is
#' excluded; \code{raw.pval}, a vector of p values for each coefficient.
#'
#' @export
#'
#' @examples
EstimateCount <- function(mat, observedCounts, filtering = FALSE,
                          pvalThreshold = 0.05, minThreshold = 0){

  #create a glm model: observedCounts = mat*r +\epsilon
  glm.model <- glm(observedCounts ~ 0 + mat,
                   family = "poisson"(link = "identity"),
                   maxit = 100, trace = TRUE)

  glm.coeffs <- glm.model$coefficients
  error.abs  <- summary(glm.model)$coefficients[, 2] # std for coeff
  raw.pval   <- summary(glm.model)$coefficients[, 4] # p-val for coeff

  if (filtering == T) {
    mod.exclude1 <- which(summary(glm.model)$coefficients[, 4] > pvalThreshold)
    mod.exclude2 <- which(abs(glm.coeffs) < minThreshold)
    mod.exclude <- sort(union(mod.exclude1, mod.exclude2))
    glm.coeffs[mod.exclude] <- 0
  }
  else {
    mod.exclude <- NULL
  }

  #redo glm in case some of the coefficients are negative
  if (sum(glm.coeffs < 0) > 0) {
    glm.coeffs <- optim.fit(observedCounts, mat, start = glm.coeffs)$par
  }

  list(glm.coeffs, error.abs, mod.exclude, raw.pval)
}


#' Accessory functions to estimate count
#' @param y, please specify
#' @param z, please specify
#' @param start, please specify
#'
#' @export
#' @return please specify
optim.fit <- function(y, z, start=rep(mean(y), ncol(z)), ...){
  deviance <- function(beta, y, z){
    mu <- z %*% beta; 2*sum(mu - y - y*log(mu/y))
  }
  grad <- function(beta, y, z){
    mu <- z %*% beta; 2 * t(1 - y/mu) %*% z
  }
  nlminb(start, deviance, grad, lower = 0, y = y, z = z,  ...)
}


#' Create empirical distribution adapted from DiTaSIC
#'
create.empdist <- function(abundance_objects, norm.fac=1) {
  abs.coeff   <- abundance_objects[[1]]
  abs.error   <- abundance_objects[[2]]
  mod.exclude <- abundance_objects[[3]]
  lambdas <- list()
  filter <- ifelse(c(1:length(abs.coeff)) %in% mod.exclude, "FALSE", "TRUE")

  a.coeff <- abs.coeff[filter == TRUE] / norm.fac # estimates of abundant taxa
  a.error <- abs.error[filter == TRUE] / sqrt(norm.fac)

  lower <- a.coeff - a.error
  upper <- a.coeff + a.error

  # 1.Step: Sample potential poisson lambdas from interval
  # [estimate +/- std.error] for every taxa
  set.seed(1442)
  N.lambda <- 100
  for (k in 1:length(a.coeff)) {
    if (lower[k] < 0) lower[k] <- 0
    lambdas[[k]] <- runif(N.lambda, min = lower[k], max = upper[k])
  }

  # 2.Step:
  # for each abundant taxa (list entry): and for each of sampled lambda's (=mean
  # read number of taxa):
  # 500(N.poi) draws from lambda-defined poisson distribution -> specific lambda
  # empirical distribution
  # Pooling all empiricial distributions of all lambdas -> overall empirical
  # ditribution for each taxa

  set.seed(1442)
  N.poi <- 500
  empdist <- lapply(lambdas, function(i){as.vector(sapply(i, function(lambd) {
    rpois(N.poi, lambd) }))})

  # All taxa data: 'empdist.all' list object (0 = not ex. Taxa, 50000 poisson
  # draws for existent taxa)
  empdist.all <- sapply(names(abs.coeff), function(x) NULL)
  empdist.all[filter == FALSE] <- 0
  empdist.all[filter == TRUE] <- empdist

  return(empdist.all)
}


#' Differential abundance analysis main function adapted from DiTASIC
#'
#' @export
DiffExpr <- function(empdist1, empdist2, norm.coeff_s1, norm.coeff_s2,
                     threshold.minreads =0, seed=1448) {

  glm.coeff.1 <- norm.coeff_s1
  glm.coeff.2 <- norm.coeff_s2
  N.draws <- 50000  # N.draws = N.lambda * N.poi (defined in <create.empdist>)

  result <- data.frame("norm.count.estimate.1" = NA,
                       "norm.count.estimate.2" = NA,
                       "diff.pval" = NA, "pval.adj" = NA, "FC.log2" = NA)

  for (k in 1:length(glm.coeff.1))   {

    if (length(empdist1[[k]]) == 1) {
      result[k, "norm.count.estimate.1"] <- NA
      if (length(empdist2[[k]]) == 1) {# case: o o
        result[k, "norm.count.estimate.2"] <- NA
      } else {# case: o x
        # normalized glm estimates
        result[k, "norm.count.estimate.2"] <- glm.coeff.2[k]

        # sig. for existence: > min.threshold input reads
        diff.emp <- empdist2[[k]] - threshold.minreads
        pval2 <- ecdf(diff.emp)(0)
        result[k, "diff.pval"] <- pval2
      }
    } else {
      if (length(empdist2[[k]]) == 1) {# case: x o
        # normalized glm estimates
        result[k, "norm.count.estimate.2"] <- NA;
        result[k, "norm.count.estimate.1"] <- glm.coeff.1[k]

        # sig. for existence: > min.threshold input reads
        diff.emp <- empdist1[[k]] - threshold.minreads
        pval1 <- ecdf(diff.emp)(0)
        result[k, "diff.pval"] <- pval1
      } else {# case: x x
        # normalized glm estimates
        result[k, "norm.count.estimate.1"] <- glm.coeff.1[k]
        result[k, "norm.count.estimate.2"] <- glm.coeff.2[k]
        result[k, "FC.log2"] <- log2(glm.coeff.2[k]/glm.coeff.1[k])

        set.seed(seed)
        emp1 <- sample(empdist1[[k]], N.draws)
        emp2 <- sample(empdist2[[k]], N.draws)
        diff.emp <- emp1 - emp2   # empirical difference distribution

        pval <- ecdf(diff.emp)(0)
        if (pval > 0.5) pval <- 1 - pval
        result[k, "diff.pval"] <- pval
      }
    }
  }

  result[, "pval.adj"] <- p.adjust(result[, "diff.pval"], "BH")
  return(result)
}



