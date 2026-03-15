#' Simulate mvSuSiE data
#'
#' Generates a simulated dataset for testing mvSuSiE.
#'
#' @param n Number of samples.
#' @param p Number of variables.
#' @param r Number of outcomes.
#' @param s Number of effect variables per outcome if >= 1; proportion
#'   of effect variables per outcome if < 1.
#' @param center_scale If \code{TRUE}, center and scale X and Y.
#' @param y_missing Fraction of entries in Y to set to \code{NA} (per
#'   sample, per outcome). If \code{NULL}, no missing values are created.
#'
#' @return A list with elements \code{X}, \code{y}, \code{y_missing},
#'   \code{d}, \code{n}, \code{p}, \code{r}, \code{V}, and \code{b}.
#'
#' @importFrom stats cov
#' @importFrom stats runif
#' @importFrom stats rnorm
#'
#' @export
mvsusie_sim1 <- function(n = 200, p = 500, r = 2, s = 4,
                         center_scale = FALSE, y_missing = NULL) {
  X <- matrix(rnorm(n * p, 0, 1), n, p)
  if (s >= 1) {
    beta <- matrix(0, p, r)
    for (i in 1:r) {
      beta[sample(1:p, s), i] <- 1
    }
  } else {
    beta <- matrix(runif(p * r) > s, p, r)
  }
  y <- X %*% beta + do.call(cbind, lapply(1:r, function(i) rnorm(n)))
  if (center_scale) {
    X <- scale(X)
    y <- t(t(y) - apply(y, 2, mean))
  }
  if (!is.null(y_missing)) {
    y2 <- y
    for (i in 1:nrow(y2)) {
      y2[i, runif(r) <= y_missing] <- NA
    }
    y_missing <- y2
  }
  scaled_prior_variance <- 0.2
  return(list(
    X = X, y = y, y_missing = y_missing, d = diag(t(X) %*% X),
    n = n, p = p, r = r, V = scaled_prior_variance * cov(y),
    b = beta
  ))
}

#' Simulated multi-trait fine-mapping data
#'
#' Simulated fine-mapping data set used to illustrate mvSuSiE.
#' Includes genotype and phenotype data for 574 samples, 1,001
#' genetic markers and 20 traits. The traits were simulated from
#' the mvSuSiE model with coefficients \code{simdata$Btrue} and
#' residual covariance \code{simdata$par$V}. Three causal genetic
#' variants at positions 255, 335 and 493 have nonzero coefficients.
#'
#' @name simdata
#'
#' @docType data
#'
#' @format \code{simdata} is a list with the following elements:
#'
#' \describe{
#'
#'   \item{raw$X}{The matrix of simulated genotypes.}
#'
#'   \item{raw$Y}{The matrix of simulated traits.}
#'
#'   \item{Btrue}{The coefficients used to simulate the data.}
#'
#'   \item{par$V}{The residual covariance matrix used to simulate the data.}
#'
#'   \item{par$U}{The collection of covariance matrices specifying the
#'     mvsusie prior.}
#'
#'   \item{par$w}{The weights associated with the covariance matrices.}
#'
#'   \item{sumstats$n}{The sample size.}
#'
#'   \item{sumstats$LD}{The LD computed from \code{raw$X}.}
#'
#'   \item{sumstats$bhat}{The least-squares effect estimates from the
#'     single-marker association tests computed using
#'     \code{\link[susieR]{univariate_regression}}. (Note that \code{X}
#'     was standardized before computing \code{bhat}.)}
#'
#'   \item{sumstats$sehat}{The standard errors of the least-squares
#'     effect estimates computed using
#'     \code{\link[susieR]{univariate_regression}}. (Note that \code{X}
#'     was standardized before computing \code{sehat}.)}}
#'
#' @keywords data
#'
NULL
