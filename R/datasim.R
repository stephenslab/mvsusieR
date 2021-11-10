#' @title Simulate mvSuSiE Data
#'
#' @description A simple simulation function to simulate some test
#'   data.
#' 
#' @param n Number of samples to simulate.
#' 
#' @param p Number of features to simulate.
#' 
#' @param r Number of conditions to simulate.
#' 
#' @param s Number of effect variables per condition if greater than
#'   1; otherwise, proportion of effect variables per condition if less
#'   than 1.
#' 
#' @param center_scale Describe what input argument "center" is for.
#'
#' @param y_missing Describe what input argument "y_missing" is for.
#'
#' @return Describe the outputs here.
#'
#' @importFrom stats cov
#' @importFrom stats runif
#' @importFrom stats rnorm
#'
#' @export
#' 
mvsusie_sim1 = function (n = 200, p = 500, r = 2, s = 4,
                         center_scale = FALSE, y_missing = NULL) {
  X = matrix(rnorm(n*p,0,1),n,p)
  if (s >= 1) {
    beta = matrix(0,p,r)
    for (i in 1:r)
      beta[sample(1:p,s),i] = 1
  } else
    beta = matrix(runif(p*r) > s,p,r)
  y = X %*% beta + do.call(cbind,lapply(1:r,function(i) rnorm(n)))
  if (center_scale) {
    X = scale(X)
    y = t(t(y) - apply(y,2,mean))
  }
  if (!is.null(y_missing)) {
    y2 = y
    for (i in 1:nrow(y2))
      y2[i,runif(r) <= y_missing] = NA
    y_missing = y2
  }
  scaled_prior_variance = 0.2
  return(list(X = X,y = y,y_missing = y_missing,d = diag(t(X) %*% X),
              n = n,p = p,r = r,V = scaled_prior_variance * cov(y),
              b = beta))
}
