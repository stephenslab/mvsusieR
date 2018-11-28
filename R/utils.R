#' @title compute value_j * weight_j / sum(value_j * weight_j)
#' @keywords internal
safe_compute_weight = function(value, weight, log = TRUE) {
    if (!log) value = log(value)
    mvalue = max(value)
    w = exp(value-mvalue)
    w_weighted = w * weight
    weighted_sum_w = sum(w_weighted)
    return(list(alpha = as.vector(w_weighted / weighted_sum_w), log_total = log(weighted_sum_w) + mvalue))
}

#' @title SuSiE model extractor
#' @keywords internal
report_susie_model = function(d, m) {
    s = list(
        alpha = t(m$pip),
        mu = t(do.call(cbind, m$posterior_b1)),
        mu2 = t(do.call(cbind, m$posterior_b2)),
        KL = m$kl,
        lbf = m$lbf,
        sigma2 = m$residual_variance,
        V = m$prior_variance,
        elbo = m$get_objective(dump=TRUE),
        niter = m$get_niter(),
        fitted = d$fitted,
        coef = m$coef(d),
        null_index = -9
        )
    s$intercept = s$coef[1]
    class(s) = 'susie'
    return(s)
}

#' @title A null progressbar, because currently `progressbar_enabled` feature does not work for `progress_bar`
#' @importFrom R6 R6Class
#' @keywords internal
null_progress_bar = R6Class('null_progress_bar', public = list(tick = function(...) {}))

#' @title check if all elements are the same in matrix of J by R, J >> R
#' @keywords internal
is_mat_common = function(mat) {
  all((t(mat) - mat[1,]) == 0)
}

#' @title A simple simulation function to simulate some test data
#' @param n number of samples
#' @param p number of features
#' @param r number of conditions
#' @param s percentage of signals
#' @export
mmbr_sim2 = function(n=200,p=500,r=2,s=0.9) {
  X = matrix(rnorm(n*p,0,1),n,p)
  beta = matrix(runif(p*r)>s, p, r)
  y = X %*% beta + do.call(cbind, lapply(1:r, function(i) rnorm(n)))
  X = scale(X)
  y = y - apply(y,2,mean)
  scaled_prior_variance = 0.2
  return(list(X=X,y=y, d=diag(t(X)%*%X), n=n,p=p,r=r,V=scaled_prior_variance * cov(y), b=beta))
}

#' @title A simple simulation function to simulate some test data
#' @param n number of samples
#' @param p number of features
#' @param r number of conditions
#' @param m number of signals per condition
#' @export
mmbr_sim1 = function(n=200,p=500,r=2,m=4) {
  X = matrix(rnorm(n*p,0,1),n,p)
  beta = matrix(0, p, r)
  for (i in 1:r) beta[sample(1:p,1), i] = 1
  y = X %*% beta + do.call(cbind, lapply(1:r, function(i) rnorm(n)))
  X = scale(X)
  y = y - apply(y,2,mean)
  scaled_prior_variance = 0.2
  return(list(X=X,y=y, d=diag(t(X)%*%X), n=n,p=p,r=r,V=scaled_prior_variance * cov(y), b=beta))
}