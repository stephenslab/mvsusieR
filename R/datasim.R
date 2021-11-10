#' @title A simple simulation function to simulate some test data
#' @param n number of samples
#' @param p number of features
#' @param r number of conditions
#' @param s number of effect variables per condition if greater than 1; otherwise percentage of effect variables per condition
#' @param center_scale FALSE by default
#' @export
mvsusie_sim1 = function(n=200,p=500,r=2,s=4,center_scale=FALSE,y_missing=NULL) {
  X = matrix(rnorm(n*p,0,1),n,p)
  if (s>=1) {
    beta = matrix(0, p, r)
    for (i in 1:r) beta[sample(1:p,s), i] = 1
  } else {
    beta = matrix(runif(p*r)>s, p, r)
  }
  y = X %*% beta + do.call(cbind, lapply(1:r, function(i) rnorm(n)))
  if (center_scale) {
    X = scale(X)
    y = t(t(y) - apply(y,2,mean))
  }
  if (!is.null(y_missing)) {
    y2 = y
    for (i in 1:nrow(y2)) {
      y2[i,runif(r) <= y_missing] = NA
    }
    y_missing = y2
  }
  scaled_prior_variance = 0.2
  return(list(X=X,y=y,y_missing=y_missing,d=diag(t(X)%*%X), n=n,p=p,r=r,V=scaled_prior_variance * cov(y),b=beta))
}
