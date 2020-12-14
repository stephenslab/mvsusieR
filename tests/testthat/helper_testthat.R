create_sparsity_mat = function(sparsity, n, p){
  nonzero = round(n*p*(1-sparsity))
  nonzero.idx = sample(n*p, nonzero)
  mat = numeric(n*p)
  mat[nonzero.idx] = 1
  mat = matrix(mat, nrow=n, ncol=p)
  return(mat)
}

simulate_univariate = function(n=100, p=200, sparse=F, summary = F) {
  set.seed(1)
  beta = rep(0,p)
  beta[1:4] = 1
  if (sparse) {
    X = create_sparsity_mat(0.99,n,p)
    X.sparse = as(X,'dgCMatrix')
  } else {
    X = matrix(rnorm(n*p,3,4),n,p)
    X.sparse = NA
  }
  y = c(X %*% beta + rnorm(n))
  L = 10
  residual_variance = 1
  scaled_prior_variance = 0.2
  X = susieR:::set_X_attributes(X)
  y = y - mean(y)
  if(summary){
    sumstats = susieR:::univariate_regression(X,y)
    z = sumstats$betahat/sumstats$sebetahat
    R = cor(X)
  }
  s = list(alpha=matrix(1/p,nrow=L,ncol=p),
           mu=matrix(2,nrow=L,ncol=p),
           mu2=matrix(3,nrow=L,ncol=p),
           Xr=rep(5,n), KL=rep(1.2,L),
           sigma2=residual_variance,
           V=scaled_prior_variance * as.numeric(var(y)))
  if (sparse) {
    # FIXME: sparse data not supported
    data = NA
  } else if(summary){
    data = RSSData$new(z, R, eigenR=NULL, tol=1e-08)
    data$set_residual_variance(residual_variance)
  }else {
    data = DenseData$new(X,y)
    data$standardize(TRUE,TRUE)
    data$set_residual_variance(residual_variance)
  }
  if(summary){
    return(list(X=X, X.sparse=X.sparse, z = z, R = R, s=s, d=data, y=y, n=n, p=p, V=s$V, b=beta, L=L))
  }else{
    return(list(X=X, X.sparse=X.sparse, s=s, d=data, y=y, n=n, p=p, V=s$V, b=beta, L=L))
  }

}

simulate_multivariate = function(n=100,p=100,r=2,center_scale=TRUE,y_missing=0) {
  set.seed(1)
  res = mmbr_sim1(n,p,r,4,center_scale=center_scale,y_missing=y_missing)
  res$L = 10
  return(res)
}

compute_cov_diag <- function(Y){
  covar <- diag(apply(Y, 2, var, na.rm=T))
  return(covar)
}

expect_susieR_equal = function(A, BA, estimate_prior_variance = FALSE, estimate_residual_variance = FALSE, tol = 1E-8) {
  expect_equal(A$alpha, BA$alpha, tolerance = tol)
  expect_equal(A$lbf, BA$lbf, tolerance = tol)
  if (!is.na(A$KL) && !is.na(BA$KL)) expect_equal(A$KL, BA$KL, tolerance = tol)
  expect_equal(A$alpha * A$mu, BA$b1, tolerance = tol)
  expect_equal(A$alpha * A$mu2, BA$b2, tolerance = tol)
  if (!is.na(A$elbo) && !is.na(BA$elbo)) expect_equal(A$elbo, BA$elbo, tolerance = tol)
  expect_equal(as.vector(A$fitted), as.vector(BA$fitted), tolerance = tol)
  expect_equal(coef(A), BA$coef, tolerance = tol)
  if (estimate_residual_variance) expect_equal(A$sigma2, BA$sigma2, tolerance = tol)
  if (estimate_prior_variance) expect_equal(A$V, BA$V, tolerance = tol)
}

expect_susie_equal = function(A, B, estimate_prior_variance = FALSE, estimate_residual_variance = FALSE, tol = 1E-8) {
  expect_equal(A$alpha, B$alpha, tolerance = tol)
  if (!is.na(A$elbo) && !is.na(B$elbo)) expect_equal(A$lbf, B$lbf, tolerance = tol)
  expect_equal(A$b1, B$b1, tolerance = tol)
  expect_equal(A$b2, B$b2, tolerance = tol)
  expect_equal(A$fitted, B$fitted, tolerance = tol)
  expect_equal(A$coef, B$coef, tolerance = tol)
  if (!is.na(A$KL) && !is.na(B$KL)) expect_equal(A$KL, B$KL, tolerance = tol)
  if (!is.na(A$elbo) && !is.na(B$elbo)) expect_equal(A$elbo, B$elbo, tolerance = tol)
  if (estimate_residual_variance) expect_equal(A$sigma2, B$sigma2, tolerance = tol)
  if (estimate_prior_variance) expect_equal(A$V, B$V, tolerance = tol)
}