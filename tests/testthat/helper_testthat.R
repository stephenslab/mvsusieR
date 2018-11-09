create_sparsity_mat = function(sparsity, n, p){
  nonzero = round(n*p*(1-sparsity))
  nonzero.idx = sample(n*p, nonzero)
  mat = numeric(n*p)
  mat[nonzero.idx] = 1
  mat = matrix(mat, nrow=n, ncol=p)
  return(mat)
}

simulate_univariate = function(n=100, p=200, sparse=F) {
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
  X = susieR:::safe_colScale(X)
  y = y - mean(y)
  s = list(alpha=matrix(1/p,nrow=L,ncol=p),
           mu=matrix(2,nrow=L,ncol=p),
           mu2=matrix(3,nrow=L,ncol=p),
           Xr=rep(5,n), KL=rep(1.2,L),
           sigma2=residual_variance,
           V=scaled_prior_variance * as.numeric(var(y)))
  if (sparse) {
    # FIXME: sparse data not supported
    data = NA
  } else {
    data = DenseData$new(X,y,TRUE,TRUE)
  }
  return(list(X=X, X.sparse=X.sparse, s=s, d=data, y=y, n=n, p=p, V=s$V, b=beta, L=L))
}

simulate_multivariate = function(n=100,p=100,r=2) {
  set.seed(1)
  X = matrix(rnorm(n*p,3,4),n,p)
  beta = matrix(0, p, r)
  for (i in 1:r) beta[1+(i-1)*4 : 4*i, r] = 1
  y = X %*% beta + do.call(cbind, lapply(1:r, function(i) rnorm(n)))
  X = susieR:::safe_colScale(X)
  y = y - apply(y,2,mean)
  scaled_prior_variance = 0.2
  L = 10
  data = DenseData$new(X,y,TRUE,TRUE)
  return(list(X=X,y=y,n=n,p=p,r=r,V=scaled_prior_variance * cov(y),b=beta,L=L))
}

expect_susieR_equal = function(A, BA, estimate_prior_variance = FALSE, estimate_residual_variance = FALSE, tol = 1E-10) {
  expect_equal(A$alpha, BA$alpha, tolerance = tol)
  # do not compare lbf when not using the same convergence check -- because it is rather sensitive
  if (!is.null(A$elbo) && !is.null(BA$elbo)) expect_equal(A$lbf, BA$lbf, tolerance = tol)
  if (!is.null(A$KL) && !is.null(BA$KL)) expect_equal(A$KL, BA$KL, tolerance = tol)
  expect_equal(A$alpha * A$mu, BA$mu, tolerance = tol)
  expect_equal(A$alpha * A$mu2, BA$mu2, tolerance = tol)
  if (!is.null(A$elbo) && !is.null(BA$elbo)) expect_equal(A$elbo[-1], BA$elbo, tolerance = tol)
  expect_equal(A$fitted, BA$fitted, tolerance = tol)
  expect_equal(coef(A), BA$coef, tolerance = tol)
  if (estimate_residual_variance) expect_equal(A$sigma2, BA$sigma2, tolerance = tol)
  if (estimate_prior_variance) expect_equal(A$V, BA$V, tolerance = tol)
}

expect_susie_equal = function(A, B, estimate_prior_variance = FALSE, estimate_residual_variance = FALSE, tol = 1E-10) {
  expect_equal(A$alpha, B$alpha, tolerance = tol)
  if (!is.null(A$elbo) && !is.null(B$elbo)) expect_equal(A$lbf, B$lbf, tolerance = tol)
  expect_equal(A$mu, B$mu, tolerance = tol)
  expect_equal(A$mu2, B$mu2, tolerance = tol)
  expect_equal(A$fitted, B$fitted, tolerance = tol)
  expect_equal(A$coef, B$coef, tolerance = tol)
  if (!is.null(A$KL) && !is.null(B$KL)) expect_equal(A$KL, B$KL, tolerance = tol)
  if (!is.null(A$elbo) && !is.null(B$elbo)) expect_equal(A$elbo, B$elbo, tolerance = tol)
  if (estimate_residual_variance) expect_equal(A$sigma2, B$sigma2, tolerance = tol)
  if (estimate_prior_variance) expect_equal(A$V, B$V, tolerance = tol)
}
