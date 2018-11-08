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

expect_susie_equal = function(A, BA, estimate_prior_variance = FALSE, estimate_residual_variance = FALSE) {
  expect_equal(A$alpha, BA$alpha)
  expect_equal(A$lbf, BA$lbf)
  expect_equal(A$KL, BA$KL)
  expect_equal(A$alpha * A$mu, BA$mu)
  expect_equal(A$alpha * A$mu2, BA$mu2)
  expect_equal(A$elbo[-1], BA$elbo)
  expect_equal(A$fitted, BA$fitted)
  expect_equal(coef(A), BA$coef)
  if (estimate_residual_variance) expect_equal(A$sigma2, BA$sigma2)
  if (estimate_prior_variance) expect_equal(A$V, BA$V)
}
