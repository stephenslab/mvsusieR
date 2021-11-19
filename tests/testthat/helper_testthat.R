# @title sets three attributes for matrix X
# @param X an n by p data matrix that can be either a trend filtering
#   matrix or a regular dense/sparse matrix
# @param center boolean indicating centered by column means or not
# @param scale boolean indicating scaled by column standard deviations or not
# @return X with three attributes e.g. `attr(X, 'scaled:center') is a
# p vector of column means of X if center=TRUE, a p vector of zeros
# otherwise. 'attr(X, 'scaled:scale') is a p vector of columan standard
# deviations of X if scale=TRUE, a p vector of 1s otherwise. 'attr(X,
# 'd') is a p vector of column sums of X.standardized^2,' where
# X.standardized is the matrix X centered by attr(X, 'scaled:center')
# and scaled by attr(X, 'scaled:scale').
#
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
set_X_attributes = function (X, center = TRUE, scale = TRUE) {
    
  # if X is a trend filtering matrix
  if (!is.null(attr(X,"matrix.type"))) {
    order = attr(X,"order")
    n = ncol(X)
    
    # Set three attributes for X.
    attr(X,"scaled:center") = compute_tf_cm(order,n)
    attr(X,"scaled:scale") = compute_tf_csd(order,n)
    attr(X,"d") = compute_tf_d(order,n,attr(X,"scaled:center"),
                               attr(X,"scaled:scale"),scale,center)
    if (!center)
      attr(X,"scaled:center") = rep(0,n)
    if (!scale)
      attr(X,"scaled:scale") = rep(1,n)
  } else {
      
    # If X is either a dense or sparse ordinary matrix.
    # Get column means.
    cm = colMeans(X,na.rm = TRUE)
    
    # Get column standard deviations.
    csd = susieR:::compute_colSds(X)
    
    # Set sd = 1 when the column has variance 0.
    csd[csd == 0] = 1
    if (!center)
      cm = rep(0,length = length(cm))
    if (!scale) 
      csd = rep(1,length = length(cm))
    X.std = (t(X) - cm)/csd
    
    # Set three attributes for X.
    attr(X,"d") = rowSums(X.std * X.std)
    attr(X,"scaled:center") = cm
    attr(X,"scaled:scale") = csd
  }
  return(X)
}

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
  X = set_X_attributes(X)
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
    data = SSData$new(R, z, 1, 2, NULL, NULL)
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
  res = mvsusie_sim1(n,p,r,4,center_scale=center_scale,y_missing=y_missing)
  res$L = 10
  return(res)
}

compute_cov_diag <- function(Y){
  covar <- diag(apply(Y, 2, var, na.rm=T))
  return(covar)
}

expect_susieR_equal = function(A, BA, estimate_prior_variance = FALSE, estimate_residual_variance = FALSE, 
                               tol = 1E-8, rss = FALSE) {
  expect_equal(A$alpha, BA$alpha, tolerance = tol)
  expect_equal(A$lbf, BA$lbf, tolerance = tol)
  if (!is.na(A$KL) && !is.na(BA$KL)) expect_equal(A$KL, BA$KL, tolerance = tol)
  expect_equal(A$alpha * A$mu, BA$b1, tolerance = tol)
  expect_equal(A$alpha * A$mu2, BA$b2, tolerance = tol)
  if (!is.na(A$elbo) && !is.na(BA$elbo)) expect_equal(A$elbo, BA$elbo, tolerance = tol)
  if (rss) {
    expect_equal(as.vector(A$Rr), as.vector(BA$fitted), tolerance = tol)
    expect_equal(coef(A)[-1], BA$coef[-1], tolerance = tol)
  } else {
    expect_equal(as.vector(A$fitted), as.vector(BA$fitted), tolerance = tol)
    expect_equal(coef(A), BA$coef, tolerance = tol)
  }
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
