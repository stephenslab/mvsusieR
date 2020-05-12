context("Test MASH regression")

test_that("Degenerated mash regression is identical to univariate BR", with(simulate_multivariate(r=1), {
    # Run univariate BMR
    prior_var = V[1,1]
    residual_var = as.numeric(var(y))
    data = DenseData$new(X,y)
    data$standardize(TRUE,TRUE)
    data$set_residual_variance(residual_var)
    A = BayesianSimpleRegression$new(ncol(X), prior_var)
    A$fit(data, save_summary_stats = T)
    # Run MASH BMR
    # null_weight = 1 - 1 / ncol(X)
    null_weight = 0
    mash_init = MashInitializer$new(list(V), 1, 1 - null_weight, null_weight)
    residual_covar = cov(y)
    data$set_residual_variance(residual_covar)
    B = MashRegression$new(ncol(X), mash_init)
    B$fit(data, save_summary_stats = T)
    # compare result
    expect_equal(A$bhat, B$bhat)
    expect_equal(A$sbhat, B$sbhat)
    expect_equal(A$posterior_b1, B$posterior_b1)
    expect_equal(A$posterior_b2, B$posterior_b2)
    expect_equal(A$lbf, B$lbf)
}))

test_that("Single component mash regression is identical to multivariate BR", with(simulate_multivariate(r=3), {
    # Run multivariate regression
    residual_var = cov(y)
    data = DenseData$new(X,y)
    data$standardize(TRUE,TRUE)
    data$set_residual_variance(residual_var)
    A = BayesianMultivariateRegression$new(ncol(X), V)
    A$fit(data, save_summary_stats = T)
    # Run MASH regression EE model
    mash_init = MashInitializer$new(list(V), 1)
    B = MashRegression$new(ncol(X), mash_init)
    B$fit(data, save_summary_stats = T)
    # compare result
    expect_equal(A$bhat, B$bhat)
    expect_equal(A$sbhat, B$sbhat)
    expect_equal(A$posterior_b1, B$posterior_b1)
    # expect_equal cannot propoerly compare posterior_b2 a 3D array
    expect_equal(as.vector(A$posterior_b2), as.vector(B$posterior_b2))
    expect_equal(A$lbf, B$lbf)
    # mix-up X and don't scale it such that sbhat will be different
    X = matrix(runif(ncol(X) * nrow(X)), nrow(X), ncol(X))
    data = DenseData$new(X,y)
    data$set_residual_variance(residual_var)
    A = BayesianMultivariateRegression$new(ncol(X), V)
    A$fit(data, save_summary_stats = T)
    # Run MASH regression EE model
    mash_init = MashInitializer$new(list(V), 1)
    B = MashRegression$new(ncol(X), mash_init)
    B$fit(data, save_summary_stats = T)
    # compare result
    expect_equal(A$bhat, B$bhat)
    expect_equal(A$sbhat, B$sbhat)
    expect_equal(A$posterior_b1, B$posterior_b1)
    # expect_equal cannot properly compare posterior_b2 a 3D array
    expect_equal(as.vector(A$posterior_b2), as.vector(B$posterior_b2))
    expect_equal(A$lbf, B$lbf)
}))

test_that("Mash regression + precomputed cov is identical to not precompute", with(simulate_multivariate(r=2), {
    # Run univariate BMR
    prior_var = V[1,1]
    residual_covar = cov(y)
    data = DenseData$new(X,y)
    data$standardize(TRUE,TRUE)
    data$set_residual_variance(residual_covar)
    null_weight = 0

    A_init = MashInitializer$new(list(V), 1, 1 - null_weight, null_weight)
    A = MashRegression$new(ncol(X), A_init)
    A$fit(data, save_summary_stats = T)
    B_init = MashInitializer$new(list(V), 1, 1 - null_weight, null_weight)
    data$set_residual_variance(residual_covar, precompute_covariances = T)
    B_init$precompute_cov_matrices(data, algorithm = 'cpp')
    B = MashRegression$new(ncol(X), B_init)
    B$fit(data, save_summary_stats = T)
    # compare result
    expect_equal(A$bhat, B$bhat)
    expect_equal(A$posterior_b1, B$posterior_b1)
    expect_equal(A$posterior_b2, B$posterior_b2)
    expect_equal(A$lbf, B$lbf)
}))

test_that("Degenerated mash regression is identical to univariate BR RSS", with(simulate_multivariate(r=1), {
  # Run univariate BMR
  prior_var = V[1,1]
  residual_var = as.numeric(var(y))
  ss = susieR:::univariate_regression(X, y)
  z = ss$betahat/ss$sebetahat
  R = cor(X)
  data = RSSData$new(z,R,1e-08)
  data$set_residual_variance(residual_var)
  A = BayesianSimpleRegression$new(ncol(X), prior_var)
  A$fit(data, save_summary_stats = T)
  # Run MASH BMR
  # null_weight = 1 - 1 / ncol(X)
  null_weight = 0
  mash_init = MashInitializer$new(list(V), 1, 1 - null_weight, null_weight)
  residual_covar = cov(y)
  data$set_residual_variance(residual_covar)
  B = MashRegression$new(ncol(X), mash_init)
  B$fit(data, save_summary_stats = T)
  # compare result
  expect_equal(A$bhat, B$bhat)
  expect_equal(A$sbhat, B$sbhat)
  expect_equal(A$posterior_b1, B$posterior_b1)
  expect_equal(A$posterior_b2, B$posterior_b2)
  expect_equal(A$lbf, B$lbf)
}))

test_that("Single component mash regression is identical to multivariate BR RSS", with(simulate_multivariate(r=3), {
  # Run multivariate regression
  residual_var = cov(y)
  z = sapply(1:ncol(y), function(j){
    ss = susieR:::univariate_regression(X, y[,j])
    ss$betahat/ss$sebetahat
    })
  R = cor(X)
  data = RSSData$new(z,R, 1e-08)
  data$set_residual_variance(residual_var)
  A = BayesianMultivariateRegression$new(ncol(X), V)
  A$fit(data, save_summary_stats = T)
  # Run MASH regression EE model
  mash_init = MashInitializer$new(list(V), 1)
  B = MashRegression$new(ncol(X), mash_init)
  B$fit(data, save_summary_stats = T)
  # compare result
  expect_equal(A$bhat, B$bhat)
  expect_equal(A$sbhat, B$sbhat)
  expect_equal(A$posterior_b1, B$posterior_b1)
  # expect_equal cannot propoerly compare posterior_b2 a 3D array
  expect_equal(as.vector(A$posterior_b2), as.vector(B$posterior_b2))
  expect_equal(A$lbf, B$lbf)
  # mix-up X and don't scale it such that sbhat will be different
  X = matrix(runif(ncol(X) * nrow(X)), nrow(X), ncol(X))
  z = sapply(1:ncol(y), function(j){
    ss = susieR:::univariate_regression(X, y[,j])
    ss$betahat/ss$sebetahat
  })
  R = cor(X)
  data = RSSData$new(z,R,1e-08)
  data$set_residual_variance(residual_var)
  A = BayesianMultivariateRegression$new(ncol(X), V)
  A$fit(data, save_summary_stats = T)
  # Run MASH regression EE model
  mash_init = MashInitializer$new(list(V), 1)
  B = MashRegression$new(ncol(X), mash_init)
  B$fit(data, save_summary_stats = T)
  # compare result
  expect_equal(A$bhat, B$bhat)
  expect_equal(A$sbhat, B$sbhat)
  expect_equal(A$posterior_b1, B$posterior_b1)
  # expect_equal cannot properly compare posterior_b2 a 3D array
  expect_equal(as.vector(A$posterior_b2), as.vector(B$posterior_b2))
  expect_equal(A$lbf, B$lbf)
}))

test_that("Mash regression + precomputed cov is identical to not precompute (RSS)", with(simulate_multivariate(r=2), {
  # Run univariate BMR
  prior_var = V[1,1]
  residual_covar = cov(y)
  z = sapply(1:ncol(y), function(j){
    ss = susieR:::univariate_regression(X, y[,j])
    ss$betahat/ss$sebetahat
  })
  R = cor(X)
  data = RSSData$new(z,R,1e-08)
  data$set_residual_variance(residual_covar)
  null_weight = 0

  A_init = MashInitializer$new(list(V), 1, 1 - null_weight, null_weight)
  A = MashRegression$new(ncol(X), A_init)
  A$fit(data, save_summary_stats = T)
  B_init = MashInitializer$new(list(V), 1, 1 - null_weight, null_weight)
  data$set_residual_variance(residual_covar, precompute_covariances = T)
  B_init$precompute_cov_matrices(data, algorithm = 'cpp')
  B = MashRegression$new(ncol(X), B_init)
  B$fit(data, save_summary_stats = T)
  # compare result
  expect_equal(A$bhat, B$bhat)
  expect_equal(A$posterior_b1, B$posterior_b1)
  expect_equal(A$posterior_b2, B$posterior_b2)
  expect_equal(A$lbf, B$lbf)
}))
