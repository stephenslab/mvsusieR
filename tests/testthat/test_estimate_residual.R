context("Test Estimate Residual Variance")

# NOTE: Tests using R6 SingleEffectModel/BayesianSimpleRegression/
# BayesianMultivariateRegression/SuSiE/DenseData have been removed because
# R6 classes have been deleted. Those tests are superseded by
# test_r1_susieR_identity.R which compares the S3 path against susieR for R=1.

test_that("estimated residual variance: matrix prior vs mash_prior agree", with(simulate_multivariate(r = 1), {
  # Compare S3 matrix prior vs S3 mash_prior (both go through S3 path)
  null_weight = 0
  mash_init = create_mash_prior(Ulist = list(V), grid = 1, null_weight = null_weight)
  CA_s3 = mvsusie(X, y, L = L, prior_variance = V,
                  residual_variance = cov(y),
                  estimate_residual_variance = T,
                  estimate_prior_variance = FALSE,
                  )
  DA = mvsusie(X, y, L = L, prior_variance = mash_init,
               residual_variance = cov(y),
               estimate_residual_variance = T,
               estimate_prior_variance = FALSE,
               )
  expect_susie_equal(CA_s3, DA, F, T)
}))
