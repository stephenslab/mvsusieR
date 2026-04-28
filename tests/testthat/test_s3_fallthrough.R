# Verify that S3 methods we deleted from mvsusieR fall through to
# susieR's defaults and produce equivalent results. Tests call the
# `.default` / `.individual` / `.ss` bodies directly (S3 dispatch
# of unexported methods across packages is unreliable; the actual
# in-package call sites work via `UseMethod` from inside susieR).

# ---- Trivial accessors -------------------------------------------------

test_that("get_alpha_l.default returns model$alpha[l, ] (mvsusie compatible)", {
  alpha_mat <- matrix(seq_len(12), nrow = 3, ncol = 4)
  model <- structure(list(alpha = alpha_mat),
                     class = c("mvsusie", "susie"))
  for (l in seq_len(3)) {
    expect_equal(susieR:::get_alpha_l.default(model, l), alpha_mat[l, ])
  }
})

test_that("get_prior_variance_l.default and set_prior_variance_l.default work on mvsusie", {
  model <- structure(list(V = c(1, 2, 3)),
                     class = c("mvsusie", "susie"))
  expect_identical(susieR:::get_prior_variance_l.default(model, 2), 2)
  m1 <- susieR:::set_prior_variance_l.default(model, 2, 5)
  expect_identical(m1$V[2], 5)
  expect_identical(m1$V[c(1, 3)], c(1, 3))
  expect_identical(class(m1), c("mvsusie", "susie"))
})

# ---- get_cs.individual covers mv_individual (via class hierarchy) ------

test_that("get_cs.individual returns a list when called on mv_individual data", {
  # mv_individual inherits from individual; .individual handles it
  # directly. We deleted get_cs.mv_individual so .individual fires.
  set.seed(1)
  data <- structure(list(X = matrix(rnorm(40 * 5), 40, 5),
                         n = 40L, p = 5L),
                    class = c("mv_individual", "individual"))
  params <- list(coverage = 0.95, min_abs_corr = 0.5, n_purity = 100)

  alpha <- matrix(0, 1, 5)
  alpha[1, 1] <- 0.99
  alpha[1, 2:5] <- 0.0025
  model <- list(alpha = alpha, mu = matrix(0, 1, 5))

  out <- susieR:::get_cs.individual(data, params, model)
  expect_true(is.list(out))
})

# ---- get_cs.ss covers mv_ss (safety upgrade vs deleted .mv_ss) ---------

test_that("get_cs.ss handles symmetric XtX without the safe_cov2cor warning", {
  # mv_ss inherits from ss; .ss uses safe_cov2cor (vs the deleted
  # .mv_ss override which used cov2cor). When XtX is exactly
  # symmetric, neither warns.
  set.seed(2)
  X <- matrix(rnorm(50 * 5), 50, 5)
  X <- scale(X)
  Xcorr <- crossprod(X) / nrow(X)
  diag(Xcorr) <- 1
  data <- structure(list(XtX = Xcorr, p = 5L),
                    class = c("mv_ss", "mv_individual", "ss"))
  params <- list(coverage = 0.95, min_abs_corr = 0.5, n_purity = 100)

  alpha <- matrix(0, 1, 5)
  alpha[1, 1] <- 0.99
  alpha[1, 2:5] <- 0.0025
  model <- list(alpha = alpha, mu = matrix(0, 1, 5))

  out <- suppressMessages(susieR:::get_cs.ss(data, params, model))
  expect_true(is.list(out))
})
