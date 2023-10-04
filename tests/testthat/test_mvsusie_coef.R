context("fit_glmpca_pois")

test_that("coef() gives the write coefficients with 1 condition",{

  # Simulate a 500 x 100 data set with 1 response.
  set.seed(1)
  n <- 500
  p <- 100
  maf <- c(c(0.5,0.2,0.1,0.05),0.05 + 0.45*runif(96))
  X   <- (runif(n*p) < maf) +
         (runif(n*p) < maf)
  X   <- matrix(as.double(X),n,p,byrow = TRUE)
  rownames(X) <- paste0("s",1:n)
  colnames(X) <- paste0("p",1:p)
  b <- rep(0,p)
  b[1:4] <- 3
  Y <- -1 + X %*% b + rnorm(n)

  # Fit an mvsusie model to the data.
  fit <- mvsusie(X,Y,L = 10,standardize = TRUE)

  # The estimated coefficients (including the intercept) should
  # closely match the coefficients used to simulate the data.
  expect_gt(cor(coef(fit),c(-1,b)),0.999)
  #
  # print(round(head(cbind(b,coef(fit)[-1]),n = 6),digits = 4))
  # plot(b,coef(fit)[-1],pch = 20,xlab = "true coef",ylab = "estimated coef")
  # abline(a = 0,b = 1,pch = 20,lty = "dotted",col = "magenta")

  # Run the same test again, but now without standardizing X.
  fit <- mvsusie(X,Y,L = 10,standardize = FALSE)
  expect_gt(cor(coef(fit),c(-1,b)),0.999)
})

test_that("coef() gives the write coefficients with 3 conditions",{

  # Simulate a 500 x 100 data set with 3 outcomes.
  set.seed(1)
  n <- 500
  p <- 100
  r <- 3
  n1 <- 3
  maf <- c(c(0.5,0.2,0.1,0.05),0.05 + 0.45*runif(96))
  X   <- (runif(n*p) < maf) +
         (runif(n*p) < maf)
  X   <- matrix(as.double(X),n,p,byrow = TRUE)
  rownames(X) <- paste0("s",1:n)
  colnames(X) <- paste0("p",1:p)
  b1 <- rep(0,p)
  b2 <- rep(0,p)
  b3 <- rep(0,p)
  b1[1:n1] <- 3*rnorm(n1)
  b2[1:n1] <- 3*rnorm(n1)
  b3[1:n1] <- 3*rnorm(n1)
  y1 <- -1 + X %*% b1 + rnorm(n)
  y2 <- +1 + X %*% b2 + rnorm(n)
  y3 <- -2 + X %*% b3 + rnorm(n)
  b  <- cbind(b1,b2,b3)
  Y  <- cbind(y1,y2,y3)

  # Fit an mvsusie model to the data.
  prior <- create_mixture_prior(R = 3)
  fit <- mvsusie(X,Y,prior_variance = prior,standardize = TRUE)

  # The estimated coefficients (including the intercept) should
  # closely match the coefficients used to simulate the data.
  expect_gt(cor(as.vector(b[1:4,]),as.vector(coef(fit)[2:5,])),0.999)
  #
  # print(coef(fit)[1,])
  # plot(b,coef(fit)[-1,],pch = 20,xlab = "true coef",ylab = "estimated coef")
  # abline(a = 0,b = 1,pch = 20,lty = "dotted",col = "magenta")

  # Run the same test again, but now without standardizing X.
  fit <- mvsusie(X,Y,L = 10,prior_variance = prior,standardize = FALSE)
  expect_gt(cor(as.vector(b[1:4,]),as.vector(coef(fit)[2:5,])),0.999)
})
