context("mvsusie_predict")

test_that("predict() gives accurate estimates in 1 condition",{

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

  # The mvsusie model predictions should be close to the true values
  # (an RMSE close to 1).
  Yest <- predict(fit,X)
  rmse <- sqrt(mean((Y - Yest)^2))
  expect_lt(rmse,1.1)
  # 
  # plot(Y,Yest,pch = 20,xlab = "true Y",ylab = "estimated Y")
  # abline(a = 0,b = 1,pch = 20,lty = "dotted",col = "magenta")

  # Note that the "fitted" output currently seems to be incorrect:
  #
  #   print(sqrt(mean((Y - fit$fitted)^2)))
  #   plot(Y,fit$fitted,pch = 20,xlab = "true Y",ylab = "estimated Y")
  #   abline(a = 0,b = 1,pch = 20,lty = "dotted",col = "magenta")
  #
})

test_that("predict() gives accurate estimates in 3 conditions",{

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

  # The mvsusie model predictions should be close to the true values
  # (an RMSE close to 1).
  Yest <- predict(fit,X)
  rmse <- sqrt(mean((Y - Yest)^2))
  expect_lt(rmse,1.1)
  # 
  # plot(Y,Yest,pch = 20,xlab = "true Y",ylab = "estimated Y")
  # abline(a = 0,b = 1,pch = 20,lty = "dotted",col = "magenta")

  # Note that the "fitted" output currently seems to be incorrect:
  #
  #   print(sqrt(mean((Y - fit$fitted)^2)))
  #   plot(Y,fit$fitted,pch = 20,xlab = "true Y",ylab = "estimated Y")
  #   abline(a = 0,b = 1,pch = 20,lty = "dotted",col = "magenta")
  #
})
