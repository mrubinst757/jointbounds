test_that("influence_functions returns n x 6 matrix", {
  n <- 40
  Y <- rbinom(n, 1, 0.4)
  A <- rbinom(n, 1, 0.5)
  C <- rbinom(n, 1, 0.2)
  e <- runif(n, 0.2, 0.8)
  pi0 <- runif(n, 0.05, 0.3)
  pi1 <- runif(n, 0.05, 0.3)
  mu0 <- runif(n, 0.1, 0.5)
  mu1 <- runif(n, 0.2, 0.6)
  phi <- jointbounds:::influence_functions(Y, A, C, e, pi0, pi1, mu0, mu1)
  expect_equal(dim(phi), c(n, 6))
  expect_equal(colnames(phi), c("phi_1_0", "phi_1_1", "phi_2_0", "phi_2_1", "phi_3_0", "phi_3_1"))
  expect_true(all(is.finite(phi)))
})

test_that("influence_functions handles NA in Y", {
  n <- 20
  Y <- rep(NA_real_, n)
  Y[1:10] <- rbinom(10, 1, 0.5)
  A <- c(rep(0, 10), rep(1, 10))
  C <- c(rep(0, 10), rep(1, 10))
  e <- rep(0.5, n)
  pi0 <- rep(0.2, n)
  pi1 <- rep(0.2, n)
  mu0 <- rep(0.4, n)
  mu1 <- rep(0.5, n)
  phi <- jointbounds:::influence_functions(Y, A, C, e, pi0, pi1, mu0, mu1)
  expect_equal(dim(phi), c(n, 6))
  expect_true(all(is.finite(phi)))
})

test_that("theta_hat returns length-6 vector", {
  phi <- matrix(rnorm(60), 10, 6)
  th <- jointbounds:::theta_hat(phi)
  expect_length(th, 6)
  expect_equal(th, colMeans(phi))
})

test_that("phi_var returns 6x6 positive semi-definite matrix", {
  phi <- matrix(rnorm(120), 20, 6)
  V <- jointbounds:::phi_var(phi)
  expect_equal(dim(V), c(6, 6))
  expect_true(all(eigen(V, only.values = TRUE)$values >= -1e-10))
})
