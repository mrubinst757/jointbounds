# Regression tests for fixes to the binary-outcome (formerly marbounds) code.

legacy_fit <- function(n = 400, seed = 11, ...) {
  suppressPackageStartupMessages(library(SuperLearner))
  set.seed(seed)
  X <- runif(n, -2, 2)
  A <- rbinom(n, 1, plogis(0.5 * X))
  C <- rbinom(n, 1, 0.2 + 0.1 * A)
  Y <- ifelse(C == 1, NA_real_, rbinom(n, 1, plogis(0.8 * A + 0.8 * X)))
  dat <- data.frame(Y = Y, A = A, C = C, X = X)
  mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X", sl_lib = "SL.glm",
             seed = 1, ...)
}

test_that("V = 1 nuisances are fitted rather than marginal means", {
  skip_if_not_installed("SuperLearner")
  fit <- legacy_fit(estimand = "ate", assumption = "general")
  expect_gt(stats::sd(fit$nuisance$e), 0)
  expect_gt(stats::sd(fit$nuisance$mu1), 0)
})

test_that("pooled outcome model lets mu0 and mu1 differ", {
  skip_if_not_installed("SuperLearner")
  dat <- make_test_data(n = 300)
  nuis <- suppressWarnings(estimate_nuisance(
    as.matrix(dat$X), dat$A, dat$C, dat$Y, V = 2, stratify_mu = FALSE,
    seed = 1
  ))
  expect_gt(max(abs(nuis$mu1 - nuis$mu0)), 0)
})

test_that("smooth bounded_risk bounds converge to the indicator bounds", {
  skip_if_not_installed("SuperLearner")
  fit <- legacy_fit(estimand = "ate", assumption = "bounded_risk",
                    delta_0u = 1, delta_1u = 1, tau_0 = 1.5)
  nuis <- fit$nuisance
  smooth <- compute_bounds(
    fit$phi, "ate", "bounded_risk", delta_0u = 1, delta_1u = 1,
    tau_0 = 1.5, tau_1 = 1.5, mu0 = nuis$mu0, mu1 = nuis$mu1,
    pi0 = nuis$pi0, pi1 = nuis$pi1, smooth_approximation = TRUE,
    epsilon = 1e-4
  )
  expect_equal(smooth$lower, fit$result$lower, tolerance = 1e-3)
  expect_equal(smooth$upper, fit$result$upper, tolerance = 1e-3)
})

test_that("grid delta and tau shorthands are applied", {
  skip_if_not_installed("SuperLearner")
  fit <- legacy_fit(estimand = "ate", assumption = "bounded_delta",
                    param_grid = list(delta = c(0.2, 1)), B = 50)
  lower <- fit$lower_grid$grid$estimate
  expect_length(unique(round(lower, 10)), 2L)

  br <- legacy_fit(estimand = "ate", assumption = "bounded_risk",
                   param_grid = list(tau = c(1.2, 3)), B = 50)
  expect_length(unique(round(br$upper_grid$grid$estimate, 10)), 2L)

  scalar <- legacy_fit(estimand = "ate", assumption = "bounded_risk",
                       delta_0u = 1, delta_1u = 1, tau_0 = 3)
  grid <- legacy_fit(estimand = "ate", assumption = "bounded_risk",
                     param_grid = list(tau_0 = 3), B = 50)
  expect_equal(grid$lower_grid$grid$estimate, scalar$result$lower)
  expect_equal(grid$upper_grid$grid$estimate, scalar$result$upper)
})

test_that("multiplier_bootstrap_grid bounded_risk matches mar_bounds", {
  skip_if_not_installed("SuperLearner")
  fit <- legacy_fit(estimand = "ate", assumption = "bounded_risk",
                    delta_0u = 1, delta_1u = 1, tau_0 = 2)
  set.seed(3)
  up <- multiplier_bootstrap_grid(fit, data.frame(tau_0 = 2), "upper",
                                  "bounded_risk", B = 50)
  lo <- multiplier_bootstrap_grid(fit, data.frame(tau_0 = 2), "lower",
                                  "bounded_risk", B = 50)
  expect_equal(up$grid$estimate, fit$result$upper)
  expect_equal(lo$grid$estimate, fit$result$lower)
})
