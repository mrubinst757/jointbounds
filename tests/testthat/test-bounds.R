test_that("compute_bounds general ATE returns expected structure", {
  phi <- matrix(runif(600, 0, 1), 100, 6)
  out <- jointbounds:::compute_bounds(phi, estimand = "ate", assumption = "general")
  expect_named(out, c("naive", "lower", "upper", "se_lower", "se_upper"))
  expect_true(all(is.finite(c(out$lower, out$upper, out$se_lower, out$se_upper))))
  expect_true(out$se_lower >= 0)
  expect_true(out$se_upper >= 0)
})

test_that("compute_bounds bounded_delta ATE works with delta < 1", {
  phi <- matrix(runif(600, 0, 1), 100, 6)
  out <- jointbounds:::compute_bounds(phi, estimand = "ate", assumption = "bounded_delta",
                                    delta_0u = 0.8, delta_1u = 0.8)
  expect_named(out, c("naive", "lower", "upper", "se_lower", "se_upper"))
  expect_true(all(is.finite(unlist(out))))
})

test_that("compute_bounds monotonicity_pos returns valid bounds", {
  phi <- matrix(runif(600, 0, 1), 100, 6)
  out <- jointbounds:::compute_bounds(phi, estimand = "ate", assumption = "monotonicity_pos",
                                    delta_0u = 0.7, delta_1u = 0.7)
  expect_true(all(is.finite(unlist(out))))
  expect_named(out, c("naive", "lower", "upper", "se_lower", "se_upper"))
})

test_that("compute_bounds point_ate returns estimate and se", {
  phi <- matrix(runif(600, 0, 1), 100, 6)
  out <- jointbounds:::compute_bounds(phi, estimand = "ate", assumption = "point_ate",
                                    delta_0 = 0.5, delta_1 = 0.5, tau = 2)
  expect_named(out, c("estimate", "se", "naive"))
  expect_length(out$estimate, 1)
  expect_true(is.finite(out$estimate))
  expect_true(out$se >= 0)
})

test_that("compute_bounds psi1 general returns bounds", {
  phi <- matrix(runif(600, 0, 1), 100, 6)
  out <- jointbounds:::compute_bounds(phi, estimand = "psi1", assumption = "general")
  expect_named(out, c("naive", "lower", "upper", "se_lower", "se_upper"))
  expect_true(all(is.finite(c(out$lower, out$upper, out$se_lower, out$se_upper))))
})

test_that("compute_bounds point_psi2 errors without pi0 and mu1", {
  phi <- matrix(runif(60, 0, 1), 10, 6)
  expect_error(
    jointbounds:::compute_bounds(phi, estimand = "psi2", assumption = "point_psi2", delta_0 = 0.5),
    "pi0 and mu1"
  )
})

test_that("compute_bounds point_psi2 works with pi0 and mu1", {
  phi <- matrix(runif(60, 0, 1), 10, 6)
  pi0 <- runif(10, 0.1, 0.3)
  mu1 <- runif(10, 0.2, 0.6)
  out <- jointbounds:::compute_bounds(phi, estimand = "psi2", assumption = "point_psi2",
                                    delta_0 = 0.5, pi0 = pi0, mu1 = mu1)
  expect_named(out, c("estimate", "se", "naive"))
  expect_true(is.finite(out$estimate))
})

test_that("zero informative-missingness radius collapses Psi2 bounds", {
  set.seed(18)
  phi <- matrix(stats::rnorm(600), 100, 6)
  mu0 <- stats::runif(100, -.2, .4)
  mu1 <- stats::runif(100, 0, .7)
  pi0 <- stats::runif(100, .05, .4)
  for (smooth in c(FALSE, TRUE)) {
    out <- jointbounds:::compute_bounds(phi, estimand = "psi2",
      assumption = "bounded_delta", delta_0u = 0,
      mu0 = mu0, mu1 = mu1, pi0 = pi0,
      smooth_approximation = smooth, epsilon = .01)
    expect_equal(out$lower, out$naive)
    expect_equal(out$upper, out$naive)
  }
})

test_that("coefficient vectors have length 6", {
  expect_length(jointbounds:::coef_general_lower_ate(), 6)
  expect_length(jointbounds:::coef_general_upper_ate(), 6)
  expect_length(jointbounds:::coef_delta_lower_ate(0.5, 0.5), 6)
  expect_length(jointbounds:::coef_mono_pos_lower_ate(0.5), 6)
  expect_length(jointbounds:::coef_point_ate(0.5, 0.5, 2), 6)
})
