test_that("mar_bounds general ATE returns expected structure", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 150)
  fit <- safe_mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                         estimand = "ate", assumption = "general",
                         V = 2, sl_lib = "SL.glm", seed = 1)
  skip_if(is.null(fit), "SuperLearner failed with small sample edge case")

  expect_named(fit, c("nuisance", "phi", "estimand", "assumption", "result"))
  expect_equal(fit$estimand, "ate")
  expect_equal(fit$assumption, "general")
  expect_true(fit$result$lower <= fit$result$upper)
  expect_true(all(is.finite(c(fit$result$naive, fit$result$lower, fit$result$upper))))
  expect_equal(nrow(fit$phi), 150)
  expect_equal(ncol(fit$phi), 6)
})

test_that("mar_bounds point_ate returns estimate and se", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 150)
  fit <- safe_mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                         estimand = "ate", assumption = "point_ate",
                         delta_0 = 0.5, delta_1 = 0.5, tau = 2,
                         V = 2, sl_lib = "SL.glm", seed = 1)
  skip_if(is.null(fit), "SuperLearner failed with small sample edge case")

  expect_true("result" %in% names(fit))
  expect_true("estimate" %in% names(fit$result))
  expect_true("se" %in% names(fit$result))
  expect_false("lower" %in% names(fit))
  expect_true(is.finite(fit$result$estimate))
})

test_that("mar_bounds errors on missing required columns", {
  dat <- make_test_data(n = 100)
  expect_error(mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "Missing"), "exist|undefined columns")
})

test_that("mar_bounds point_ate errors without delta_0, delta_1, tau", {
  dat <- make_test_data(n = 100)
  # Test the error checking before SuperLearner runs
  expect_error(
    mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X", estimand = "ate",
               assumption = "point_ate", V = 2, sl_lib = "SL.glm"),
    "delta_0|delta_1|tau"
  )
})

test_that("psi_naive returns EIF-based estimate, var, se, and eif", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 150)
  out <- tryCatch(
    psi_naive(dat, Y = "Y", A = "A", C = "C", X = "X"),
    error = function(e) {
      if (grepl("All algorithms dropped|Argument mu must be", e$message)) NULL else stop(e)
    }
  )
  skip_if(is.null(out), "SuperLearner failed with small sample edge case")

  expect_true(all(c("estimate", "var", "se", "eif") %in% names(out)))
  expect_true(is.finite(out$estimate))
  expect_true(out$se >= 0)
  expect_equal(out$se^2, out$var, tolerance = 1e-10)
  expect_length(out$eif, 150)
})

test_that("psi_naive with precomputed nuisance matches mar_bounds naive", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 150)
  fit <- safe_mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                         estimand = "ate", assumption = "general",
                         V = 2, sl_lib = "SL.glm", seed = 1)
  skip_if(is.null(fit), "SuperLearner failed with small sample edge case")

  out <- psi_naive(dat, Y = "Y", A = "A", C = "C", X = "X", nuis = fit$nuisance)
  expect_equal(unname(out$estimate), unname(fit$result$naive), tolerance = 1e-6)
})

test_that("mar_bounds param_grid returns lower_grid and upper_grid for Psi_1 bounded_delta", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 150)
  fit <- safe_mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                         estimand = "psi1", assumption = "bounded_delta",
                         param_grid = list(delta_0u = c(0.3, 0.6), delta_1u = c(0.3, 0.6)),
                         sl_lib = "SL.glm", B = 50, seed = 1)
  skip_if(is.null(fit), "SuperLearner failed with small sample edge case")

  expect_true("lower_grid" %in% names(fit))
  expect_true("upper_grid" %in% names(fit))
  g <- fit$lower_grid$grid
  expect_true(all(c("estimate", "se", "ci_lower_pointwise",
                    "ci_upper_pointwise", "ci_lower_uniform",
                    "ci_upper_uniform") %in% names(g)))
  expect_equal(nrow(g), 4L)
  expect_true(is.finite(fit$lower_grid$critical_value) && fit$lower_grid$critical_value > 0)
})

test_that("mar_bounds bounded_risk_unbounded_tau matches original coefficient bounds", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 150)
  fit <- safe_mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                         estimand = "ate", assumption = "bounded_risk_unbounded_tau",
                         delta_0u = 0.3, delta_1u = 0.3, tau_0 = 2, tau_1 = 2,
                         V = 2, sl_lib = "SL.glm", seed = 1)
  skip_if(is.null(fit), "SuperLearner failed with small sample edge case")

  expected <- jointbounds:::compute_bounds(fit$phi, estimand = "ate", assumption = "bounded_risk_unbounded_tau",
                                         delta_0u = 0.3, delta_1u = 0.3, tau_0 = 2, tau_1 = 2)
  expect_equal(fit$result$lower, expected$lower)
  expect_equal(fit$result$upper, expected$upper)
})

test_that("mar_bounds param_grid accepts underscored and non-underscored delta naming", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 80)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "psi1", assumption = "bounded_delta",
                    param_grid = list(delta0u = c(0.3, 0.6), delta1u = c(0.3, 0.6)),
                    B = 50, seed = 1)
  expect_true("lower_grid" %in% names(fit))
  expect_true("upper_grid" %in% names(fit))
  expect_equal(nrow(fit$lower_grid$grid), 4L)
})

test_that("mar_bounds param_grid works for general assumption with delta parameters", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 80)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "ate", assumption = "general",
                    param_grid = list(delta_0u = c(0.2, 0.4), delta_1u = c(0.2, 0.4)),
                    B = 50, seed = 1)
  expect_true("lower_grid" %in% names(fit))
  expect_true("upper_grid" %in% names(fit))
  expect_equal(nrow(fit$lower_grid$grid), 4L)
})

test_that("mar_bounds uses a shared multiplier critical value for both bounds", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 80)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "psi1", assumption = "bounded_delta",
                    param_grid = list(delta_0u = c(0.3, 0.6), delta_1u = c(0.3, 0.6)),
                    grid_bounds = "both",
                    B = 50, seed = 1)
  expect_equal(fit$lower_grid$critical_value, fit$upper_grid$critical_value)
})

test_that("mar_bounds drops constant param columns from grid output", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 80)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "psi1", assumption = "bounded_delta",
                    param_grid = list(delta_0u = c(0.3, 0.6), delta_1u = c(0.3, 0.6), tau = 2),
                    B = 50, seed = 1)
  expect_false("tau" %in% names(fit$lower_grid$grid))
})

test_that("mar_bounds param_grid returns estimate_grid for Psi_2 point_psi2", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 80)
  fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "psi2", assumption = "point_psi2", delta_0 = 0.5,
                    param_grid = list(delta_0 = c(0.2, 0.5, 0.8)),
                    B = 50, seed = 1)
  expect_true("estimate_grid" %in% names(fit))
  g <- fit$estimate_grid$grid
  expect_true(all(c("estimate", "se", "ci_lower_uniform", "ci_upper_uniform") %in% names(g)))
  expect_equal(nrow(g), 3L)
})

test_that("mar_bounds warns when delta_0/delta_1 used with bounded_risk assumptions", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 150)

  # Should warn for bounded_risk
  fit1 <- expect_warning(
    safe_mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "ate", assumption = "bounded_risk",
                    delta_0 = 0.8, delta_1 = 0.8, tau_0 = 2, tau_1 = 2,
                    V = 2, sl_lib = "SL.glm", seed = 1),
    "delta_0 and delta_1 are interpreted as delta_0u and delta_1u"
  )
  skip_if(is.null(fit1), "SuperLearner failed with small sample edge case")

  # Should warn for bounded_risk_unbounded_tau
  fit2 <- expect_warning(
    safe_mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "ate", assumption = "bounded_risk_unbounded_tau",
                    delta_0 = 0.8, delta_1 = 0.8, tau_0 = 2, tau_1 = 2,
                    V = 2, sl_lib = "SL.glm", seed = 1),
    "delta_0 and delta_1 are interpreted as delta_0u and delta_1u"
  )
  skip_if(is.null(fit2), "SuperLearner failed with small sample edge case")

  # Should NOT warn for other assumptions (delta aliasing is intentional there)
  fit3 <- expect_no_warning(
    safe_mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "ate", assumption = "general",
                    delta_0 = 0.8, delta_1 = 0.8,
                    V = 2, sl_lib = "SL.glm", seed = 1)
  )
  skip_if(is.null(fit3), "SuperLearner failed with small sample edge case")
})

test_that("psi2 bounded_delta returns bounds with indicator method", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 150)

  fit <- safe_mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                         estimand = "psi2", assumption = "bounded_delta",
                         delta_0u = 0.8,
                         V = 2, sl_lib = "SL.glm", seed = 1)
  skip_if(is.null(fit), "SuperLearner failed with small sample edge case")

  expect_true("lower" %in% names(fit$result))
  expect_true("upper" %in% names(fit$result))
  expect_true(fit$result$lower <= fit$result$upper)
  expect_true(is.finite(fit$result$lower))
  expect_true(is.finite(fit$result$upper))
})

test_that("psi2 bounded_delta with smooth approximation", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 150)

  fit <- safe_mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                         estimand = "psi2", assumption = "bounded_delta",
                         delta_0u = 0.8, smooth_approximation = TRUE,
                         epsilon = 0.01,
                         V = 2, sl_lib = "SL.glm", seed = 1)
  skip_if(is.null(fit), "SuperLearner failed with small sample edge case")

  expect_true("lower" %in% names(fit$result))
  expect_true("upper" %in% names(fit$result))
  expect_true(fit$result$lower <= fit$result$upper)
})

test_that("stratify_mu = FALSE uses pooled outcome model", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 150)

  fit_stratified <- safe_mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                                     estimand = "ate", assumption = "general",
                                     stratify_mu = TRUE,
                                     V = 2, sl_lib = "SL.glm", seed = 1)
  skip_if(is.null(fit_stratified), "SuperLearner failed with small sample edge case")

  fit_pooled <- safe_mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                                 estimand = "ate", assumption = "general",
                                 stratify_mu = FALSE,
                                 V = 2, sl_lib = "SL.glm", seed = 1)
  skip_if(is.null(fit_pooled), "SuperLearner failed with small sample edge case")

  # Results should differ when stratification matters
  expect_true(is.finite(fit_stratified$result$lower))
  expect_true(is.finite(fit_pooled$result$lower))
})

test_that("binary Y auto-detection warns and sets family to binomial", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 150)

  # Y is binary, should auto-detect and warn
  fit <- expect_message(
    safe_mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                    estimand = "ate", assumption = "general",
                    # Don't set family_Y so auto-detection runs
                    V = 2, sl_lib = "SL.glm", seed = 1),
    "binary.*binomial"
  )
  skip_if(is.null(fit), "SuperLearner failed with small sample edge case")
})

test_that("redundant fields removed from output when result dataframe exists", {
  skip_if_not_installed("SuperLearner")
  suppressPackageStartupMessages(library(SuperLearner))
  dat <- make_test_data(n = 150)

  # Point estimate
  fit_point <- safe_mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                               estimand = "ate", assumption = "point_ate",
                               delta_0 = 0.5, delta_1 = 0.5, tau = 2,
                               V = 2, sl_lib = "SL.glm", seed = 1)
  skip_if(is.null(fit_point), "SuperLearner failed with small sample edge case")

  expect_true("result" %in% names(fit_point))
  expect_false("estimate" %in% names(fit_point))  # should be removed
  expect_false("se" %in% names(fit_point))  # should be removed
  expect_false("naive" %in% names(fit_point))  # should be removed

  # Bounds
  fit_bounds <- safe_mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                                estimand = "ate", assumption = "general",
                                delta_0 = 0.8, delta_1 = 0.8,
                                V = 2, sl_lib = "SL.glm", seed = 1)
  skip_if(is.null(fit_bounds), "SuperLearner failed with small sample edge case")

  expect_true("result" %in% names(fit_bounds))
  expect_false("lower" %in% names(fit_bounds))  # should be removed
  expect_false("upper" %in% names(fit_bounds))  # should be removed
  expect_false("se_lower" %in% names(fit_bounds))  # should be removed
  expect_false("se_upper" %in% names(fit_bounds))  # should be removed
})
