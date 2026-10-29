make_paper_band_fixture <- function() {
  v <- rep(c(-.01, .01), 20)
  grid <- data.frame(delta_M = 0, delta_K = c(0, 1))
  fits <- lapply(c(0, 1), function(s) list(
    ate = data.frame(endpoint = c("lower", "upper"), estimator = "eif",
      estimate = c(1 - 2*s, 1 + 2*s)),
    scores = list(eif_lower = v, eif_upper = -v)))
  surface <- structure(list(grid = grid, fits = fits, bound_method = "cs"),
    class = "jointbounds_l2_sensitivity_band")
  bench <- structure(list(results = data.frame(benchmark_id = 1,
    label = "omit", estimator = "aipw", full_estimate = 1, reduced_estimate = .5),
    full_scores = list(aipw = v), target_scores = list(`1_aipw` = -v)),
    class = "jointbounds_l2_benchmarks")
  list(surface = surface, bench = bench)
}

test_that("frontier confidence inversion uses the actual grid by default", {
  f <- make_paper_band_fixture()
  b <- l2_calibration_band(f$surface, f$bench, B = 29)
  expect_equal(b$frontier$estimate, 1)
  expect_equal(b$frontier$simultaneous_low, 1)
  expect_equal(b$frontier$simultaneous_high, 1)
  expect_false(b$interpolate)
  expect_equal(b$compatibility$estimate, c(FALSE, TRUE))
  expect_true(all(!b$compatibility$simultaneous_inner | b$compatibility$estimate))
  expect_true(all(!b$compatibility$estimate | b$compatibility$simultaneous_outer))
  expect_warning(smooth <- l2_calibration_band(f$surface, f$bench, B = 29, interpolate = TRUE), "coverage between evaluated grid points")
  expect_equal(smooth$frontier$estimate, .25)
  expect_match(smooth$inference_note, "descriptive")
  f$surface$fits[[2]]$ate$estimate <- c(.8, 1.2)
  missing <- l2_calibration_band(f$surface, f$bench, B = 29)
  expect_true(is.na(missing$frontier$estimate))
})

test_that("fixed-target tipping bands need no estimated benchmark", {
  f <- make_paper_band_fixture()
  b <- l2_tipping_band(f$surface, B = 29)
  expect_equal(b$frontier$estimate, 1)
  expect_warning(smooth <- l2_tipping_band(f$surface, B = 29, interpolate = TRUE), "coverage between evaluated grid points")
  expect_equal(smooth$frontier$estimate, .5)
  expect_equal(b$target_type, "fixed")
  expect_equal(unique(b$gaps$target), 0)
  expected <- f$surface$fits[[1]]$scores$eif_lower
  expect_equal(b$scores[, 1], expected - mean(expected))
  expect_equal(unique(b$gaps$endpoint), c("lower", "upper"))
  expect_error(l2_tipping_band(f$surface, calibration = "signed"), "fixed target")
})

test_that("tipping requires both endpoints for positive and negative centers", {
  d <- data.frame(delta_M = 0, delta_K = c(0, 1), method = "cs",
    lower = c(-2, -3), upper = c(-1, 1))
  f <- l2_tipping_frontier(d)
  expect_equal(f$delta_K, .5)
  expect_equal(f$upper, 0)
  expect_equal(l2_tipping_frontier(d, endpoint = "lower")$delta_K, 0)
})

test_that("joint endpoint covariance retains singular cases and correct scaling", {
  d <- simulate_l2_data(120, seed = 812)
  for (fun in list(l2_cs_crossfit, l2_sharp_crossfit)) {
    z <- fun(d, "Y", "A", "C", c("X1", "X2"), delta_M = 0,
      delta_K = 0, folds = 3, nuisance_method = "glm")
    est <- if (inherits(z, "jointbounds_l2_cs_cf")) "eif" else "dr"
    sc <- cbind(z$scores[[paste0(est, "_lower")]], z$scores[[paste0(est, "_upper")]])
    expect_equal(unname(z$influence_covariance[[est]]), unname(cov(sc)))
    expect_equal(unname(z$estimate_covariance[[est]]), unname(cov(sc) / nrow(d)))
    expect_equal(det(z$influence_covariance[[est]]), 0, tolerance = 1e-12)
    expect_equal(unname(diag(z$estimate_covariance[[est]])),
      z$ate$se[z$ate$estimator == est]^2)
  }
})

test_that("direct convex frontier matches an analytic finite-sample answer", {
  skip_if_not_installed("CVXR")
  d <- data.frame(Y = c(0, 0, -1, 1), A = c(0, 0, 1, 1), C = 0)
  nuis <- list(e0 = rep(.5, 4), e1 = rep(.5, 4), rho0 = rep(1, 4),
    rho1 = rep(1, 4), mu0 = rep(0, 4), mu1 = rep(0, 4))
  for (target in c(-.2, .2)) {
    z <- l2_frontier_direct(d, "Y", "A", "C", character(), delta_M = 0,
      target = target, nuisance = nuis, solver = "ECOS")
    expect_equal(z$radius, .4, tolerance = 1e-4)
    expect_true(z$feasible)
    expect_lt(z$target_error, 1e-5)
    expect_equal(z$arms[[2]]$missingness_cost, 0, tolerance = 1e-8)
  }
  no <- l2_frontier_direct(d, "Y", "A", "C", character(), delta_M = 0,
    target = 2, nuisance = nuis, solver = "ECOS")
  expect_equal(no$status, "infeasible")
  expect_equal(no$radius, Inf)
  expect_false(no$feasible)
})
test_that("CS reports nonregular zero moments instead of silently certifying them", {
  d <- simulate_l2_data(120, seed = 912)
  d$Y[d$C == 0] <- 0
  z <- l2_cs_crossfit(d, "Y", "A", "C", c("X1", "X2"),
    delta_M = 0, delta_K = .2, folds = 3, nuisance_method = "glm")
  flags <- z$moment_diagnostics
  expect_true(all(!flags$positive_required_moments[flags$estimator == "eif"]))
})