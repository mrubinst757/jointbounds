test_that("continuous simulation has finite empirical moments", {
  dat <- simulate_l2_data(500, seed = 11)
  expect_true(is.numeric(dat$Y))
  expect_true(all(is.na(dat$Y) == (dat$C == 1)))
  expect_true(all(dat$C[dat$U_I == 1] == 1))
  expect_true(is.finite(mean((dat$Yfull - mean(dat$Yfull))^4)))
  expect_true(all(c("p_A1_given_X", "p_UI0_given_YX", "p_UI1_given_YX",
                    "p_C0_given_YX", "p_C1_given_YX", "p_C0_given_X",
                    "p_C1_given_X", "C0", "C1", "UI0", "UI1") %in% names(dat)))
})

test_that("known DGP simulation returns oracle Monte Carlo CS bounds", {
  dat <- simulate_l2_data(500, seed = 17)
  fit <- l2_bounds(dat, "Y", "A", "C", c("X1", "X2"),
                   delta = .5, delta_R = .2, delta_A = .1,
                   method = "cs", folds = 2, nuisance_method = "glm",
                   return_true_bounds = TRUE)
  expect_true(all(c("arm", "ate") %in% names(fit$true_bounds)))
  expect_true(fit$true_bounds$ate$lower <= fit$true_bounds$ate$upper)
})

test_that("CS bounds collapse at zero radii and widen monotonically", {
  dat <- simulate_l2_data(600, seed = 12)
  fit0 <- l2_bounds(dat, "Y", "A", "C", c("X1", "X2"),
                    delta = 0.5, delta_R = 0, delta_A = 0,
                    method = "cs", folds = 3, seed = 3)
  fit1 <- l2_bounds(dat, "Y", "A", "C", c("X1", "X2"),
                    delta = 0.5, delta_R = 0.25, delta_A = 0.2,
                    method = "cs", folds = 3, seed = 3)
  expect_equal(fit0$ate$lower, fit0$ate$upper, tolerance = 1e-10)
  expect_lte(fit1$ate$lower, fit0$ate$lower)
  expect_gte(fit1$ate$upper, fit0$ate$upper)
  expect_true(fit1$moment_check$finite_fourth_empirical_moment)
})

test_that("net-M and separated CS models return valid arm and ATE intervals", {
  dat <- simulate_l2_data(500, seed = 13)
  sep <- l2_bounds(dat, "Y", "A", "C", c("X1", "X2"),
                   delta = c(0.4, 0.6), delta_R = c(0.2, 0.25),
                   delta_A = 0.15, model = "separated",
                   method = "cs", folds = 3)
  net <- l2_bounds(dat, "Y", "A", "C", c("X1", "X2"),
                   delta_M = c(0.1, 0.15), delta_A = 0.15,
                   model = "net", method = "cs", folds = 3)
  expect_true(all(sep$arm$lower <= sep$arm$upper))
  expect_true(all(net$arm$lower <= net$arm$upper))
  expect_true(sep$ate$lower <= sep$ate$upper)
  expect_true(net$ate$lower <= net$ate$upper)
})

test_that("sharp sieve dual returns diagnostics", {
  dat <- simulate_l2_data(250, seed = 14)
  fit <- l2_bounds(dat, "Y", "A", "C", c("X1", "X2"),
                   delta = 0.5, delta_R = 0.2, delta_A = 0.15,
                   method = "sharp", folds = 2, nuisance_method = "glm",
                   control = list(maxit = 300),
                   diagnostic_control = list(warn = FALSE))
  expect_equal(fit$ate$method, "sharp_sieve")
  expect_true(is.finite(fit$ate$lower))
  expect_true(is.finite(fit$ate$upper))
  expect_true(all(vapply(fit$diagnostics, function(x) {
    is.list(x$sharp) && is.finite(x$sharp$upper$value)
  }, logical(1))))
  expect_true(all(vapply(fit$diagnostics, function(x) {
    all(c("primal_value", "primal_dual_gap",
          "maximum_normalization_error", "maximum_budget_violation") %in%
        names(x$sharp$upper))
  }, logical(1))))
  expect_equal(nrow(fit$diagnostic_table), 4L)
  expect_equal(fit$diagnostic_summary$n_endpoint_fits, 4L)
  expect_true(all(c("lambda_missingness", "lambda_K",
    "missingness_budget_ratio", "confounding_budget_ratio",
    "boundary_M_fraction", "boundary_K_fraction",
    "relative_minimum_hessian_eigenvalue") %in%
    names(fit$diagnostic_table)))
})

test_that("sharp dual enforces a zero missingness budget exactly", {
  dat <- simulate_l2_data(180, seed = 141)
  nuis <- l2_nuisance_glm(dat, "Y", "A", "C", c("X1", "X2"), 2, 1)
  B <- l2_basis(dat, c("X1", "X2"), NULL)
  r <- l2_reference_quantities(dat$Y, dat$A, dat$C, 1, nuis)
  fit <- l2_dual_arm(r, B, .5 * r$pi, 0, .1, 0,
                     "separated", 1, list(maxit = 100))
  expect_true(fit$fixed_missingness_tilt)
  expect_equal(fit$divergence_R, 0)
  ev <- l2_dual_evaluate(dat$Y, dat$A, dat$C, 1, nuis, B, .5,
    fit$coefficients, 1)
  expect_equal(ev$M[ev$obs], rep(1, sum(ev$obs)))
})

test_that("sensitivity grids and tipping frontiers work", {
  dat <- simulate_l2_data(350, seed = 15)
  grid <- expand.grid(delta_R = c(0, 0.15), delta_A = c(0, 0.1))
  out <- l2_sensitivity_grid(dat, "Y", "A", "C", c("X1", "X2"),
                             grid = grid, delta = 0.5,
                             method = "cs", folds = 2)
  expect_equal(nrow(out$results), nrow(grid))
  front <- l2_tipping_frontier(out)
  expect_true(all(c("delta_R", "delta_A", "lower") %in% names(front)))
  expect_s3_class(front, "marbounds_l2_frontier")
})

test_that("sharp sensitivity grids report diagnostics for every radius pair", {
  dat <- simulate_l2_data(220, seed = 151)
  grid <- data.frame(delta_M = c(.1, .2), delta_K = c(.1, .15))
  out <- l2_sensitivity_grid(dat, "Y", "A", "C", c("X1", "X2"),
    grid = grid, model = "net", method = "sharp", folds = 2,
    nuisance_method = "glm",
    diagnostic_control = list(warn = FALSE, compute_hessian = FALSE))
  expect_equal(nrow(out$diagnostic_summary), nrow(grid))
  expect_equal(nrow(out$diagnostic_table), 4L * nrow(grid))
  expect_equal(sort(unique(out$diagnostic_table$grid_id)), seq_len(nrow(grid)))
  expect_true(all(c("diagnostic_ok", "n_diagnostic_flags") %in%
                  names(out$results)))
})

test_that("diagnostic warnings are consolidated and identify grid pairs", {
  summary <- data.frame(diagnostic_available = c(TRUE, TRUE),
    diagnostic_ok = c(FALSE, TRUE), flagged_codes = c("primal_dual_gap", ""),
    grid_id = 1:2)
  expect_warning(l2_warn_diagnostics(summary, "Sensitivity-grid"),
    "Affected grid_id values: 1")
})

test_that("tipping frontier requires an actual threshold crossing", {
  d <- data.frame(delta_R = 0, delta_A = c(0, .5), method = "cs",
                  lower = c(.4, .2), upper = c(.6, .8))
  no_tip <- l2_tipping_frontier(d, threshold = 0)
  expect_false(no_tip$tipped_on_grid)
  expect_true(is.na(no_tip$delta_A))
  d$lower <- c(.4, -.2)
  tip <- l2_tipping_frontier(d, threshold = 0)
  expect_true(tip$tipped_on_grid)
  expect_true(tip$interpolated)
  expect_equal(tip$delta_A, 1 / 3, tolerance = 1e-10)
})

test_that("benchmark tilts calibrate to the advertised L2 radii", {
  R <- c(0.8, 1, 1.2)
  K <- c(0.9, 1, 1.1)
  M <- c(0.95, 1, 1.05)
  out <- l2_calibrate(R, K, M)
  expect_equal(unname(out["delta_R"]), sqrt(mean((R - 1)^2)))
  expect_equal(unname(out["delta_M"]), sqrt(mean((M - 1)^2)))
  expect_equal(unname(out["delta_A"]), sqrt(mean(M * (K - 1)^2)))
  expect_equal(unname(out["delta_K"]), unname(out["delta_A"]))
  expect_equal(attr(out, "benchmark"), "user_supplied")
})

test_that("net sensitivity aliases and frontier defaults use delta_M and delta_K", {
  d <- data.frame(delta_M = c(0, 0, .2, .2),
                  delta_K = c(0, .2, 0, .2), method = "cs",
                  lower = c(.3, .1, .1, -.1),
                  upper = c(.5, .7, .7, .9))
  tip <- l2_tipping_frontier(d, threshold = 0)
  expect_true(all(c("delta_M", "delta_K") %in% names(tip)))
  expect_equal(attr(tip, "radius_x"), "delta_M")
  expect_equal(attr(tip, "radius_y"), "delta_K")
  cal <- l2_calibration_frontier(transform(d, estimate_reference = .4),
    benchmark = .3, type = "discrepancy")
  expect_true(all(c("delta_M", "minimum_delta_K") %in% names(cal)))
})

test_that("observed-data calibration enforces identification claims", {
  M <- K <- c(0.9, 1, 1.1)
  R <- c(0.8, 1, 1.2)
  net <- l2_calibrate(M = M, benchmark = "observed_data")
  expect_equal(attr(net, "interpretation"), "reduced-form combined distortion")
  expect_error(l2_calibrate(K = K, benchmark = "observed_data"), "requires")
  expect_error(l2_calibrate(R = R, K = K, benchmark = "observed_data"),
               "not jointly identified")
  cond <- l2_calibrate(K = K, benchmark = "observed_data",
                       conditional_on = "MAR")
  expect_equal(attr(cond, "conditional_on"), "MAR")
})

test_that("calibration frontier inverts joint sensitivity bounds", {
  grid <- expand.grid(delta_R = c(0, 0.2), delta_A = c(0, 0.1, 0.2))
  grid$method <- "sharp_sieve"
  grid$estimate_reference <- 1
  radius <- grid$delta_R + grid$delta_A
  grid$lower <- 1 - radius
  grid$upper <- 1 + radius
  out <- l2_calibration_frontier(grid, benchmark = 0.15,
                                 type = "discrepancy",
                                 radius_x = "delta_R")
  expect_equal(out$minimum_delta_A[out$delta_R == 0], 0.2)
  expect_equal(out$minimum_delta_A[out$delta_R == 0.2], 0)
  expect_true(all(out$compatible_on_grid))
})

test_that("bootstrap returns endpoint replicates", {
  dat <- simulate_l2_data(220, seed = 16)
  out <- l2_bootstrap(dat, "Y", "A", "C", c("X1", "X2"),
                      delta = 0.5, delta_R = 0.1, delta_A = 0.1,
                      method = "cs", folds = 2, B = 3, seed = 4)
  expect_equal(dim(out$replicates), c(3, 2))
  expect_true(out$successful >= 2)
  expect_equal(dim(out$ci), c(2, 2))
})

test_that("cross-fitted sharp estimator returns plug-in and DR endpoints", {
  dat <- simulate_l2_data(300, seed = 31)
  fit <- l2_sharp_crossfit(dat, "Y", "A", "C", c("X1", "X2"),
                           delta = .5, delta_R = .15, delta_A = .1,
                           folds = 2, nuisance_method = "glm",
                           control = list(maxit = 150),
                           diagnostic_control = list(warn = FALSE))
  expect_true(all(c("plugin", "dr", "true_eif") %in% fit$ate$estimator))
  expect_true(all(c("lower", "upper") %in% fit$ate$endpoint))
  expect_true(all(is.finite(fit$ate$estimate)))
  expect_true(all(is.finite(fit$ate$se)))
  expect_equal(length(fit$scores), 6L)
  expect_equal(nrow(fit$diagnostic_table), 12L)
  expect_equal(fit$diagnostic_summary$n_endpoint_fits, 12L)
  expect_equal(nrow(fit$crossfit_roles), 3L)
  expect_true(all(fit$crossfit_roles$evaluation_fold !=
                  fit$crossfit_roles$regression_fold))
  expect_true(all(nzchar(fit$crossfit_roles$fit_folds)))
})

test_that("single-mechanism affine extremizers recover the CS endpoints", {
  y <- c(-1, -.5, .5, 1)
  r <- list(obs = rep(TRUE, 4), y = y, e = rep(.5, 4),
    q = rep(.5, 4), rho = rep(.8, 4), w = rep(1, 4))
  basis <- matrix(1, 4, 1)
  cfg <- l2_diagnostic_control(list(warn = FALSE), compute_hessian = FALSE)

  dm <- .1; vm <- mean(y^2)
  lo_m <- l2_dual_arm(r, basis, rep(0, 4), 0, 0, dm,
    "net", -1, list(maxit = 500), cfg)
  hi_m <- l2_dual_arm(r, basis, rep(0, 4), 0, 0, dm,
    "net", 1, list(maxit = 500), cfg)
  expect_equal(c(-lo_m$value, hi_m$value),
    c(-dm * sqrt(vm), dm * sqrt(vm)), tolerance = 1e-5)

  dk <- .1; vk <- mean(.5^2 * y^2)
  lo_k <- l2_dual_arm(r, basis, rep(0, 4), 0, dk, 0,
    "net", -1, list(maxit = 500), cfg)
  hi_k <- l2_dual_arm(r, basis, rep(0, 4), 0, dk, 0,
    "net", 1, list(maxit = 500), cfg)
  expect_equal(c(-lo_k$value, hi_k$value),
    c(-dk * sqrt(vk), dk * sqrt(vk)), tolerance = 1e-5)
})

test_that("dual derivative evaluator includes an active boundary multiplier", {
  dat <- simulate_l2_data(120, seed = 32)
  nuis <- l2_nuisance_glm(dat, "Y", "A", "C", c("X1", "X2"), 2, 1)
  B <- l2_basis(dat, c("X1", "X2"), NULL)
  r <- l2_reference_quantities(dat$Y, dat$A, dat$C, 1, nuis)
  dual <- l2_dual_arm(r, B, .5 * r$pi, .15, .1, .15,
                      "separated", 1, list(maxit = 100))
  ev <- l2_dual_evaluate(dat$Y, dat$A, dat$C, 1, nuis, B, .5,
                         dual$coefficients, 1, boundary_tol = 1e-6)
  expect_true(all(ev$eta >= 0))
  expect_true(all(is.finite(ev$d_rho[ev$obs])))
})

test_that("analytic dual-loss derivatives agree with finite differences", {
  y <- .2; A <- 1; C <- 0; B <- matrix(1, 1, 1)
  nuisance <- list(e1 = .6, rho1 = .7)
  par <- c(log(2), log(2), 0, 0); h <- 1e-6
  base <- l2_dual_evaluate(y, A, C, 1, nuisance, B, .5, par, 1)
  ep <- nuisance; em <- nuisance
  ep$e1 <- ep$e1 + h; em$e1 <- em$e1 - h
  num_e <- (l2_dual_evaluate(y, A, C, 1, ep, B, .5, par, 1)$H -
    l2_dual_evaluate(y, A, C, 1, em, B, .5, par, 1)$H) / (2 * h)
  rp <- nuisance; rm <- nuisance
  rp$rho1 <- rp$rho1 + h; rm$rho1 <- rm$rho1 - h
  num_rho <- (l2_dual_evaluate(y, A, C, 1, rp, B, .5, par, 1)$H -
    l2_dual_evaluate(y, A, C, 1, rm, B, .5, par, 1)$H) / (2 * h)
  expect_false(base$active)
  expect_equal(base$d_e, num_e, tolerance = 1e-5)
  expect_equal(base$d_rho, num_rho, tolerance = 1e-5)
})

test_that("oracle nuisance error mode is finite and records its settings", {
  dat <- simulate_l2_data(260, seed = 33)
  fit <- l2_sharp_crossfit(dat, "Y", "A", "C", c("X1", "X2"),
    delta = .5, delta_R = .15, delta_A = .1, folds = 2,
    nuisance_source = "oracle_error", nuisance_error_rate = .25,
    nuisance_error_scale = c(e = 1, rho = .8, regression = .5),
    nuisance_error_sign = c(e = -1, rho = 1, regression = -1),
    nuisance_method = "glm", control = list(maxit = 100))
  expect_true(all(is.finite(fit$ate$estimate)))
  expect_equal(fit$parameters$nuisance_source, "oracle_error")
  expect_equal(fit$parameters$nuisance_error_rate, .25)
  expect_true("true_eif" %in% fit$ate$estimator)
  expect_equal(fit$parameters$folds, 1L)
  expect_match(fit$parameters$sample_splitting, "none")
  expect_true(all(vapply(fit$diagnostics, function(endpoint_folds)
    endpoint_folds[[1L]]$dual_nuisance == "controlled-error estimate",
    logical(1))))
})

test_that("controlled errors perturb outcome and second-stage regressions", {
  dat <- simulate_l2_data(180, seed = 331)
  use <- rep(TRUE, nrow(dat))
  exact <- l2_split_nuisance_oracle_error(dat, c("X1", "X2"), use, use,
    rate = 0, scale = c(e = 0, rho = 0, mu = 0, regression = 0),
    sign = c(e = 1, rho = 1, mu = 1, regression = 1))
  perturbed <- l2_split_nuisance_oracle_error(dat, c("X1", "X2"), use, use,
    rate = .25, scale = c(e = 0, rho = 0, mu = 1, regression = 1),
    sign = c(e = 1, rho = 1, mu = -1, regression = 1))
  expect_equal(exact$test$mu0, dat$mu_Y0_given_XC0)
  expect_false(isTRUE(all.equal(perturbed$test$mu0,
    dat$mu_Y0_given_XC0)))
  expect_gt(stats::sd(perturbed$test$err_H), 0)
  expect_gt(stats::sd(perturbed$test$err_nu4), 0)
})

test_that("oracle sensitivity calibration reports the minimal DGP parameters", {
  dat <- simulate_l2_data(500, seed = 34)
  op <- l2_oracle_sensitivity_parameters(dat, nodes = 40)
  expect_equal(op$delta, c(arm0 = 1, arm1 = 1))
  expect_true(all(op$delta_A > 0))
  expect_equal(op$delta_K, op$delta_A)
  expect_true(all(is.finite(op$delta_M) & op$delta_M > 0))
  expect_true(all(is.finite(op$delta_R) & op$delta_R > 0))
  expect_true(all(op$empirical_delta <= op$delta))
  dat0 <- simulate_l2_data(500, seed = 34, confounding_strength = 0)
  op0 <- l2_oracle_sensitivity_parameters(dat0, nodes = 40)
  expect_equal(op0$delta_A, c(arm0 = 0, arm1 = 0), tolerance = 1e-10)
})

test_that("paper simulation design defaults to net M/K and controlled errors", {
  d <- l2_simulation_design("main")
  expect_true(all(d$model == "net"))
  expect_true("joint_primary" %in% d$scenario_id)
  expect_false(any(d$scenario_id %in%
    c("joint_same", "joint_opposed", "nonlinear")))
  expect_true(all(d$nuisance_simulation == "oracle_error"))
})

test_that("simulation regimes switch MAR and NUC violations independently", {
  getp <- function(v) l2_oracle_sensitivity_parameters(
    simulate_l2_data(400, seed = 35, violation = v), nodes = 40)
  mar <- getp("mar"); nuc <- getp("nuc"); both <- getp("both"); none <- getp("none")
  expect_true(all(mar$delta_R > 0) && all(mar$delta_A == 0))
  expect_true(all(nuc$delta_R == 0) && all(nuc$delta_A > 0))
  expect_true(all(both$delta_R > 0) && all(both$delta_A > 0))
  expect_true(all(none$delta_R == 0) && all(none$delta_A == 0))
})

test_that("bounded and linear links give finite and infinite Linf radii", {
  db <- simulate_l2_data(400, seed = 36, violation = "both",
                         density_ratio_link = "bounded")
  dl <- simulate_l2_data(400, seed = 36, violation = "both",
                         density_ratio_link = "linear")
  pb <- l2_oracle_sensitivity_parameters(db, nodes = 40)
  pl <- l2_oracle_sensitivity_parameters(dl, nodes = 40)
  expect_true(all(is.finite(c(pb$delta_R_inf, pb$delta_A_inf))))
  expect_true(all(is.infinite(c(pl$delta_R_inf, pl$delta_A_inf))))
  bi <- l2_oracle_linf_outer_bounds(db, pb, nodes = 40)
  expect_true(bi$finite && all(is.finite(bi$ate)))
})

test_that("sharp Linf LP respects the outer interval", {
  skip_if_not_installed("lpSolve")
  d <- simulate_l2_data(250, seed = 37, violation = "both",
                        density_ratio_link = "bounded")
  p <- l2_oracle_sensitivity_parameters(d, nodes = 40)
  outer <- l2_oracle_linf_outer_bounds(d, p, nodes = 40)
  sharp <- l2_oracle_linf_sharp_bounds(d, p, max_n = 250)
  expect_true(sharp$finite)
  expect_gte(sharp$ate[["lower"]], outer$ate[["lower"]] - 1e-7)
  expect_lte(sharp$ate[["upper"]], outer$ate[["upper"]] + 1e-7)
})

test_that("CS crossfit returns plug-in, EIF, and oracle-EIF endpoints", {
  d <- simulate_l2_data(280, seed = 38, violation = "both")
  fit <- l2_cs_crossfit(d, "Y", "A", "C", c("X1", "X2"),
    delta = .8, delta_R = .3, delta_A = .2, folds = 2,
    nuisance_source = "oracle_error", nuisance_method = "glm")
  expect_true(all(c("plugin", "eif", "true_eif") %in% fit$ate$estimator))
  expect_true(all(is.finite(fit$ate$estimate)))
  expect_true(all(is.finite(fit$ate$se)))
  expect_equal(fit$parameters$folds, 1L)
  expect_match(fit$parameters$sample_splitting, "none")
})

test_that("cross-fitted CS and sharp estimators support the net-M model", {
  d <- simulate_l2_data(220, seed = 381, violation = "both")
  cs <- l2_cs_crossfit(d, "Y", "A", "C", c("X1", "X2"),
    delta_M = .15, delta_A = .1, model = "net", folds = 2,
    nuisance_source = "oracle_error", nuisance_method = "glm")
  sharp <- l2_sharp_crossfit(d, "Y", "A", "C", c("X1", "X2"),
    delta_M = .15, delta_A = .1, model = "net", folds = 2,
    nuisance_source = "oracle_error", nuisance_method = "glm",
    control = list(maxit = 100))
  expect_true(all(is.finite(cs$ate$estimate)))
  expect_true(all(is.finite(sharp$ate$estimate)))
  expect_equal(cs$parameters$model, "net")
  expect_equal(sharp$parameters$model, "net")
})

test_that("continuous multiplier bands use all endpoint scores jointly", {
  make_fit <- function(shift) list(
    ate = data.frame(endpoint = c("lower", "upper"), estimator = "eif",
      estimate = c(-.2, .3) + shift),
    scores = list(eif_lower = c(-1, 0, 1, 0),
                  eif_upper = c(0, -1, 0, 1)))
  out <- l2_multiplier_band(list(first = make_fit(0), second = make_fit(.1)),
    B = 49, conf_level = .9, seed = 7)
  expect_equal(nrow(out$results), 4L)
  expect_equal(nrow(out$outward), 2L)
  expect_true(is.finite(out$critical_value) && out$critical_value > 0)
  expect_true(all(out$outward$outward_lower <= out$outward$lower))
  expect_true(all(out$outward$outward_upper >= out$outward$upper))
})

test_that("sensitivity-band wrapper retains grid parameters", {
  d <- simulate_l2_data(160, seed = 384, violation = "both")
  g <- data.frame(delta_R = c(.1, .2), delta_A = c(.1, .1))
  out <- l2_sensitivity_band(d, "Y", "A", "C", c("X1", "X2"), g,
    delta = .5, bound_method = "cs", folds = 2,
    nuisance_method = "glm", B = 9, seed = 2)
  expect_equal(nrow(out$outward), 2L)
  expect_equal(out$outward$delta_R, g$delta_R)
  expect_equal(out$outward$delta_A, g$delta_A)
})

test_that("L-infinity sharp lower envelope respects nonnegative R", {
  skip_if_not_installed("lpSolve")
  d <- simulate_l2_data(180, seed = 382, violation = "both",
                        density_ratio_link = "bounded")
  fit <- l2_linf_bounds(d, "Y", "A", "C", c("X1", "X2"),
    delta = .5, delta_R_inf = 2, delta_A_inf = .2,
    method = "sharp", folds = 2, nuisance_method = "glm", max_n = 100)
  expect_true(all(vapply(fit$diagnostics,
    function(z) isTRUE(z$mixture_envelope_respected), logical(1))))
})

test_that("capped L2 support and alternative binary bounds are available", {
  g <- c(0, 1, 2, 3)
  opt <- l2_capped_support(g, .4)
  expect_equal(opt$norm, .4, tolerance = 1e-7)
  expect_true(all(opt$q >= 0 & opt$q <= 1))
  d <- simulate_l2_data(240, seed = 383, violation = "mar")
  d$Yfull <- as.numeric(d$Yfull > 0)
  d$Y <- d$Yfull; d$Y[d$C == 1] <- NA_real_
  fit <- l2_prevalence_bounds(d, "Y", "A", "C", c("X1", "X2"),
    estimand = "composite", kappa = .25, folds = 2,
    nuisance_method = "glm")
  expect_true(all(is.finite(fit$bounds$estimate)))
  expect_true(fit$support$arm0$value >= 0)
  expect_true(fit$support$arm1$value >= 0)
})

test_that("simultaneous AIPW benchmarks cap subsets within p_z", {
  d <- simulate_l2_data(300, seed = 39, violation = "both")
  fit <- l2_aipw_benchmarks(d, "Y", "A", "C", c("X1", "X2"),
    p_z = 0:1, max_per_pz = 1, estimator = "both", variance = TRUE,
    folds = 2, nuisance_method = "glm", seed = 4)
  expect_equal(nrow(fit$results), 4L)
  expect_equal(as.integer(table(fit$results$p_z)), c(2L, 2L))
  expect_true(all(is.finite(fit$results$delta_bench)))
  expect_true(all(is.finite(fit$results$se)))
  expect_true(any(abs(fit$results$delta_bench) > 1e-8))
  expect_true(any(fit$results$se > 1e-8))
  expect_equal(sort(unique(fit$results$estimator)), c("aipw", "plugin"))
  band <- l2_benchmark_band(fit, B = 29, conf_level = .9, seed = 8)
  expect_equal(nrow(band$results), 2L)
  expect_true(all(band$results$simultaneous_low <=
                  band$results$delta_bench))
  expect_true(all(band$results$simultaneous_high >=
                  band$results$delta_bench))
  expect_equal(length(fit$target_scores), nrow(fit$results))
  expect_equal(length(fit$full_scores$aipw), nrow(d))
})

test_that("calibration frontier bands combine endpoint and benchmark scores", {
  n <- 8L
  grid <- expand.grid(delta_M = c(0, .1), delta_K = c(0, .1, .2))
  base_score <- rep(c(-1, 1), length.out = n)
  fits <- lapply(seq_len(nrow(grid)), function(i) {
    width <- grid$delta_M[i] + grid$delta_K[i]
    list(ate = data.frame(endpoint = c("lower", "upper"),
      estimator = "eif", estimate = c(.5 - width, .5 + width)),
      scores = list(eif_lower = base_score * (1 + width),
                    eif_upper = -base_score * (1 + width)))
  })
  surface <- structure(list(fits = fits, grid = grid, bound_method = "cs"),
    class = c("marbounds_l2_sensitivity_band", "list"))
  target_score <- rep(c(-.4, .4), length.out = n)
  benchmarks <- structure(list(
    results = data.frame(benchmark_id = 1L, label = "omit X1",
      estimator = "aipw", full_estimate = .5, reduced_estimate = .35),
    target_scores = list(`1_aipw` = target_score),
    scores = list(`1_aipw` = target_score - base_score),
    full_scores = list(aipw = base_score)),
    class = "marbounds_l2_benchmarks")
  out <- l2_calibration_band(surface, benchmarks, B = 29L,
    conf_level = .9, seed = 3L)
  expect_equal(nrow(out$frontier), 2L)
  expect_equal(out$frontier$endpoint, rep("lower", 2L))
  expect_true(all(is.finite(out$frontier$estimate)))
  expect_true(is.finite(out$critical_value) && out$critical_value > 0)
  expect_equal(nrow(out$scores), n)
  expect_equal(ncol(out$scores), nrow(grid))
  expect_true(all(c("pointwise_low", "pointwise_high",
    "simultaneous_low", "simultaneous_high") %in%
    names(out$frontier)))
})

test_that("simulation DGP exposes prespecified stress factors", {
  d <- simulate_l2_data(180, seed = 391, violation = "both",
    informative_strength = .8, informative_intercept = -3,
    treatment_scale = 1.5, covariate_complexity = "nonlinear",
    informative_concentration = "high_variance",
    outcome_distribution = "t", outcome_df = 12)
  expect_equal(unique(d$outcome_distribution), "t")
  expect_equal(unique(d$covariate_complexity), "nonlinear")
  expect_true(all(d$informative_outcome_coef[d$informative_outcome_coef != 0] == .8))
  expect_true(all(is.finite(d$p_A1_given_X)))
  expect_equal(unique(d$true_ate), .7)
  expect_error(simulate_l2_data(50, outcome_distribution = "t", outcome_df = 8),
    "exceed eight")
})

test_that("paper simulation registry is complete and fast profile is a subset", {
  main <- l2_simulation_design("main")
  fast <- l2_simulation_design("fast")
  expect_true(all(c("joint_n500", "joint_n2000", "slow_nuisance",
    "heavy_tail", "unbounded_ratio") %in% main$scenario_id))
  expect_true(all(fast$scenario_id %in% main$scenario_id))
  expect_equal(anyDuplicated(main$scenario_id), 0L)
})

test_that("multi-benchmark frontiers retain benchmark labels", {
  grid <- expand.grid(delta_R = c(0, .2), delta_A = c(0, .2),
                      method = "cs", stringsAsFactors = FALSE)
  grid$lower <- -grid$delta_R - grid$delta_A
  grid$upper <- grid$delta_R + grid$delta_A
  grid$estimate_reference <- 0
  b <- data.frame(benchmark_id = 1:2, label = c("Z1", "Z2"), p_z = c(1, 2),
                  estimator = "aipw", delta_bench = c(.1, -.1))
  f <- l2_benchmark_frontiers(grid, b)
  expect_s3_class(f, "marbounds_l2_benchmark_frontiers")
  expect_equal(sort(unique(f$label)), c("Z1", "Z2"))
  expect_true(any(f$compatible_on_grid))
  eq <- l2_equal_radius(f)
  expect_equal(nrow(eq), 2L)
  expect_equal(eq$equal_radius, c(.2, .2))
  expect_equal(unique(f$equal_radius), .2)
  if (requireNamespace("ggplot2", quietly = TRUE))
    expect_s3_class(plot(f), "ggplot")
})

test_that("benchmark inversion uses its paired full-X reference", {
  grid <- expand.grid(delta_R = c(0, .2), delta_A = c(0, .2),
                      method = "cs", stringsAsFactors = FALSE)
  grid$lower <- .5 - grid$delta_R - grid$delta_A
  grid$upper <- .5 + grid$delta_R + grid$delta_A
  grid$estimate_reference <- .5
  b <- data.frame(benchmark_id = 1, label = "Z", p_z = 1,
    estimator = "aipw", full_estimate = .6, reduced_estimate = .7,
    delta_bench = .1)
  f <- l2_benchmark_frontiers(grid, b)
  expect_equal(unique(f$reference), .6)
  expect_equal(unique(f$target), .7)
})

test_that("sensitivity grid forwards nuisance estimation options", {
  d <- simulate_l2_data(180, seed = 40)
  g <- data.frame(delta_R = c(0, .1), delta_A = c(0, .1))
  s <- l2_sensitivity_grid(d, "Y", "A", "C", c("X1", "X2"), g,
    delta = .5, method = "cs", folds = 2, nuisance_method = "glm")
  expect_equal(nrow(s$results), 2L)
  expect_true(all(is.finite(s$results$lower)))
})

test_that("overall equal radius compares benchmark and tipping strength", {
  grid <- expand.grid(delta_R = c(0, .2), delta_A = c(0, .2),
                      method = "cs", stringsAsFactors = FALSE)
  grid$lower <- .3 - grid$delta_R - grid$delta_A
  grid$upper <- .3 + grid$delta_R + grid$delta_A
  grid$estimate_reference <- .3
  b <- data.frame(benchmark_id = 1, label = "Z1", p_z = 1, p_x = 2,
                  n_omitted = 1, omitted = "X2",
                  estimator = "aipw", delta_bench = -.1,
                  se = .025,
                  conf_low = -.15, conf_high = -.05,
                  full_estimate = .3, full_se = .05,
                  full_conf_low = .2, full_conf_high = .4)
  f <- l2_benchmark_frontiers(grid, b)
  overall <- l2_overall_equal_radius(grid)
  expect_equal(overall$overall_equal_radius, .2)
  cmp <- l2_equal_radius_comparison(grid, f)
  expect_equal(cmp$strength_ratio, 1)
  cmp_u <- l2_equal_radius_comparison(grid, f, uncertainty = TRUE,
                                      conf_level = .90)
  expect_equal(attr(cmp_u, "conf_level"), .90)
  expect_true(all(c("benchmark_radius_low", "benchmark_radius_high",
    "overall_radius_low", "overall_radius_high") %in% names(cmp_u)))
  if (requireNamespace("ggplot2", quietly = TRUE))
    expect_s3_class(plot(cmp_u), "ggplot")
})
