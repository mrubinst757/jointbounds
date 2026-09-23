#' Coverage simulation for pointwise and uniform CIs (mar_bounds param_grid inference)
#'
#' Uses a known DGP to compute true bounds (via large-MC theta), then checks
#' empirical coverage of pointwise and uniform CIs over replicates.
#' Run from package root: source("inst/scripts/coverage_simulation.R")
#' or: Rscript inst/scripts/coverage_simulation.R (after installing the package).

suppressPackageStartupMessages({
  if (requireNamespace("SuperLearner", quietly = TRUE)) library(SuperLearner)
})
devtools::load_all(".", quiet = TRUE)

# --- DGP: simple parametric (C = 1 means missing) ---
dgp <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- runif(n, -1.5, 1.5)
  e_x <- plogis(0.4 * X)
  A <- rbinom(n, 1, e_x)
  pi0_x <- 0.1 + 0.15 * plogis(X)
  pi1_x <- 0.1 + 0.2 * plogis(X)
  C <- ifelse(A == 0, rbinom(n, 1, pi0_x), rbinom(n, 1, pi1_x))
  mu0_x <- plogis(0.2 * X)
  mu1_x <- plogis(0.3 + 0.2 * X)
  Y <- ifelse(C == 1, NA_real_,
              ifelse(A == 0, rbinom(n, 1, mu0_x), rbinom(n, 1, mu1_x)))
  list(
    data = data.frame(Y = Y, A = A, C = C, X = X),
    e = e_x,
    pi0 = pi0_x,
    pi1 = pi1_x,
    mu0 = mu0_x,
    mu1 = mu1_x
  )
}

# True theta = E[phi(O; eta)] approximated by large sample
compute_true_theta <- function(N = 200000, seed = 1) {
  sim <- dgp(N, seed = seed)
  phi <- jointbounds:::influence_functions(
    sim$data$Y, sim$data$A, sim$data$C,
    sim$e, sim$pi0, sim$pi1, sim$mu0, sim$mu1
  )
  colMeans(phi)
}

# True bound at grid point (e.g. monotonicity_pos lower: b = coef_mono_pos_lower_ate(delta_0u))
true_bounds_at_grid <- function(theta, assumption, bound_type, param_grid) {
  grid_df <- expand.grid(param_grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  n_grid <- nrow(grid_df)
  vals <- numeric(n_grid)
  for (j in seq_len(n_grid)) {
    row <- as.list(grid_df[j, , drop = FALSE])
    delta_0u <- row$delta_0u %||% 1
    delta_1u <- row$delta_1u %||% 1
    if (assumption == "monotonicity_pos") {
      b <- if (bound_type == "lower") jointbounds:::coef_mono_pos_lower_ate(delta_0u)
      else jointbounds:::coef_mono_pos_upper_ate(delta_1u)
    } else if (assumption == "bounded_delta") {
      b <- if (bound_type == "lower") jointbounds:::coef_delta_lower_ate(delta_0u, delta_1u)
      else jointbounds:::coef_delta_upper_ate(delta_0u, delta_1u)
    } else {
      stop("assumption not implemented for true_bounds_at_grid")
    }
    vals[j] <- sum(b * theta)
  }
  list(grid = grid_df, true_bounds = vals)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- Main: run coverage simulation ---
run_coverage_simulation <- function(n = 800,
                                    R = 200,
                                    alpha = 0.05,
                                    assumption = "monotonicity_pos",
                                    bound_type = "lower",
                                    param_grid = list(delta_0u = c(0.3, 0.5, 0.7),
                                                      delta_1u = c(0.3, 0.5, 0.7)),
                                    B_boot = 500,
                                    seed_theta = 1,
                                    seed_repl = 100) {
  message("Computing true theta (large MC)...")
  theta_true <- compute_true_theta(N = 200000, seed = seed_theta)
  true_spec <- true_bounds_at_grid(theta_true, assumption, bound_type, param_grid)
  true_bounds <- true_spec$true_bounds
  n_grid <- length(true_bounds)

  pointwise_covered <- matrix(NA, R, n_grid)
  uniform_covered <- logical(R)

  for (r in seq_len(R)) {
    if (r %% 25 == 0) message("Replicate ", r, "/", R)
    set.seed(seed_repl + r)
    sim <- dgp(n)
    fit <- tryCatch(
      jointbounds::mar_bounds(
        sim$data, Y = "Y", A = "A", C = "C", X = "X",
        estimand = "ate", assumption = assumption,
        delta_0u = 1, delta_1u = 1,
        param_grid = param_grid,
        grid_bounds = if (bound_type == "lower") "lower" else "upper",
        alpha = alpha, B = B_boot,
        V = 3, sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm", sl_lib_outcome = "SL.glm"
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      pointwise_covered[r, ] <- NA
      uniform_covered[r] <- NA
      next
    }
    # Extract grid from appropriate bound
    grid_obj <- if (bound_type == "lower") fit$lower_grid else fit$upper_grid
    if (is.null(grid_obj)) {
      pointwise_covered[r, ] <- NA
      uniform_covered[r] <- NA
      next
    }
    g <- grid_obj$grid
    for (j in seq_len(n_grid)) {
      pointwise_covered[r, j] <- (true_bounds[j] >= g$ci_lower_pointwise[j] &
                                   true_bounds[j] <= g$ci_upper_pointwise[j])
    }
    uniform_covered[r] <- all(true_bounds >= g$ci_lower_uniform &
                               true_bounds <= g$ci_upper_uniform)
  }

  # Drop failed replicates
  ok <- !is.na(uniform_covered)
  R_ok <- sum(ok)
  pointwise_cov_j <- colMeans(pointwise_covered[ok, , drop = FALSE], na.rm = TRUE)
  pointwise_cov_avg <- mean(pointwise_cov_j)
  uniform_cov <- mean(uniform_covered[ok], na.rm = TRUE)

  list(
    n = n,
    R = R,
    R_ok = R_ok,
    alpha = alpha,
    assumption = assumption,
    bound_type = bound_type,
    n_grid = n_grid,
    pointwise_coverage_by_point = pointwise_cov_j,
    pointwise_coverage_avg = pointwise_cov_avg,
    uniform_coverage = uniform_cov,
    true_bounds = true_bounds,
    param_grid = param_grid
  )
}

# Run with small grid for a quick check (few minutes)
message("Running coverage simulation (n=800, R=200, 9 grid points)...")
results <- run_coverage_simulation(
  n = 800,
  R = 200,
  alpha = 0.05,
  assumption = "monotonicity_pos",
  bound_type = "lower",
  param_grid = list(delta_0u = c(0.3, 0.5, 0.7), delta_1u = c(0.3, 0.5, 0.7)),
  B_boot = 500
)

message("\n--- Results ---")
message("Replicates used: ", results$R_ok, " / ", results$R)
message("Pointwise coverage (average over grid points): ", round(results$pointwise_coverage_avg, 3))
message("Pointwise coverage by grid point: ", paste(round(results$pointwise_coverage_by_point, 3), collapse = ", "))
message("Uniform (simultaneous) coverage: ", round(results$uniform_coverage, 3))
message("Nominal level: ", 1 - results$alpha)
