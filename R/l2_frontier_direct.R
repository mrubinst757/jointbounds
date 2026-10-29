#' Direct convex optimization of a sharp net-model frontier
#'
#' Solves the paper's finite-basis empirical frontier program without a
#' sensitivity grid. For fixed missingness budgets, minimize the common squared
#' confounding radius subject to an ATE target equality. This is a plug-in
#' computation, not an influence-function estimator or confidence interval.
#'
#' @inheritParams l2_bounds
#' @param target A finite ATE target. Zero gives tipping. Supply a signed or
#'   toward-zero calibration target on the estimand scale for calibration.
#' @param nuisance Optional list of full-length vectors `e0`, `e1`, `rho0`,
#'   `rho1`, `mu0`, `mu1`, for matching an existing fitted empirical program.
#'   If NULL, nuisances are estimated with the usual cross-fitting helpers.
#' @param solver Optional installed CVXR conic solver name.
#' @param solver_control Named list of additional arguments to `CVXR::solve()`.
#' @param feasibility_tol Tolerance for the reported constraint residual check.
#' @return A list with `radius` (common delta_K), `squared_radius`, solver
#'   `status`, arm-specific tilts and budgets, target and balance residuals,
#'   and `feasible`. Infeasible programs return an infinite radius; solver
#'   failures return NA. An inaccurate or infeasible numerical candidate is
#'   flagged and must not be treated as a certified frontier point.
#' @details Uses the net-M model and the same inverse-probability weighted
#'   empirical measure and basis constraints as `l2_bounds(method = "sharp")`.
#'   The optional CVXR package is required. Conditional normalization is relaxed
#'   to finite-basis balance, so this computes an empirical sieve frontier.
#' @export
l2_frontier_direct <- function(data, Y, A, C, X, delta_M = 0, target = 0,
    basis = NULL, folds = 5L, nuisance_method = c("SuperLearner", "glm"),
    sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm", sl_lib_outcome = "SL.glm",
    seed = 1L, nuisance = NULL, solver = NULL, solver_control = list(),
    feasibility_tol = 1e-5) {
  if (!requireNamespace("CVXR", quietly = TRUE))
    stop("Direct frontier optimization requires the suggested package CVXR")
  nuisance_method <- match.arg(nuisance_method)
  delta_M <- l2_as_pair(delta_M, "delta_M")
  if (length(target) != 1L || !is.finite(target)) stop("target must be finite")
  if (length(feasibility_tol) != 1L || !is.finite(feasibility_tol) ||
      feasibility_tol <= 0) stop("feasibility_tol must be positive")
  if (!is.list(solver_control) || (length(solver_control) &&
      (is.null(names(solver_control)) || any(!nzchar(names(solver_control))))))
    stop("solver_control must be a named list")
  if (any(names(solver_control) %in% c("a", "object", "solver")))
    stop("Use solver rather than overriding object or solver in solver_control")
  inp <- l2_validate_inputs(data, Y, A, C, X, c(1, 1), c(0, 0),
                            c(0, 0), delta_M)
  if (is.null(nuisance)) nuisance <- if (nuisance_method == "SuperLearner")
    l2_nuisance_sl(data, Y, A, C, X, folds, seed, sl_lib_prop,
                  sl_lib_miss, sl_lib_outcome) else
    l2_nuisance_glm(data, Y, A, C, X, folds, seed)
  n <- nrow(data)
  for (nm in c("e0", "e1", "rho0", "rho1", "mu0", "mu1")) {
    if (length(nuisance[[nm]]) != n || any(!is.finite(nuisance[[nm]])))
      stop("nuisance$", nm, " must be a finite vector of length nrow(data)")
    if (!startsWith(nm, "mu") && any(nuisance[[nm]] <= 0 | nuisance[[nm]] > 1))
      stop("Fitted exposure and response probabilities must be in (0, 1]")
  }
  BB <- l2_basis(data, X, basis)
  if (any(!is.finite(BB))) stop("basis must be finite")
  t <- CVXR::Variable(1L, nonneg = TRUE)
  constraints <- list(); arms <- vector("list", 2L); effects <- list()
  for (a in 0:1) {
    r <- l2_reference_quantities(inp$y, inp$A, inp$C, a, nuisance)
    use <- which(r$obs)
    if (!length(use)) stop("Each arm must contain observed outcomes")
    w <- r$w[use] / n; Bmat <- BB[use, , drop = FALSE]
    m <- CVXR::Variable(length(use)); l <- CVXR::Variable(length(use))
    constraints <- c(constraints, list(m >= r$rho[use], l >= 0,
      t(Bmat) %*% CVXR::multiply(w, m - 1) == 0,
      t(Bmat) %*% CVXR::multiply(w, l - 1) == 0))
    if (delta_M[a + 1L] == 0) constraints <- c(constraints, list(m == 1))
    else constraints <- c(constraints, list(
      CVXR::sum_entries(CVXR::multiply(w, CVXR::square(m - 1))) <= delta_M[a + 1L]^2))
    perspective <- lapply(seq_along(use), function(i)
      w[i] * CVXR::quad_over_lin(l[i] - m[i], m[i]))
    constraints <- c(constraints, list(Reduce(`+`, perspective) <= t))
    effects[[a + 1L]] <- CVXR::sum_entries(CVXR::multiply(w * r$y[use],
      CVXR::multiply(r$e[use], m) + CVXR::multiply(r$q[use], l)))
    arms[[a + 1L]] <- list(m = m, l = l, w = w, B = Bmat, r = r, use = use)
  }
  constraints <- c(constraints, list(effects[[2L]] - effects[[1L]] == target))
  problem <- CVXR::Problem(CVXR::Minimize(t), constraints)
  if (!CVXR::is_dcp(problem)) stop("Internal frontier program is not convex")
  args <- c(list(a = problem), if (!is.null(solver)) list(solver = solver),
            solver_control)
  fit <- do.call(CVXR::solve, args)
  status <- as.character(fit$status)
  success <- status %in% c("optimal", "optimal_inaccurate")
  out <- list(radius = if (status == "infeasible") Inf else NA_real_,
    squared_radius = if (status == "infeasible") Inf else NA_real_,
    status = status, feasible = FALSE, delta_M = delta_M, target = target,
    model = "net", method = "direct empirical sieve frontier", call = match.call(),
    inference_note = "Plug-in optimization only; no frontier confidence interval is provided.")
  if (!success) return(out)
  tv <- as.numeric(fit$getValue(t))
  tilts <- lapply(seq_along(arms), function(j) {
    z <- arms[[j]]; m <- as.numeric(fit$getValue(z$m)); l <- as.numeric(fit$getValue(z$l))
    r <- z$r; use <- z$use
    list(arm = j - 1L, indices = use, M = m, L = l,
      missingness_cost = sum(z$w * (m - 1)^2),
      confounding_cost = if (any(m <= 0)) Inf else sum(z$w * (l - m)^2 / m),
      balance_error = max(abs(c(crossprod(z$B, z$w * (m - 1)),
                                  crossprod(z$B, z$w * (l - 1))))),
      pointwise_violation = max(0, r$rho[use] - m, -l),
      value = sum(z$w * r$y[use] * (r$e[use] * m + r$q[use] * l)))
  })
  target_error <- abs(tilts[[2L]]$value - tilts[[1L]]$value - target)
  violation <- max(0, -tv, target_error, vapply(seq_along(tilts), function(j) {
    z <- tilts[[j]]
    max(z$balance_error, z$pointwise_violation,
        z$missingness_cost - delta_M[j]^2, z$confounding_cost - tv)
  }, numeric(1)))
  out$radius <- sqrt(max(0, tv)); out$squared_radius <- max(0, tv)
  out$arms <- tilts; out$target_error <- target_error
  out$maximum_constraint_violation <- violation
  out$feasible <- is.finite(violation) && violation <= feasibility_tol
  if (!out$feasible || status != "optimal")
    warning("Direct frontier candidate requires review: status = ", status,
            ", maximum constraint violation = ", signif(violation, 4))
  out
}
