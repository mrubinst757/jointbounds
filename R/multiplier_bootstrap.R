#' Multiplier bootstrap for simultaneous inference over a parameter grid
#'
#' For a grid of sensitivity parameters, each defining a linear functional b' theta
#' (e.g. a bound or point estimate), computes point estimates and then uses the
#' multiplier bootstrap to construct simultaneous confidence intervals or bands.
#' Let psi_i = b' phi_i be the influence function for the functional. The bootstrap
#' draws G_1,...,G_n iid (e.g. Rademacher or N(0,1)) and forms T* = sqrt(n) * mean(G_i * psi_i).
#' Repeating B times approximates the distribution of sqrt(n)(thetahat - theta).
#'
#' @param mar_result Output from \code{mar_bounds()} containing \code{phi} and optionally \code{nuisance}.
#' @param param_grid Data.frame or list; each row (or combination) gives a parameter set. Columns must include those needed for the coefficient vector (e.g. delta_0u, delta_1u, tau, etc.). Alternatively, a list of vectors with names \code{delta_0u}, \code{delta_1u}, etc., which will be expanded via \code{expand.grid}.
#' @param bound_spec Character: "lower", "upper", or "point". For "lower"/"upper" the grid rows define the sensitivity parameters for that bound (e.g. delta_0u, delta_1u for monotonicity_pos). For "point" the grid gives (delta_0, delta_1, tau) for point_ate.
#' @param assumption Same as in \code{mar_bounds}: "bounded_delta", "monotonicity_pos", "monotonicity_neg", "bounded_risk", "point_ate".
#' @param B Number of multiplier bootstrap replications (default 1000).
#' @param alpha Significance level for simultaneous CIs (default 0.05).
#' @param multiplier "rademacher" (default) or "gaussian". Multiplier draws use
#'   the current RNG state; call \code{set.seed()} beforehand for reproducibility.
#' @return List with \code{grid} (data.frame with columns from param_grid plus \code{estimate}, \code{se}, \code{ci_lower}, \code{ci_upper}, \code{simultaneous_lower}, \code{simultaneous_upper}), \code{boot_replicates} (B x nrow(grid) matrix of bootstrap T* values), and \code{critical_value} (estimated 1-alpha quantile of max over grid of |T*|).
#' @export
multiplier_bootstrap_grid <- function(mar_result,
                                      param_grid,
                                      bound_spec = c("lower", "upper", "point"),
                                      assumption = c("bounded_delta", "monotonicity_pos", "monotonicity_neg", "bounded_risk", "point_ate"),
                                      B = 1000,
                                      alpha = 0.05,
                                      multiplier = c("rademacher", "gaussian")) {
  bound_spec <- match.arg(bound_spec)
  assumption <- match.arg(assumption)
  multiplier <- match.arg(multiplier)
  phi <- mar_result$phi
  if (is.null(phi)) stop("mar_result must contain phi.")
  n <- nrow(phi)

  if (is.list(param_grid) && !is.data.frame(param_grid)) {
    param_grid <- expand.grid(param_grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  }
  grid_df <- as.data.frame(param_grid)
  n_grid <- nrow(grid_df)

  # Coefficient vector b for each row
  get_b <- function(row) {
    delta_0u <- row$delta_0u %||% 1
    delta_1u <- row$delta_1u %||% 1
    delta_0l <- row$delta_0l %||% 0
    delta_1l <- row$delta_1l %||% 0
    tau_0 <- row$tau_0 %||% row$tau %||% 2
    tau_1 <- row$tau_1 %||% row$tau %||% tau_0
    tau <- row$tau %||% 2
    delta_0 <- row$delta_0
    delta_1 <- row$delta_1
    if (assumption == "point_ate" && bound_spec == "point") {
      return(coef_point_ate(delta_0 %||% 0.5, delta_1 %||% 0.5, tau))
    }
    if (assumption == "monotonicity_pos") {
      if (bound_spec == "lower") return(coef_mono_pos_lower_ate(delta_0u))
      return(coef_mono_pos_upper_ate(delta_1u))
    }
    if (assumption == "monotonicity_neg") {
      if (bound_spec == "lower") return(coef_mono_neg_lower_ate(delta_1u))
      return(coef_mono_neg_upper_ate(delta_0u))
    }
    if (assumption == "bounded_delta") {
      if (bound_spec == "lower") return(coef_delta_lower_ate(delta_0u, delta_1u))
      return(coef_delta_upper_ate(delta_0u, delta_1u))
    }
    coef_point_ate(delta_0 %||% 0.5, delta_1 %||% 0.5, tau)
  }

  # Point estimates and IF values (n x n_grid)
  estimates <- numeric(n_grid)
  psi_mat <- matrix(NA, n, n_grid)
  if (assumption == "bounded_risk") {
    if (bound_spec == "point") stop("bounded_risk supports bound_spec = 'lower' or 'upper'.")
    nuis <- mar_result$nuisance
    if (is.null(nuis$mu0) || is.null(nuis$mu1)) {
      stop("For bounded_risk, mar_result must contain nuisance$mu0 and nuisance$mu1.")
    }
  }
  for (j in seq_len(n_grid)) {
    row <- as.list(grid_df[j, , drop = FALSE])
    if (assumption == "bounded_risk") {
      # Same masked min{1 - mu, (tau - 1) mu} score as mar_bounds()
      tau_0 <- row$tau_0 %||% row$tau %||% 2
      tau_1 <- row$tau_1 %||% row$tau %||% tau_0
      score <- bounded_risk_scores(
        phi, nuis$mu0, nuis$mu1, row$delta_0u %||% 1, row$delta_1u %||% 1,
        tau_0, tau_1, bound_spec
      )
    } else {
      score <- drop(phi %*% get_b(row))
    }
    estimates[j] <- mean(score)
    psi_mat[, j] <- score - estimates[j]  # centered IF
  }

  # Multiplier bootstrap
  if (multiplier == "rademacher") {
    G_fun <- function() sample(c(-1, 1), n, replace = TRUE)
  } else {
    G_fun <- function() stats::rnorm(n)
  }
  boot_mat <- matrix(NA, B, n_grid)
  for (b in seq_len(B)) {
    G <- G_fun()
    for (j in seq_len(n_grid)) {
      boot_mat[b, j] <- sqrt(n) * mean(G * psi_mat[, j])
    }
  }

  # Pointwise SE and CI
  se <- sqrt(apply(psi_mat, 2, function(x) sum(x^2) / n) / n)
  z <- stats::qnorm(1 - alpha / 2)
  ci_lower <- estimates - z * se
  ci_upper <- estimates + z * se

  # Simultaneous: multiplier bootstrap for uniform bands. Let
  #   T*_j = sqrt(n) * mean(G * psi_j),
  # and form C as the (1-alpha) quantile of max_j |T*_j| / sigma_j, where
  # sigma_j^2 = Var(psi_j). Then uniform bands are estimate +/- C * SE_pointwise.
  max_abs <- apply(boot_mat, 1, function(row) max(abs(row)))
  sigma_j <- sqrt(apply(psi_mat, 2, function(x) sum(x^2) / n))
  max_std <- apply(boot_mat, 1, function(row) max(abs(row) / sigma_j))
  critical_value <- stats::quantile(max_std, probs = 1 - alpha)
  simultaneous_lower <- estimates - critical_value * se
  simultaneous_upper <- estimates + critical_value * se

  grid_df$estimate <- estimates
  grid_df$se <- se
  grid_df$ci_lower <- ci_lower
  grid_df$ci_upper <- ci_upper
  grid_df$simultaneous_lower <- simultaneous_lower
  grid_df$simultaneous_upper <- simultaneous_upper

  list(
    grid = grid_df,
    boot_replicates = boot_mat,
    critical_value = critical_value,
    alpha = alpha
  )
}
