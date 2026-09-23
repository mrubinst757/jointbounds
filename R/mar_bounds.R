#' Bounds and point estimates for causal effects under mixed missingness
#'
#' Implements the methods from the paper: bounds on ATE, composite ATE (Psi_1), and
#' separable direct effect (Psi_2), with user-specified estimand, assumptions, and
#' sensitivity parameters. Nuisance functions are estimated via SuperLearner with
#' cross-fitting.
#'
#' @param data A data.frame containing the analysis data.
#' @param Y Character; column name of the outcome (numeric; use NA when missing).
#' @param A Character; column name of treatment (0/1).
#' @param C Character; column name of missingness indicator (1 = missing, 0 = observed).
#' @param X Character vector; column names of covariates.
#' @param estimand One of "ate" (average treatment effect), "psi1" (composite ATE), "psi2" (separable direct effect).
#' @param assumption One of "general", "bounded_delta", "monotonicity_pos", "monotonicity_neg", "bounded_risk", "bounded_risk_unbounded_tau", "point_ate", "point_psi1", "point_psi2".
#'   \describe{
#'     \item{general}{No extra assumptions; widest bounds (ATE/Psi_1).}
#'     \item{bounded_delta}{Proportion of informative missingness bounded by delta_0u, delta_1u (and delta_0l, delta_1l for Psi_1).}
#'     \item{monotonicity_pos}{Informative missingness increases outcome risk; use with delta_0u, delta_1u.}
#'     \item{monotonicity_neg}{Informative missingness decreases outcome risk.}
#'     \item{bounded_risk}{Outcome risk ratio in informatively missing bounded by tau_0, tau_1, with an upper bound that uses the minimum of 1-mu and mu*(tau-1) (as in the paper).}
#'     \item{bounded_risk_unbounded_tau}{Original bounded-risk formula without applying the minimum-based tau constraint (provides backward-compatible behavior).}
#'     \item{point_ate}{Point identification under known tau and delta_0, delta_1 (Assumption known-risk + known-missingness).}
#'     \item{point_psi1}{Point identification of Psi_1 under known delta_0, delta_1.}
#'     \item{point_psi2}{Point identification of Psi_2 under known delta_0.}
#'   }
#' @param delta_0u,delta_1u Upper bounds on proportion of missingness that is informative (in [0,1]). Default 1.
#' @param delta_0l,delta_1l Lower bounds (for Psi_1 bounds). Default 0.
#' @param delta_0,delta_1 Point values for point identification (when assumption is point_*).
#' @param tau_0,tau_1 Bounded risk parameters (>= 1) for "bounded_risk". Single \code{tau} for "point_ate" (same in both arms).
#' @param tau Single sensitivity parameter for "point_ate" (risk ratio in both arms).
#' @param V Number of cross-fitting folds for nuisance estimation (default 2; uses 1 when
#'   \code{SL.glm} is the only library for all nuisance functions).
#' @param sl_lib A single SuperLearner library specification used for all nuisance models
#'   by default (propensity, missingness, outcome); defaults to \code{"SL.glm"}.
#' @param sl_lib_prop,sl_lib_miss,sl_lib_outcome Optional SuperLearner libraries for
#'   propensity, missingness, and outcome models. If \code{NULL}, they fall back to
#'   \code{sl_lib}.
#' @param stratify_mu Logical; if TRUE (default), estimate separate outcome models for A=0 and A=1. If FALSE, estimate a single pooled model E(Y|X,C=0).
#' @param family_Y Character; family for outcome model ("gaussian" or "binomial"). If NULL (default), automatically detects: "binomial" for binary Y (only 0/1 values), otherwise "gaussian". Set explicitly to override auto-detection.
#' @param smooth_approximation Logical; for Psi_2 bounds under bounded_delta, use smooth approximation (TRUE) or indicator function method (FALSE, default).
#' @param epsilon Smoothing parameter for Psi_2 smooth approximation (default 0.001). Only used when smooth_approximation = TRUE.
#' @param seed Optional seed for cross-fitting splits.
#' @param param_grid Optional list or data.frame describing a grid of sensitivity
#'   parameters for which to compute bounds and/or point estimates and confidence
#'   bands. For \code{estimand = "ate"}: \code{delta_0u}, \code{delta_1u}, \code{tau_0},
#'   \code{tau_1}, \code{tau}, \code{delta_0}, \code{delta_1} (by assumption). For
#'   \code{estimand = "psi1"} (bounded_delta): \code{delta_0u}, \code{delta_1u},
#'   \code{delta_0l}, \code{delta_1l}; (point_psi1): \code{delta_0}, \code{delta_1}.
#'   For \code{estimand = "psi2"} (point_psi2): \code{delta_0}.
#'   Column names may also be provided without underscores (e.g., \code{delta0u}).
#' @param grid_bounds Which bounds to compute over the grid: \code{"both"} (default),
#'   \code{"lower"}, or \code{"upper"}. For point-identified grids (point_ate, point_psi1,
#'   point_psi2) this is ignored and a single \code{estimate_grid} is returned.
#' @param alpha Significance level for pointwise and uniform (multiplier bootstrap)
#'   confidence intervals over the grid (default 0.05).
#' @param B Number of multiplier bootstrap replications for uniform bands over the grid
#'   (default 1000).
#' @param multiplier Multiplier distribution for the bootstrap: \code{"rademacher"}
#'   (default) or \code{"gaussian"}.
#' @return A list with components: \code{naive} (naive estimate assuming MAR), \code{lower}, \code{upper} (for bounds), or \code{estimate} (for point ID), with \code{se_*} or \code{se} and optional \code{nuisance}, \code{phi} for downstream use.
#'   \code{result} contains a convenient dataframe summary: for scalar results, it includes the estimates and confidence intervals; for grid results (\code{param_grid} provided), it contains the combined grid output in long format with \code{bound_type} column for bounds.
#'   When \code{param_grid} is provided, the list also contains \code{lower_grid}, \code{upper_grid}, and/or \code{estimate_grid} as separate list objects with pointwise and uniform inference.
#' @export
mar_bounds <- function(data,
                       Y,
                       A,
                       C,
                       X,
                       estimand = c("ate", "psi1", "psi2"),
                       assumption = c("general", "bounded_delta", "monotonicity_pos", "monotonicity_neg", "bounded_risk", "bounded_risk_unbounded_tau", "point_ate", "point_psi1", "point_psi2"),
                       delta_0u = 1,
                       delta_1u = 1,
                       delta_0l = 0,
                       delta_1l = 0,
                       delta_0 = NULL,
                       delta_1 = NULL,
                       tau_0 = NULL,
                       tau_1 = NULL,
                       tau = NULL,
                       V = 2,
                       sl_lib = "SL.glm",
                       sl_lib_prop = NULL,
                       sl_lib_miss = NULL,
                       sl_lib_outcome = NULL,
                       stratify_mu = TRUE,
                       family_Y = NULL,
                       smooth_approximation = FALSE,
                       epsilon = 0.001,
                       seed = NULL,
                       param_grid = NULL,
                       grid_bounds = c("both", "lower", "upper"),
                       alpha = 0.05,
                       B = 1000,
                       multiplier = c("rademacher", "gaussian")) {
  estimand <- match.arg(estimand)
  assumption <- match.arg(assumption)
  grid_bounds <- match.arg(grid_bounds)
  multiplier <- match.arg(multiplier)

  # Alias delta parameters: if delta_0/delta_1 provided for bounds assumptions,
  # use them as delta_0u/delta_1u
  if (assumption %in% c("general", "bounded_delta", "monotonicity_pos", "monotonicity_neg",
                        "bounded_risk", "bounded_risk_unbounded_tau")) {
    if (!is.null(delta_0) || !is.null(delta_1)) {
      # Add warning for bounded_risk assumptions
      if (assumption %in% c("bounded_risk", "bounded_risk_unbounded_tau")) {
        warning("For ", assumption, ", delta_0 and delta_1 are interpreted as delta_0u and delta_1u (upper bounds on proportion of informative missingness).")
      }
      if (!is.null(delta_0)) delta_0u <- delta_0
      if (!is.null(delta_1)) delta_1u <- delta_1
    }
  }

  # Resolve SuperLearner libraries: default to sl_lib when not provided
  if (is.null(sl_lib_prop)) sl_lib_prop <- sl_lib
  if (is.null(sl_lib_miss)) sl_lib_miss <- sl_lib
  if (is.null(sl_lib_outcome)) sl_lib_outcome <- sl_lib

  prep <- prepare_data(data, Y, A, C, X)
  Y_vec <- prep$Y_vec
  A_vec <- prep$A_vec
  C_vec <- prep$C_vec
  X_mat <- prep$X_mat
  n <- prep$n

  # Warn about small sample sizes
  if (n < 100) {
    warning("Small sample size (n = ", n, "). SuperLearner estimation may be unstable with n < 100. Consider increasing sample size for more reliable results.")
  }

  # Auto-detect binary outcome if family_Y not specified
  if (is.null(family_Y)) {
    Y_obs <- Y_vec[!is.na(Y_vec)]
    if (length(Y_obs) > 0) {
      is_binary <- all(Y_obs %in% c(0, 1))
      if (is_binary) {
        message("Y appears to be binary (only 0/1 values). Automatically setting family_Y = 'binomial'. To override, explicitly set family_Y = 'gaussian'.")
        family_Y <- "binomial"
      } else {
        family_Y <- "gaussian"
      }
    } else {
      family_Y <- "gaussian"
    }
  }

  # Default sample splitting: V = 2, but if only SL.glm is used everywhere and V was
  # not explicitly set by the user, use V = 1.
  V_cf <- V
  if (missing(V)) {
    simple_libs <- length(sl_lib_prop) == 1L && length(sl_lib_miss) == 1L &&
      length(sl_lib_outcome) == 1L &&
      identical(sl_lib_prop, "SL.glm") &&
      identical(sl_lib_miss, "SL.glm") &&
      identical(sl_lib_outcome, "SL.glm")
    if (simple_libs) {
      V_cf <- 1L
    }
  }

  nuis <- estimate_nuisance(
    X = X_mat, A = A_vec, C = C_vec, Y = Y_vec,
    V = V_cf,
    sl_lib_prop = sl_lib_prop,
    sl_lib_miss = sl_lib_miss,
    sl_lib_outcome = sl_lib_outcome,
    stratify_mu = stratify_mu,
    family_Y = family_Y,
    seed = seed
  )

  phi <- influence_functions(
    Y = Y_vec, A = A_vec, C = C_vec,
    e = nuis$e, pi0 = nuis$pi0, pi1 = nuis$pi1,
    mu0 = nuis$mu0, mu1 = nuis$mu1
  )

  if (assumption == "point_ate") {
    if (is.null(param_grid) && (is.null(delta_0) || is.null(delta_1) || is.null(tau))) {
      stop("For point_ate specify delta_0, delta_1, and tau when param_grid is not provided.")
    }
  }
  if (assumption == "point_psi1") {
    if (is.null(param_grid) && (is.null(delta_0) || is.null(delta_1))) {
      stop("For point_psi1 specify delta_0 and delta_1 when param_grid is not provided.")
    }
  }
  if (assumption == "point_psi2") {
    if (is.null(param_grid) && is.null(delta_0)) stop("For point_psi2 specify delta_0 when param_grid is not provided.")
  }
  if (assumption == "bounded_risk") {
    if (is.null(tau_0)) tau_0 <- 2
    if (is.null(tau_1)) tau_1 <- tau_0
  }

  if (assumption == "point_psi2") {
    out <- compute_bounds(
      phi, estimand = estimand, assumption = assumption,
      delta_0 = delta_0, delta_1 = delta_1,
      pi0 = nuis$pi0, mu1 = nuis$mu1,
      smooth_approximation = smooth_approximation,
      epsilon = epsilon
    )
  } else if (assumption == "bounded_risk") {
    out <- compute_bounds(
      phi,
      estimand = estimand,
      assumption = assumption,
      delta_0u = delta_0u, delta_1u = delta_1u,
      delta_0l = delta_0l, delta_1l = delta_1l,
      delta_0 = delta_0, delta_1 = delta_1,
      tau_0 = tau_0, tau_1 = tau_1,
      tau = tau,
      mu0 = nuis$mu0, mu1 = nuis$mu1,
      pi0 = nuis$pi0, pi1 = nuis$pi1,
      smooth_approximation = smooth_approximation,
      epsilon = epsilon
    )
  } else {
    out <- compute_bounds(
      phi,
      estimand = estimand,
      assumption = assumption,
      delta_0u = delta_0u, delta_1u = delta_1u,
      delta_0l = delta_0l, delta_1l = delta_1l,
      delta_0 = delta_0, delta_1 = delta_1,
      tau_0 = tau_0, tau_1 = tau_1,
      tau = tau,
      mu0 = nuis$mu0, mu1 = nuis$mu1,
      pi0 = nuis$pi0, pi1 = nuis$pi1,
      smooth_approximation = smooth_approximation,
      epsilon = epsilon
    )
  }

  out$nuisance <- nuis
  out$phi <- phi
  out$estimand <- estimand
  out$assumption <- assumption

  # Optional: simultaneous bounds/estimates over a parameter grid (ATE, Psi_1, Psi_2)
  grid_ok_ate <- estimand == "ate" &&
    assumption %in% c("general", "bounded_delta", "bounded_risk_unbounded_tau", "monotonicity_pos", "monotonicity_neg", "bounded_risk", "point_ate")
  grid_ok_psi1 <- estimand == "psi1" && assumption %in% c("bounded_delta", "point_psi1")
  grid_ok_psi2 <- estimand == "psi2" && assumption == "point_psi2"
  if (!is.null(param_grid) && (grid_ok_ate || grid_ok_psi1 || grid_ok_psi2)) {
    grid_res <- mar_bounds_grid_bands(
      phi = phi,
      estimand = estimand,
      assumption = assumption,
      param_grid = param_grid,
      grid_bounds = grid_bounds,
      alpha = alpha,
      B = B,
      multiplier = multiplier,
      pi0 = nuis$pi0,
      pi1 = nuis$pi1,
      mu0 = nuis$mu0,
      mu1 = nuis$mu1
    )
    if (!is.null(grid_res$lower_grid)) out$lower_grid <- grid_res$lower_grid
    if (!is.null(grid_res$upper_grid)) out$upper_grid <- grid_res$upper_grid
    if (!is.null(grid_res$estimate_grid)) out$estimate_grid <- grid_res$estimate_grid

    # For grid results, put the combined output in 'result'
    if (!is.null(grid_res$estimate_grid)) {
      # Point estimate grid
      out$result <- grid_res$estimate_grid$grid
    } else if (!is.null(grid_res$lower_grid) && !is.null(grid_res$upper_grid)) {
      # Bounds grid: combine into long format
      lower_df <- grid_res$lower_grid$grid
      lower_df$bound_type <- "lower"
      upper_df <- grid_res$upper_grid$grid
      upper_df$bound_type <- "upper"
      out$result <- rbind(lower_df, upper_df)
    } else if (!is.null(grid_res$lower_grid)) {
      out$result <- grid_res$lower_grid$grid
      out$result$bound_type <- "lower"
    } else if (!is.null(grid_res$upper_grid)) {
      out$result <- grid_res$upper_grid$grid
      out$result$bound_type <- "upper"
    }
  } else {
    # Create summary dataframe for scalar results (no grid)
    z <- stats::qnorm(1 - alpha / 2)
    if (!is.null(out$estimate)) {
      # Point estimate
      out$result <- data.frame(
        naive = out$naive,
        estimate = out$estimate,
        se = out$se,
        ci_lower = out$estimate - z * out$se,
        ci_upper = out$estimate + z * out$se
      )
      # Clean up redundant fields now in result dataframe
      out$naive <- NULL
      out$estimate <- NULL
      out$se <- NULL
    } else if (!is.null(out$lower) && !is.null(out$upper)) {
      # Bounds
      out$result <- data.frame(
        naive = out$naive,
        lower = out$lower,
        upper = out$upper,
        se_lower = out$se_lower,
        se_upper = out$se_upper,
        ci_lower_lower = out$lower - z * out$se_lower,
        ci_lower_upper = out$lower + z * out$se_lower,
        ci_upper_lower = out$upper - z * out$se_upper,
        ci_upper_upper = out$upper + z * out$se_upper
      )
      # Clean up redundant fields now in result dataframe
      out$naive <- NULL
      out$lower <- NULL
      out$upper <- NULL
      out$se_lower <- NULL
      out$se_upper <- NULL
    }
  }

  out
}

# Internal helper: compute pointwise and uniform (multiplier bootstrap) bands over a grid
mar_bounds_grid_bands <- function(phi,
                                  estimand,
                                  assumption,
                                  param_grid,
                                  grid_bounds = c("both", "lower", "upper"),
                                  alpha = 0.05,
                                  B = 1000,
                                  multiplier = c("rademacher", "gaussian"),
                                  pi0 = NULL,
                                  pi1 = NULL,
                                  mu0 = NULL,
                                  mu1 = NULL) {
  grid_bounds <- match.arg(grid_bounds)
  multiplier <- match.arg(multiplier)

  if (is.list(param_grid) && !is.data.frame(param_grid)) {
    grid_df <- expand.grid(param_grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  } else {
    grid_df <- as.data.frame(param_grid)
  }
  n_grid <- nrow(grid_df)
  if (n_grid < 1L) stop("param_grid must define at least one grid point.")

  # Support alternate naming conventions for grid columns (e.g., delta0u vs delta_0u)
  get_grid_param <- function(row, names, default = NULL) {
    for (nm in names) {
      if (!is.null(row[[nm]])) return(row[[nm]])
    }
    default
  }

  n <- nrow(phi)
  z <- stats::qnorm(1 - alpha / 2)
  if (multiplier == "rademacher") {
    G_fun <- function() sample(c(-1, 1), n, replace = TRUE)
  } else {
    G_fun <- function() stats::rnorm(n)
  }

  # Build grid output from estimates + psi_mat (n x n_grid centered influence)
  compute_critical_value <- function(psi_mats) {
    # psi_mats: list of matrices (n x n_grid)
    sigma_list <- lapply(psi_mats, function(psi) sqrt(apply(psi, 2L, function(x) sum(x^2) / n)))
    max_std <- numeric(B)
    for (bb in seq_len(B)) {
      G <- G_fun()
      max_val <- -Inf
      for (k in seq_along(psi_mats)) {
        psi <- psi_mats[[k]]
        sigma_j <- sigma_list[[k]]
        vals <- sqrt(n) * colMeans(G * psi) / sigma_j
        max_val <- max(max_val, max(abs(vals), na.rm = TRUE))
      }
      max_std[bb] <- max_val
    }
    as.numeric(stats::quantile(max_std, probs = 1 - alpha))
  }

  make_grid_result <- function(estimates, psi_mat, grid_df, critical_value = NULL) {
    se <- sqrt(apply(psi_mat, 2L, function(x) sum(x^2) / n) / n)
    ci_lower_pt <- estimates - z * se
    ci_upper_pt <- estimates + z * se
    if (is.null(critical_value)) {
      critical_value <- compute_critical_value(list(psi_mat))
    }
    ci_lower_unif <- estimates - critical_value * se
    ci_upper_unif <- estimates + critical_value * se

    grid_out <- grid_df
    grid_out$estimate <- estimates
    grid_out$se <- se
    grid_out$ci_lower_pointwise <- ci_lower_pt
    grid_out$ci_upper_pointwise <- ci_upper_pt
    grid_out$ci_lower_uniform <- ci_lower_unif
    grid_out$ci_upper_uniform <- ci_upper_unif

    # Drop parameter columns that do not vary (i.e., the bound does not depend on them)
    metric_cols <- c("estimate", "se", "ci_lower_pointwise", "ci_upper_pointwise",
                     "ci_lower_uniform", "ci_upper_uniform")
    param_cols <- setdiff(names(grid_out), metric_cols)
    if (length(param_cols) > 0) {
      const_params <- vapply(grid_out[param_cols], function(x) {
        is.numeric(x) && length(unique(x)) == 1
      }, logical(1))
      keep_params <- param_cols[!const_params]
      if (length(keep_params) < length(param_cols)) {
        grid_out <- grid_out[c(keep_params, metric_cols)]
      }
    }

    list(grid = grid_out, critical_value = critical_value, alpha = alpha, multiplier = multiplier)
  }

  # ----- ATE -----
  # Resolve per-arm ATE parameters from a grid row. Arm-specific columns take
  # precedence, then delta_0/delta_1, then the delta/tau shorthands; numeric
  # defaults are applied only after the shorthands (tau_1 defaults to tau_0,
  # matching the scalar bounded_risk path).
  resolve_ate_params <- function(row) {
    delta_0u <- get_grid_param(row, c("delta_0u", "delta0u"))
    delta_1u <- get_grid_param(row, c("delta_1u", "delta1u"))
    delta <- get_grid_param(row, c("delta"))
    tau_0 <- get_grid_param(row, c("tau_0", "tau0"))
    tau_1 <- get_grid_param(row, c("tau_1", "tau1"))
    tau <- get_grid_param(row, c("tau"))
    delta_0 <- get_grid_param(row, c("delta_0", "delta0"))
    delta_1 <- get_grid_param(row, c("delta_1", "delta1"))

    if (is.null(delta_0u)) delta_0u <- delta_0 %||% delta %||% 1
    if (is.null(delta_1u)) delta_1u <- delta_1 %||% delta %||% 1
    if (is.null(tau_0)) tau_0 <- tau
    if (is.null(tau_1)) tau_1 <- tau
    if (assumption == "bounded_risk") {
      if (is.null(tau_0)) tau_0 <- 2
      if (is.null(tau_1)) tau_1 <- tau_0
    }
    list(delta_0u = delta_0u, delta_1u = delta_1u, tau_0 = tau_0, tau_1 = tau_1,
         tau = tau, delta_0 = delta_0, delta_1 = delta_1)
  }

  get_b_ate <- function(row, bound_type) {
    p <- resolve_ate_params(row)
    delta_0u <- p$delta_0u
    delta_1u <- p$delta_1u
    tau_0 <- p$tau_0
    tau_1 <- p$tau_1
    tau <- p$tau
    delta_0 <- p$delta_0
    delta_1 <- p$delta_1

    if (assumption == "point_ate" && !is.null(delta_0) && !is.null(delta_1)) {
      return(coef_point_ate(delta_0, delta_1, tau %||% 2))
    }
    if (assumption == "monotonicity_pos") {
      if (bound_type == "lower") return(coef_mono_pos_lower_ate(delta_0u))
      return(coef_mono_pos_upper_ate(delta_1u))
    }
    if (assumption == "monotonicity_neg") {
      if (bound_type == "lower") return(coef_mono_neg_lower_ate(delta_1u))
      return(coef_mono_neg_upper_ate(delta_0u))
    }
    if (assumption == "bounded_delta") {
      if (bound_type == "lower") return(coef_delta_lower_ate(delta_0u, delta_1u))
      return(coef_delta_upper_ate(delta_0u, delta_1u))
    }
    if (assumption == "bounded_risk_unbounded_tau") {
      if (bound_type == "lower") return(coef_bounded_risk_lower_ate(delta_0u, delta_1u, tau_0, tau_1))
      return(coef_bounded_risk_upper_ate(delta_0u, delta_1u, tau_0, tau_1))
    }
    if (assumption == "bounded_risk") {
      if (bound_type == "lower") return(coef_bounded_risk_lower_ate(delta_0u, delta_1u, tau_0, tau_1))
      return(coef_bounded_risk_upper_ate(delta_0u, delta_1u, tau_0, tau_1))
    }
    if (assumption == "general") {
      # Corollary 1: when delta parameters are specified, general bounds shrink
      if (!is.null(delta_0u) && !is.null(delta_1u) && (delta_0u != 1 || delta_1u != 1)) {
        if (bound_type == "lower") return(coef_delta_lower_ate(delta_0u, delta_1u))
        return(coef_delta_upper_ate(delta_0u, delta_1u))
      }
      if (bound_type == "lower") return(coef_general_lower_ate())
      return(coef_general_upper_ate())
    }
    coef_point_ate(delta_0 %||% 0.5, delta_1 %||% 0.5, tau %||% 2)
  }

  # ----- Psi_1 -----
  get_b_psi1 <- function(row, bound_type) {
    delta_0u <- get_grid_param(row, c("delta_0u", "delta0u"))
    delta_1u <- get_grid_param(row, c("delta_1u", "delta1u"))
    delta_0l <- get_grid_param(row, c("delta_0l", "delta0l"), 0)
    delta_1l <- get_grid_param(row, c("delta_1l", "delta1l"), 0)
    delta <- get_grid_param(row, c("delta"))
    delta_0 <- get_grid_param(row, c("delta_0", "delta0"))
    delta_1 <- get_grid_param(row, c("delta_1", "delta1"))

    # Arm-specific columns first, then delta_0/delta_1, then the delta shorthand
    if (is.null(delta_0u)) delta_0u <- delta_0 %||% delta %||% 1
    if (is.null(delta_1u)) delta_1u <- delta_1 %||% delta %||% 1

    if (assumption == "point_psi1" && !is.null(delta_0) && !is.null(delta_1)) {
      return(coef_point_psi1(delta_0, delta_1))
    }
    if (assumption == "bounded_delta") {
      if (bound_type == "lower") return(coef_lower_psi1(delta_1l, delta_0u))
      return(coef_upper_psi1(delta_1u, delta_0l))
    }
    coef_point_psi1(delta_0 %||% 0.5, delta_1 %||% 0.5)
  }

  res <- list(lower_grid = NULL, upper_grid = NULL, estimate_grid = NULL)

  if (estimand == "ate") {
    if (assumption == "point_ate") {
      estimates <- numeric(n_grid)
      psi_mat <- matrix(NA_real_, nrow = n, ncol = n_grid)
      for (j in seq_len(n_grid)) {
        row <- as.list(grid_df[j, , drop = FALSE])
        b <- get_b_ate(row, "lower")
        val <- as.numeric(mean(phi %*% b))
        estimates[j] <- val
        psi_mat[, j] <- drop(phi %*% b) - val
      }
      res$estimate_grid <- make_grid_result(estimates, psi_mat, grid_df)
      return(res)
    }

    if (assumption == "bounded_risk") {
      if (is.null(mu0) || is.null(mu1) || is.null(pi0) || is.null(pi1)) stop("For bounded_risk grid bounds provide mu0, mu1, pi0, and pi1.")

      if (grid_bounds %in% c("both", "lower")) {
        lower_estimates <- numeric(n_grid)
        lower_psi <- matrix(NA_real_, nrow = n, ncol = n_grid)
        for (j in seq_len(n_grid)) {
          p <- resolve_ate_params(as.list(grid_df[j, , drop = FALSE]))
          lower_vec <- bounded_risk_scores(phi, mu0, mu1, p$delta_0u, p$delta_1u,
                                           p$tau_0, p$tau_1, "lower")

          val <- mean(lower_vec)
          lower_estimates[j] <- val
          lower_psi[, j] <- lower_vec - val
        }
      }
      if (grid_bounds %in% c("both", "upper")) {
        upper_estimates <- numeric(n_grid)
        upper_psi <- matrix(NA_real_, nrow = n, ncol = n_grid)
        for (j in seq_len(n_grid)) {
          p <- resolve_ate_params(as.list(grid_df[j, , drop = FALSE]))
          upper_vec <- bounded_risk_scores(phi, mu0, mu1, p$delta_0u, p$delta_1u,
                                           p$tau_0, p$tau_1, "upper")

          val <- mean(upper_vec)
          upper_estimates[j] <- val
          upper_psi[, j] <- upper_vec - val
        }
      }

      if (grid_bounds == "both") {
        shared_cv <- compute_critical_value(list(lower_psi, upper_psi))
        res$lower_grid <- make_grid_result(lower_estimates, lower_psi, grid_df, critical_value = shared_cv)
        res$upper_grid <- make_grid_result(upper_estimates, upper_psi, grid_df, critical_value = shared_cv)
        return(res)
      }
      if (grid_bounds == "lower") {
        res$lower_grid <- make_grid_result(lower_estimates, lower_psi, grid_df)
        return(res)
      }
      if (grid_bounds == "upper") {
        res$upper_grid <- make_grid_result(upper_estimates, upper_psi, grid_df)
        return(res)
      }
      return(res)
    }

    if (grid_bounds %in% c("both", "lower")) {
      lower_estimates <- numeric(n_grid)
      lower_psi <- matrix(NA_real_, nrow = n, ncol = n_grid)
      for (j in seq_len(n_grid)) {
        row <- as.list(grid_df[j, , drop = FALSE])
        b <- get_b_ate(row, "lower")
        val <- as.numeric(mean(phi %*% b))
        lower_estimates[j] <- val
        lower_psi[, j] <- drop(phi %*% b) - val
      }
    }
    if (grid_bounds %in% c("both", "upper")) {
      upper_estimates <- numeric(n_grid)
      upper_psi <- matrix(NA_real_, nrow = n, ncol = n_grid)
      for (j in seq_len(n_grid)) {
        row <- as.list(grid_df[j, , drop = FALSE])
        b <- get_b_ate(row, "upper")
        val <- as.numeric(mean(phi %*% b))
        upper_estimates[j] <- val
        upper_psi[, j] <- drop(phi %*% b) - val
      }
    }

    if (grid_bounds == "both") {
      shared_cv <- compute_critical_value(list(lower_psi, upper_psi))
      res$lower_grid <- make_grid_result(lower_estimates, lower_psi, grid_df, critical_value = shared_cv)
      res$upper_grid <- make_grid_result(upper_estimates, upper_psi, grid_df, critical_value = shared_cv)
      return(res)
    }

    if (grid_bounds == "lower") {
      res$lower_grid <- make_grid_result(lower_estimates, lower_psi, grid_df)
      return(res)
    }
    if (grid_bounds == "upper") {
      res$upper_grid <- make_grid_result(upper_estimates, upper_psi, grid_df)
      return(res)
    }
    return(res)
  }

  if (estimand == "psi1") {
    if (assumption == "point_psi1") {
      estimates <- numeric(n_grid)
      psi_mat <- matrix(NA_real_, nrow = n, ncol = n_grid)
      for (j in seq_len(n_grid)) {
        row <- as.list(grid_df[j, , drop = FALSE])
        b <- get_b_psi1(row, "lower")
        val <- as.numeric(mean(phi %*% b))
        estimates[j] <- val
        psi_mat[, j] <- drop(phi %*% b) - val
      }
      res$estimate_grid <- make_grid_result(estimates, psi_mat, grid_df)
      return(res)
    }

    if (assumption == "bounded_delta") {
      if (grid_bounds %in% c("both", "lower")) {
        lower_estimates <- numeric(n_grid)
        lower_psi <- matrix(NA_real_, nrow = n, ncol = n_grid)
        for (j in seq_len(n_grid)) {
          row <- as.list(grid_df[j, , drop = FALSE])
          b <- get_b_psi1(row, "lower")
          val <- as.numeric(mean(phi %*% b))
          lower_estimates[j] <- val
          lower_psi[, j] <- drop(phi %*% b) - val
        }
      }
      if (grid_bounds %in% c("both", "upper")) {
        upper_estimates <- numeric(n_grid)
        upper_psi <- matrix(NA_real_, nrow = n, ncol = n_grid)
        for (j in seq_len(n_grid)) {
          row <- as.list(grid_df[j, , drop = FALSE])
          b <- get_b_psi1(row, "upper")
          val <- as.numeric(mean(phi %*% b))
          upper_estimates[j] <- val
          upper_psi[, j] <- drop(phi %*% b) - val
        }
      }

      if (grid_bounds == "both") {
        shared_cv <- compute_critical_value(list(lower_psi, upper_psi))
        res$lower_grid <- make_grid_result(lower_estimates, lower_psi, grid_df, critical_value = shared_cv)
        res$upper_grid <- make_grid_result(upper_estimates, upper_psi, grid_df, critical_value = shared_cv)
        return(res)
      }
      if (grid_bounds == "lower") {
        res$lower_grid <- make_grid_result(lower_estimates, lower_psi, grid_df)
        return(res)
      }
      if (grid_bounds == "upper") {
        res$upper_grid <- make_grid_result(upper_estimates, upper_psi, grid_df)
        return(res)
      }
    }
    return(res)
  }

  if (estimand == "psi2" && assumption == "point_psi2") {
    if (is.null(pi0) || is.null(mu1)) stop("For psi2 point_psi2 grid, pi0 and mu1 are required.")
    phi_3_01 <- phi[, 2L] * pi0 + phi[, 3L] * mu1 - mu1 * pi0
    estimates <- numeric(n_grid)
    psi_mat <- matrix(NA_real_, nrow = n, ncol = n_grid)
    for (j in seq_len(n_grid)) {
      row <- as.list(grid_df[j, , drop = FALSE])
      d0 <- get_grid_param(row, c("delta_0", "delta0"), 0.5)
      psi2_vec <- phi[, 2L] - phi[, 1L] - d0 * (phi_3_01 - phi[, 5L])
      val <- mean(psi2_vec)
      estimates[j] <- val
      psi_mat[, j] <- psi2_vec - val
    }
    res$estimate_grid <- make_grid_result(estimates, psi_mat, grid_df)
    return(res)
  }

  stop("Unsupported estimand/assumption combination for param_grid.")
}

#' Naive estimate of the ATE (Psi_tilde) assuming MAR
#'
#' Returns the influence-function-based estimate and variance using the EIF for
#' the naive estimand Psi_tilde = E[mu_1(X) - mu_0(X)], where mu_a(X) = E[Y | X, A = a, C = 0]
#' (outcome mean among observed). Package convention: C = 0 denotes observed (non-missing),
#' C = 1 denotes missing. The EIF is
#' I(A=1,C=0)/(P(C=0|X,A=1)*P(A=1|X)) * (Y - mu_1(X)) - I(A=0,C=0)/(P(C=0|X,A=0)*P(A=0|X)) * (Y - mu_0(X))
#' + mu_1(X) - mu_0(X), with variance estimated by (1/n) * sample variance of the EIF.
#'
#' @param data Data.frame.
#' @param Y,A,C,X Column names as in \code{mar_bounds}.
#' @param nuis Optional pre-computed nuisance list (e, pi0, pi1, mu0, mu1).
#' @return List with \code{estimate} (EIF-based point estimate), \code{var} (variance estimate
#'   (1/n)*Var(EIF)), \code{se} (sqrt(var)), and optionally \code{eif} (length-n vector of EIF values).
#' @export
psi_naive <- function(data, Y, A, C, X, nuis = NULL) {
  prep <- prepare_data(data, Y, A, C, X)
  if (is.null(nuis)) {
    nuis <- estimate_nuisance(
      X = prep$X_mat, A = prep$A_vec, C = prep$C_vec, Y = prep$Y_vec,
      V = 5, sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm", sl_lib_outcome = "SL.glm"
    )
  }
  Y_ <- replace(prep$Y_vec, is.na(prep$Y_vec), 0)
  A <- prep$A_vec
  C <- prep$C_vec
  e <- clip_probs(nuis$e)
  e0 <- clip_probs(1 - e)
  # P(C=0|X,A=a) = 1 - pi_a (observed given X, A=a)
  p_obs_1 <- clip_probs(1 - nuis$pi1)
  p_obs_0 <- clip_probs(1 - nuis$pi0)
  mu1 <- nuis$mu1
  mu0 <- nuis$mu0
  # EIF: I(A=1,C=0)/(P(C=0|X,A=1)*P(A=1|X))*(Y - mu_1) - I(A=0,C=0)/(P(C=0|X,A=0)*P(A=0|X))*(Y - mu_0) + mu_1 - mu_0
  term1 <- (A == 1 & C == 0) / (p_obs_1 * e) * (Y_ - mu1)
  term2 <- (A == 0 & C == 0) / (p_obs_0 * e0) * (Y_ - mu0)
  eif <- term1 - term2 + (mu1 - mu0)
  n <- prep$n
  estimate <- mean(eif)
  # Variance of psi_hat estimated by (1/n) * Var(EIF); R's var() uses (n-1) denominator
  var_est <- if (n > 1) as.numeric(stats::var(eif)) / n else 0
  list(estimate = estimate, var = var_est, se = sqrt(var_est), eif = eif)
}
