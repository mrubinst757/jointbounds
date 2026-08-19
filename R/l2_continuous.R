#' Continuous-outcome L2 sensitivity bounds
#'
#' Implements the continuous-outcome sensitivity models in the companion paper.
#' The observed-data nuisance functions are estimated by cross-fitted
#' SuperLearner or generalized linear models. Closed-form Cauchy--Schwarz outer bounds and empirical sharp
#' dual bounds are available.  The latter enforce conditional normalization on
#' a user-chosen linear basis in `X`, and therefore become exact for the empirical
#' sieve problem defined by that basis.
#'
#' @param data A data frame.
#' @param Y,A,C Column names. `A` and `C` are binary; `C = 1` means missing.
#' @param X Character vector of baseline covariate names.
#' @param delta Length-two prevalence-envelope parameters for arms 0 and 1.
#' @param delta_R,delta_A Length-two L2 radii for informative missingness and
#'   unmeasured confounding. `delta_A` is retained for backward compatibility.
#' @param delta_K Preferred alias for the confounding density-ratio radius.
#' @param delta_M Length-two L2 radii for the reduced-form net-missingness model.
#' @param model Either `"separated"` for the prevalence/severity model or
#'   `"net"` for the reduced-form model on `M - 1`.
#' @param method One of `"cs"`, `"sharp"`, or `"both"`.
#' @param folds Number of cross-fitting folds for nuisance estimation.
#' @param nuisance_method Either `"SuperLearner"` (default) or `"glm"`.
#' @param sl_lib_prop,sl_lib_miss,sl_lib_outcome SuperLearner libraries for the
#'   treatment, observation, and continuous-outcome regressions.
#' @param basis Optional formula used for conditional-normalization multipliers.
#' @param seed Random seed used to construct folds.
#' @param control Optional list passed to the empirical dual optimizer.
#' @param keep_nuisance Return row-level nuisance estimates.
#' @param return_true_bounds If `TRUE`, also compute the requested oracle Monte
#'   Carlo CS and/or sharp-sieve bounds using DGP columns returned by
#'   `simulate_l2_data()`. This is only for simulations in which the complete
#'   DGP is known.
#'
#' @return An object of class `marbounds_l2` containing arm-specific and ATE
#'   bounds, observed reference means, diagnostics, and sensitivity parameters.
#' @export
l2_bounds <- function(data, Y, A, C, X,
                      delta = c(1, 1),
                      delta_R = c(0, 0),
                      delta_A = c(0, 0),
                      delta_M = delta_R,
                      model = c("net", "separated"),
                      method = c("both", "cs", "sharp"),
                      folds = 5L,
                      nuisance_method = c("SuperLearner", "glm"),
                      sl_lib_prop = "SL.glm",
                      sl_lib_miss = "SL.glm",
                      sl_lib_outcome = "SL.glm",
                      basis = NULL,
                      seed = 1L,
                      control = list(),
                      keep_nuisance = FALSE,
                      return_true_bounds = FALSE,
                      delta_K = NULL) {
  model <- match.arg(model)
  method <- match.arg(method)
  nuisance_method <- match.arg(nuisance_method)
  delta <- l2_as_pair(delta, "delta", 1)
  delta_R <- l2_as_pair(delta_R, "delta_R")
  if (!is.null(delta_K)) delta_A <- delta_K
  delta_A <- l2_as_pair(delta_A, "delta_K")
  delta_M <- l2_as_pair(delta_M, "delta_M")
  inp <- l2_validate_inputs(data, Y, A, C, X, delta, delta_R, delta_A, delta_M)
  nuisance <- if (nuisance_method == "SuperLearner") {
    l2_nuisance_sl(data, Y, A, C, X, folds, seed, sl_lib_prop,
                   sl_lib_miss, sl_lib_outcome)
  } else l2_nuisance_glm(data, Y, A, C, X, folds = folds, seed = seed)
  basis_mat <- l2_basis(data, X, basis)

  arm <- vector("list", 2L)
  for (a in 0:1) {
    arm[[a + 1L]] <- l2_arm_bounds(
      y = inp$y, A = inp$A, C = inp$C, a = a,
      nuisance = nuisance, basis = basis_mat,
      delta = delta[a + 1L], delta_R = delta_R[a + 1L],
      delta_A = delta_A[a + 1L], delta_M = delta_M[a + 1L],
      model = model, method = method, control = control
    )
  }

  arm_df <- do.call(rbind, lapply(arm, function(z) z$summary))
  rownames(arm_df) <- NULL
  ate <- data.frame(
    method = unique(arm_df$method),
    lower = NA_real_, upper = NA_real_, estimate_reference = arm[[2]]$psi - arm[[1]]$psi
  )
  for (j in seq_len(nrow(ate))) {
    meth <- ate$method[j]
    z0 <- arm_df[arm_df$arm == 0 & arm_df$method == meth, ]
    z1 <- arm_df[arm_df$arm == 1 & arm_df$method == meth, ]
    ate$lower[j] <- z1$lower - z0$upper
    ate$upper[j] <- z1$upper - z0$lower
  }

  out <- list(
    arm = arm_df,
    ate = ate,
    diagnostics = lapply(arm, `[[`, "diagnostics"),
    parameters = list(delta = delta, delta_R = delta_R,
                      delta_A = delta_A, delta_K = delta_A, delta_M = delta_M,
                      model = model, method = method),
    nuisance_method = nuisance_method,
    moment_check = l2_moment_check(inp$y[inp$C == 0]),
    call = match.call()
  )
  if (keep_nuisance) out$nuisance <- nuisance
  if (isTRUE(return_true_bounds)) {
    true_list <- list()
    if (method %in% c("cs", "both"))
      true_list$cs <- l2_oracle_cs(data, delta, delta_R, delta_A, delta_M, model)
    if (method %in% c("sharp", "both"))
      true_list$sharp <- l2_oracle_sharp_reference(data, delta, delta_R,
        delta_A, delta_M = delta_M, model = model, basis = basis,
        control = control, return_details = TRUE)
    out$true_bounds <- if (length(true_list) == 1L) true_list[[1L]] else true_list
    out$true_bounds_note <- paste(
      "Oracle Monte Carlo bounds use complete potential outcomes and true",
      "DGP probabilities from simulate_l2_data(); sharp values retain the",
      "specified empirical normalization sieve and are unavailable for real data.")
  }
  class(out) <- "marbounds_l2"
  out
}

l2_nuisance_sl <- function(data, Y, A, C, X, folds, seed,
                           sl_lib_prop, sl_lib_miss, sl_lib_outcome,
                           family_Y = "gaussian") {
  if (!requireNamespace("SuperLearner", quietly = TRUE))
    stop("nuisance_method = 'SuperLearner' requires the SuperLearner package")
  z <- estimate_nuisance(
    X = data[, X, drop = FALSE], A = data[[A]], C = data[[C]], Y = data[[Y]],
    V = folds, sl_lib_prop = sl_lib_prop, sl_lib_miss = sl_lib_miss,
    sl_lib_outcome = sl_lib_outcome, family_Y = family_Y, seed = seed
  )
  list(e0 = clip_probs(1 - z$e, 1e-4), e1 = clip_probs(z$e, 1e-4),
       rho0 = clip_probs(1 - z$pi0, 1e-4),
       rho1 = clip_probs(1 - z$pi1, 1e-4),
       mu0 = z$mu0, mu1 = z$mu1, fold_id = z$fold_id)
}

#' @export
print.marbounds_l2 <- function(x, ...) {
  cat("Continuous-outcome L2 sensitivity bounds\n")
  cat("Model:", x$parameters$model, "\n\n")
  print(x$ate, row.names = FALSE)
  invisible(x)
}

#' Continuous-outcome L-infinity sensitivity bounds
#'
#' Computes the Holder outer comparison and/or the sharp empirical
#' conditional-normalization sieve bound under symmetric likelihood-ratio boxes
#' for informative missingness and confounding.  The outcome itself need not be
#' bounded; finite bounds require finite likelihood-ratio radii and an
#' integrable outcome under the respondent law.
#'
#' @inheritParams l2_bounds
#' @param delta_R_inf,delta_A_inf Length-two L-infinity radii for the
#'   missingness and confounding likelihood ratios.
#' @param max_n Maximum number of respondent observations per arm in the sharp
#'   linear program.  The outer bound always uses the full sample.
#' @return An object with arm-specific and ATE bounds and sharp-LP diagnostics.
#' @export
l2_linf_bounds <- function(data, Y, A, C, X, delta = c(1, 1),
                           delta_R_inf = c(0, 0),
                           delta_A_inf = c(0, 0),
                           method = c("both", "outer", "sharp"), folds = 5L,
                           nuisance_method = c("SuperLearner", "glm"),
                           sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm",
                           sl_lib_outcome = "SL.glm", basis = NULL,
                           max_n = 1000L, seed = 1L) {
  method <- match.arg(method)
  nuisance_method <- match.arg(nuisance_method)
  delta <- l2_as_pair(delta, "delta", 1)
  delta_R_inf <- l2_as_pair(delta_R_inf, "delta_R_inf")
  delta_A_inf <- l2_as_pair(delta_A_inf, "delta_A_inf")
  inp <- l2_validate_inputs(data, Y, A, C, X, delta, delta_R_inf,
                            delta_A_inf, delta_R_inf)
  max_n <- as.integer(max_n)
  if (length(max_n) != 1L || !is.finite(max_n) || max_n < 20L)
    stop("max_n must be an integer of at least 20")
  nuisance <- if (nuisance_method == "SuperLearner") {
    l2_nuisance_sl(data, Y, A, C, X, folds, seed, sl_lib_prop,
                   sl_lib_miss, sl_lib_outcome)
  } else l2_nuisance_glm(data, Y, A, C, X, folds = folds, seed = seed)
  B <- l2_basis(data, X, basis)
  arm_rows <- list(); diagnostics <- vector("list", 2L); h <- 0L
  for (a in 0:1) {
    r <- l2_reference_quantities(inp$y, inp$A, inp$C, a, nuisance)
    b <- delta[a + 1L] * r$pi
    zabs <- abs(r$z)
    rad <- mean(r$w * zabs *
      (b * delta_R_inf[a + 1L] + r$q *
       (1 + b * delta_R_inf[a + 1L]) * delta_A_inf[a + 1L]))
    if (method %in% c("outer", "both")) {
      h <- h + 1L
      arm_rows[[h]] <- data.frame(arm = a, method = "linf_outer",
        lower = r$psi - rad, upper = r$psi + rad, reference = r$psi)
    }
    if (method %in% c("sharp", "both")) {
      lp <- l2_linf_arm_lp(r, B, b, delta_R_inf[a + 1L],
                           delta_A_inf[a + 1L], max_n)
      h <- h + 1L
      arm_rows[[h]] <- data.frame(arm = a, method = "linf_sharp_sieve",
        lower = lp$lower, upper = lp$upper, reference = r$psi)
      diagnostics[[a + 1L]] <- lp$diagnostics
    }
  }
  arm <- do.call(rbind, arm_rows); rownames(arm) <- NULL
  ate <- do.call(rbind, lapply(unique(arm$method), function(meth) {
    z0 <- arm[arm$arm == 0 & arm$method == meth, , drop = FALSE]
    z1 <- arm[arm$arm == 1 & arm$method == meth, , drop = FALSE]
    data.frame(method = meth, lower = z1$lower - z0$upper,
      upper = z1$upper - z0$lower,
      estimate_reference = z1$reference - z0$reference)
  }))
  out <- list(arm = arm, ate = ate, diagnostics = diagnostics,
    parameters = list(delta = delta, delta_R_inf = delta_R_inf,
      delta_A_inf = delta_A_inf, method = method), call = match.call())
  class(out) <- c("marbounds_linf", "list")
  out
}

l2_linf_arm_lp <- function(r, B, b, delta_R_inf, delta_A_inf, max_n) {
  if (!requireNamespace("lpSolve", quietly = TRUE))
    stop("Sharp L-infinity bounds require the suggested package lpSolve")
  all_use <- which(r$obs)
  use <- all_use
  if (length(use) > max_n)
    use <- use[unique(round(seq(1, length(use), length.out = max_n)))]
  if (!length(use)) stop("No observed outcomes in arm")
  ww <- r$w[use] / length(r$w)
  if (length(use) < length(all_use)) {
    target_mass <- sum(r$w[all_use]) / length(r$w)
    ww <- ww * target_mass / sum(ww)
  }
  BB <- B[use, , drop = FALSE]
  y <- r$y[use]; e <- r$e[use]; q <- r$q[use]; bb <- b[use]
  m <- length(use); Z <- matrix(0, m, m); I <- diag(m)
  ## R is nonnegative in addition to satisfying |R-1| <= delta_R_inf.
  ## Therefore M cannot fall below the mixture envelope 1-b even when the
  ## symmetric radius is greater than one.
  mlo <- pmax(1 - bb, 1 - bb * delta_R_inf)
  mup <- 1 + bb * delta_R_inf
  klo <- max(0, 1 - delta_A_inf); kup <- 1 + delta_A_inf
  eqM <- cbind(t(BB * ww), matrix(0, ncol(BB), m))
  eqL <- cbind(matrix(0, ncol(BB), m), t(BB * ww))
  rhsM <- colSums(BB * ww)
  Acon <- rbind(eqM, eqL, cbind(I, Z), cbind(I, Z),
                cbind(-kup * I, I), cbind(-klo * I, I))
  dirs <- c(rep("=", 2 * ncol(BB)), rep("<=", m), rep(">=", m),
            rep("<=", m), rep(">=", m))
  rhs <- c(rhsM, rhsM, mup, mlo, rep(0, 2 * m))
  objective <- c(ww * y * e, ww * y * q)
  lo <- lpSolve::lp("min", objective, Acon, dirs, rhs)
  hi <- lpSolve::lp("max", objective, Acon, dirs, rhs)
  if (lo$status != 0 || hi$status != 0)
    stop("Sharp L-infinity LP failed (statuses ", lo$status, ", ",
         hi$status, ")")
  list(lower = lo$objval, upper = hi$objval,
       diagnostics = list(n_respondents = m, lower_status = lo$status,
         upper_status = hi$status, minimum_M_allowed = min(mlo),
         respondent_subsampled = length(use) < length(all_use),
         mixture_envelope_respected = all(mlo >= 1 - bb - 1e-12)))
}

l2_validate_inputs <- function(data, Y, A, C, X, delta, delta_R, delta_A, delta_M) {
  if (!is.data.frame(data)) stop("data must be a data.frame")
  cols <- c(Y, A, C, X)
  if (!all(cols %in% names(data))) stop("All Y, A, C, and X columns must exist")
  aa <- data[[A]]
  cc <- data[[C]]
  yy <- data[[Y]]
  if (anyNA(aa) || !all(aa %in% 0:1)) stop("A must contain only 0 and 1")
  if (anyNA(cc) || !all(cc %in% 0:1)) stop("C must contain only 0 and 1")
  if (anyNA(yy[cc == 0])) stop("Observed outcomes (C = 0) cannot be missing")
  if (!is.numeric(yy)) stop("Y must be numeric")
  list(y = yy, A = aa, C = cc)
}

l2_as_pair <- function(z, nm, upper = Inf) {
  if (length(z) == 1L) z <- rep(z, 2L)
  if (length(z) != 2L || any(!is.finite(z)) || any(z < 0) || any(z > upper))
    stop(nm, " must be one nonnegative value or a length-two vector")
  as.numeric(z)
}

l2_fit_predict <- function(formula, family, train, test, fallback) {
  fit <- try(stats::glm(formula, data = train, family = family), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(fallback, nrow(test)))
  pred <- try(stats::predict.glm(fit, newdata = test, type = "response"), silent = TRUE)
  if (inherits(pred, "try-error") || any(!is.finite(pred))) rep(fallback, nrow(test)) else as.numeric(pred)
}

l2_nuisance_glm <- function(data, Y, A, C, X, folds = 5L, seed = 1L,
                            outcome_family = stats::gaussian()) {
  n <- nrow(data)
  folds <- max(2L, min(as.integer(folds), n))
  set.seed(seed)
  fold <- sample(rep(seq_len(folds), length.out = n))
  rhs <- if (length(X)) paste(X, collapse = " + ") else "1"
  f_e <- stats::as.formula(paste(A, "~", rhs))
  e1 <- rho0 <- rho1 <- mu0 <- mu1 <- rep(NA_real_, n)
  for (v in seq_len(folds)) {
    tr <- data[fold != v, , drop = FALSE]
    te <- data[fold == v, , drop = FALSE]
    id <- which(fold == v)
    e1[id] <- l2_fit_predict(f_e, stats::binomial(), tr, te, mean(tr[[A]]))
    for (a in 0:1) {
      tra <- tr[tr[[A]] == a, , drop = FALSE]
      f_r <- stats::as.formula(paste("I(1 -", C, ") ~", rhs))
      fallback_r <- mean(1 - tra[[C]])
      if (!is.finite(fallback_r)) fallback_r <- mean(1 - tr[[C]])
      if (!is.finite(fallback_r)) fallback_r <- .5
      rr <- l2_fit_predict(f_r, stats::binomial(), tra, te,
                           fallback_r)
      tr_y <- tra[tra[[C]] == 0, , drop = FALSE]
      f_y <- stats::as.formula(paste(Y, "~", rhs))
      fallback_y <- mean(tr_y[[Y]], na.rm = TRUE)
      if (!is.finite(fallback_y)) fallback_y <- mean(tr[[Y]], na.rm = TRUE)
      if (!is.finite(fallback_y)) fallback_y <- 0
      yy <- l2_fit_predict(f_y, outcome_family, tr_y, te,
                           fallback_y)
      if (a == 0) { rho0[id] <- rr; mu0[id] <- yy }
      else { rho1[id] <- rr; mu1[id] <- yy }
    }
  }
  list(e0 = clip_probs(1 - e1, 1e-4), e1 = clip_probs(e1, 1e-4),
       rho0 = clip_probs(rho0, 1e-4), rho1 = clip_probs(rho1, 1e-4),
       mu0 = mu0, mu1 = mu1, fold_id = fold)
}

l2_basis <- function(data, X, basis) {
  if (is.null(basis)) {
    if (!length(X)) return(matrix(1, nrow(data), 1L, dimnames = list(NULL, "(Intercept)")))
    basis <- stats::as.formula(paste("~", paste(X, collapse = " + ")))
  }
  B <- stats::model.matrix(basis, data = data)
  for (j in seq_len(ncol(B))) {
    if (colnames(B)[j] != "(Intercept)") {
      s <- stats::sd(B[, j])
      if (is.finite(s) && s > 0) B[, j] <- (B[, j] - mean(B[, j])) / s
    }
  }
  B
}

l2_reference_quantities <- function(y, A, C, a, nuisance) {
  e <- nuisance[[paste0("e", a)]]
  rho <- nuisance[[paste0("rho", a)]]
  mu <- nuisance[[paste0("mu", a)]]
  obs <- A == a & C == 0
  y0 <- y
  y0[!obs] <- 0
  w <- as.numeric(obs) / (e * rho)
  z <- y0 - mu
  z[!obs] <- 0
  psi <- mean(mu + w * (y0 - mu))
  list(e = e, q = 1 - e, rho = rho, pi = 1 - rho,
       mu = mu, obs = obs, w = w, y = y0, z = z, psi = psi)
}

l2_arm_bounds <- function(y, A, C, a, nuisance, basis, delta, delta_R,
                          delta_A, delta_M, model, method, control) {
  r <- l2_reference_quantities(y, A, C, a, nuisance)
  b <- delta * r$pi
  VR <- mean(r$w * b^2 * r$z^2)
  VA <- mean(r$w * r$q^2 * r$z^2)
  WRA <- mean(r$w * r$q^4 * b^2 * r$z^4)
  VM <- mean(r$w * r$z^2)
  WMA <- mean(r$w * r$q^4 * r$z^4)
  radius <- if (model == "separated") {
    delta_R * sqrt(max(VR, 0)) +
      delta_A * sqrt(max(VA + delta_R * sqrt(max(WRA, 0)), 0))
  } else {
    delta_M * sqrt(max(VM, 0)) +
      delta_A * sqrt(max(VA + delta_M * sqrt(max(WMA, 0)), 0))
  }
  rows <- list()
  diag <- list(moments = c(VR = VR, VA = VA, WRA = WRA, VM = VM, WMA = WMA))
  if (method %in% c("cs", "both")) {
    rows[[length(rows) + 1L]] <- data.frame(
      arm = a, method = "cs", lower = r$psi - radius,
      upper = r$psi + radius, reference = r$psi
    )
  }
  if (method %in% c("sharp", "both")) {
    lo <- l2_dual_arm(r, basis, b, delta_R, delta_A, delta_M,
                      model, sigma = -1, control = control)
    hi <- l2_dual_arm(r, basis, b, delta_R, delta_A, delta_M,
                      model, sigma = 1, control = control)
    rows[[length(rows) + 1L]] <- data.frame(
      arm = a, method = "sharp_sieve", lower = -lo$value,
      upper = hi$value, reference = r$psi
    )
    diag$sharp <- list(lower = lo, upper = hi)
  }
  list(summary = do.call(rbind, rows), diagnostics = diag, psi = r$psi)
}

l2_dual_arm <- function(r, B, b, delta_R, delta_A, delta_M,
                        model, sigma, control) {
  use <- which(r$obs)
  if (!length(use)) stop("No observed outcomes in arm")
  y <- r$y[use]; e <- r$e[use]; q <- r$q[use]
  rho <- r$rho[use]; bb <- b[use]; ww <- r$w[use] / length(r$w)
  BB <- B[use, , drop = FALSE]
  p <- ncol(BB)
  eps <- 1e-8
  radR <- if (model == "separated") delta_R else delta_M
  fixed_M <- radR == 0 || (model == "separated" && all(bb <= eps))
  if (fixed_M) {
    ## A zero missingness budget (or a zero prevalence envelope everywhere)
    ## identifies M=1 exactly.  Handling this as an equality avoids asking the
    ## numerical optimizer to approximate the equality by sending lambda_R to
    ## infinity.
    if (delta_A == 0) {
      primal <- sum(ww * sigma * y)
      return(list(value = primal, convergence = 0L,
        message = "Both divergence components are fixed at their null values",
        lambda_R = 0, lambda_A = 0, divergence_R = 0,
        divergence_A = 0, normalization_M = rep(0, p),
        normalization_L = rep(0, p),
        coefficients = c(Inf, Inf, rep(0, 2 * p)),
        primal_value = primal, primal_dual_gap = 0,
        maximum_normalization_error = 0, maximum_budget_violation = 0,
        fixed_missingness_tilt = TRUE, no_confounding = TRUE))
    }
    objectiveA <- function(par, details = FALSE) {
      lamA <- exp(par[1])
      beta <- drop(BB %*% par[1 + seq_len(p)])
      L <- pmax(0, 1 + (sigma * y * q + beta) / (2 * lamA))
      divA <- (L - 1)^2
      H <- sigma * y * (e + q * L) + beta * (L - 1) - lamA * divA
      val <- lamA * delta_A^2 + sum(ww * H)
      if (!details) return(val)
      normL <- colSums(BB * (ww * (L - 1)))
      primal <- sum(ww * sigma * y * (e + q * L))
      list(value = val, L = L, divA = sum(ww * divA),
        normL = normL, lamA = lamA, primal = primal)
    }
    initA <- c(log(max(delta_A, .1)), rep(0, p))
    ctlA <- utils::modifyList(list(maxit = 1000, reltol = 1e-9), control)
    fitA <- stats::optim(initA, objectiveA, method = "BFGS", control = ctlA)
    detA <- objectiveA(fitA$par, details = TRUE)
    max_norm <- max(abs(detA$normL), 0)
    max_budget <- max(0, detA$divA - delta_A^2)
    return(list(value = detA$value, convergence = fitA$convergence,
      message = fitA$message, lambda_R = 0, lambda_A = detA$lamA,
      divergence_R = 0, divergence_A = detA$divA,
      normalization_M = rep(0, p), normalization_L = detA$normL,
      coefficients = c(Inf, fitA$par[1], rep(0, p), fitA$par[-1]),
      primal_value = detA$primal,
      primal_dual_gap = detA$value - detA$primal,
      maximum_normalization_error = max_norm,
      maximum_budget_violation = max_budget,
      fixed_missingness_tilt = TRUE))
  }
  if (delta_A == 0) {
    objective0 <- function(par, details = FALSE) {
      lamR <- exp(par[1])
      alpha <- drop(BB %*% par[1 + seq_len(p)])
      if (model == "separated") {
        M0 <- 1 + bb^2 / (2 * lamR) * (sigma * y + alpha)
        M <- pmax(1 - bb, M0)
        divR <- (M - 1)^2 / pmax(bb^2, eps)
        radR <- delta_R
      } else {
        M0 <- 1 + (sigma * y + alpha) / (2 * lamR)
        M <- pmax(rho, M0)
        divR <- (M - 1)^2
        radR <- delta_M
      }
      H <- sigma * y * M + alpha * (M - 1) - lamR * divR
      val <- lamR * radR^2 + sum(ww * H)
      if (!details) return(val)
       normM <- colSums(BB * (ww * (M - 1)))
       primal <- sum(ww * sigma * y * M)
       list(value = val, M = M, divR = sum(ww * divR),
            normM = normM, lamR = lamR, primal = primal)
    }
    init0 <- c(log(max(delta_R, delta_M, .1)), rep(0, p))
    ctl0 <- utils::modifyList(list(maxit = 1000, reltol = 1e-9), control)
    fit0 <- stats::optim(init0, objective0, method = "BFGS", control = ctl0)
    det0 <- objective0(fit0$par, details = TRUE)
    full_par <- c(fit0$par[1], Inf, fit0$par[-1], rep(0, p))
    max_norm <- max(abs(det0$normM), 0)
    max_budget <- max(0, det0$divR - radR^2)
    return(list(value = det0$value, convergence = fit0$convergence,
      message = fit0$message, lambda_R = det0$lamR, lambda_A = Inf,
      divergence_R = det0$divR, divergence_A = 0,
      normalization_M = det0$normM, normalization_L = det0$normM,
      coefficients = full_par, no_confounding = TRUE,
      primal_value = det0$primal,
      primal_dual_gap = det0$value - det0$primal,
      maximum_normalization_error = max_norm,
      maximum_budget_violation = max_budget))
  }
  objective <- function(par, details = FALSE) {
    lamR <- exp(par[1]); lamA <- exp(par[2])
    alpha <- drop(BB %*% par[2 + seq_len(p)])
    beta <- drop(BB %*% par[2 + p + seq_len(p)])
    K <- pmax(0, 1 + (sigma * y * q + beta) / (2 * lamA))
    if (model == "separated") {
      M0 <- 1 + bb^2 / (2 * lamR) *
        (sigma * y * e + alpha - lamA * (1 - K^2))
      M <- pmax(1 - bb, M0)
      divR <- (M - 1)^2 / pmax(bb^2, eps)
      radR <- delta_R
    } else {
      M0 <- 1 + 1 / (2 * lamR) *
        (sigma * y * e + alpha - lamA * (1 - K^2))
      M <- pmax(rho, M0)
      divR <- (M - 1)^2
      radR <- delta_M
    }
    L <- M * K
    divA <- (L - M)^2 / pmax(M, eps)
    H <- sigma * y * (e * M + q * L) + alpha * (M - 1) +
      beta * (L - 1) - lamR * divR - lamA * divA
    val <- lamR * radR^2 + lamA * delta_A^2 + sum(ww * H)
    if (!details) return(val)
    normM <- colSums(BB * (ww * (M - 1)))
    normL <- colSums(BB * (ww * (L - 1)))
    primal <- sum(ww * sigma * y * (e * M + q * L))
    list(value = val, M = M, L = L, K = K, divR = sum(ww * divR),
         divA = sum(ww * divA), normM = normM,
         normL = normL, lamR = lamR, lamA = lamA, primal = primal)
  }
  init <- c(log(max(delta_R, delta_M, 0.1)), log(max(delta_A, 0.1)), rep(0, 2 * p))
  ctl <- utils::modifyList(list(maxit = 1000, reltol = 1e-9), control)
  fit <- stats::optim(init, objective, method = "BFGS", control = ctl)
  det <- objective(fit$par, details = TRUE)
  max_norm <- max(abs(c(det$normM, det$normL)), 0)
  max_budget <- max(0, det$divR - radR^2,
                    det$divA - delta_A^2)
  list(value = det$value, convergence = fit$convergence, message = fit$message,
       lambda_R = det$lamR, lambda_A = det$lamA,
       divergence_R = det$divR, divergence_A = det$divA,
       normalization_M = det$normM, normalization_L = det$normL,
       coefficients = fit$par, primal_value = det$primal,
       primal_dual_gap = det$value - det$primal,
       maximum_normalization_error = max_norm,
       maximum_budget_violation = max_budget)
}

l2_moment_check <- function(y) {
  y <- y[is.finite(y)]
  if (!length(y)) return(list(finite_fourth_empirical_moment = FALSE, fourth_moment = NA_real_))
  m4 <- mean((y - mean(y))^4)
  list(finite_fourth_empirical_moment = is.finite(m4), fourth_moment = m4,
       note = "A finite sample cannot verify a population moment assumption; inspect tail behavior substantively.")
}

#' Simulate a continuous-outcome mixed-missingness dataset
#'
#' Errors are Gaussian by default, so all polynomial moment restrictions in the
#' paper hold. A standardized Student t distribution with more than eight
#' degrees of freedom is available for finite-eighth-moment tail stress tests.
#' Informative missingness is generated by a latent type whose probability is
#' bounded by the total missingness probability and whose missingness depends on
#' the potential outcome.
#'
#' @param n Sample size.
#' @param seed Random seed.
#' @param violation Which identifying assumptions the DGP violates: `"mar"`,
#'   `"nuc"`, `"both"`, or `"none"`.
#' @param density_ratio_link `"bounded"` uses tanh latent effects so oracle
#'   density-ratio L-infinity radii are finite; `"linear"` retains unbounded
#'   logistic tilts and therefore infinite L-infinity radii.
#' @param confounding_strength Coefficient linking treatment to the shared
#'   latent potential-outcome disturbance. Zero recovers conditional
#'   exchangeability; the default generates nonzero oracle confounding radii.
#' @param informative_strength Coefficient linking informative missingness to
#'   the shared outcome disturbance. Its sign controls whether the missingness
#'   and treatment tilts move with or against one another.
#' @param informative_intercept Baseline log odds of informative missingness.
#' @param noninformative_intercept Baseline log odds of noninformative
#'   missingness.
#' @param treatment_scale Multiplier on the covariate part of the treatment
#'   linear predictor. Values above one create poorer overlap.
#' @param covariate_complexity Either `"linear"` or `"nonlinear"`; the latter
#'   adds quadratic and sinusoidal terms omitted by a main-effects GLM.
#' @param informative_concentration Either `"diffuse"` or `"high_variance"`.
#'   The latter applies the outcome-dependent missingness tilt only when
#'   `abs(X1) > 1`, concentrating its L2 budget in a smaller covariate region.
#' @param outcome_distribution `"normal"` or standardized `"t"`.
#' @param outcome_df Degrees of freedom for the t distribution. It must exceed
#'   eight so the sufficient moment conditions for the endpoint EIFs hold.
#' @param outcome_scale Marginal standard deviation of the outcome disturbance.
#' @return A data frame with observed data, complete potential outcomes and
#'   potential latent indicators, and true conditional DGP quantities.
#' @export
simulate_l2_data <- function(n = 1000L, seed = 1L,
                             violation = c("both", "mar", "nuc", "none"),
                             density_ratio_link = c("bounded", "linear"),
                             confounding_strength = 0.6,
                             informative_strength = 0.45,
                             informative_intercept = -2.4,
                             noninformative_intercept = -1.8,
                             treatment_scale = 1,
                             covariate_complexity = c("linear", "nonlinear"),
                             informative_concentration = c("diffuse", "high_variance"),
                             outcome_distribution = c("normal", "t"),
                             outcome_df = 12,
                             outcome_scale = 1) {
  violation <- match.arg(violation)
  density_ratio_link <- match.arg(density_ratio_link)
  covariate_complexity <- match.arg(covariate_complexity)
  informative_concentration <- match.arg(informative_concentration)
  outcome_distribution <- match.arg(outcome_distribution)
  scalar_finite <- function(x, nm) if (length(x) != 1L || !is.finite(x))
    stop(nm, " must be one finite number")
  scalar_finite(confounding_strength, "confounding_strength")
  scalar_finite(informative_strength, "informative_strength")
  scalar_finite(informative_intercept, "informative_intercept")
  scalar_finite(noninformative_intercept, "noninformative_intercept")
  scalar_finite(treatment_scale, "treatment_scale")
  scalar_finite(outcome_scale, "outcome_scale")
  if (treatment_scale <= 0) stop("treatment_scale must be positive")
  if (outcome_scale <= 0) stop("outcome_scale must be positive")
  if (outcome_distribution == "t" &&
      (length(outcome_df) != 1L || !is.finite(outcome_df) || outcome_df <= 8))
    stop("outcome_df must exceed eight for the t outcome distribution")
  violate_mar <- violation %in% c("mar", "both")
  violate_nuc <- violation %in% c("nuc", "both")
  gamma_A <- if (violate_nuc) confounding_strength else 0
  latent <- function(z) if (density_ratio_link == "bounded") tanh(z) else z
  set.seed(seed)
  X1 <- stats::rnorm(n)
  X2 <- stats::runif(n, -1, 1)
  nonlinear_A <- if (covariate_complexity == "nonlinear")
    0.3 * (X1^2 - 1) + 0.25 * sin(pi * X2) else 0
  lp_A <- 0.2 + treatment_scale * (0.5 * X1 - 0.35 * X2 + nonlinear_A)
  mu_Y0 <- 0.3 + 0.6 * X1 - 0.25 * X2
  if (covariate_complexity == "nonlinear")
    mu_Y0 <- mu_Y0 + 0.35 * (X1^2 - 1) + 0.25 * sin(pi * X2)
  mu_Y1 <- mu_Y0 + 0.7 + 0.2 * X1
  draw_error <- function(nn) if (outcome_distribution == "normal")
    stats::rnorm(nn, sd = outcome_scale) else
    stats::rt(nn, df = outcome_df) * sqrt((outcome_df - 2) / outcome_df) * outcome_scale
  error_quantiles <- function(nn) {
    pp <- (seq_len(nn) - 0.5) / nn
    if (outcome_distribution == "normal") stats::qnorm(pp) * outcome_scale else
      stats::qt(pp, df = outcome_df) * sqrt((outcome_df - 2) / outcome_df) * outcome_scale
  }
  eps0 <- draw_error(n)
  eps1 <- eps0
  Y0 <- mu_Y0 + eps0
  Y1 <- mu_Y1 + eps1
  qA <- error_quantiles(200L)
  e1 <- rowMeans(stats::plogis(outer(lp_A, gamma_A * latent(qA), "+")))
  p_A1_given_YX <- stats::plogis(lp_A + gamma_A * latent(eps0))
  A <- stats::rbinom(n, 1, p_A1_given_YX)
  Yfull <- ifelse(A == 1, Y1, Y0)
  informative_region <- if (informative_concentration == "high_variance")
    as.numeric(abs(X1) > 1) else rep(1, n)
  informative_coef <- informative_strength * informative_region
  p_UI0 <- if (violate_mar) stats::plogis(informative_intercept + 0.25 * X1 +
    informative_coef * latent(eps0)) else
    rep(0, n)
  p_UI1 <- if (violate_mar) stats::plogis(informative_intercept + 0.35 +
    0.25 * X1 + informative_coef * latent(eps1)) else
    rep(0, n)
  p_N0 <- stats::plogis(noninformative_intercept - 0.2 * X2)
  p_N1 <- stats::plogis(noninformative_intercept + 0.25 - 0.2 * X2)
  p_C0 <- p_UI0 + (1 - p_UI0) * p_N0
  p_C1 <- p_UI1 + (1 - p_UI1) * p_N1
  u_ui <- stats::runif(n); u_n <- stats::runif(n)
  UI0 <- as.integer(u_ui < p_UI0); UI1 <- as.integer(u_ui < p_UI1)
  C0 <- as.integer(UI0 == 1 | u_n < p_N0)
  C1 <- as.integer(UI1 == 1 | u_n < p_N1)
  U_I <- ifelse(A == 1, UI1, UI0)
  C <- ifelse(A == 1, C1, C0)
  ## Deterministic quantile quadrature gives P(C(a)=1|X) and E(Y(a)|X,C(a)=0).
  q <- error_quantiles(200L)
  integrate_arm <- function(mu, a, pN) {
    eta <- outer(mu, q, "+")
    pui <- if (violate_mar)
      stats::plogis(outer(informative_intercept + 0.35 * a + 0.25 * X1,
                          rep(1, length(q))) +
        outer(informative_coef, latent(q))) else
      matrix(0, nrow = n, ncol = length(q))
    pA1 <- stats::plogis(outer(lp_A, gamma_A * latent(q), "+"))
    ga <- if (a == 1) pA1 else 1 - pA1
    ea <- if (a == 1) e1 else 1 - e1
    po <- (1 - pui) * (1 - pN)
    joint_o <- rowMeans(ga * po)
    rho <- joint_o / ea
    list(rho = rho, muO = rowMeans(eta * ga * po) / joint_o,
         piI = rowMeans(ga * pui) / ea)
  }
  o0 <- integrate_arm(mu_Y0, 0, p_N0)
  o1 <- integrate_arm(mu_Y1, 1, p_N1)
  Y <- Yfull
  Y[C == 1] <- NA_real_
  data.frame(Y = Y, A = A, C = C, X1 = X1, X2 = X2,
             Y0 = Y0, Y1 = Y1, Yfull = Yfull,
             mu_Y0 = mu_Y0, mu_Y1 = mu_Y1,
             sd_Y0 = outcome_scale, sd_Y1 = outcome_scale,
             f_Y0_given_X = if (outcome_distribution == "normal")
               stats::dnorm(Y0, mu_Y0, outcome_scale) else
               stats::dt((Y0 - mu_Y0) /
                 (sqrt((outcome_df - 2) / outcome_df) * outcome_scale),
                 df = outcome_df) /
                 (sqrt((outcome_df - 2) / outcome_df) * outcome_scale),
             f_Y1_given_X = if (outcome_distribution == "normal")
               stats::dnorm(Y1, mu_Y1, outcome_scale) else
               stats::dt((Y1 - mu_Y1) /
                 (sqrt((outcome_df - 2) / outcome_df) * outcome_scale),
                 df = outcome_df) /
                 (sqrt((outcome_df - 2) / outcome_df) * outcome_scale),
             p_A1_given_X = e1,
             p_A1_given_YX = p_A1_given_YX,
             treatment_lp_X = lp_A,
             confounding_strength = gamma_A,
             informative_strength = informative_strength,
             informative_intercept = informative_intercept,
             informative_outcome_coef = informative_coef,
             noninformative_intercept = noninformative_intercept,
             treatment_scale = treatment_scale,
             covariate_complexity = covariate_complexity,
             informative_concentration = informative_concentration,
             outcome_distribution = outcome_distribution,
             outcome_df = outcome_df,
             outcome_scale = outcome_scale,
             density_ratio_link = density_ratio_link,
             violation = violation,
             informative_missingness = as.integer(violate_mar),
             UI0 = UI0, UI1 = UI1, U_I = U_I,
             p_UI0_given_YX = p_UI0, p_UI1_given_YX = p_UI1,
             C0 = C0, C1 = C1,
             p_C0_given_YX = p_C0, p_C1_given_YX = p_C1,
             p_C0_given_X = 1 - o0$rho, p_C1_given_X = 1 - o1$rho,
             p_UI0_given_AX = o0$piI, p_UI1_given_AX = o1$piI,
             mu_Y0_given_XC0 = o0$muO, mu_Y1_given_XC0 = o1$muO,
             true_ate = 0.7)
}

l2_dgp_scalar <- function(data, name, default) {
  if (name %in% names(data)) data[[name]][1L] else default
}

l2_dgp_error_quantiles <- function(data, nodes) {
  pp <- (seq_len(nodes) - 0.5) / nodes
  distribution <- as.character(l2_dgp_scalar(data, "outcome_distribution", "normal"))
  scale <- as.numeric(l2_dgp_scalar(data, "outcome_scale", 1))
  if (distribution == "normal") return(stats::qnorm(pp) * scale)
  df <- as.numeric(l2_dgp_scalar(data, "outcome_df", 12))
  stats::qt(pp, df = df) * sqrt((df - 2) / df) * scale
}

l2_oracle_cs <- function(data, delta, delta_R, delta_A, delta_M, model) {
  needed <- c("Y0", "Y1", "p_A1_given_X", "p_C0_given_YX",
              "p_C1_given_YX", "p_C0_given_X", "p_C1_given_X",
              "mu_Y0_given_XC0", "mu_Y1_given_XC0")
  if (!all(needed %in% names(data)))
    stop("return_true_bounds = TRUE requires DGP columns from simulate_l2_data()")
  ans <- vector("list", 2L)
  for (a in 0:1) {
    ya <- data[[paste0("Y", a)]]
    e <- if (a == 1) data$p_A1_given_X else 1 - data$p_A1_given_X
    q <- 1 - e
    pc_y <- data[[paste0("p_C", a, "_given_YX")]]
    rho <- 1 - data[[paste0("p_C", a, "_given_X")]]
    mu <- data[[paste0("mu_Y", a, "_given_XC0")]]
    gF <- if (a == 1) data$p_A1_given_YX else 1 - data$p_A1_given_YX
    wO <- (gF / e) * (1 - pc_y) / rho
    z <- ya - mu
    b <- delta[a + 1L] * (1 - rho)
    VR <- mean(wO * b^2 * z^2)
    VA <- mean(wO * q^2 * z^2)
    WRA <- mean(wO * q^4 * b^2 * z^4)
    VM <- mean(wO * z^2)
    WMA <- mean(wO * q^4 * z^4)
    rad <- if (model == "separated")
      delta_R[a + 1L] * sqrt(VR) +
        delta_A[a + 1L] * sqrt(VA + delta_R[a + 1L] * sqrt(WRA))
    else delta_M[a + 1L] * sqrt(VM) +
      delta_A[a + 1L] * sqrt(VA + delta_M[a + 1L] * sqrt(WMA))
    psi <- mean(mu)
    ans[[a + 1L]] <- data.frame(arm = a, method = "oracle_cs_monte_carlo",
      lower = psi - rad, upper = psi + rad, reference = psi)
  }
  arm <- do.call(rbind, ans)
  ate <- data.frame(method = "oracle_cs_monte_carlo",
    lower = arm$lower[arm$arm == 1] - arm$upper[arm$arm == 0],
    upper = arm$upper[arm$arm == 1] - arm$lower[arm$arm == 0],
    estimate_reference = arm$reference[arm$arm == 1] - arm$reference[arm$arm == 0])
  list(arm = arm, ate = ate)
}

#' Evaluate continuous-outcome bounds over a sensitivity grid
#'
#' @inheritParams l2_bounds
#' @param grid Data frame containing `delta_R` and `delta_A`, and optionally
#'   `delta`, `delta_M`. Scalar columns are applied to both arms; arm-specific
#'   columns may be supplied with suffixes `_0` and `_1`.
#' @return An object of class `marbounds_l2_grid` with one ATE row per method
#'   and grid point.
#' @export
l2_sensitivity_grid <- function(data, Y, A, C, X, grid,
                                delta = c(1, 1), model = c("net", "separated"),
                                method = c("both", "cs", "sharp"), folds = 5L,
                                nuisance_method = c("SuperLearner", "glm"),
                                sl_lib_prop = "SL.glm",
                                sl_lib_miss = "SL.glm",
                                sl_lib_outcome = "SL.glm",
                                basis = NULL, seed = 1L, control = list()) {
  model <- match.arg(model); method <- match.arg(method)
  nuisance_method <- match.arg(nuisance_method)
  if (!is.data.frame(grid) || !nrow(grid)) stop("grid must be a nonempty data frame")
  pair <- function(row, nm, default) {
    if (all(paste0(nm, c("_0", "_1")) %in% names(row)))
      return(as.numeric(row[1, paste0(nm, c("_0", "_1"))]))
    if (nm %in% names(row)) return(rep(as.numeric(row[[nm]]), 2L))
    default
  }
  ans <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    rr <- pair(grid[i, , drop = FALSE], "delta_R", c(0, 0))
    aa <- pair(grid[i, , drop = FALSE], "delta_K",
      pair(grid[i, , drop = FALSE], "delta_A", c(0, 0)))
    dd <- pair(grid[i, , drop = FALSE], "delta", delta)
    mm <- pair(grid[i, , drop = FALSE], "delta_M", rr)
    fit <- l2_bounds(data, Y, A, C, X, delta = dd, delta_R = rr,
                     delta_A = aa, delta_M = mm, model = model,
                     method = method, folds = folds, basis = basis,
                     nuisance_method = nuisance_method,
                     sl_lib_prop = sl_lib_prop, sl_lib_miss = sl_lib_miss,
                     sl_lib_outcome = sl_lib_outcome,
                     seed = seed, control = control)
    z <- fit$ate
    z$grid_id <- i
    ans[[i]] <- cbind(grid[rep(i, nrow(z)), , drop = FALSE], z)
  }
  out <- list(results = do.call(rbind, ans), call = match.call())
  class(out) <- "marbounds_l2_grid"
  out
}

#' Extract a tipping frontier from an L2 sensitivity grid
#'
#' @param x Output from `l2_sensitivity_grid()` or its results data frame.
#' @param threshold Null threshold, usually zero.
#' @param endpoint Endpoint used for tipping; defaults to `"lower"`.
#' @param plot Draw the frontier before returning it.
#' @return A data frame containing the first threshold crossing for each value
#'   of the missingness radius and method. Crossings are linearly interpolated between grid
#'   points; groups that do not tip are retained with `tipped_on_grid = FALSE`.
#' @export
l2_tipping_frontier <- function(x, threshold = 0,
                                endpoint = c("lower", "upper"), plot = FALSE) {
  endpoint <- match.arg(endpoint)
  d <- if (inherits(x, "marbounds_l2_grid")) x$results else x
  radius_x <- if ("delta_M" %in% names(d)) "delta_M" else "delta_R"
  radius_y <- if ("delta_K" %in% names(d)) "delta_K" else "delta_A"
  if (!is.data.frame(d) || !all(c(radius_x, endpoint, "method") %in% names(d)))
    stop("x must contain a missingness radius, method, and the selected endpoint")
  if (!radius_y %in% names(d))
    stop("x must contain delta_K (or legacy delta_A) to trace a frontier")
  split_d <- split(d, interaction(d[[radius_x]], d$method, drop = TRUE))
  out <- lapply(split_d, function(z) {
    z <- z[order(z[[radius_y]]), , drop = FALSE]
    hit <- if (endpoint == "lower") z[[endpoint]] <= threshold else
      z[[endpoint]] >= threshold
    hit[!is.finite(hit)] <- FALSE
    if (!any(hit)) {
      ans <- z[nrow(z), , drop = FALSE]
      ans[[radius_y]] <- NA_real_; ans[[endpoint]] <- NA_real_
      ans$tipped_on_grid <- FALSE; ans$interpolated <- FALSE
      return(ans)
    }
    j <- which(hit)[1L]; ans <- z[j, , drop = FALSE]
    ans$tipped_on_grid <- TRUE; ans$interpolated <- FALSE
    if (j > 1L && is.finite(z[[endpoint]][j - 1L]) &&
        z[[endpoint]][j] != z[[endpoint]][j - 1L]) {
      frac <- (threshold - z[[endpoint]][j - 1L]) /
        (z[[endpoint]][j] - z[[endpoint]][j - 1L])
      if (is.finite(frac) && frac >= 0 && frac <= 1) {
        ans[[radius_y]] <- z[[radius_y]][j - 1L] + frac *
          (z[[radius_y]][j] - z[[radius_y]][j - 1L])
        ans[[endpoint]] <- threshold; ans$interpolated <- TRUE
      }
    }
    ans
  })
  out <- do.call(rbind, out)
  class(out) <- c("marbounds_l2_frontier", "data.frame")
  attr(out, "endpoint") <- endpoint
  attr(out, "threshold") <- threshold
  attr(out, "radius_x") <- radius_x
  attr(out, "radius_y") <- radius_y
  if (isTRUE(plot)) graphics::plot(out)
  out
}

#' Plot an L2 tipping frontier
#' @param x Output from `l2_tipping_frontier()`.
#' @param ... Additional graphical arguments passed to `plot()`.
#' @return The frontier, invisibly.
#' @export
plot.marbounds_l2_frontier <- function(x, ...) {
  radius_x <- attr(x, "radius_x"); radius_y <- attr(x, "radius_y")
  if (is.null(radius_x)) radius_x <- if ("delta_M" %in% names(x)) "delta_M" else "delta_R"
  if (is.null(radius_y)) radius_y <- if ("delta_K" %in% names(x)) "delta_K" else "delta_A"
  if (!all(c(radius_x, radius_y) %in% names(x)))
    stop("Frontier must contain its two sensitivity-radius columns")
  x <- x[is.finite(x[[radius_x]]) & is.finite(x[[radius_y]]), , drop = FALSE]
  if (!nrow(x)) stop("The threshold was not reached on the evaluated grid")
  methods <- unique(x$method)
  cols <- seq_along(methods)
  graphics::plot(range(x[[radius_x]], finite = TRUE),
                 range(x[[radius_y]], finite = TRUE), type = "n",
                 xlab = radius_x, ylab = radius_y, ...)
  for (j in seq_along(methods)) {
    z <- x[x$method == methods[j], , drop = FALSE]
    z <- z[order(z[[radius_x]]), , drop = FALSE]
    graphics::lines(z[[radius_x]], z[[radius_y]], type = "b", col = cols[j], pch = 19)
  }
  graphics::legend("topright", legend = methods, col = cols, lty = 1, pch = 19,
                   bty = "n")
  invisible(x)
}

#' Invert a sensitivity grid to obtain a joint calibration frontier
#'
#' A grid point is calibration-compatible when its bound contains the observed
#' benchmark target. For each value of the first radius, this function returns
#' the smallest compatible value of the second radius. Thus the result describes
#' jointly sufficient sensitivity budgets, not an estimate of their true
#' decomposition.
#'
#' @param x Output from `l2_sensitivity_grid()` or its results data frame.
#' @param benchmark Scalar observed benchmark, interpreted according to `type`.
#' @param type Use `"target"` when `benchmark` is on the estimand scale, or
#'   `"discrepancy"` when it is a change from `reference`.
#' @param reference Reference value for a discrepancy. By default, uses the
#'   grid's `estimate_reference` column.
#' @param radius_x,radius_y Names of the two sensitivity-radius columns. The
#'   first defaults to `delta_M` when available and otherwise `delta_R`; the
#'   second defaults to `delta_K` when available and otherwise legacy
#'   `delta_A`.
#' @return A data frame giving the minimum compatible `radius_y` at each
#'   `radius_x` and method. CS inversions are labeled as outer approximations.
#' @export
l2_calibration_frontier <- function(x, benchmark,
                                    type = c("target", "discrepancy"),
                                    reference = NULL, radius_x = NULL,
                                    radius_y = NULL) {
  type <- match.arg(type)
  d <- if (inherits(x, "marbounds_l2_grid")) x$results else x
  if (!is.data.frame(d) || !all(c("lower", "upper", "method") %in% names(d)))
    stop("x must be a sensitivity-grid result with lower, upper, and method columns")
  if (length(benchmark) != 1L || !is.finite(benchmark))
    stop("benchmark must be one finite number")
  if (is.null(radius_x)) {
    radius_x <- if ("delta_M" %in% names(d)) "delta_M" else "delta_R"
  }
  if (is.null(radius_y)) {
    radius_y <- if ("delta_K" %in% names(d)) "delta_K" else "delta_A"
  }
  if (!all(c(radius_x, radius_y) %in% names(d)))
    stop("radius_x and radius_y must name columns in the grid")
  if (is.null(reference)) {
    if (!"estimate_reference" %in% names(d))
      stop("Supply reference or include estimate_reference in the grid")
    refs <- unique(d$estimate_reference[is.finite(d$estimate_reference)])
    if (length(refs) != 1L)
      stop("The grid must have one reference value; otherwise supply reference")
    reference <- refs
  }
  if (length(reference) != 1L || !is.finite(reference))
    stop("reference must be one finite number")
  target <- if (type == "discrepancy") reference + benchmark else benchmark
  d$compatible <- d$lower <= target & target <= d$upper
  groups <- split(d, interaction(d[[radius_x]], d$method, drop = TRUE))
  ans <- lapply(groups, function(z) {
    ok <- z[z$compatible, , drop = FALSE]
    j <- if (nrow(ok)) which.min(ok[[radius_y]]) else NA_integer_
    data.frame(
      radius_x = z[[radius_x]][1L],
      minimum_radius_y = if (is.na(j)) NA_real_ else ok[[radius_y]][j],
      method = z$method[1L], compatible_on_grid = nrow(ok) > 0,
      benchmark = benchmark, reference = reference, target = target,
      approximation = if (z$method[1L] == "cs")
        "outer approximation from CS bounds" else
        "inversion of empirical sharp sieve bounds",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, ans)
  names(out)[names(out) == "radius_x"] <- radius_x
  names(out)[names(out) == "minimum_radius_y"] <- paste0("minimum_", radius_y)
  rownames(out) <- NULL
  out[order(out$method, out[[radius_x]]), , drop = FALSE]
}

#' Calibrate L2 sensitivity radii from benchmark tilts
#'
#' @param R,K,M Numeric vectors representing benchmark likelihood ratios.
#'   At least one must be supplied. `M` is the reduced-form, jointly
#'   calibratable distortion. If `K` is supplied, `M` also supplies its
#'   reference-to-factual weighting and defaults to one.
#' @param weights Optional nonnegative sampling weights.
#' @param benchmark Either `"user_supplied"` for elicited or externally
#'   validated tilts, or `"observed_data"` for a benchmark estimated from the
#'   same observed data used in the analysis.
#' @param conditional_on Identifying condition for an observed-data,
#'   mechanism-specific benchmark: `"exchangeability"` for `R` or `"MAR"`
#'   for `K`. It is unnecessary for a reduced-form `M` benchmark. Separate
#'   observed-data `R` and `K` tilts cannot be calibrated jointly by this
#'   function.
#' @return Named sensitivity radii corresponding to the supplied tilts, with
#'   attributes recording their calibration interpretation.
#' @export
l2_calibrate <- function(R = NULL, K = NULL, M = NULL, weights = NULL,
                         benchmark = c("user_supplied", "observed_data"),
                         conditional_on = NULL) {
  benchmark <- match.arg(benchmark)
  if (!is.null(conditional_on)) {
    conditional_on <- match.arg(conditional_on, c("exchangeability", "MAR"))
  }
  if (benchmark == "observed_data") {
    if (!is.null(R) && !is.null(K)) {
      stop("Separate R and K benchmarks are not jointly identified from the observed data; calibrate M, or analyze one mechanism conditionally on the other.")
    }
    if (!is.null(R) && !identical(conditional_on, "exchangeability")) {
      stop("An observed-data R benchmark requires conditional_on = 'exchangeability'.")
    }
    if (!is.null(K) && !identical(conditional_on, "MAR")) {
      stop("An observed-data K benchmark requires conditional_on = 'MAR'.")
    }
  }
  n <- max(length(R), length(K), length(M), length(weights), 0L)
  if (n == 0L) stop("Supply at least one of R, K, or M")
  recycle <- function(z, default = NULL) {
    if (is.null(z)) return(if (is.null(default)) NULL else rep(default, n))
    if (length(z) == 1L) z <- rep(z, n)
    if (length(z) != n || any(!is.finite(z))) stop("Tilt inputs must have a common finite length")
    as.numeric(z)
  }
  R <- recycle(R)
  K <- recycle(K)
  M_in <- recycle(M)
  w <- recycle(weights, 1)
  if (any(w < 0) || sum(w) <= 0) stop("weights must be nonnegative with positive sum")
  w <- w / sum(w)
  out <- numeric(0)
  if (!is.null(R)) {
    if (any(R < 0)) stop("R must be nonnegative")
    out["delta_R"] <- sqrt(sum(w * (R - 1)^2))
  }
  if (!is.null(M_in)) {
    if (any(M_in < 0)) stop("M must be nonnegative")
    out["delta_M"] <- sqrt(sum(w * (M_in - 1)^2))
  }
  if (!is.null(K)) {
    if (any(K < 0)) stop("K must be nonnegative")
    M_weight <- if (is.null(M_in)) rep(1, n) else M_in
    delta_K <- sqrt(sum(w * M_weight * (K - 1)^2))
    out["delta_K"] <- delta_K
    out["delta_A"] <- delta_K
  }
  attr(out, "benchmark") <- benchmark
  attr(out, "conditional_on") <- conditional_on
  attr(out, "interpretation") <- if (benchmark == "observed_data" &&
                                      is.null(R) && is.null(K)) {
    "reduced-form combined distortion"
  } else if (benchmark == "observed_data") {
    paste("mechanism-specific calibration conditional on", conditional_on)
  } else {
    "user-supplied or externally validated tilt"
  }
  out
}

#' Bootstrap continuous-outcome L2 bounds
#'
#' @inheritParams l2_bounds
#' @param B Number of nonparametric bootstrap samples.
#' @param conf_level Confidence level for percentile intervals.
#' @return A list containing the original fit, bootstrap replicates, standard
#'   errors, and percentile confidence intervals for each ATE endpoint.
#' @export
l2_bootstrap <- function(data, Y, A, C, X,
                         delta = c(1, 1), delta_R = c(0, 0),
                         delta_A = c(0, 0), delta_M = delta_R,
                         model = c("net", "separated"),
                         method = c("cs", "sharp"), folds = 5L,
                         nuisance_method = c("SuperLearner", "glm"),
                         sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm",
                         sl_lib_outcome = "SL.glm",
                         basis = NULL, seed = 1L, control = list(),
                         B = 200L, conf_level = 0.95,
                         delta_K = NULL) {
  model <- match.arg(model); method <- match.arg(method)
  nuisance_method <- match.arg(nuisance_method)
  if (!is.null(delta_K)) delta_A <- delta_K
  B <- as.integer(B)
  if (B < 2L) stop("B must be at least 2")
  if (!is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("conf_level must lie strictly between zero and one")
  original <- l2_bounds(data, Y, A, C, X, delta = delta, delta_R = delta_R,
                        delta_K = delta_A, delta_M = delta_M, model = model,
                        method = method, folds = folds, basis = basis,
                        nuisance_method = nuisance_method,
                        sl_lib_prop = sl_lib_prop, sl_lib_miss = sl_lib_miss,
                        sl_lib_outcome = sl_lib_outcome,
                        seed = seed, control = control)
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, B)
  reps <- matrix(NA_real_, B, 2L, dimnames = list(NULL, c("lower", "upper")))
  for (b in seq_len(B)) {
    id <- sample.int(nrow(data), nrow(data), replace = TRUE)
    fit <- try(l2_bounds(data[id, , drop = FALSE], Y, A, C, X,
                         delta = delta, delta_R = delta_R, delta_K = delta_A,
                          delta_M = delta_M, model = model, method = method,
                          folds = folds, nuisance_method = nuisance_method,
                          sl_lib_prop = sl_lib_prop, sl_lib_miss = sl_lib_miss,
                          sl_lib_outcome = sl_lib_outcome,
                          basis = basis, seed = seeds[b],
                         control = control), silent = TRUE)
    if (!inherits(fit, "try-error")) reps[b, ] <- unlist(fit$ate[1, c("lower", "upper")])
  }
  ok <- stats::complete.cases(reps)
  if (sum(ok) < 2L) stop("Fewer than two bootstrap fits succeeded")
  if (sum(ok) < max(2L, ceiling(B / 2)))
    warning("Fewer than half of bootstrap fits converged")
  alpha <- (1 - conf_level) / 2
  ci <- apply(reps[ok, , drop = FALSE], 2L, stats::quantile,
              probs = c(alpha, 1 - alpha), na.rm = TRUE)
  list(original = original, replicates = reps, successful = sum(ok),
       se = apply(reps, 2L, stats::sd, na.rm = TRUE), ci = ci,
       conf_level = conf_level, call = match.call())
}
