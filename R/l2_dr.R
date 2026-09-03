#' Cross-fitted plug-in and doubly robust sharp L2 endpoint estimators
#'
#' Fits the empirical net-distortion sieve dual by default (or the optional
#' separated prevalence--severity dual) using rotating three-stage
#' cross-fitting. Structural nuisances and the dual optimizer, pseudo-outcome
#' regressions, and final score evaluation use disjoint observations in every
#' rotation. The plug-in estimator uses inverse respondent weights.
#' The DR estimator augments the held-out dual loss by its respondent-law
#' conditional mean and adds the structural propensity and response-probability
#' chain-rule corrections described in the companion paper.
#'
#' @inheritParams l2_bounds
#' @param conf_level Wald confidence level.
#' @param folds Number of cross-fitting folds when nuisances are fitted. At
#'   least three folds are then used so the three stages remain disjoint.
#'   Controlled-error simulations directly construct nuisance estimates and
#'   ignore this argument.
#' @param boundary_tol Numerical tolerance for declaring the mixture lower
#'   envelope active.
#' @param nuisance_source Either `"estimated"`, which fits nuisance models with
#'   cross-fitting, or `"oracle_error"`, which directly constructs estimated
#'   nuisances by perturbing known DGP functions. The latter fits no regressions,
#'   uses no sample splitting, and is intended only for simulation.
#' @param nuisance_error_rate Exponent alpha in the simulated nuisance-error
#'   magnitude n^(-alpha).
#' @param nuisance_error_scale Scalar or named vector with entries `e`, `rho`,
#'   `mu`, and `regression` multiplying errors in the structural probabilities,
#'   outcome regression, and conditional dual-loss/derivative regressions.
#'   For backward compatibility, a three-entry vector omitting `mu` uses the
#'   `regression` magnitude for `mu`.
#' @param nuisance_error_sign Named vector giving independently configurable
#'   perturbation directions for `e`, `rho`, and `regression`.
#' @return A list with arm and ATE endpoint estimates, standard errors,
#'   confidence intervals, fold diagnostics, and row-level influence scores.
#' @export
l2_sharp_crossfit <- function(data, Y, A, C, X,
                              delta = c(1, 1), delta_R = c(.2, .2),
                              delta_A = c(.2, .2), delta_M = delta_R,
                              model = c("net", "separated"), folds = 5L,
                              basis = NULL,
                              nuisance_method = c("SuperLearner", "glm"),
                              sl_lib_prop = "SL.glm",
                              sl_lib_miss = "SL.glm",
                              sl_lib_outcome = "SL.glm",
                              seed = 1L, control = list(),
                              diagnostic_control = list(),
                              conf_level = .95, boundary_tol = 1e-7,
                              nuisance_source = c("estimated", "oracle_error"),
                              nuisance_error_rate = .25,
                              nuisance_error_scale = 1,
                              nuisance_error_sign = c(e = 1, rho = -1,
                                                      regression = 1),
                              delta_K = NULL) {
  nuisance_method <- match.arg(nuisance_method)
  nuisance_source <- match.arg(nuisance_source)
  model <- match.arg(model)
  diagnostic_control <- l2_diagnostic_control(diagnostic_control,
                                               compute_hessian = FALSE)
  if (length(conf_level) != 1L || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1)
    stop("conf_level must be a scalar strictly between zero and one")
  delta <- l2_as_pair(delta, "delta", 1)
  delta_R <- l2_as_pair(delta_R, "delta_R")
  if (!is.null(delta_K)) delta_A <- delta_K
  delta_A <- l2_as_pair(delta_A, "delta_K")
  delta_M <- l2_as_pair(delta_M, "delta_M")
  l2_validate_inputs(data, Y, A, C, X, delta, delta_R, delta_A, delta_M)
  n <- nrow(data)
  if (n < 3L) stop("l2_sharp_crossfit() requires at least three observations")
  controlled_error <- nuisance_source == "oracle_error"
  requested_folds <- as.integer(folds)
  folds <- if (controlled_error) 1L else
    max(3L, min(requested_folds, n))
  set.seed(seed)
  fold_id <- if (controlled_error) rep(1L, n) else
    sample(rep(seq_len(folds), length.out = n))
  zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)
  arm_ep <- vector("list", 4L); scores <- vector("list", 4L)
  diags <- list()
  has_oracle <- all(c("p_A1_given_X", "p_C0_given_X", "p_C1_given_X",
    "mu_Y0_given_XC0", "mu_Y1_given_XC0", "mu_Y0", "mu_Y1", "X1",
    "treatment_lp_X", "confounding_strength") %in% names(data))

  for (a in 0:1) for (sigma in c(-1, 1)) {
    endpoint <- if (sigma == -1) "lower" else "upper"
    score_plugin <- score_dr <- score_true <- rep(NA_real_, n)
    fold_constant <- fold_constant_true <- numeric(folds)
    fold_diag <- vector("list", folds)
    for (v in seq_len(folds)) {
      ## Rotate three disjoint roles.  Fold v evaluates the final score, the
      ## preceding fold fits the pseudo-outcome regressions, and all remaining
      ## folds fit the structural nuisances and dual optimizer.
      regression_fold <- if (controlled_error) 1L else
        if (v == 1L) folds else v - 1L
      if (controlled_error) {
        ## Controlled-error simulations directly construct every nuisance
        ## estimate from its DGP value plus a prespecified error. There are no
        ## learned regressions and hence no cross-fitting roles.
        te <- reg <- tr <- rep(TRUE, n)
      } else {
        te <- fold_id == v
        reg <- fold_id == regression_fold
        tr <- !(te | reg)
      }
      nonfit <- if (controlled_error) rep(TRUE, n) else !tr
      nonfit_id <- which(nonfit)
      reg_pos <- match(which(reg), nonfit_id)
      te_pos <- match(which(te), nonfit_id)
      split_nonfit <- function(z) list(
        train = z$train,
        regression = lapply(z$test, `[`, reg_pos),
        test = lapply(z$test, `[`, te_pos)
      )
      ns <- if (nuisance_source == "oracle_error")
        l2_split_nuisance_oracle_error(data, X, tr, nonfit,
          nuisance_error_rate, nuisance_error_scale, nuisance_error_sign)
      else l2_split_nuisance(data, Y, A, C, X, tr, nonfit,
        nuisance_method, sl_lib_prop, sl_lib_miss, sl_lib_outcome)
      ns <- split_nonfit(ns)
      bp0 <- l2_basis_split(data, X, basis, tr, nonfit)
      bp <- list(train = bp0$train,
        regression = bp0$test[reg_pos, , drop = FALSE],
        test = bp0$test[te_pos, , drop = FALSE])
      if (has_oracle) {
        ns0 <- l2_split_nuisance_oracle_error(data, X, tr, nonfit, 0,
          c(e = 0, rho = 0, regression = 0),
          c(e = 1, rho = 1, regression = 1))
        ns0 <- split_nonfit(ns0)
        r0 <- l2_reference_quantities(data[[Y]][tr], data[[A]][tr],
          data[[C]][tr], a, ns0$train)
        dual0 <- l2_dual_arm(r0, bp$train,
          delta[a + 1L] * r0$pi, delta_R[a + 1L], delta_A[a + 1L],
          delta_M[a + 1L], model, sigma, control,
          diagnostic_control = diagnostic_control)
      }
      dual_nuisance <- ns$train
      rtr <- l2_reference_quantities(data[[Y]][tr], data[[A]][tr],
                                    data[[C]][tr], a, dual_nuisance)
      btr <- delta[a + 1L] * rtr$pi
      dual <- l2_dual_arm(rtr, bp$train, btr, delta_R[a + 1L],
                          delta_A[a + 1L], delta_M[a + 1L], model,
                          sigma, control,
                          diagnostic_control = diagnostic_control)
      par <- dual$coefficients
      regression_eval <- l2_dual_evaluate(
        y = data[[Y]][reg], A = data[[A]][reg], C = data[[C]][reg], a = a,
         nuisance = ns$regression, B = bp$regression,
         delta = delta[a + 1L],
         par = par, sigma = sigma, boundary_tol = boundary_tol,
         model = model, delta_M = delta_M[a + 1L])
      test_eval <- l2_dual_evaluate(
        y = data[[Y]][te], A = data[[A]][te], C = data[[C]][te], a = a,
         nuisance = ns$test, B = bp$test, delta = delta[a + 1L],
         par = par, sigma = sigma, boundary_tol = boundary_tol,
         model = model, delta_M = delta_M[a + 1L])

      obs_reg <- data[[A]][reg] == a & data[[C]][reg] == 0
      if (nuisance_source == "oracle_error") {
        cm_err <- l2_oracle_dual_regressions(data[te, , drop = FALSE], a,
          bp$test, delta[a + 1L], par, sigma, boundary_tol,
          nuisance_override = ns$test, model = model,
          delta_M = delta_M[a + 1L])
        mH <- cm_err$mH; de <- cm_err$d_e; drho <- cm_err$d_rho
        esc <- l2_nuisance_error_scale(nuisance_error_scale)
        esign <- l2_nuisance_error_sign(nuisance_error_sign)
        emag <- sum(reg)^(-nuisance_error_rate) * esc[["regression"]] *
          esign[["regression"]]
        add_error <- function(pred, direction) {
          pred + emag * direction
        }
        mH <- add_error(mH, ns$test$err_H)
        de <- add_error(de, ns$test$err_de)
        drho <- add_error(drho, ns$test$err_drho)
      } else {
        mH <- l2_pseudo_regression(data[reg, , drop = FALSE], X,
          regression_eval$H[obs_reg], obs_reg, data[te, , drop = FALSE],
          nuisance_method, sl_lib_outcome)
        de <- l2_pseudo_regression(data[reg, , drop = FALSE], X,
          regression_eval$d_e[obs_reg], obs_reg, data[te, , drop = FALSE],
          nuisance_method, sl_lib_outcome)
        drho <- l2_pseudo_regression(data[reg, , drop = FALSE], X,
          regression_eval$d_rho[obs_reg], obs_reg, data[te, , drop = FALSE],
          nuisance_method, sl_lib_outcome)
      }
      e <- ns$test[[paste0("e", a)]]
      rho <- ns$test[[paste0("rho", a)]]
      S <- as.numeric(data[[A]][te] == a & data[[C]][te] == 0)
      Hobs <- test_eval$H; Hobs[!is.finite(Hobs)] <- 0
      plugin <- S / (e * rho) * Hobs
      dr <- mH + S / (e * rho) * (Hobs - mH) +
        de * (as.numeric(data[[A]][te] == a) - e) +
        drho * as.numeric(data[[A]][te] == a) / e *
          (as.numeric(data[[C]][te] == 0) - rho)
      score_plugin[te] <- plugin
      score_dr[te] <- dr
      if (has_oracle) {
        ev0 <- l2_dual_evaluate(data[[Y]][te], data[[A]][te], data[[C]][te],
          a, ns0$test, bp$test, delta[a + 1L], dual0$coefficients, sigma,
          boundary_tol, model = model, delta_M = delta_M[a + 1L])
        cm0 <- l2_oracle_dual_regressions(data[te, , drop = FALSE], a,
          bp$test, delta[a + 1L], dual0$coefficients, sigma, boundary_tol,
          model = model, delta_M = delta_M[a + 1L])
        e0 <- ns0$test[[paste0("e", a)]]
        rho0 <- ns0$test[[paste0("rho", a)]]
        H0 <- ev0$H; H0[!is.finite(H0)] <- 0
        score_true[te] <- cm0$mH + S / (e0 * rho0) * (H0 - cm0$mH) +
          cm0$d_e * (as.numeric(data[[A]][te] == a) - e0) +
          cm0$d_rho * as.numeric(data[[A]][te] == a) / e0 *
            (as.numeric(data[[C]][te] == 0) - rho0)
        radR <- if (model == "separated") delta_R[a + 1L] else delta_M[a + 1L]
        fold_constant_true[v] <- dual0$lambda_R * radR^2 +
          if (delta_A[a + 1L] == 0) 0 else
            dual0$lambda_A * delta_A[a + 1L]^2
      }
      radR <- if (model == "separated") delta_R[a + 1L] else delta_M[a + 1L]
      fold_constant[v] <- dual$lambda_R * radR^2 +
        if (delta_A[a + 1L] == 0) 0 else
          dual$lambda_A * delta_A[a + 1L]^2
      fold_diag[[v]] <- list(dual = dual,
        dual_nuisance = if (controlled_error)
          "controlled-error estimate" else "estimated",
        active_fraction = mean(test_eval$active[test_eval$obs]),
        arm = a, endpoint = endpoint,
        fold = v, regression_fold = regression_fold,
        fit_folds = paste(setdiff(seq_len(folds), c(v, regression_fold)),
                          collapse = ","))
    }
    constants <- fold_constant[fold_id]
    raw_plugin <- constants + score_plugin
    raw_dr <- constants + score_dr
    raw_true <- fold_constant_true[fold_id] + score_true
    sign_out <- if (sigma == -1) -1 else 1
    endpoint <- if (sigma == -1) "lower" else "upper"
    idx <- (a * 2L) + if (sigma == -1) 1L else 2L
    arm_ep[[idx]] <- rbind(
      data.frame(arm = a, endpoint = endpoint, estimator = "plugin",
                 estimate = sign_out * mean(raw_plugin),
                 se = stats::sd(sign_out * raw_plugin) / sqrt(n)),
      data.frame(arm = a, endpoint = endpoint, estimator = "dr",
                 estimate = sign_out * mean(raw_dr),
                 se = stats::sd(sign_out * raw_dr) / sqrt(n)),
      if (has_oracle) data.frame(arm = a, endpoint = endpoint,
        estimator = "true_eif", estimate = sign_out * mean(raw_true),
        se = stats::sd(sign_out * raw_true) / sqrt(n)))
    scores[[idx]] <- list(plugin = sign_out * (raw_plugin - mean(raw_plugin)),
                          dr = sign_out * (raw_dr - mean(raw_dr)))
    if (has_oracle) scores[[idx]]$true_eif <-
      sign_out * (raw_true - mean(raw_true))
    diags[[paste0("a", a, "_", endpoint)]] <- fold_diag
  }
  arm <- do.call(rbind, arm_ep); rownames(arm) <- NULL
  arm$conf_low <- arm$estimate - zcrit * arm$se
  arm$conf_high <- arm$estimate + zcrit * arm$se
  ate <- list(); ate_scores <- list(); j <- 0L
  estimators <- c("plugin", "dr", if (has_oracle) "true_eif")
  for (est in estimators) for (ep in c("lower", "upper")) {
    j <- j + 1L
    if (ep == "lower") { i1 <- 1L; i0 <- 2L } else { i1 <- 2L; i0 <- 1L }
    z1 <- arm[arm$arm == 1 & arm$endpoint == c("lower", "upper")[i1] & arm$estimator == est, ]
    z0 <- arm[arm$arm == 0 & arm$endpoint == c("lower", "upper")[i0] & arm$estimator == est, ]
    sc <- scores[[(1L * 2L) + i1]][[est]] - scores[[(0L * 2L) + i0]][[est]]
    val <- z1$estimate - z0$estimate; se <- stats::sd(sc) / sqrt(n)
    ate[[j]] <- data.frame(endpoint = ep, estimator = est, estimate = val,
                           se = se, conf_low = val - zcrit * se,
                           conf_high = val + zcrit * se)
    ate_scores[[paste(est, ep, sep = "_")]] <- sc
  }
  ate_df <- do.call(rbind, ate)
  outward <- do.call(rbind, lapply(unique(ate_df$estimator), function(est) {
    lo <- ate_df[ate_df$estimator == est & ate_df$endpoint == "lower", ]
    hi <- ate_df[ate_df$estimator == est & ate_df$endpoint == "upper", ]
    data.frame(estimator = est, lower = lo$conf_low, upper = hi$conf_high)
  }))
  diagnostic_table <- do.call(rbind, lapply(diags, function(endpoint_folds)
    do.call(rbind, lapply(endpoint_folds, function(z) {
      row <- l2_dual_diagnostic_row(z$dual, z$arm, z$endpoint)
      row$fold <- z$fold
      row$regression_fold <- z$regression_fold
      row$fit_folds <- z$fit_folds
      row$dual_nuisance <- z$dual_nuisance
      row$active_fraction <- z$active_fraction
      row
    }))))
  rownames(diagnostic_table) <- NULL
  diagnostic_summary <- l2_diagnostic_summary(diagnostic_table)
  out <- list(arm = arm, ate = ate_df, identified_set_ci = outward,
              scores = ate_scores,
              diagnostics = diags, fold_id = fold_id,
              crossfit_roles = if (controlled_error) data.frame(
                evaluation_fold = NA_integer_, regression_fold = NA_integer_,
                fit_folds = NA_character_,
                scheme = "not used for controlled nuisance errors",
                stringsAsFactors = FALSE) else data.frame(
                  evaluation_fold = seq_len(folds),
                  regression_fold = c(folds, seq_len(folds - 1L)),
                  fit_folds = vapply(seq_len(folds), function(v) {
                    rv <- if (v == 1L) folds else v - 1L
                    paste(setdiff(seq_len(folds), c(v, rv)), collapse = ",")
                  }, character(1)), scheme = "rotating three-stage",
                  stringsAsFactors = FALSE),
              diagnostic_table = diagnostic_table,
              diagnostic_summary = diagnostic_summary,
              diagnostic_note = l2_diagnostic_note(),
              parameters = list(delta = delta, delta_R = delta_R,
                                delta_A = delta_A, delta_K = delta_A,
                                delta_M = delta_M,
                                model = model, folds = folds,
                                requested_folds = requested_folds,
                                sample_splitting = if (controlled_error)
                                  "none: nuisances directly perturbed" else
                                  "rotating three-stage cross-fitting",
                                nuisance_source = nuisance_source,
                                nuisance_error_rate = nuisance_error_rate,
                                nuisance_error_scale = nuisance_error_scale,
                                nuisance_error_sign = nuisance_error_sign),
              inference_note = paste("DR intervals use the augmented endpoint",
                "score under the paper's regularity and product-rate conditions;",
                "the three fitting stages are separated by cyclic rotation;",
                "plug-in intervals condition on fitted nuisances."),
              call = match.call())
  class(out) <- "marbounds_l2_sharp_cf"
  if (isTRUE(diagnostic_control$warn))
    l2_warn_diagnostics(diagnostic_summary, "Cross-fitted sharp-bound")
  out
}

l2_basis_split <- function(data, X, basis, tr, te) {
  if (is.null(basis)) basis <- if (length(X))
    stats::as.formula(paste("~", paste(X, collapse = " + "))) else ~ 1
  B <- stats::model.matrix(basis, data = data)
  for (j in seq_len(ncol(B))) if (colnames(B)[j] != "(Intercept)") {
    cen <- mean(B[tr, j]); sc <- stats::sd(B[tr, j]); if (!is.finite(sc) || sc == 0) sc <- 1
    B[, j] <- (B[, j] - cen) / sc
  }
  list(train = B[tr, , drop = FALSE], test = B[te, , drop = FALSE])
}

l2_learner_predict <- function(y, xtrain, xnew, family, method, library) {
  fallback <- mean(y, na.rm = TRUE)
  if (!is.finite(fallback)) fallback <- 0
  if (length(y) < 5L || (family == "binomial" && length(unique(y)) < 2L))
    return(rep(fallback, nrow(xnew)))
  if (method == "SuperLearner") {
    fit <- try(SuperLearner::SuperLearner(Y = y, X = as.data.frame(xtrain),
      newX = as.data.frame(xnew), family = family, SL.library = library,
      cvControl = list(V = min(5L, max(2L, floor(length(y) / 20))))), silent = TRUE)
    if (!inherits(fit, "try-error")) return(as.numeric(fit$SL.predict))
  }
  dd <- data.frame(.y = y, xtrain, check.names = FALSE)
  nd <- data.frame(xnew, check.names = FALSE)
  fit <- try(stats::glm(.y ~ ., data = dd, family = family), silent = TRUE)
  if (inherits(fit, "try-error")) rep(fallback, nrow(xnew)) else
    as.numeric(stats::predict.glm(fit, newdata = nd, type = "response"))
}

l2_split_nuisance <- function(data, Y, A, C, X, tr, te, method,
                              lib_e, lib_r, lib_y) {
  xt <- data[tr, X, drop = FALSE]; xa <- data[, X, drop = FALSE]
  e1 <- l2_learner_predict(data[[A]][tr], xt, xa, "binomial", method, lib_e)
  out <- list(e0 = clip_probs(1 - e1, 1e-4), e1 = clip_probs(e1, 1e-4))
  for (a in 0:1) {
    ia <- tr & data[[A]] == a
    rho <- l2_learner_predict(as.numeric(data[[C]][ia] == 0),
      data[ia, X, drop = FALSE], xa, "binomial", method, lib_r)
    iy <- ia & data[[C]] == 0
    mu <- l2_learner_predict(data[[Y]][iy], data[iy, X, drop = FALSE], xa,
                             "gaussian", method, lib_y)
    out[[paste0("rho", a)]] <- clip_probs(rho, 1e-4)
    out[[paste0("mu", a)]] <- mu
  }
  list(train = lapply(out, `[`, tr), test = lapply(out, `[`, te))
}

l2_split_nuisance_oracle_error <- function(data, X, tr, te, rate, scale, sign) {
  needed <- c("p_A1_given_X", "p_C0_given_X", "p_C1_given_X",
              "mu_Y0_given_XC0", "mu_Y1_given_XC0")
  if (!all(needed %in% names(data)))
    stop("nuisance_source = 'oracle_error' requires data from simulate_l2_data()")
  if (length(rate) != 1L || !is.finite(rate) || rate < 0)
    stop("nuisance_error_rate must be one nonnegative number")
  scale <- l2_nuisance_error_scale(scale)
  sign <- l2_nuisance_error_sign(sign)
  xx <- as.matrix(data[, X, drop = FALSE]); storage.mode(xx) <- "double"
  G <- cbind(1, xx, xx^2)
  if (ncol(xx) > 1L) G <- cbind(G, xx[, 1L] * xx[, 2L])
  smooth_error <- function(index) {
    ## Keep the error direction fixed across Monte Carlo replications so that
    ## first-order plug-in bias does not cancel by construction.
    coef <- sin(seq_len(ncol(G)) * (index + .5)) +
      .5 * cos(seq_len(ncol(G)) * (index + 1.5))
    g <- drop(G %*% coef)
    s <- stats::sd(g[tr]); if (!is.finite(s) || s == 0) s <- 1
    (g - mean(g[tr])) / s
  }
  mag <- sum(tr)^(-rate)
  logit <- function(p) stats::qlogis(clip_probs(p, 1e-5))
  e1 <- stats::plogis(logit(data$p_A1_given_X) +
    sign[["e"]] * scale[["e"]] * mag * smooth_error(1))
  out <- list(e0 = clip_probs(1 - e1, 1e-4), e1 = clip_probs(e1, 1e-4))
  for (a in 0:1) {
    rho0 <- 1 - data[[paste0("p_C", a, "_given_X")]]
    rho <- stats::plogis(logit(rho0) +
      sign[["rho"]] * scale[["rho"]] * mag * smooth_error(2 + a))
    mu <- data[[paste0("mu_Y", a, "_given_XC0")]] +
      sign[["mu"]] * scale[["mu"]] * mag * smooth_error(4 + a)
    out[[paste0("rho", a)]] <- clip_probs(rho, 1e-4)
    out[[paste0("mu", a)]] <- mu
  }
  out$err_H <- smooth_error(7)
  out$err_de <- smooth_error(8)
  out$err_drho <- smooth_error(9)
  for (k in 1:4) out[[paste0("err_nu", k)]] <- smooth_error(9 + k)
  list(train = lapply(out, `[`, tr), test = lapply(out, `[`, te))
}

l2_nuisance_error_scale <- function(scale) {
  target <- c("e", "rho", "mu", "regression")
  if (length(scale) == 1L && is.null(names(scale))) {
    scale <- rep(scale, 4L); names(scale) <- target
  } else if (is.null(names(scale))) {
    old_target <- c("e", "rho", "regression")
    names(scale) <- if (length(scale) == 3L) old_target else
      target[seq_along(scale)]
  }
  if (!"mu" %in% names(scale) && "regression" %in% names(scale))
    scale[["mu"]] <- scale[["regression"]]
  if (!all(target %in% names(scale)) || any(!is.finite(scale[target])) ||
      any(scale[target] < 0))
    stop("nuisance_error_scale must be nonnegative with entries e, rho, mu, and regression")
  scale[target]
}

l2_nuisance_error_sign <- function(sign) {
  target <- c("e", "rho", "mu", "regression")
  if (length(sign) == 1L && is.null(names(sign))) {
    sign <- rep(sign, 4L); names(sign) <- target
  } else if (is.null(names(sign))) {
    old_target <- c("e", "rho", "regression")
    names(sign) <- if (length(sign) == 3L) old_target else
      target[seq_along(sign)]
  }
  if (!"mu" %in% names(sign) && "regression" %in% names(sign))
    sign[["mu"]] <- sign[["regression"]]
  if (!all(target %in% names(sign)) || any(!is.finite(sign[target])) ||
      any(!sign[target] %in% c(-1, 1)))
    stop("nuisance_error_sign must have entries e, rho, mu, and regression equal to -1 or 1")
  sign[target]
}

l2_oracle_dual_regressions <- function(data, a, B, delta, par, sigma,
                                       boundary_tol, nodes = 100L,
                                       nuisance_override = NULL,
                                       model = c("net", "separated"),
                                       delta_M = 0) {
  model <- match.arg(model)
  ## Under simulate_l2_data(), observation given Y,X
  ## is proportional to 1 - p_UI(a|Y,X). The noninformative observation factor
  ## depends only on X and therefore cancels from the respondent-law weights.
  q <- l2_dgp_error_quantiles(data, nodes)
  mu <- data[[paste0("mu_Y", a)]]
  e <- if (is.null(nuisance_override)) {
    if (a == 1) data$p_A1_given_X else 1 - data$p_A1_given_X
  } else nuisance_override[[paste0("e", a)]]
  rho <- if (is.null(nuisance_override))
    1 - data[[paste0("p_C", a, "_given_X")]] else
    nuisance_override[[paste0("rho", a)]]
  nr <- nrow(data); ii <- rep(seq_len(nr), times = nodes)
  yy <- as.vector(outer(mu, q, "+"))
  link <- if ("density_ratio_link" %in% names(data))
    as.character(data$density_ratio_link[1L]) else "linear"
  latent <- function(z) if (link == "bounded") tanh(z) else z
  has_mar <- if ("informative_missingness" %in% names(data))
    data$informative_missingness[1L] == 1 else TRUE
  informative_intercept <- l2_dgp_scalar(data, "informative_intercept", -2.4)
  informative_coef <- if ("informative_outcome_coef" %in% names(data))
    data$informative_outcome_coef[ii] else rep(.45, length(ii))
  pui <- if (has_mar)
    stats::plogis(informative_intercept + .35 * a + .25 * data$X1[ii] +
                   informative_coef * latent(yy - mu[ii])) else
    rep(0, length(yy))
  pA1 <- stats::plogis(data$treatment_lp_X[ii] +
    data$confounding_strength[ii] * latent(yy - mu[ii]))
  ga <- if (a == 1) pA1 else 1 - pA1
  w <- ga * (1 - pui)
  nuisance <- list()
  nuisance[[paste0("e", a)]] <- e[ii]
  nuisance[[paste0("rho", a)]] <- rho[ii]
  ev <- l2_dual_evaluate(yy, rep(a, length(yy)), rep(0, length(yy)),
    a, nuisance, B[ii, , drop = FALSE], delta, par, sigma, boundary_tol,
    model = model, delta_M = delta_M)
  W <- matrix(w, nrow = nr, ncol = nodes)
  den <- rowSums(W)
  weighted <- function(z) rowSums(W * matrix(z, nr, nodes)) / den
  list(mH = weighted(ev$H), d_e = weighted(ev$d_e),
       d_rho = weighted(ev$d_rho))
}

l2_pseudo_regression <- function(train, X, pseudo, use, test, method, library) {
  ## pseudo is already restricted to use rows.
  xtr <- train[use, X, drop = FALSE]
  l2_learner_predict(pseudo, xtr, test[, X, drop = FALSE], "gaussian", method, library)
}

l2_dual_evaluate <- function(y, A, C, a, nuisance, B, delta, par, sigma,
                             boundary_tol = 1e-7,
                             model = c("separated", "net"), delta_M = 0) {
  model <- match.arg(model)
  n <- length(A); p <- ncol(B); eps <- 1e-8
  e <- nuisance[[paste0("e", a)]]; rho <- nuisance[[paste0("rho", a)]]
  b <- delta * (1 - rho); q <- 1 - e
  lamR <- exp(par[1]); lamA <- exp(par[2])
  alpha <- drop(B %*% par[2 + seq_len(p)])
  beta <- drop(B %*% par[2 + p + seq_len(p)])
  yy <- y; yy[!is.finite(yy)] <- 0
  no_A <- is.infinite(lamA)
  fixed_M <- is.infinite(lamR)
  K <- if (no_A) rep(1, n) else
    pmax(0, 1 + (sigma * yy * q + beta) / (2 * lamA))
  Fbase <- if (no_A) sigma * yy + alpha else
    sigma * yy * e + alpha - lamA * (1 - K^2)
  if (fixed_M) {
    lower <- rep(-Inf, n); M0 <- M <- rep(1, n)
  } else if (model == "separated") {
    M0 <- 1 + b^2 / (2 * lamR) * Fbase
    lower <- 1 - b; M <- pmax(lower, M0)
  } else {
    M0 <- 1 + Fbase / (2 * lamR)
    lower <- rho; M <- pmax(lower, M0)
  }
  L <- M * K
  divR <- if (fixed_M) rep(0, n) else if (model == "separated")
    (M - 1)^2 / pmax(b^2, eps) else (M - 1)^2
  divA <- (L - M)^2 / pmax(M, eps)
  H <- sigma * yy * (e * M + q * L) + alpha * (M - 1) +
    beta * (L - 1) - (if (fixed_M) 0 else lamR * divR) -
    (if (no_A) 0 else lamA * divA)
  penalty_gradient <- if (fixed_M) rep(0, n) else if (model == "separated")
    2 * lamR * (M - 1) / pmax(b^2, eps) else 2 * lamR * (M - 1)
  Fm <- Fbase - penalty_gradient
  active <- !fixed_M & M <= lower + boundary_tol
  eta <- ifelse(active, pmax(0, -Fm), 0)
  d_e <- sigma * yy * (M - L)
  d_rho <- if (fixed_M) rep(0, n) else if (model == "separated")
    -2 * lamR * delta * (M - 1)^2 / pmax(b^3, eps) - delta * eta else
    -eta
  obs <- A == a & C == 0
  H[!obs] <- d_e[!obs] <- d_rho[!obs] <- NA_real_
  list(H = H, d_e = d_e, d_rho = d_rho, M = M, L = L,
       eta = eta, active = active, obs = obs)
}

#' Simulate plug-in versus DR sharp-bound performance
#'
#' @param B Number of Monte Carlo replications.
#' @param n Sample size per replication.
#' @param confounding_strength Passed to `simulate_l2_data()`; zero recovers
#'   conditional exchangeability.
#' @param violation Passed to `simulate_l2_data()`; identifies which of MAR and
#'   NUC are violated.
#' @param density_ratio_link Passed to `simulate_l2_data()`.
#' @param truth Optional named vector `c(lower, upper)` of sharp ATE endpoints.
#' @param truth_n Size of the Monte Carlo reference sample when `truth` is NULL.
#' @param nuisance_simulation Use fitted nuisance models or fast, controlled
#'   perturbations of the simulator's oracle nuisance functions.
#' @param parameter_evaluation `"auto"` uses the minimal oracle parameters when
#'   no sensitivity parameters are supplied, `"oracle"` always uses them, and
#'   `"specified"` uses the supplied/default arguments. Oracle parameters are
#'   returned in all cases.
#' @param linf_method Which L-infinity comparison to compute: the Holder outer
#'   bound, the sharp empirical sieve LP, or both.
#' @param linf_max_n Maximum calibration-sample size used by the sharp LP.
#' @param dgp_args Named list of additional arguments passed to
#'   `simulate_l2_data()`, such as informative-missingness strength, overlap,
#'   covariate complexity, or the outcome distribution.
#' @param ... Arguments passed to `l2_sharp_crossfit()`.
#' @return Replication-level results, endpoint bias/coverage summaries, the
#'   reference endpoints, and bound coverage of the simulated true ATE.
#' @export
l2_compare_dr_plugin <- function(B = 200L, n = 1000L,
                                 violation = c("both", "mar", "nuc", "none"),
                                 density_ratio_link = c("bounded", "linear"),
                                 confounding_strength = .6, truth = NULL,
                                 truth_n = 20000L, seed = 1L,
                                 conf_level = .95,
                                 nuisance_simulation = c("oracle_error", "estimated"),
                                 parameter_evaluation = c("auto", "oracle", "specified"),
                                 linf_method = c("outer", "sharp", "both"),
                                 l2_method = c("sharp", "cs", "both"),
                                 linf_max_n = 1000L,
                                 dgp_args = list(),
                                 ...) {
  B <- as.integer(B); if (B < 2L) stop("B must be at least 2")
  violation <- match.arg(violation)
  density_ratio_link <- match.arg(density_ratio_link)
  nuisance_simulation <- match.arg(nuisance_simulation)
  parameter_evaluation <- match.arg(parameter_evaluation)
  linf_method <- match.arg(linf_method)
  l2_method <- match.arg(l2_method)
  if (!is.list(dgp_args) || (length(dgp_args) &&
      (is.null(names(dgp_args)) || any(!nzchar(names(dgp_args))))))
    stop("dgp_args must be a named list")
  forbidden_dgp <- intersect(names(dgp_args), c("n", "seed", "violation",
    "density_ratio_link", "confounding_strength"))
  if (length(forbidden_dgp)) stop("Specify ", paste(forbidden_dgp, collapse = ", "),
    " through the corresponding l2_compare_dr_plugin argument, not dgp_args")
  truth_supplied <- !is.null(truth)
  dots <- list(...)
  if (nuisance_simulation == "estimated" &&
      l2_method %in% c("sharp", "both") && !is.null(dots$folds) &&
      is.finite(dots$folds) && dots$folds < 3L) {
    warning("The sharp one-step estimator requires three disjoint roles; folds was promoted to 3")
    dots$folds <- 3L
  }
  if (is.null(dots$diagnostic_control))
    dots$diagnostic_control <- list(warn = FALSE, compute_hessian = FALSE)
  if (!is.null(dots$delta_K)) {
    dots$delta_A <- dots$delta_K
    dots$delta_K <- NULL
  }
  sensitivity_given <- any(c("delta", "delta_R", "delta_A", "delta_M") %in%
    names(dots))
  if ("nuisance_source" %in% names(dots))
    stop("Use nuisance_simulation, rather than nuisance_source, in this function")
  dots$nuisance_source <- if (nuisance_simulation == "oracle_error")
    "oracle_error" else "estimated"
  if (nuisance_simulation == "oracle_error" &&
      !"nuisance_method" %in% names(dots)) dots$nuisance_method <- "glm"
  merge_args <- function(base, extra) {
    if (!length(extra)) return(base)
    if (is.null(names(extra)) || any(!nzchar(names(extra))))
      stop("All arguments supplied through ... must be named")
    base[names(extra)] <- NULL
    c(base, extra)
  }
  calibration_data <- do.call(simulate_l2_data, c(list(n = truth_n,
    seed = seed + 100000L, violation = violation,
    density_ratio_link = density_ratio_link,
    confounding_strength = confounding_strength), dgp_args))
  oracle_parameters <- l2_oracle_sensitivity_parameters(calibration_data)
  linf_bounds <- l2_oracle_linf_outer_bounds(calibration_data,
    oracle_parameters)
  linf_sharp_bounds <- if (linf_method %in% c("sharp", "both"))
    l2_oracle_linf_sharp_bounds(calibration_data, oracle_parameters,
      basis = dots$basis, max_n = linf_max_n) else NULL
  use_oracle <- parameter_evaluation == "oracle" ||
    (parameter_evaluation == "auto" && !sensitivity_given)
  if (use_oracle) {
    dots$delta <- oracle_parameters$delta
    dots$delta_R <- oracle_parameters$delta_R
    dots$delta_A <- oracle_parameters$delta_A
    dots$delta_M <- oracle_parameters$delta_M
  }
  if (is.null(dots$delta)) dots$delta <- c(1, 1)
  if (is.null(dots$delta_R)) dots$delta_R <- c(.2, .2)
  if (is.null(dots$delta_A)) dots$delta_A <- c(.2, .2)
  if (is.null(dots$delta_M)) dots$delta_M <- dots$delta_R
  if (is.null(dots$model)) dots$model <- "net"
  if (is.null(truth)) {
    truth <- l2_oracle_sharp_reference(calibration_data,
      delta = if (is.null(dots$delta)) c(1, 1) else dots$delta,
      delta_R = if (is.null(dots$delta_R)) c(.2, .2) else dots$delta_R,
      delta_A = if (is.null(dots$delta_A)) c(.2, .2) else dots$delta_A,
      delta_M = dots$delta_M, model = dots$model,
      basis = dots$basis, control = if (is.null(dots$control)) list() else dots$control)
  }
  if (length(truth) != 2L || any(!is.finite(truth)))
    stop("truth must contain finite lower and upper endpoints")
  names(truth) <- c("lower", "upper")
  cs_ref <- l2_oracle_cs(calibration_data, dots$delta, dots$delta_R,
    dots$delta_A, dots$delta_M, dots$model)$ate[1L, c("lower", "upper")]
  cs_truth <- c(lower = cs_ref$lower, upper = cs_ref$upper)
  reps <- list(); fit_errors <- character();
  z <- stats::qnorm(1 - (1 - conf_level) / 2); k <- 0L
  for (b in seq_len(B)) {
    dat <- do.call(simulate_l2_data, c(list(n = n, seed = seed + b,
      violation = violation, density_ratio_link = density_ratio_link,
      confounding_strength = confounding_strength), dgp_args))
    args <- merge_args(list(data = dat, Y = "Y", A = "A", C = "C",
                   X = c("X1", "X2"), seed = seed + b,
                   conf_level = conf_level), dots)
    fits <- list()
    if (l2_method %in% c("sharp", "both"))
      fits$sharp <- try(do.call(l2_sharp_crossfit, args), silent = TRUE)
    if (l2_method %in% c("cs", "both")) {
      cs_args <- args
      cs_args[c("basis", "control", "boundary_tol", "diagnostic_control")] <- NULL
      fits$cs <- try(do.call(l2_cs_crossfit, cs_args), silent = TRUE)
    }
    bad <- vapply(fits, inherits, logical(1), "try-error")
    if (any(bad)) { fit_errors <- c(fit_errors, as.character(fits[bad][[1L]])); next }
    for (meth in names(fits)) for (j in seq_len(nrow(fits[[meth]]$ate))) {
      k <- k + 1L; zz <- fits[[meth]]$ate[j, ]
      tv <- if (meth == "cs") cs_truth[[zz$endpoint]] else truth[[zz$endpoint]]
      opt <- if (meth == "sharp") {
        duals <- unlist(lapply(fits[[meth]]$diagnostics,
          function(ep) lapply(ep, `[[`, "dual")), recursive = FALSE)
        list(converged = all(vapply(duals, function(x) x$convergence == 0,
          logical(1))),
          diagnostic_ok = all(vapply(duals, function(x) x$diagnostic_ok,
            logical(1))),
          diagnostic_codes = paste(unique(unlist(lapply(duals,
            `[[`, "diagnostic_codes"))), collapse = ";"),
          gap = max(vapply(duals, function(x) abs(x$primal_dual_gap),
            numeric(1)), na.rm = TRUE),
          norm = max(vapply(duals, function(x) x$maximum_normalization_error,
            numeric(1)), na.rm = TRUE),
          budget = max(vapply(duals, function(x) x$maximum_budget_violation,
            numeric(1)), na.rm = TRUE),
          active = max(vapply(unlist(fits[[meth]]$diagnostics,
            recursive = FALSE), function(x) x$active_fraction,
            numeric(1)), na.rm = TRUE))
      } else list(converged = NA, diagnostic_ok = NA,
                  diagnostic_codes = "", gap = NA, norm = NA,
                  budget = NA, active = NA)
      reps[[k]] <- data.frame(replication = b, bound_method = meth,
        estimator = zz$estimator,
        endpoint = zz$endpoint, estimate = zz$estimate, se = zz$se,
        truth = tv, covered = zz$conf_low <= tv & tv <= zz$conf_high,
        true_ate = dat$true_ate[1], optimizer_converged = opt$converged,
        diagnostic_ok = opt$diagnostic_ok,
        diagnostic_codes = opt$diagnostic_codes,
        maximum_primal_dual_gap = opt$gap,
        maximum_normalization_error = opt$norm,
        maximum_budget_violation = opt$budget,
        maximum_active_fraction = opt$active,
        stringsAsFactors = FALSE)
    }
  }
  d <- do.call(rbind, reps)
  if (is.null(d) || !nrow(d)) stop("All simulation fits failed. First error: ",
    if (length(fit_errors)) fit_errors[1L] else "unknown error")
  spl <- split(d, interaction(d$bound_method, d$estimator, d$endpoint, drop = TRUE))
  summary <- do.call(rbind, lapply(spl, function(x) data.frame(
    bound_method = x$bound_method[1], estimator = x$estimator[1],
    endpoint = x$endpoint[1], n_success = nrow(x),
    bias = mean(x$estimate - x$truth), rmse = sqrt(mean((x$estimate - x$truth)^2)),
    empirical_sd = stats::sd(x$estimate), mean_se = mean(x$se),
    coverage = mean(x$covered),
    diagnostic_flag_rate = if (all(is.na(x$diagnostic_ok))) NA_real_ else
      mean(!x$diagnostic_ok, na.rm = TRUE))))
  wide <- stats::reshape(d[, c("replication", "bound_method", "estimator", "endpoint",
                        "estimate", "se", "true_ate")],
                  idvar = c("replication", "bound_method", "estimator", "true_ate"),
                  timevar = "endpoint", direction = "wide")
  point_contains <- wide$estimate.lower <= wide$true_ate &
    wide$true_ate <= wide$estimate.upper
  wald_contains <- wide$estimate.lower - z * wide$se.lower <= wide$true_ate &
    wide$true_ate <= wide$estimate.upper + z * wide$se.upper
  bound_coverage <- stats::aggregate(cbind(point_contains, wald_contains),
    list(bound_method = wide$bound_method, estimator = wide$estimator), mean)
  names(bound_coverage)[3:4] <- c("true_ate_bound_coverage",
    "true_ate_outward_wald_coverage")
  bound_comparison <- data.frame(model = paste("L2", dots$model, "sharp sieve"),
    lower = truth[["lower"]], upper = truth[["upper"]])
  if (l2_method %in% c("cs", "both")) bound_comparison <- rbind(
    bound_comparison, data.frame(model = paste("L2", dots$model,
      "Cauchy-Schwarz outer"),
      lower = cs_truth[["lower"]], upper = cs_truth[["upper"]]))
  if (linf_method %in% c("outer", "both")) bound_comparison <- rbind(
    bound_comparison, data.frame(model = "L-infinity/L1 Holder outer",
      lower = linf_bounds$ate[["lower"]], upper = linf_bounds$ate[["upper"]]))
  if (!is.null(linf_sharp_bounds)) bound_comparison <- rbind(
    bound_comparison, data.frame(model = "L-infinity sharp sieve",
      lower = linf_sharp_bounds$ate[["lower"]],
      upper = linf_sharp_bounds$ate[["upper"]]))
  comparison_model <- if (dots$model == "net") "separated" else "net"
  comparison_truth <- l2_oracle_sharp_reference(calibration_data,
    delta = oracle_parameters$delta,
    delta_R = oracle_parameters$delta_R,
    delta_A = oracle_parameters$delta_A,
    delta_M = oracle_parameters$delta_M,
    model = comparison_model, basis = dots$basis,
    control = if (is.null(dots$control)) list() else dots$control)
  comparison_cs <- l2_oracle_cs(calibration_data,
    oracle_parameters$delta, oracle_parameters$delta_R,
    oracle_parameters$delta_A, oracle_parameters$delta_M,
    comparison_model)$ate[1L, c("lower", "upper")]
  bound_comparison <- rbind(bound_comparison,
    data.frame(model = paste("L2", comparison_model,
      "sharp sieve (parameterization comparison)"),
      lower = comparison_truth[["lower"]],
      upper = comparison_truth[["upper"]]),
    data.frame(model = paste("L2", comparison_model,
      "Cauchy-Schwarz outer (parameterization comparison)"),
      lower = comparison_cs$lower, upper = comparison_cs$upper))
  bound_comparison$width <- bound_comparison$upper - bound_comparison$lower
  list(replicates = d, summary = summary, bound_coverage = bound_coverage,
       n_requested = B, n_success = length(unique(d$replication)),
       failure_rate = 1 - length(unique(d$replication)) / B,
       fit_errors = unique(fit_errors),
       truth = truth, truth_note = if (truth_supplied) "user supplied" else
         "Monte Carlo oracle empirical-sieve reference",
       cs_truth = cs_truth,
       oracle_parameters = oracle_parameters,
       linf_bounds = linf_bounds,
       linf_sharp_bounds = linf_sharp_bounds,
       bound_comparison = bound_comparison,
       evaluated_parameters = list(delta = dots$delta, delta_R = dots$delta_R,
                                    delta_A = dots$delta_A,
                                    delta_K = dots$delta_A,
                                    delta_M = dots$delta_M,
                                    model = dots$model),
       parameter_evaluation = if (use_oracle) "oracle" else "specified",
       confounding_strength = confounding_strength,
       effective_confounding_strength = calibration_data$confounding_strength[1L],
       dgp_args = dgp_args,
       violation = violation,
       density_ratio_link = density_ratio_link,
       nuisance_simulation = nuisance_simulation, call = match.call())
}

#' Minimal oracle sensitivity parameters for the continuous DGP
#'
#' Computes the smallest arm-specific prevalence envelope and informative-
#' missingness L2 radius under the known law generated by `simulate_l2_data()`.
#' Treatment depends on the shared potential-outcome disturbance, producing a
#' nonzero oracle confounding radius unless `confounding_strength = 0`.
#' Integrals over the Gaussian potential-outcome law use deterministic normal
#' quadrature.
#'
#' @param data Data returned by `simulate_l2_data()`.
#' @param nodes Number of normal quadrature nodes.
#' @return A list containing arm-specific `delta`, `delta_R`, `delta_A`, and
#'   covariate-level prevalence ratios.
#' @export
l2_oracle_sensitivity_parameters <- function(data, nodes = 200L) {
  needed <- c("mu_Y0", "mu_Y1", "X1", "p_C0_given_X", "p_C1_given_X",
              "treatment_lp_X", "confounding_strength")
  if (!all(needed %in% names(data)))
    stop("Oracle sensitivity parameters require data from simulate_l2_data()")
  nodes <- as.integer(nodes)
  if (length(nodes) != 1L || !is.finite(nodes) || nodes < 20L)
    stop("nodes must be an integer of at least 20")
  q <- l2_dgp_error_quantiles(data, nodes)
  link <- if ("density_ratio_link" %in% names(data))
    as.character(data$density_ratio_link[1L]) else "linear"
  latent <- function(z) if (link == "bounded") tanh(z) else z
  has_mar <- if ("informative_missingness" %in% names(data))
    data$informative_missingness[1L] == 1 else TRUE
  ## With informative missingness, the Gaussian tail makes the population
  ## essential supremum one. Under MAR, the informative prevalence is zero.
  delta <- if (has_mar) c(1, 1) else c(0, 0)
  delta_R <- delta_M <- delta_A <- delta_R_inf <- delta_A_inf <-
    empirical_delta <- numeric(2L)
  ratios <- vector("list", 2L)
  for (a in 0:1) {
    mu <- data[[paste0("mu_Y", a)]]
    yy <- outer(mu, q, "+")
    informative_intercept <- l2_dgp_scalar(data, "informative_intercept", -2.4)
    informative_coef <- if ("informative_outcome_coef" %in% names(data))
      data$informative_outcome_coef else rep(.45, nrow(data))
    pui <- if (has_mar)
      stats::plogis(outer(informative_intercept + .35 * a + .25 * data$X1,
                          rep(1, length(q))) +
        outer(informative_coef, latent(q))) else
      matrix(0, nrow = nrow(data), ncol = nodes)
    pA1 <- stats::plogis(outer(data$treatment_lp_X,
      data$confounding_strength[1L] * latent(q), "+"))
    gF <- if (a == 1) pA1 else 1 - pA1
    gC <- 1 - gF
    eF <- if (a == 1) data$p_A1_given_X else 1 - data$p_A1_given_X
    eC <- 1 - eF
    pi_star <- rowMeans(gF * pui) / eF
    pi_total <- data[[paste0("p_C", a, "_given_X")]]
    ratio <- ifelse(pi_total > 0, pi_star / pi_total, 0)
    empirical_delta[a + 1L] <- min(1, max(ratio))
    ## R = f_I/f_O. Conditional respondent weights are proportional to 1-pui.
    joint_I <- rowMeans(gF * pui)
    joint_O0 <- rowMeans(gF * (1 - pui))
    if (has_mar) {
      R <- sweep(pui / (1 - pui), 1L, joint_O0 / joint_I, "*")
      wO <- gF * (1 - pui)
      chi_x <- rowSums(wO * (R - 1)^2) / rowSums(wO)
      delta_R[a + 1L] <- sqrt(mean(chi_x))
      M <- sweep(1 / (1 - pui), 1L, joint_O0 / eF, "*")
      chi_M_x <- rowSums(wO * (M - 1)^2) / rowSums(wO)
      delta_M[a + 1L] <- sqrt(mean(chi_M_x))
      if (link == "bounded") {
        pui_ext <- stats::plogis(
          outer(informative_intercept + .35 * a + .25 * data$X1,
                rep(1, 2L)) + outer(informative_coef, c(-1, 1)))
        R_ext <- sweep(pui_ext / (1 - pui_ext), 1L,
                       joint_O0 / joint_I, "*")
        delta_R_inf[a + 1L] <- max(abs(R_ext - 1))
      } else delta_R_inf[a + 1L] <- Inf
    } else {
      delta_R[a + 1L] <- 0
      delta_M[a + 1L] <- 0
      delta_R_inf[a + 1L] <- 0
    }
    K <- (gC / eC) / (gF / eF)
    chi_A_x <- rowMeans((gF / eF) * (K - 1)^2)
    delta_A[a + 1L] <- sqrt(mean(chi_A_x))
    if (data$confounding_strength[1L] == 0) delta_A_inf[a + 1L] <- 0 else
      if (link == "bounded") {
        pA_ext <- stats::plogis(outer(data$treatment_lp_X,
          data$confounding_strength[1L] * c(-1, 1), "+"))
        gf_ext <- if (a == 1) pA_ext else 1 - pA_ext
        gc_ext <- 1 - gf_ext
        K_ext <- (gc_ext / eC) / (gf_ext / eF)
        delta_A_inf[a + 1L] <- max(abs(K_ext - 1))
      } else delta_A_inf[a + 1L] <- Inf
    ratios[[a + 1L]] <- ratio
  }
  names(delta) <- names(delta_R) <- names(delta_M) <- names(delta_A) <-
    names(delta_R_inf) <- names(delta_A_inf) <-
    names(empirical_delta) <- c("arm0", "arm1")
  list(delta = delta, delta_R = delta_R, delta_M = delta_M,
       delta_A = delta_A, delta_K = delta_A,
       delta_R_inf = delta_R_inf, delta_A_inf = delta_A_inf,
       empirical_delta = empirical_delta,
       prevalence_ratio = ratios,
       note = paste(if (has_mar)
         paste("The population prevalence envelope is one because X1 has",
               "unbounded Gaussian support.") else
         "The prevalence envelope and informative-missingness radius are zero under MAR.",
         "L2 radii use Monte Carlo over the supplied X distribution with",
         nodes, "deterministic quantile-quadrature nodes per X.",
         "delta_K is the reported name for the backward-compatible delta_A argument.",
         "Set confounding_strength=0 to recover delta_K=0."))
}

l2_oracle_linf_outer_bounds <- function(data, parameters, nodes = 200L) {
  rR <- parameters$delta_R_inf; rA <- parameters$delta_A_inf
  if (any(!is.finite(c(rR, rA)))) return(list(
    arm = NULL, ate = c(lower = -Inf, upper = Inf), finite = FALSE,
    method = "L-infinity/L1 Holder outer bound",
    note = "The linear-link DGP has infinite oracle L-infinity density-ratio radii."))
  q <- l2_dgp_error_quantiles(data, nodes)
  link <- as.character(data$density_ratio_link[1L])
  latent <- function(z) if (link == "bounded") tanh(z) else z
  arm <- matrix(NA_real_, 2L, 2L, dimnames = list(c("arm0", "arm1"),
    c("lower", "upper")))
  for (a in 0:1) {
    mu <- data[[paste0("mu_Y", a)]]
    yy <- outer(mu, q, "+")
    informative_intercept <- l2_dgp_scalar(data, "informative_intercept", -2.4)
    informative_coef <- if ("informative_outcome_coef" %in% names(data))
      data$informative_outcome_coef else rep(.45, nrow(data))
    pui <- if (data$informative_missingness[1L] == 1)
      stats::plogis(outer(informative_intercept + .35 * a + .25 * data$X1,
                          rep(1, length(q))) +
        outer(informative_coef, latent(q))) else
      matrix(0, nrow(data), nodes)
    pA1 <- stats::plogis(outer(data$treatment_lp_X,
      data$confounding_strength[1L] * latent(q), "+"))
    gF <- if (a == 1) pA1 else 1 - pA1
    wO <- gF * (1 - pui)
    center <- data[[paste0("mu_Y", a, "_given_XC0")]]
    absY <- rowSums(wO * abs(yy - center)) / rowSums(wO)
    e <- if (a == 1) data$p_A1_given_X else 1 - data$p_A1_given_X
    b <- parameters$delta[a + 1L] * data[[paste0("p_C", a, "_given_X")]]
    rad <- mean(absY * (b * rR[a + 1L] +
      (1 - e) * (1 + b * rR[a + 1L]) * rA[a + 1L]))
    psi <- mean(data[[paste0("mu_Y", a, "_given_XC0")]])
    arm[a + 1L, ] <- c(psi - rad, psi + rad)
  }
  list(arm = arm,
    ate = c(lower = arm[2L, "lower"] - arm[1L, "upper"],
            upper = arm[2L, "upper"] - arm[1L, "lower"]),
    finite = TRUE, method = "L-infinity/L1 Holder outer bound",
    note = "Valid outer comparison; conditional normalization can make the sharp L-infinity interval narrower.")
}

l2_oracle_linf_sharp_bounds <- function(data, parameters, basis = NULL,
                                        max_n = 1000L) {
  if (!requireNamespace("lpSolve", quietly = TRUE))
    stop("linf_method = 'sharp' requires the lpSolve package")
  if (any(!is.finite(c(parameters$delta_R_inf,
                       parameters$delta_A_inf))))
    return(list(arm = NULL, ate = c(lower = -Inf, upper = Inf),
      finite = FALSE, method = "L-infinity sharp empirical sieve",
      note = "The DGP has infinite oracle L-infinity radii."))
  max_n <- as.integer(max_n)
  if (!is.finite(max_n) || max_n < 100L) stop("linf_max_n must be at least 100")
  if (nrow(data) > max_n) {
    ii <- unique(round(seq(1, nrow(data), length.out = max_n)))
    data <- data[ii, , drop = FALSE]
  }
  B <- l2_basis(data, c("X1", "X2"), basis)
  nuisance <- list(
    e1 = data$p_A1_given_X, e0 = 1 - data$p_A1_given_X,
    rho0 = 1 - data$p_C0_given_X, rho1 = 1 - data$p_C1_given_X,
    mu0 = data$mu_Y0_given_XC0, mu1 = data$mu_Y1_given_XC0)
  arm <- matrix(NA_real_, 2L, 2L, dimnames = list(c("arm0", "arm1"),
    c("lower", "upper")))
  diagnostics <- vector("list", 2L)
  for (a in 0:1) {
    r <- l2_reference_quantities(data$Y, data$A, data$C, a, nuisance)
    use <- which(r$obs); m <- length(use)
    if (!m) stop("No respondents in arm ", a)
    ww <- r$w[use] / nrow(data); BB <- B[use, , drop = FALSE]
    y <- r$y[use]; e <- r$e[use]; q <- r$q[use]
    b <- parameters$delta[a + 1L] * r$pi[use]
    rR <- parameters$delta_R_inf[a + 1L]
    rA <- parameters$delta_A_inf[a + 1L]
    ## Nonnegativity of R additionally requires M >= 1-b.  This constraint is
    ## active whenever the symmetric L-infinity radius exceeds one.
    mlo <- pmax(1 - b, 1 - b * rR); mup <- 1 + b * rR
    klo <- pmax(0, 1 - rA); kup <- 1 + rA
    Z <- matrix(0, m, m); I <- diag(m)
    eqM <- cbind(t(BB * ww), matrix(0, ncol(BB), m))
    eqL <- cbind(matrix(0, ncol(BB), m), t(BB * ww))
    rhsM <- colSums(BB * ww)
    Acon <- rbind(eqM, eqL,
      cbind(I, Z), cbind(I, Z),
      cbind(-kup * I, I), cbind(-klo * I, I))
    dirs <- c(rep("=", 2 * ncol(BB)), rep("<=", m), rep(">=", m),
              rep("<=", m), rep(">=", m))
    rhs <- c(rhsM, rhsM, mup, mlo, rep(0, 2 * m))
    objective <- c(ww * y * e, ww * y * q)
    lo <- lpSolve::lp("min", objective, Acon, dirs, rhs)
    hi <- lpSolve::lp("max", objective, Acon, dirs, rhs)
    if (lo$status != 0 || hi$status != 0)
      stop("Sharp L-infinity LP failed in arm ", a,
           " (statuses ", lo$status, ", ", hi$status, ")")
    arm[a + 1L, ] <- c(lo$objval, hi$objval)
    diagnostics[[a + 1L]] <- list(n_respondents = m,
      lower_status = lo$status, upper_status = hi$status)
  }
  list(arm = arm,
    ate = c(lower = arm[2L, "lower"] - arm[1L, "upper"],
            upper = arm[2L, "upper"] - arm[1L, "lower"]),
    finite = TRUE, method = "L-infinity sharp empirical sieve",
    diagnostics = diagnostics,
    note = paste("Sharp for the empirical conditional-normalization sieve;",
                 "calibration sample capped at", max_n, "rows."))
}

l2_oracle_sharp_reference <- function(data, delta, delta_R, delta_A,
                                      delta_M = delta_R,
                                      model = c("net", "separated"),
                                      basis = NULL, control = list(),
                                      return_details = FALSE) {
  needed <- c("p_A1_given_X", "p_C0_given_X", "p_C1_given_X",
    "mu_Y0_given_XC0", "mu_Y1_given_XC0", "Y", "A", "C", "X1", "X2")
  if (!all(needed %in% names(data)))
    stop("Oracle sharp bounds require DGP columns from simulate_l2_data()")
  model <- match.arg(model)
  delta <- l2_as_pair(delta, "delta", 1)
  delta_R <- l2_as_pair(delta_R, "delta_R")
  delta_A <- l2_as_pair(delta_A, "delta_A")
  delta_M <- l2_as_pair(delta_M, "delta_M")
  B <- l2_basis(data, c("X1", "X2"), basis)
  arm <- vector("list", 2L)
  nuisance <- list(
    e1 = data$p_A1_given_X, e0 = 1 - data$p_A1_given_X,
    rho0 = 1 - data$p_C0_given_X, rho1 = 1 - data$p_C1_given_X,
    mu0 = data$mu_Y0_given_XC0, mu1 = data$mu_Y1_given_XC0)
  for (a in 0:1) {
    r <- l2_reference_quantities(data$Y, data$A, data$C, a, nuisance)
    b <- delta[a + 1L] * r$pi
    lo <- l2_dual_arm(r, B, b, delta_R[a + 1L], delta_A[a + 1L],
      delta_M[a + 1L], model, -1, control)
    hi <- l2_dual_arm(r, B, b, delta_R[a + 1L], delta_A[a + 1L],
      delta_M[a + 1L], model, 1, control)
    arm[[a + 1L]] <- c(lower = -lo$value, upper = hi$value)
  }
  ate <- c(lower = arm[[2L]][["lower"]] - arm[[1L]][["upper"]],
    upper = arm[[2L]][["upper"]] - arm[[1L]][["lower"]])
  if (!isTRUE(return_details)) return(ate)
  arm_df <- data.frame(arm = 0:1, method = "oracle_sharp_sieve",
    lower = vapply(arm, `[[`, numeric(1), "lower"),
    upper = vapply(arm, `[[`, numeric(1), "upper"))
  list(arm = arm_df, ate = data.frame(method = "oracle_sharp_sieve",
    lower = ate[["lower"]], upper = ate[["upper"]]))
}
