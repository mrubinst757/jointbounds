#' Cross-fitted doubly robust ATE under MAR and no unmeasured confounding
#'
#' Estimates each potential-outcome mean with the cross-fitted augmented
#' inverse-probability weighted score
#' \deqn{\mu_a(X) + \frac{I(A=a,C=0)}{e_a(X)\rho_a(X)}
#'       \{Y-\mu_a(X)\},}
#' where `C = 1` denotes a missing outcome. The ATE is the difference between
#' the two arm means. Under MAR, consistency, positivity, and no unmeasured
#' confounding, the estimator is consistent if either both outcome regressions
#' are consistently estimated or both joint observation scores
#' `e_a(X) * rho_a(X)` are consistently estimated.
#'
#' @param data A data frame.
#' @param Y,A,C Column names. `A` and `C` must be binary and `C = 1`
#'   denotes a missing outcome.
#' @param X Character vector of baseline covariate names. Use `character()` for
#'   an intercept-only analysis.
#' @param folds Number of cross-fitting folds.
#' @param nuisance_method Either `"SuperLearner"` or `"glm"`.
#' @param sl_lib_prop,sl_lib_miss,sl_lib_outcome SuperLearner libraries for the
#'   treatment, observation, and outcome regressions.
#' @param outcome_family Either `"gaussian"` or `"binomial"`.
#' @param seed Random seed used to construct folds.
#' @param trim Lower and upper truncation level for estimated treatment and
#'   response probabilities.
#' @param conf_level Wald confidence level.
#' @param diagnostic_control Optional list with entries `warn`, `max_weight`,
#'   `min_ess`, and `min_ess_fraction` controlling positivity diagnostics and
#'   the consolidated warning.
#' @param keep_nuisance Return row-level cross-fitted nuisance estimates.
#'
#' @return An object of class `jointbounds_ate`. Its `ate` component contains the
#'   cross-fitted AIPW estimate, standard error, and Wald interval; `arm`
#'   contains the corresponding potential-outcome means. Centered influence
#'   scores and positivity diagnostics are also returned.
#' @export
generate_ate <- function(data, Y, A, C, X,
                         folds = 5L,
                         nuisance_method = c("SuperLearner", "glm"),
                         sl_lib_prop = "SL.glm",
                         sl_lib_miss = "SL.glm",
                         sl_lib_outcome = "SL.glm",
                         outcome_family = c("gaussian", "binomial"),
                         seed = 1L,
                         trim = 1e-4,
                         conf_level = .95,
                         diagnostic_control = list(),
                         keep_nuisance = FALSE) {
  nuisance_method <- match.arg(nuisance_method)
  outcome_family <- match.arg(outcome_family)
  if (!is.data.frame(data)) stop("data must be a data.frame")
  if (!is.character(X)) stop("X must be a character vector")
  if (!all(c(Y, A, C, X) %in% names(data)))
    stop("All Y, A, C, and X columns must exist")
  if (length(Y) != 1L || length(A) != 1L || length(C) != 1L)
    stop("Y, A, and C must each name one column")
  aa <- data[[A]]; cc <- data[[C]]; yy <- data[[Y]]
  if (anyNA(aa) || !all(aa %in% 0:1)) stop("A must contain only 0 and 1")
  if (anyNA(cc) || !all(cc %in% 0:1)) stop("C must contain only 0 and 1")
  if (!is.numeric(yy)) stop("Y must be numeric")
  if (any(!is.finite(yy[cc == 0])))
    stop("Observed outcomes (C = 0) must be finite")
  n <- nrow(data)
  folds <- as.integer(folds)
  if (length(folds) != 1L || !is.finite(folds) || folds < 2L || folds > n)
    stop("folds must be an integer between 2 and nrow(data)")
  if (length(trim) != 1L || !is.finite(trim) || trim <= 0 || trim >= .5)
    stop("trim must lie strictly between zero and one half")
  if (length(conf_level) != 1L || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1)
    stop("conf_level must lie strictly between zero and one")
  dc <- ate_diagnostic_control(diagnostic_control)

  fit_data <- data
  fit_X <- X
  if (!length(fit_X) && nuisance_method == "SuperLearner") {
    intercept_name <- ".jointbounds_intercept"
    while (intercept_name %in% names(fit_data))
      intercept_name <- paste0(intercept_name, "_")
    fit_data[[intercept_name]] <- 0
    fit_X <- intercept_name
  }
  nuisance <- if (nuisance_method == "SuperLearner") {
    l2_nuisance_sl(fit_data, Y, A, C, fit_X, folds, seed,
      sl_lib_prop, sl_lib_miss, sl_lib_outcome,
      family_Y = outcome_family)
  } else {
    family_object <- if (outcome_family == "gaussian")
      stats::gaussian() else stats::binomial()
    l2_nuisance_glm(fit_data, Y, A, C, fit_X, folds, seed,
                    outcome_family = family_object)
  }
  for (nm in c("e0", "e1", "rho0", "rho1"))
    nuisance[[nm]] <- clip_probs(nuisance[[nm]], trim)
  needed <- c("e0", "e1", "rho0", "rho1", "mu0", "mu1")
  if (any(!vapply(nuisance[needed], function(z)
    length(z) == n && all(is.finite(z)), logical(1))))
    stop("Nuisance estimation returned nonfinite or incorrectly sized predictions")

  y0 <- yy
  y0[!is.finite(y0)] <- 0
  zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)
  arm_scores <- vector("list", 2L)
  arm_rows <- vector("list", 2L)
  diagnostic_rows <- vector("list", 2L)
  for (a in 0:1) {
    e <- nuisance[[paste0("e", a)]]
    rho <- nuisance[[paste0("rho", a)]]
    mu <- nuisance[[paste0("mu", a)]]
    observed <- as.numeric(aa == a & cc == 0)
    weight <- observed / (e * rho)
    score <- mu + weight * (y0 - mu)
    estimate <- mean(score)
    influence <- score - estimate
    se <- stats::sd(influence) / sqrt(n)
    arm_scores[[a + 1L]] <- influence
    arm_rows[[a + 1L]] <- data.frame(
      arm = a, estimator = "cross-fitted AIPW", estimate = estimate,
      plugin = mean(mu), se = se,
      conf_low = estimate - zcrit * se,
      conf_high = estimate + zcrit * se
    )
    positive_weight <- weight[observed == 1]
    ess <- if (length(positive_weight) && sum(positive_weight^2) > 0)
      sum(positive_weight)^2 / sum(positive_weight^2) else 0
    codes <- character()
    if (sum(observed) < max(10L, 2L * folds))
      codes <- c(codes, "few_observed_outcomes")
    if (min(c(e, rho)) <= trim * (1 + 1e-8) ||
        max(c(e, rho)) >= 1 - trim * (1 + 1e-8))
      codes <- c(codes, "probability_trimming")
    if (length(positive_weight) && max(positive_weight) > dc$max_weight)
      codes <- c(codes, "large_inverse_probability_weight")
    if (ess < max(dc$min_ess, dc$min_ess_fraction * n))
      codes <- c(codes, "low_weighted_effective_sample_size")
    diagnostic_rows[[a + 1L]] <- data.frame(
      arm = a, diagnostic_ok = !length(codes),
      diagnostic_codes = paste(unique(codes), collapse = ";"),
      n_observed_outcomes = sum(observed),
      min_treatment_probability = min(e),
      max_treatment_probability = max(e),
      min_response_probability = min(rho),
      max_response_probability = max(rho),
      maximum_inverse_probability_weight = if (length(positive_weight))
        max(positive_weight) else Inf,
      weighted_effective_sample_size = ess,
      stringsAsFactors = FALSE
    )
  }
  arm <- do.call(rbind, arm_rows)
  ate_influence <- arm_scores[[2L]] - arm_scores[[1L]]
  ate_estimate <- arm$estimate[arm$arm == 1] - arm$estimate[arm$arm == 0]
  ate_se <- stats::sd(ate_influence) / sqrt(n)
  ate <- data.frame(
    estimand = "ATE", estimator = "cross-fitted AIPW",
    estimate = ate_estimate,
    plugin = arm$plugin[arm$arm == 1] - arm$plugin[arm$arm == 0],
    se = ate_se, conf_low = ate_estimate - zcrit * ate_se,
    conf_high = ate_estimate + zcrit * ate_se
  )
  diagnostic_table <- do.call(rbind, diagnostic_rows)
  rownames(arm) <- NULL
  rownames(diagnostic_table) <- NULL
  bad <- !diagnostic_table$diagnostic_ok
  diagnostic_summary <- data.frame(
    diagnostic_ok = !any(bad), n_arms_flagged = sum(bad),
    flagged_codes = paste(unique(unlist(strsplit(
      diagnostic_table$diagnostic_codes[bad], ";", fixed = TRUE))),
      collapse = ";"), stringsAsFactors = FALSE
  )
  out <- list(
    ate = ate, arm = arm,
    influence_function = list(
      arm0 = arm_scores[[1L]], arm1 = arm_scores[[2L]], ate = ate_influence),
    diagnostic_table = diagnostic_table,
    diagnostic_summary = diagnostic_summary,
    assumptions = paste("Identification requires consistency, positivity,",
      "MAR, and no unmeasured confounding. The AIPW estimator is consistent",
      "if the outcome regressions are correct or if the arm-specific joint",
      "treatment-response scores e_a(X) rho_a(X) are correct."),
    parameters = list(folds = folds, nuisance_method = nuisance_method,
      outcome_family = outcome_family, trim = trim,
      conf_level = conf_level),
    call = match.call()
  )
  if (isTRUE(keep_nuisance)) out$nuisance <- nuisance[needed]
  class(out) <- "jointbounds_ate"
  if (isTRUE(dc$warn) && any(bad)) warning(sprintf(
    "ATE diagnostics flagged %d arm%s. Inspect $diagnostic_summary and $diagnostic_table. Codes: %s",
    sum(bad), if (sum(bad) == 1L) "" else "s",
    diagnostic_summary$flagged_codes), call. = FALSE)
  out
}

ate_diagnostic_control <- function(x = list()) {
  if (!is.list(x)) stop("diagnostic_control must be a list")
  defaults <- list(warn = TRUE, max_weight = 100,
                   min_ess = 20, min_ess_fraction = .02)
  unknown <- setdiff(names(x), names(defaults))
  if (length(unknown)) stop("Unknown diagnostic_control entr",
    if (length(unknown) == 1L) "y: " else "ies: ",
    paste(unknown, collapse = ", "))
  out <- utils::modifyList(defaults, x)
  if (!is.logical(out$warn) || length(out$warn) != 1L || is.na(out$warn))
    stop("diagnostic_control$warn must be TRUE or FALSE")
  for (nm in c("max_weight", "min_ess"))
    if (length(out[[nm]]) != 1L || !is.finite(out[[nm]]) || out[[nm]] <= 0)
      stop("diagnostic_control$", nm, " must be one positive number")
  if (length(out$min_ess_fraction) != 1L ||
      !is.finite(out$min_ess_fraction) || out$min_ess_fraction < 0 ||
      out$min_ess_fraction > 1)
    stop("diagnostic_control$min_ess_fraction must lie between zero and one")
  out
}

#' @export
print.jointbounds_ate <- function(x, ...) {
  cat("Cross-fitted AIPW estimate under MAR and no unmeasured confounding\n\n")
  print(x$ate, row.names = FALSE)
  cat("\nDiagnostics:", if (isTRUE(x$diagnostic_summary$diagnostic_ok))
    "all checks passed" else paste(x$diagnostic_summary$n_arms_flagged,
      "arms flagged"), "\n")
  invisible(x)
}
