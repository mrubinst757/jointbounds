#' MAR-only L2 prevalence bounds for binary alternative estimands
#'
#' Implements the capped L2 support-function bounds for the binary composite
#' effect and separable direct effect.  These bounds maintain conditional
#' exchangeability and allow informative missingness through
#' `q_a(X) = pi_a^star(X) / pi_a(X)`, with `0 <= q_a <= 1` and
#' `||q_a||_2 <= kappa_a`.
#'
#' @inheritParams l2_bounds
#' @param estimand Either `"composite"` or `"sde"`.
#' @param kappa One value or arm-specific length-two L2 prevalence radii.  The
#'   SDE uses the arm-zero radius.
#' @param conf_level Wald confidence level.
#' @param keep_nuisance Retain nuisance estimates and empirical support-function
#'   optimizers.
#' @return Bounds, standard errors, Wald intervals, influence scores, and
#'   support-function diagnostics.
#' @export
l2_prevalence_bounds <- function(data, Y, A, C, X,
                                 estimand = c("composite", "sde"),
                                 kappa = c(0, 0), folds = 5L,
                                 nuisance_method = c("SuperLearner", "glm"),
                                 sl_lib_prop = "SL.glm",
                                 sl_lib_miss = "SL.glm",
                                 sl_lib_outcome = "SL.glm", seed = 1L,
                                 conf_level = .95, keep_nuisance = FALSE) {
  estimand <- match.arg(estimand)
  nuisance_method <- match.arg(nuisance_method)
  kappa <- l2_as_pair(kappa, "kappa", 1)
  inp <- l2_validate_inputs(data, Y, A, C, X, c(1, 1), kappa, kappa, kappa)
  observed_y <- inp$y[inp$C == 0]
  if (!all(observed_y %in% c(0, 1)))
    stop("l2_prevalence_bounds requires a binary outcome coded 0/1")
  if (length(conf_level) != 1L || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1)
    stop("conf_level must lie strictly between zero and one")
  nuisance <- if (nuisance_method == "SuperLearner") {
    l2_nuisance_sl(data, Y, A, C, X, folds, seed, sl_lib_prop,
      sl_lib_miss, sl_lib_outcome, family_Y = "binomial")
  } else l2_nuisance_glm(data, Y, A, C, X, folds, seed,
                         outcome_family = stats::binomial())
  for (a in 0:1) nuisance[[paste0("mu", a)]] <-
    pmin(1, pmax(0, nuisance[[paste0("mu", a)]]))
  n <- nrow(data); y <- inp$y; y[!is.finite(y)] <- 0
  Rind <- as.numeric(inp$C == 0)
  arm_score <- vector("list", 2L)
  for (a in 0:1) {
    e <- nuisance[[paste0("e", a)]]
    rho <- nuisance[[paste0("rho", a)]]
    mu <- nuisance[[paste0("mu", a)]]
    S <- as.numeric(inp$A == a) * Rind
    arm_score[[a + 1L]] <- mu + S / (e * rho) * (y - mu)
  }
  reference_score <- arm_score[[2L]] - arm_score[[1L]]
  reference <- mean(reference_score)

  support <- list(); endpoint_score <- list()
  if (estimand == "composite") {
    for (a in 0:1) {
      e <- nuisance[[paste0("e", a)]]
      rho <- nuisance[[paste0("rho", a)]]
      mu <- nuisance[[paste0("mu", a)]]
      pi <- 1 - rho; Ia <- as.numeric(inp$A == a)
      S <- Ia * Rind
      g <- pi * (1 - mu)
      opt <- l2_capped_support(g, kappa[a + 1L])
      q <- opt$q
      fixed_q_score <- q * pi * (1 - mu) -
        q * pi * S / (e * rho) * (y - mu) -
        q * (1 - mu) * Ia / e * (Rind - rho)
      opt$score <- fixed_q_score - opt$lambda *
        (q^2 - kappa[a + 1L]^2)
      support[[paste0("arm", a)]] <- opt
    }
    endpoint_score$lower <- reference_score - support$arm0$score
    endpoint_score$upper <- reference_score + support$arm1$score
  } else {
    e0 <- nuisance$e0; e1 <- nuisance$e1
    rho0 <- nuisance$rho0; rho1 <- nuisance$rho1
    mu0 <- nuisance$mu0; mu1 <- nuisance$mu1
    I0 <- as.numeric(inp$A == 0); I1 <- as.numeric(inp$A == 1)
    S0 <- I0 * Rind; S1 <- I1 * Rind
    pi0 <- 1 - rho0; d <- mu1 - mu0
    sde_support_score <- function(sign) {
      ## sign=1 uses (mu1-mu0)_+; sign=-1 uses (mu0-mu1)_+.
      contrast <- sign * d
      positive <- as.numeric(contrast > 0)
      g <- pi0 * pmax(contrast, 0)
      opt <- l2_capped_support(g, kappa[1L]); q <- opt$q
      outcome_score <- sign * positive *
        (S1 / (e1 * rho1) * (y - mu1) -
         S0 / (e0 * rho0) * (y - mu0))
      fixed_q_score <- q * g + q * pi0 * outcome_score -
        q * pmax(contrast, 0) * I0 / e0 * (Rind - rho0)
      opt$score <- fixed_q_score - opt$lambda * (q^2 - kappa[1L]^2)
      opt
    }
    support$positive <- sde_support_score(1)
    support$negative <- sde_support_score(-1)
    endpoint_score$lower <- reference_score - support$positive$score
    endpoint_score$upper <- reference_score + support$negative$score
  }
  zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)
  bounds <- do.call(rbind, lapply(c("lower", "upper"), function(ep) {
    sc <- endpoint_score[[ep]]; est <- mean(sc); se <- stats::sd(sc) / sqrt(n)
    data.frame(endpoint = ep, estimate = est, se = se,
      conf_low = est - zcrit * se, conf_high = est + zcrit * se)
  }))
  scores <- lapply(endpoint_score, function(z) z - mean(z))
  out <- list(bounds = bounds, reference = reference, scores = scores,
    support = lapply(support, function(z) z[setdiff(names(z), c("q", "score"))]),
    parameters = list(estimand = estimand, kappa = kappa), call = match.call())
  if (isTRUE(keep_nuisance)) {
    out$nuisance <- nuisance
    out$optimizers <- lapply(support, `[[`, "q")
  }
  class(out) <- "jointbounds_l2_prevalence"
  out
}

l2_capped_support <- function(g, kappa, tol = 1e-10) {
  if (any(!is.finite(g)) || any(g < 0))
    stop("The capped support function requires a finite nonnegative loading")
  if (kappa <= tol || !any(g > 0))
    return(list(value = 0, q = rep(0, length(g)), lambda = 0,
      norm = 0, constraint_active = kappa <= tol, cap_active = FALSE))
  q_full <- as.numeric(g > 0)
  if (mean(q_full^2) <= kappa^2 + tol)
    return(list(value = mean(g * q_full), q = q_full, lambda = 0,
      norm = sqrt(mean(q_full^2)), constraint_active = FALSE,
      cap_active = TRUE))
  norm_equation <- function(c) mean(pmin(1, c * g)^2) - kappa^2
  upper <- 1 / max(g)
  while (norm_equation(upper) < 0) upper <- upper * 2
  cstar <- stats::uniroot(norm_equation, c(0, upper), tol = tol)$root
  q <- pmin(1, cstar * g)
  list(value = mean(g * q), q = q, lambda = 1 / (2 * cstar),
    norm = sqrt(mean(q^2)), constraint_active = TRUE,
    cap_active = any(q >= 1 - sqrt(tol)), multiplier_c = cstar)
}
