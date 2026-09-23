#' Cross-fitted plug-in and EIF estimators of the L2 CS outer bounds
#'
#' @inheritParams l2_sharp_crossfit
#' When `nuisance_source = "oracle_error"`, treatment, response, outcome-
#' regression, and conditional residual-moment nuisances are directly
#' perturbed. No regressions or cross-fitting are used in that simulation mode.
#'
#' @return Arm and ATE CS endpoints with plug-in and EIF standard errors. For
#'   data returned by `simulate_l2_data()`, an additional `true_eif` diagnostic
#'   evaluates the EIF one-step estimator using the known nuisance functions.
#' @export
l2_cs_crossfit <- function(data, Y, A, C, X, delta = c(1, 1),
                           delta_R = c(.2, .2), delta_A = c(.2, .2),
                           delta_M = delta_R,
                           model = c("net", "separated"),
                           folds = 5L,
                           nuisance_method = c("SuperLearner", "glm"),
                           sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm",
                           sl_lib_outcome = "SL.glm", seed = 1L,
                           conf_level = .95,
                           nuisance_source = c("estimated", "oracle_error"),
                           nuisance_error_rate = .25,
                           nuisance_error_scale = 1,
                           nuisance_error_sign = c(e = 1, rho = -1,
                                                   regression = 1),
                           delta_K = NULL) {
  nuisance_method <- match.arg(nuisance_method)
  nuisance_source <- match.arg(nuisance_source)
  model <- match.arg(model)
  if (length(conf_level) != 1L || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1)
    stop("conf_level must lie strictly between zero and one")
  delta <- l2_as_pair(delta, "delta", 1)
  delta_R <- l2_as_pair(delta_R, "delta_R")
  if (!is.null(delta_K)) delta_A <- delta_K
  delta_A <- l2_as_pair(delta_A, "delta_K")
  delta_M <- l2_as_pair(delta_M, "delta_M")
  l2_validate_inputs(data, Y, A, C, X, delta, delta_R, delta_A, delta_M)
  n <- nrow(data)
  controlled_error <- nuisance_source == "oracle_error"
  requested_folds <- as.integer(folds)
  folds <- if (controlled_error) 1L else
    max(2L, min(requested_folds, n))
  set.seed(seed)
  fold_id <- if (controlled_error) rep(1L, n) else
    sample(rep(seq_len(folds), length.out = n))
  zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)
  has_oracle <- all(c("p_A1_given_X", "p_C0_given_X", "p_C1_given_X",
    "mu_Y0_given_XC0", "mu_Y1_given_XC0", "mu_Y0", "mu_Y1", "X1",
    "treatment_lp_X", "confounding_strength") %in% names(data))
  estimators <- c("plugin", "eif", if (has_oracle) "true_eif")
  arm_rows <- vector("list", 4L * length(estimators))
  arm_scores <- vector("list", 4L * length(estimators))
  for (a in 0:1) {
    sc <- lapply(c("psi", "VR", "VA", "W"), function(x)
      stats::setNames(lapply(estimators, function(est) rep(NA_real_, n)), estimators))
    names(sc) <- c("psi", "VR", "VA", "W")
    for (v in seq_len(folds)) {
      if (controlled_error) {
        ## No regressions are fitted in this simulation mode: all nuisance
        ## estimates are directly constructed as DGP functions plus controlled
        ## errors, so sample splitting is neither needed nor used.
        te <- tr <- rep(TRUE, n)
      } else {
        te <- fold_id == v; tr <- !te
      }
      ns <- if (nuisance_source == "oracle_error")
        l2_split_nuisance_oracle_error(data, X, tr, te,
          nuisance_error_rate, nuisance_error_scale, nuisance_error_sign) else
        l2_split_nuisance(data, Y, A, C, X, tr, te, nuisance_method,
          sl_lib_prop, sl_lib_miss, sl_lib_outcome)
      if (has_oracle) ns0 <- l2_split_nuisance_oracle_error(data, X, tr, te,
        0, c(e = 0, rho = 0, regression = 0),
        c(e = 1, rho = 1, regression = 1))
      e <- ns$test[[paste0("e", a)]]; rho <- ns$test[[paste0("rho", a)]]
      mu <- ns$test[[paste0("mu", a)]]
      yy <- data[[Y]][te]; y0 <- yy; y0[!is.finite(y0)] <- 0
      z <- y0 - mu; Rind <- as.numeric(data[[C]][te] == 0)
      Ia <- as.numeric(data[[A]][te] == a); S <- Ia * Rind; g <- e * rho
      if (nuisance_source == "oracle_error") {
        nu <- l2_oracle_cs_moments(data[te, , drop = FALSE], a, mu)
        esc <- l2_nuisance_error_scale(nuisance_error_scale)
        esign <- l2_nuisance_error_sign(nuisance_error_sign)
        emag <- n^(-nuisance_error_rate) * esc[["regression"]] *
          esign[["regression"]]
        for (k in 1:4) nu[[paste0("nu", k)]] <-
          nu[[paste0("nu", k)]] +
          emag * ns$test[[paste0("err_nu", k)]]
        nu$nu2 <- pmax(nu$nu2, 1e-8)
        nu$nu4 <- pmax(nu$nu4, 1e-8)
      } else {
        iao <- data[[A]][tr] == a & data[[C]][tr] == 0
        mutr <- ns$train[[paste0("mu", a)]]
        ztr <- data[[Y]][tr] - mutr
        nu <- list(nu1 = rep(0, sum(te)))
        for (k in 2:4) nu[[paste0("nu", k)]] <- l2_pseudo_regression(
          data[tr, , drop = FALSE], X, ztr[iao]^k, iao,
          data[te, , drop = FALSE], nuisance_method, sl_lib_outcome)
      }
      b <- delta[a + 1L] * (1 - rho); q <- 1 - e
      moment_scores <- function(k, w, we, wr, nuk, nukm1,
                                e_use = e, rho_use = rho, z_use = z) {
        g_use <- e_use * rho_use
        plugin <- S / g_use * w * z_use^k
        eif <- w * nuk + S / g_use * w *
          (z_use^k - nuk - k * nukm1 * z_use) +
          we * nuk * (Ia - e_use) +
          wr * nuk * Ia / e_use * (Rind - rho_use)
        list(plugin = plugin, eif = eif)
      }
      sc$psi$plugin[te] <- S / g * y0
      sc$psi$eif[te] <- mu + S / g * z
       tmp <- if (model == "separated")
         moment_scores(2, b^2, 0, -2 * delta[a + 1L] * b,
                       nu$nu2, nu$nu1) else
         moment_scores(2, rep(1, length(b)), 0, 0,
                       nu$nu2, nu$nu1)
      sc$VR$plugin[te] <- tmp$plugin; sc$VR$eif[te] <- tmp$eif
      tmp <- moment_scores(2, q^2, -2 * q, 0, nu$nu2, nu$nu1)
      sc$VA$plugin[te] <- tmp$plugin; sc$VA$eif[te] <- tmp$eif
       tmp <- if (model == "separated")
         moment_scores(4, q^4 * b^2, -4 * q^3 * b^2,
           -2 * delta[a + 1L] * q^4 * b, nu$nu4, nu$nu3) else
         moment_scores(4, q^4, -4 * q^3, 0, nu$nu4, nu$nu3)
      sc$W$plugin[te] <- tmp$plugin; sc$W$eif[te] <- tmp$eif
      if (has_oracle) {
        e0 <- ns0$test[[paste0("e", a)]]
        rho0 <- ns0$test[[paste0("rho", a)]]
        mu0 <- ns0$test[[paste0("mu", a)]]
        z0 <- y0 - mu0
        nu0 <- l2_oracle_cs_moments(data[te, , drop = FALSE], a, mu0,
          nodes = 200L)
        b0 <- delta[a + 1L] * (1 - rho0); q0 <- 1 - e0
        sc$psi$true_eif[te] <- mu0 + S / (e0 * rho0) * z0
        tmp0 <- if (model == "separated")
          moment_scores(2, b0^2, 0, -2 * delta[a + 1L] * b0,
            nu0$nu2, nu0$nu1, e0, rho0, z0) else
          moment_scores(2, rep(1, length(b0)), 0, 0,
            nu0$nu2, nu0$nu1, e0, rho0, z0)
        sc$VR$true_eif[te] <- tmp0$eif
        tmp0 <- moment_scores(2, q0^2, -2 * q0, 0,
          nu0$nu2, nu0$nu1, e0, rho0, z0)
        sc$VA$true_eif[te] <- tmp0$eif
        tmp0 <- if (model == "separated")
          moment_scores(4, q0^4 * b0^2, -4 * q0^3 * b0^2,
            -2 * delta[a + 1L] * q0^4 * b0, nu0$nu4, nu0$nu3,
            e0, rho0, z0) else
          moment_scores(4, q0^4, -4 * q0^3, 0, nu0$nu4, nu0$nu3,
            e0, rho0, z0)
        sc$W$true_eif[te] <- tmp0$eif
      }
    }
    for (est in estimators) {
      val <- vapply(sc, function(x) mean(x[[est]]), numeric(1))
      val[c("VR", "VA", "W")] <- pmax(val[c("VR", "VA", "W")], 0)
      miss_radius <- if (model == "separated")
        delta_R[a + 1L] else delta_M[a + 1L]
      H <- val["VA"] + miss_radius * sqrt(val["W"])
      rad <- miss_radius * sqrt(val["VR"]) +
        delta_A[a + 1L] * sqrt(H)
      ifVR <- if (val["VR"] > 0) miss_radius /
        (2 * sqrt(val["VR"])) else 0
      ifVA <- if (H > 0) delta_A[a + 1L] / (2 * sqrt(H)) else 0
      ifW <- if (H > 0 && val["W"] > 0)
        delta_A[a + 1L] * miss_radius /
          (4 * sqrt(H) * sqrt(val["W"])) else 0
      rad_if <- ifVR * (sc$VR[[est]] - val["VR"]) +
        ifVA * (sc$VA[[est]] - val["VA"]) +
        ifW * (sc$W[[est]] - val["W"])
      psi_if <- sc$psi[[est]] - val["psi"]
      for (endpoint in c("lower", "upper")) {
        sg <- if (endpoint == "lower") -1 else 1
        idx <- length(Filter(Negate(is.null), arm_rows)) + 1L
        rawif <- psi_if + sg * rad_if
        arm_rows[[idx]] <- data.frame(arm = a, endpoint = endpoint,
          estimator = est, estimate = val["psi"] + sg * rad,
          se = stats::sd(rawif) / sqrt(n))
        arm_scores[[idx]] <- rawif
      }
    }
  }
  arm <- do.call(rbind, arm_rows); rownames(arm) <- NULL
  arm$conf_low <- arm$estimate - zcrit * arm$se
  arm$conf_high <- arm$estimate + zcrit * arm$se
  ate <- list(); scores <- list(); j <- 0L
  for (est in estimators) for (ep in c("lower", "upper")) {
    j <- j + 1L; ep0 <- if (ep == "lower") "upper" else "lower"
    i1 <- which(arm$arm == 1 & arm$endpoint == ep & arm$estimator == est)
    i0 <- which(arm$arm == 0 & arm$endpoint == ep0 & arm$estimator == est)
    val <- arm$estimate[i1] - arm$estimate[i0]
    raw <- arm_scores[[i1]] - arm_scores[[i0]]; se <- stats::sd(raw) / sqrt(n)
    ate[[j]] <- data.frame(endpoint = ep, estimator = est, estimate = val,
      se = se, conf_low = val - zcrit * se, conf_high = val + zcrit * se)
    scores[[paste(est, ep, sep = "_")]] <- raw
  }
  ate_df <- do.call(rbind, ate)
  outward <- do.call(rbind, lapply(unique(ate_df$estimator), function(est) {
    lo <- ate_df[ate_df$estimator == est & ate_df$endpoint == "lower", ]
    hi <- ate_df[ate_df$estimator == est & ate_df$endpoint == "upper", ]
    data.frame(estimator = est, lower = lo$conf_low, upper = hi$conf_high)
  }))
  out <- list(arm = arm, ate = ate_df, identified_set_ci = outward,
       scores = scores,
       fold_id = fold_id,
       parameters = list(delta = delta, delta_R = delta_R,
         delta_A = delta_A, delta_K = delta_A,
         delta_M = delta_M, model = model, folds = folds,
         requested_folds = requested_folds,
         sample_splitting = if (controlled_error)
           "none: nuisances directly perturbed" else "cross-fitting"),
       inference_note = paste("EIF intervals use the orthogonal endpoint score;",
         "plug-in intervals condition on fitted nuisances and need not account",
         "for first-order nuisance-estimation error."),
       call = match.call())
  class(out) <- c("jointbounds_l2_cs_cf", "list")
  out
}

l2_oracle_cs_moments <- function(data, a, center, nodes = 100L) {
  q <- l2_dgp_error_quantiles(data, nodes)
  mu <- data[[paste0("mu_Y", a)]]; yy <- outer(mu, q, "+")
  link <- as.character(data$density_ratio_link[1L])
  latent <- function(z) if (link == "bounded") tanh(z) else z
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
  w <- (if (a == 1) pA1 else 1 - pA1) * (1 - pui)
  z <- yy - center
  ans <- list(nu1 = rowSums(w * z) / rowSums(w))
  for (k in 2:4) ans[[paste0("nu", k)]] <- rowSums(w * z^k) / rowSums(w)
  ans
}
