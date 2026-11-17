#' Estimate covariate-omission benchmarks on the L2 radius scale
#'
#' Deliberately hides one or more covariates `W = setdiff(X, Z)` and estimates
#' the missingness and confounding density-ratio radii induced by that omission,
#' under MAR and no unmeasured confounding given the full covariate set `X`.
#' These are descriptive plug-in benchmarks on the same scale as `delta_M` and
#' `delta_K`; the function does not currently provide standard errors.
#'
#' @param data A data frame.
#' @param Y,A,C Column names. `C = 1` denotes a missing outcome.
#' @param X Character vector containing the full covariate set.
#' @param p_z Optional integer vector of retained-set sizes. When supplied and
#'   `subsets` is `NULL`, candidate retained sets `Z` are constructed at each
#'   requested size and the omitted benchmark set is `W = setdiff(X, Z)`.
#' @param max_per_pz Maximum number of retained sets selected at each requested
#'   value of `p_z`.
#' @param subsets Optional list of retained covariate sets `Z`. The omitted set
#'   is `setdiff(X, Z)`. By default, each covariate is omitted one at a time and
#'   all covariates are omitted together.
#' @param subset_selection Reproducible random sampling or the lexicographic
#'   first sets when the number of combinations exceeds `max_per_pz`.
#' @param folds Number of cross-fitting folds.
#' @param seed Random seed used to construct folds.
#' @param probability_clip Clipping threshold for treatment, response, and
#'   joint-membership probability estimates.
#' @param ratio_clip Clipping threshold for probabilities from the two weighted
#'   density-ratio classifiers.
#' @param combine_arms How arm-specific radii are combined. The default `"max"`
#'   matches a sensitivity analysis using a common radius in both arms.
#' @param return_predictions Retain held-out density-ratio predictions and
#'   evaluation weights.
#' @return An object of class `jointbounds_l2_radius_benchmarks`. `arm_results`
#'   contains arm-specific radii and diagnostics; `results` contains one pair
#'   per omitted set after combining arms.
#' @details For each arm, respondents are reweighted to represent the observed
#'   respondent (`O`), full factual (`F`), and cross-treatment (`C`) outcome
#'   distributions after hiding `W`. Weighted logistic classifiers on `(Y,Z)`
#'   estimate `M = dP_F/dP_O` and `K = dP_C/dP_F` through class odds. All
#'   predictions used in the empirical L2 norms are out of fold.
#' @export
l2_radius_benchmarks <- function(data, Y, A, C, X, p_z = NULL,
                                 max_per_pz = 10L, subsets = NULL,
                                 subset_selection = c("random", "first"),
                                 folds = 5L, seed = 1L,
                                 probability_clip = 1e-3,
                                 ratio_clip = 1e-3,
                                 combine_arms = c("max", "mean"),
                                 return_predictions = FALSE) {
  combine_arms <- match.arg(combine_arms)
  subset_selection <- match.arg(subset_selection)
  if (!is.data.frame(data)) stop("data must be a data frame")
  if (!all(c(Y, A, C, X) %in% names(data)))
    stop("Y, A, C, and every X variable must exist in data")
  if (!length(X) || anyDuplicated(X)) stop("X must contain distinct covariate names")
  if (!all(data[[A]] %in% 0:1) || !all(data[[C]] %in% 0:1))
    stop("A and C must be coded 0/1")
  if (any(!is.finite(data[[Y]][data[[C]] == 0])))
    stop("Y must be finite whenever C = 0")
  if (!is.numeric(probability_clip) || length(probability_clip) != 1L ||
      probability_clip <= 0 || probability_clip >= .5)
    stop("probability_clip must lie strictly between zero and one half")
  if (!is.numeric(ratio_clip) || length(ratio_clip) != 1L ||
      ratio_clip <= 0 || ratio_clip >= .5)
    stop("ratio_clip must lie strictly between zero and one half")
  n <- nrow(data); folds <- as.integer(folds)
  if (!is.finite(folds) || folds < 2L || folds > n)
    stop("folds must be between 2 and n")
  max_per_pz <- as.integer(max_per_pz)
  if (!is.finite(max_per_pz) || length(max_per_pz) != 1L || max_per_pz < 1L)
    stop("max_per_pz must be positive")
  if (is.null(subsets)) {
    if (is.null(p_z)) {
      subsets <- lapply(X, function(w) setdiff(X, w))
      names(subsets) <- paste0("omit ", X)
      if (length(X) > 1L) subsets[["omit all"]] <- character()
    } else {
      if (!is.numeric(p_z) || any(!is.finite(p_z)) ||
          any(p_z != as.integer(p_z)) || any(p_z < 0 | p_z >= length(X)))
        stop("p_z must contain integers between 0 and length(X)-1")
      set.seed(seed)
      subsets <- list()
      for (k in unique(as.integer(p_z))) {
        cmb <- if (k == 0L) matrix(character(), nrow = 0L, ncol = 1L) else
          utils::combn(X, k)
        nc <- ncol(cmb)
        take <- seq_len(nc)
        if (nc > max_per_pz) take <- if (subset_selection == "random")
          sort(sample.int(nc, max_per_pz)) else seq_len(max_per_pz)
        for (jj in take) subsets[[length(subsets) + 1L]] <- if (k == 0L)
          character() else cmb[, jj]
      }
    }
  }
  if (!is.list(subsets) || !length(subsets) ||
      any(!vapply(subsets, function(z) all(z %in% X), logical(1))))
    stop("subsets must be a nonempty list of retained subsets of X")
  labels <- names(subsets)
  if (is.null(labels)) labels <- rep("", length(subsets))
  for (j in seq_along(subsets)) {
    subsets[[j]] <- unique(as.character(subsets[[j]]))
    if (!nzchar(labels[j])) {
      w <- setdiff(X, subsets[[j]])
      labels[j] <- if (length(w)) paste0("omit ", paste(w, collapse = "+")) else "omit none"
    }
  }
  if (any(vapply(subsets, function(z) !length(setdiff(X, z)), logical(1))))
    stop("Every retained set must omit at least one variable from X")

  set.seed(seed)
  fold_id <- sample(rep(seq_len(folds), length.out = n))
  arm_rows <- list(); pred_rows <- list(); h <- 0L; hp <- 0L
  for (j in seq_along(subsets)) {
    Z <- subsets[[j]]; W <- setdiff(X, Z)
    accum <- lapply(0:1, function(a) list(num_M = 0, den_M = 0,
      num_K = 0, den_K = 0, sw_O = 0, sw_O2 = 0,
      sw_F = 0, sw_F2 = 0, clip_M = 0, clip_K = 0, n_eval = 0))
    for (v in seq_len(folds)) {
      tr <- fold_id != v; va <- fold_id == v
      fit <- l2_rb_nuisance_fit(data[tr, , drop = FALSE], A, C, X, Z,
                                probability_clip)
      ptr <- l2_rb_nuisance_predict(fit, data[tr, , drop = FALSE],
                                    probability_clip)
      pva <- l2_rb_nuisance_predict(fit, data[va, , drop = FALSE],
                                    probability_clip)
      for (a in 0:1) {
        iva_local <- which(data[[A]][va] == a & data[[C]][va] == 0)
        tr_local <- which(data[[A]][tr] == a & data[[C]][tr] == 0)
        if (length(tr_local) < 5L || length(iva_local) < 1L)
          stop("Too few observed outcomes in arm ", a, " for fold ", v)
        wtr <- l2_rb_weights(ptr, a)
        wva <- l2_rb_weights(pva, a)
        xtr <- data[tr, c(Y, Z), drop = FALSE][tr_local, , drop = FALSE]
        xva <- data[va, c(Y, Z), drop = FALSE][iva_local, , drop = FALSE]
        fit_M <- l2_rb_ratio_fit(xtr, wtr$O[tr_local], wtr$F[tr_local])
        fit_K <- l2_rb_ratio_fit(xtr, wtr$F[tr_local], wtr$C[tr_local])
        pM <- l2_rb_ratio_predict(fit_M, xva, ratio_clip)
        pK <- l2_rb_ratio_predict(fit_K, xva, ratio_clip)
        Mhat <- pM / (1 - pM); Khat <- pK / (1 - pK)
        o <- wva$O[iva_local]; f <- wva$F[iva_local]
        z <- accum[[a + 1L]]
        z$num_M <- z$num_M + sum(o * (Mhat - 1)^2)
        z$den_M <- z$den_M + sum(o)
        z$num_K <- z$num_K + sum(f * (Khat - 1)^2)
        z$den_K <- z$den_K + sum(f)
        z$sw_O <- z$sw_O + sum(o); z$sw_O2 <- z$sw_O2 + sum(o^2)
        z$sw_F <- z$sw_F + sum(f); z$sw_F2 <- z$sw_F2 + sum(f^2)
        z$clip_M <- z$clip_M + sum(pM %in% c(ratio_clip, 1 - ratio_clip))
        z$clip_K <- z$clip_K + sum(pK %in% c(ratio_clip, 1 - ratio_clip))
        z$n_eval <- z$n_eval + length(iva_local)
        accum[[a + 1L]] <- z
        if (isTRUE(return_predictions)) {
          hp <- hp + 1L
          pred_rows[[hp]] <- data.frame(benchmark_id = j, label = labels[j],
            fold = v, row = which(va)[iva_local], arm = a,
            M = Mhat, K = Khat, weight_O = o, weight_F = f)
        }
      }
    }
    for (a in 0:1) {
      z <- accum[[a + 1L]]; h <- h + 1L
      arm_rows[[h]] <- data.frame(benchmark_id = j, label = labels[j],
        retained = paste(Z, collapse = ","), omitted = paste(W, collapse = ","),
        p_z = length(Z), n_omitted = length(W), arm = a,
        delta_M = sqrt(z$num_M / z$den_M),
        delta_K = sqrt(z$num_K / z$den_K),
        ess_O = z$sw_O^2 / z$sw_O2, ess_F = z$sw_F^2 / z$sw_F2,
        ratio_clip_fraction_M = z$clip_M / z$n_eval,
        ratio_clip_fraction_K = z$clip_K / z$n_eval,
        n_respondents = z$n_eval, stringsAsFactors = FALSE)
    }
  }
  arm_results <- do.call(rbind, arm_rows)
  combined <- do.call(rbind, lapply(split(arm_results, arm_results$benchmark_id), function(d) {
    fun <- if (combine_arms == "max") max else mean
    data.frame(benchmark_id = d$benchmark_id[1], label = d$label[1],
      retained = d$retained[1], omitted = d$omitted[1],
      p_z = d$p_z[1], n_omitted = d$n_omitted[1],
      delta_M = fun(d$delta_M),
      delta_K = fun(d$delta_K), min_ess_O = min(d$ess_O),
      min_ess_F = min(d$ess_F), max_ratio_clip_fraction =
        max(d$ratio_clip_fraction_M, d$ratio_clip_fraction_K),
      stringsAsFactors = FALSE)
  }))
  rownames(combined) <- NULL
  if (any(combined$max_ratio_clip_fraction > .05, na.rm = TRUE) ||
      any(combined$min_ess_O < 20 | combined$min_ess_F < 20, na.rm = TRUE))
    warning("Some radius benchmarks have limited overlap: inspect effective ",
      "sample sizes and ratio-clipping fractions in the returned results.",
      call. = FALSE)
  out <- list(results = combined, arm_results = arm_results,
    predictions = if (return_predictions && length(pred_rows)) do.call(rbind, pred_rows) else NULL,
    subsets = subsets, fold_id = fold_id,
    metadata = list(X = X, folds = folds, seed = seed,
      probability_clip = probability_clip, ratio_clip = ratio_clip,
      combine_arms = combine_arms, p_z = p_z,
      max_per_pz = max_per_pz, subset_selection = subset_selection,
      uncertainty = "not estimated"),
    call = match.call())
  class(out) <- "jointbounds_l2_radius_benchmarks"
  out
}

l2_rb_formula <- function(response, predictors) {
  if (!length(predictors)) stats::as.formula(paste(response, "~ 1")) else
    stats::reformulate(predictors, response)
}

l2_rb_safe_glm <- function(formula, data, weights = NULL) {
  # Omitting the weights argument when it is NULL matters: model.frame()
  # otherwise looks for a column named `weights` in some R versions and the
  # helper silently falls back to a marginal fit.
  fit <- suppressWarnings(try(if (is.null(weights)) {
    stats::glm(formula, data = data, family = stats::binomial())
  } else {
    data$.rb_glm_weight <- weights
    stats::glm(formula, data = data, weights = .rb_glm_weight,
      family = stats::binomial())
  }, silent = TRUE))
  if (inherits(fit, "try-error") || any(!is.finite(stats::coef(fit)))) NULL else fit
}

l2_rb_prob_fit <- function(data, outcome, predictors) {
  y <- data[[outcome]]
  if (length(unique(y)) < 2L) return(list(constant = mean(y)))
  fit <- l2_rb_safe_glm(l2_rb_formula(outcome, predictors), data)
  if (is.null(fit)) list(constant = mean(y)) else list(fit = fit)
}

l2_rb_prob_predict <- function(object, newdata, clip) {
  p <- if (!is.null(object$constant)) rep(object$constant, nrow(newdata)) else
    suppressWarnings(stats::predict(object$fit, newdata = newdata, type = "response"))
  clip_probs(as.numeric(p), clip)
}

l2_rb_nuisance_fit <- function(data, A, C, X, Z, clip) {
  data$.rb_R <- 1 - data[[C]]
  data$.rb_A <- data[[A]]
  out <- list(eX = l2_rb_prob_fit(data, ".rb_A", X),
              eZ = l2_rb_prob_fit(data, ".rb_A", Z), rho = list(), gamma = list(),
              A = A, C = C, X = X, Z = Z)
  for (a in 0:1) {
    da <- data[data[[A]] == a, , drop = FALSE]
    out$rho[[a + 1L]] <- l2_rb_prob_fit(da, ".rb_R", X)
    nm <- paste0(".rb_G", a); data[[nm]] <- as.numeric(data[[A]] == a & data$.rb_R == 1)
    out$gamma[[a + 1L]] <- l2_rb_prob_fit(data, nm, Z)
  }
  out
}

l2_rb_nuisance_predict <- function(fit, data, clip) {
  e1X <- l2_rb_prob_predict(fit$eX, data, clip)
  e1Z <- l2_rb_prob_predict(fit$eZ, data, clip)
  data.frame(e0X = 1 - e1X, e1X = e1X, e0Z = 1 - e1Z, e1Z = e1Z,
    rho0 = l2_rb_prob_predict(fit$rho[[1L]], data, clip),
    rho1 = l2_rb_prob_predict(fit$rho[[2L]], data, clip),
    gamma0 = l2_rb_prob_predict(fit$gamma[[1L]], data, clip),
    gamma1 = l2_rb_prob_predict(fit$gamma[[2L]], data, clip))
}

l2_rb_weights <- function(p, a) {
  eaX <- p[[paste0("e", a, "X")]]; eoX <- p[[paste0("e", 1 - a, "X")]]
  eaZ <- p[[paste0("e", a, "Z")]]; eoZ <- p[[paste0("e", 1 - a, "Z")]]
  rho <- p[[paste0("rho", a)]]; gamma <- p[[paste0("gamma", a)]]
  list(O = 1 / gamma, F = 1 / (eaZ * rho),
       C = eoX / (eoZ * eaX * rho))
}

l2_rb_ratio_fit <- function(x, w0, w1) {
  w0 <- w0 / sum(w0); w1 <- w1 / sum(w1)
  d <- rbind(transform(x, .rb_class = 0), transform(x, .rb_class = 1))
  w <- c(w0, w1); predictors <- setdiff(names(d), ".rb_class")
  fit <- l2_rb_safe_glm(l2_rb_formula(".rb_class", predictors), d, w)
  if (is.null(fit)) list(constant = .5) else list(fit = fit)
}

l2_rb_ratio_predict <- function(object, newdata, clip) {
  p <- if (!is.null(object$constant)) rep(object$constant, nrow(newdata)) else
    suppressWarnings(stats::predict(object$fit, newdata = newdata, type = "response"))
  clip_probs(as.numeric(p), clip)
}

#' Plot radius benchmarks against tipping frontiers
#'
#' @param frontier A tipping-frontier object, or a named list of objects. Each
#'   object may come from `l2_tipping_frontier()` or `l2_tipping_band()`.
#' @param benchmarks Output from `l2_radius_benchmarks()`.
#' @param label_points Label benchmark points.
#' @param uncertainty Show `"simultaneous"`, `"pointwise"`, `"both"`, or no
#'   frontier uncertainty. The default is simultaneous for band objects and
#'   none for point-estimate frontier objects.
#' @param interpolate Connect finite frontier estimates and bands between
#'   evaluated horizontal grid values. Such connections are descriptive;
#'   finite-grid coverage does not imply coverage between grid points.
#' @param diagonal Add the equal-radius reference line.
#' @param point_size Size of omission-benchmark points.
#' @param ... Reserved for future plotting options.
#' @return A ggplot object.
#' @export
l2_plot_radius_benchmarks <- function(frontier, benchmarks,
                                      uncertainty = NULL,
                                      interpolate = TRUE, diagonal = TRUE,
                                      label_points = FALSE, point_size = 3,
                                      ...) {
  if (!inherits(benchmarks, "jointbounds_l2_radius_benchmarks"))
    stop("benchmarks must be returned by l2_radius_benchmarks()")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Plotting requires the suggested package ggplot2")
  recognized <- function(z) inherits(z, "jointbounds_l2_calibration_band") ||
    inherits(z, "jointbounds_l2_frontier")
  objects <- if (recognized(frontier)) list(frontier) else frontier
  if (!is.list(objects) || !length(objects) ||
      any(!vapply(objects, recognized, logical(1))))
    stop("frontier must be a tipping-frontier object or a list of them")
  object_names <- names(objects)
  if (is.null(object_names)) object_names <- rep("", length(objects))
  rows <- vector("list", length(objects)); has_band <- logical(length(objects))
  for (j in seq_along(objects)) {
    z <- objects[[j]]
    if (inherits(z, "jointbounds_l2_calibration_band")) {
      d <- z$frontier; rx <- z$radius_x; ry <- z$radius_y
      has_band[j] <- TRUE
    } else {
      d <- as.data.frame(z); rx <- attr(z, "radius_x"); ry <- attr(z, "radius_y")
      if (is.null(rx)) rx <- if ("delta_M" %in% names(d)) "delta_M" else "delta_R"
      if (is.null(ry)) ry <- if ("delta_K" %in% names(d)) "delta_K" else "delta_A"
      d$estimate <- d[[ry]]
    }
    if (!identical(rx, "delta_M") || !identical(ry, "delta_K"))
      stop("Radius benchmarks require frontiers on delta_M and delta_K axes")
    if (!nzchar(object_names[j])) object_names[j] <- if (length(objects) == 1L)
      "Tipping frontier" else if (!is.null(z$endpoint_estimator) &&
        identical(z$endpoint_estimator, "dr")) "Sharp sieve" else if
        (!is.null(z$endpoint_estimator) && identical(z$endpoint_estimator, "eif"))
          "CS outer" else paste("Frontier", j)
    d$.bound <- object_names[j]
    rows[[j]] <- d
  }
  if (is.null(uncertainty)) uncertainty <- if (all(has_band))
    "simultaneous" else "none"
  uncertainty <- match.arg(uncertainty,
    c("simultaneous", "pointwise", "both", "none"))
  if (uncertainty != "none" && any(!has_band))
    stop("Requested uncertainty requires every object to come from l2_tipping_band()")
  f <- do.call(rbind, rows)
  f <- f[order(f$.bound, f$delta_M), , drop = FALSE]
  f$.missing <- !is.finite(f$estimate)
  f$.segment <- stats::ave(as.integer(f$.missing), f$.bound, FUN = cumsum)
  f$.line_group <- interaction(f$.bound, f$.segment, drop = TRUE)
  finite <- f[!f$.missing, , drop = FALSE]
  p <- ggplot2::ggplot(f, ggplot2::aes(x = delta_M))
  if (isTRUE(diagonal)) p <- p + ggplot2::geom_abline(slope = 1,
    intercept = 0, colour = "grey70", linewidth = .55)
  kinds <- if (uncertainty == "both") c("simultaneous", "pointwise") else
    if (uncertainty == "none") character() else uncertainty
  for (kind in kinds) {
    low <- paste0(kind, "_low"); high <- paste0(kind, "_high")
    z <- f
    z$.low <- z[[low]]; z$.high <- z[[high]]
    z$.band_missing <- !is.finite(z$.low) | !is.finite(z$.high)
    z$.band_segment <- stats::ave(as.integer(z$.band_missing), z$.bound,
      FUN = cumsum)
    z$.band_group <- interaction(z$.bound, z$.band_segment, drop = TRUE)
    z <- z[!z$.band_missing, , drop = FALSE]
    if (isTRUE(interpolate)) p <- p + ggplot2::geom_ribbon(data = z,
      ggplot2::aes(x = delta_M, ymin = .low, ymax = .high, fill = .bound,
                   group = .band_group), inherit.aes = FALSE,
      alpha = if (kind == "simultaneous") .14 else .26, colour = NA)
    else p <- p + ggplot2::geom_errorbar(data = z,
      ggplot2::aes(x = delta_M, ymin = .low, ymax = .high,
                   group = .bound), inherit.aes = FALSE, colour = "black",
      width = 0, alpha = if (kind == "simultaneous") .65 else 1)
  }
  if (isTRUE(interpolate)) p <- p + ggplot2::geom_line(data = finite,
    ggplot2::aes(x = delta_M, y = estimate, linetype = .bound,
                 group = .line_group),
    colour = "black", linewidth = .9, inherit.aes = FALSE)
  else p <- p + ggplot2::geom_point(data = finite,
    ggplot2::aes(x = delta_M, y = estimate, shape = .bound),
    colour = "black", inherit.aes = FALSE)
  b <- benchmarks$results
  p <- p + ggplot2::geom_point(data = b,
    ggplot2::aes(x = delta_M, y = delta_K, colour = label),
    inherit.aes = FALSE, size = point_size)
  if (isTRUE(label_points)) p <- p + ggplot2::geom_text(data = b,
    ggplot2::aes(x = delta_M, y = delta_K, label = label, colour = label),
    inherit.aes = FALSE, nudge_y = .02, check_overlap = TRUE,
    show.legend = FALSE, size = 3.5)
  caption <- paste0("Omission-radius points are descriptive estimates; their ",
    "uncertainty is not shown.", if (isTRUE(interpolate))
      " Lines and ribbons between evaluated grid values are descriptive." else "")
  p + ggplot2::labs(x = expression(delta[M]), y = expression(delta[K]),
    fill = "Tipping frontier", linetype = "Tipping frontier",
    shape = "Tipping frontier", colour = "Omitted covariates",
    caption = caption) + ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom", legend.box = "vertical")
}

#' Plot omission-radius benchmarks
#'
#' @param x Output from `l2_radius_benchmarks()`.
#' @param frontier A tipping-frontier object or named list of objects.
#' @param ... Additional arguments passed to `l2_plot_radius_benchmarks()`.
#' @return A ggplot object.
#' @export
plot.jointbounds_l2_radius_benchmarks <- function(x, frontier, ...) {
  l2_plot_radius_benchmarks(frontier = frontier, benchmarks = x, ...)
}
