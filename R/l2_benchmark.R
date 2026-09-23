utils::globalVariables(c(".radius_x", "estimate", "simultaneous_low",
  "simultaneous_high", "pointwise_low", "pointwise_high",
  "equal_radius", ".n_omitted_plot", "overall_radius_low",
  "overall_radius_high", "overall_equal_radius", "benchmark_radius_low",
  "benchmark_radius_high", ".radius_R", ".radius_A", "label", "method"))

#' Estimate many leave-covariate-out calibration benchmarks
#'
#' Constructs retained covariate sets Z, capped within each requested size,
#' and compares the reduced-Z and full-X ATEs. Both outcome-regression plug-in
#' and cross-fitted AIPW estimates are returned. The paired reduced-minus-full
#' discrepancy is the combined MAR/NUC benchmark used to invert a sensitivity
#' surface; it does not identify a mechanism-specific decomposition.
#'
#' @param data A data frame.
#' @param Y,A,C Column names; C=1 denotes a missing outcome.
#' @param X Character vector containing the full covariate set.
#' @param p_z Integer vector of retained-set sizes. Defaults to 0 through p_x-1.
#' @param max_per_pz Maximum number of retained sets at each p_z.
#' @param subsets Optional named or unnamed list of explicit retained sets.
#' @param subset_selection Reproducible random sampling or lexicographic first
#'   subsets when the number of combinations exceeds max_per_pz.
#' @param estimator Return AIPW, plug-in, or both.
#' @param variance Include paired influence-score standard errors and intervals.
#' @inheritParams l2_bounds
#' @return An object of class jointbounds_l2_benchmarks with a results data frame,
#'   selected subsets, row-level centered influence scores when requested, and
#'   metadata. The signed benchmark is reduced-Z minus full-X, so it can be
#'   added directly to a full-X reference estimate during frontier inversion.
#' @export
l2_aipw_benchmarks <- function(data, Y, A, C, X, p_z = NULL,
                               max_per_pz = 10L, subsets = NULL,
                               subset_selection = c("random", "first"),
                               estimator = c("both", "aipw", "plugin"),
                               variance = FALSE, conf_level = .95,
                               folds = 5L,
                               nuisance_method = c("SuperLearner", "glm"),
                               sl_lib_prop = "SL.glm",
                               sl_lib_miss = "SL.glm",
                               sl_lib_outcome = "SL.glm", seed = 1L) {
  subset_selection <- match.arg(subset_selection)
  estimator <- match.arg(estimator); nuisance_method <- match.arg(nuisance_method)
  if (!all(c(Y, A, C, X) %in% names(data))) stop("Y, A, C, and X must exist")
  px <- length(X); if (!px) stop("X must contain at least one covariate")
  max_per_pz <- as.integer(max_per_pz)
  if (!is.finite(max_per_pz) || max_per_pz < 1L) stop("max_per_pz must be positive")
  if (is.null(subsets)) {
    if (is.null(p_z)) p_z <- 0:(px - 1L)
    if (any(p_z < 0 | p_z >= px)) stop("p_z must lie between 0 and length(X)-1")
    set.seed(seed)
    subsets <- list()
    for (k in unique(as.integer(p_z))) {
      cmb <- if (k == 0L) matrix(character(), nrow = 0L, ncol = 1L) else
        utils::combn(X, k)
      nc <- ncol(cmb); take <- seq_len(nc)
      if (nc > max_per_pz) take <- if (subset_selection == "random")
        sort(sample.int(nc, max_per_pz)) else seq_len(max_per_pz)
      for (j in take) subsets[[length(subsets) + 1L]] <- if (k == 0L)
        character() else cmb[, j]
    }
  }
  if (!is.list(subsets) || !length(subsets)) stop("subsets must be a nonempty list")
  if (any(!vapply(subsets, function(z) all(z %in% X), logical(1))))
    stop("Every retained set must be a subset of X")
  labels <- names(subsets)
  if (is.null(labels)) labels <- rep("", length(subsets))
  for (j in seq_along(subsets)) if (!nzchar(labels[j])) labels[j] <-
    if (length(subsets[[j]])) paste(subsets[[j]], collapse = "+") else "(intercept)"
  full <- l2_benchmark_ate_fit(data, Y, A, C, X, folds, nuisance_method,
    sl_lib_prop, sl_lib_miss, sl_lib_outcome, seed)
  methods <- if (estimator == "both") c("plugin", "aipw") else estimator
  rows <- list(); scores <- list(); target_scores <- list(); h <- 0L
  zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)
  for (j in seq_along(subsets)) {
    red <- l2_benchmark_ate_fit(data, Y, A, C, subsets[[j]], folds,
      nuisance_method, sl_lib_prop, sl_lib_miss, sl_lib_outcome, seed)
    omitted <- setdiff(X, subsets[[j]])
    for (meth in methods) {
      h <- h + 1L; sc <- red$scores[[meth]] - full$scores[[meth]]
      est <- mean(sc); se <- if (variance) stats::sd(sc) / sqrt(nrow(data)) else NA_real_
      full_est <- mean(full$scores[[meth]])
      full_se <- if (variance) stats::sd(full$scores[[meth]]) / sqrt(nrow(data)) else NA_real_
      rows[[h]] <- data.frame(benchmark_id = j, label = labels[j],
        p_z = length(subsets[[j]]), retained = paste(subsets[[j]], collapse = ","),
        p_x = length(X), n_omitted = length(omitted),
        omitted = paste(omitted, collapse = ","), estimator = meth,
        full_estimate = full_est, full_se = full_se,
        full_conf_low = if (variance) full_est - zcrit * full_se else NA_real_,
        full_conf_high = if (variance) full_est + zcrit * full_se else NA_real_,
        reduced_estimate = mean(red$scores[[meth]]), delta_bench = est,
        se = se, conf_low = if (variance) est - zcrit * se else NA_real_,
        conf_high = if (variance) est + zcrit * se else NA_real_,
        se_type = if (meth == "aipw") "paired EIF" else
          "conditional plug-in spread; nuisance uncertainty omitted",
        stringsAsFactors = FALSE)
      if (variance) {
        key <- paste(j, meth, sep = "_")
        scores[[key]] <- sc - est
        target_scores[[key]] <- red$scores[[meth]] - mean(red$scores[[meth]])
      }
    }
  }
  out <- list(results = do.call(rbind, rows), subsets = subsets,
    scores = scores, target_scores = target_scores,
    full_scores = if (variance) lapply(full$scores,
      function(z) z - mean(z)) else list(),
    full = full$estimates,
    metadata = list(X = X, max_per_pz = max_per_pz,
      subset_selection = subset_selection, variance = variance), call = match.call())
  class(out) <- "jointbounds_l2_benchmarks"
  out
}

l2_benchmark_ate_fit <- function(data, Y, A, C, X, folds, method,
                                 lib_e, lib_r, lib_y, seed) {
  fit_data <- data; fit_X <- X
  if (!length(fit_X)) {
    fit_X <- ".l2_benchmark_intercept"
    fit_data[[fit_X]] <- 0
  }
  nuis <- if (method == "SuperLearner") l2_nuisance_sl(fit_data, Y, A, C, fit_X,
    folds, seed, lib_e, lib_r, lib_y) else
    l2_nuisance_glm(data, Y, A, C, X, folds, seed)
  y <- data[[Y]]; y0 <- y; y0[!is.finite(y0)] <- 0
  arm <- list()
  for (a in 0:1) {
    e <- nuis[[paste0("e", a)]]; rho <- nuis[[paste0("rho", a)]]
    mu <- nuis[[paste0("mu", a)]]
    S <- as.numeric(data[[A]] == a & data[[C]] == 0)
    arm[[a + 1L]] <- list(plugin = mu,
      aipw = mu + S / (e * rho) * (y0 - mu))
  }
  scores <- list(plugin = arm[[2L]]$plugin - arm[[1L]]$plugin,
                 aipw = arm[[2L]]$aipw - arm[[1L]]$aipw)
  list(scores = scores, estimates = vapply(scores, mean, numeric(1)), nuisance = nuis)
}

#' Simultaneous multiplier intervals for calibration benchmarks
#'
#' Uses the paired cross-fitted AIPW influence scores retained by
#' `l2_aipw_benchmarks(..., variance = TRUE)` to form a simultaneous interval
#' over every selected covariate-omission benchmark.
#'
#' @param benchmarks Output from `l2_aipw_benchmarks()` estimated with
#'   `variance = TRUE`.
#' @param estimator Benchmark estimator; simultaneous influence-function
#'   inference is available for `"aipw"`.
#' @param B Number of multiplier draws.
#' @param conf_level Simultaneous confidence level.
#' @param multiplier Rademacher or standard-normal multipliers.
#' @param seed Random seed.
#' @return Benchmark results augmented with simultaneous intervals, together
#'   with the multiplier critical value and draw distribution.
#' @export
l2_benchmark_band <- function(benchmarks, estimator = "aipw", B = 1000L,
                              conf_level = .95,
                              multiplier = c("rademacher", "normal"),
                              seed = 1L) {
  if (!inherits(benchmarks, "jointbounds_l2_benchmarks"))
    stop("benchmarks must be returned by l2_aipw_benchmarks()")
  if (!identical(estimator, "aipw"))
    stop("Influence-function benchmark bands require estimator = 'aipw'")
  multiplier <- match.arg(multiplier); B <- as.integer(B)
  if (length(B) != 1L || !is.finite(B) || B < 2L)
    stop("B must be an integer of at least two")
  if (length(conf_level) != 1L || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1)
    stop("conf_level must lie strictly between zero and one")
  d <- benchmarks$results[benchmarks$results$estimator == estimator, , drop = FALSE]
  keys <- paste(d$benchmark_id, estimator, sep = "_")
  score_list <- benchmarks$scores[keys]
  if (!nrow(d) || any(vapply(score_list, is.null, logical(1))))
    stop("Re-estimate benchmarks with estimator including 'aipw' and variance = TRUE")
  scores <- do.call(cbind, score_list)
  sigma <- apply(scores, 2L, stats::sd); n <- nrow(scores)
  regular <- is.finite(sigma) & sigma > 0
  if (!any(regular)) stop("All benchmark influence scores have zero variance")
  set.seed(seed); suprema <- numeric(B)
  for (b in seq_len(B)) {
    xi <- if (multiplier == "rademacher")
      sample(c(-1, 1), n, replace = TRUE) else stats::rnorm(n)
    proc <- colSums(scores[, regular, drop = FALSE] * xi) /
      (sqrt(n) * sigma[regular])
    suprema[b] <- max(abs(proc))
  }
  critical <- unname(stats::quantile(suprema, conf_level, type = 8))
  d$simultaneous_low <- d$delta_bench - critical * d$se
  d$simultaneous_high <- d$delta_bench + critical * d$se
  out <- list(results = d, critical_value = critical, suprema = suprema,
    scores = scores, estimator = estimator, conf_level = conf_level,
    multiplier = multiplier, B = B, call = match.call())
  class(out) <- "jointbounds_l2_benchmark_band"
  out
}

#' Simultaneous confidence bands for calibration frontiers
#'
#' Combines each sensitivity-endpoint influence score with the paired
#' reduced-adjustment ATE influence score, constructs pointwise and joint
#' multiplier bands for their difference, and inverts those bands over the
#' evaluated sensitivity grid.  This implements the calibration-frontier
#' inference construction in the companion paper while retaining the
#' covariance between the endpoint and benchmark estimators.
#'
#' @param surface Output from `l2_sensitivity_band(..., keep_fits = TRUE)`.
#' @param benchmarks Output from `l2_aipw_benchmarks(..., variance = TRUE)`.
#' @param benchmark_estimator Benchmark estimator. Influence-function
#'   inference currently requires `"aipw"`.
#' @param endpoint_estimator Endpoint estimator stored in `surface`; defaults
#'   to `"eif"` for CS bounds and `"dr"` for sharp bounds.
#' @param active_endpoint Use `"auto"` to select the lower or upper endpoint
#'   from the sign of the paired reduced-minus-full benchmark. A fixed
#'   `"lower"`, `"upper"`, or `"both"` may be supplied for a prespecified
#'   direction.
#' @param B Number of multiplier draws.
#' @param conf_level Confidence level for pointwise and simultaneous bands.
#' @param multiplier Rademacher or standard-normal multipliers.
#' @param seed Random seed.
#' @param radius_x,radius_y Sensitivity-radius columns. They default to
#'   `delta_M` and `delta_K` when available.
#' @return A list containing the estimated frontiers, the grid-level gap
#'   estimates and bands, the paired influence-score matrix, and the joint
#'   multiplier critical value.
#' @export
l2_calibration_band <- function(surface, benchmarks,
                                benchmark_estimator = "aipw",
                                endpoint_estimator = NULL,
                                active_endpoint = c("auto", "lower", "upper",
                                                    "both"),
                                B = 1000L, conf_level = .95,
                                multiplier = c("rademacher", "normal"),
                                seed = 1L, radius_x = NULL,
                                radius_y = NULL) {
  if (!inherits(surface, "jointbounds_l2_sensitivity_band") ||
      is.null(surface$fits) || !length(surface$fits))
    stop("surface must come from l2_sensitivity_band(..., keep_fits = TRUE)")
  if (!inherits(benchmarks, "jointbounds_l2_benchmarks"))
    stop("benchmarks must come from l2_aipw_benchmarks()")
  if (!identical(benchmark_estimator, "aipw"))
    stop("Combined influence-function inference requires benchmark_estimator = 'aipw'")
  active_endpoint <- match.arg(active_endpoint)
  multiplier <- match.arg(multiplier)
  B <- as.integer(B)
  if (length(B) != 1L || !is.finite(B) || B < 2L)
    stop("B must be an integer of at least two")
  if (length(conf_level) != 1L || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1)
    stop("conf_level must lie strictly between zero and one")
  grid <- surface$grid
  if (!is.data.frame(grid) || nrow(grid) != length(surface$fits))
    stop("surface does not retain a valid sensitivity grid and fit list")
  if (is.null(radius_x)) radius_x <- if ("delta_M" %in% names(grid))
    "delta_M" else "delta_R"
  if (is.null(radius_y)) radius_y <- if ("delta_K" %in% names(grid))
    "delta_K" else "delta_A"
  if (!all(c(radius_x, radius_y) %in% names(grid)))
    stop("The sensitivity grid must contain radius_x and radius_y")
  if (is.null(endpoint_estimator)) endpoint_estimator <-
    if (identical(surface$bound_method, "sharp")) "dr" else "eif"

  bd <- benchmarks$results[
    benchmarks$results$estimator == benchmark_estimator, , drop = FALSE]
  keys <- paste(bd$benchmark_id, benchmark_estimator, sep = "_")
  target_scores <- benchmarks$target_scores[keys]
  if (!nrow(bd) || any(vapply(target_scores, is.null, logical(1))))
    stop("Re-estimate benchmarks with estimator including 'aipw' and variance = TRUE")
  n <- unique(vapply(target_scores, length, integer(1)))
  if (length(n) != 1L) stop("Benchmark influence scores have unequal lengths")

  gap_rows <- list(); gap_scores <- list(); h <- 0L
  for (j in seq_len(nrow(bd))) {
    target <- bd$reduced_estimate[j]
    reference <- bd$full_estimate[j]
    if (!is.finite(target) || !is.finite(reference))
      stop("Benchmark target and full-X reference estimates must be finite")
    direction <- if (active_endpoint != "auto") active_endpoint else
      if (target < reference) "lower" else if (target > reference)
        "upper" else "both"
    if (direction == "both") {
      warning("Benchmark ", bd$label[j],
        " equals its full-X estimate; its calibration frontier begins at the origin")
    }
    target_score <- target_scores[[j]] - mean(target_scores[[j]])
    for (i in seq_len(nrow(grid))) {
      fit <- surface$fits[[i]]
      endpoints <- if (direction == "both") c("lower", "upper") else direction
      for (endpoint in endpoints) {
        rr <- fit$ate[fit$ate$estimator == endpoint_estimator &
                        fit$ate$endpoint == endpoint, , drop = FALSE]
        score <- fit$scores[[paste(endpoint_estimator, endpoint, sep = "_")]]
        if (nrow(rr) != 1L || is.null(score) || length(score) != n)
          stop("Surface fits and benchmark scores must use the same ordered observations")
        if (endpoint == "lower") {
          gap <- rr$estimate - target
          score <- score - target_score
        } else {
          gap <- target - rr$estimate
          score <- target_score - score
        }
        h <- h + 1L
        score <- score - mean(score)
        gap_scores[[h]] <- score
        gap_rows[[h]] <- data.frame(benchmark_id = bd$benchmark_id[j],
          label = bd$label[j], grid_id = i, endpoint = endpoint,
          direction = direction, target = target, reference = reference,
          gap = gap, se = stats::sd(score) / sqrt(n),
          stringsAsFactors = FALSE)
      }
    }
  }
  gaps <- do.call(rbind, gap_rows)
  gaps <- cbind(gaps, grid[gaps$grid_id, , drop = FALSE])
  score_matrix <- do.call(cbind, gap_scores)
  sigma <- apply(score_matrix, 2L, stats::sd)
  regular <- is.finite(sigma) & sigma > 0
  if (!any(regular)) stop("All combined calibration-gap scores have zero variance")
  set.seed(seed)
  suprema <- numeric(B)
  for (b in seq_len(B)) {
    xi <- if (multiplier == "rademacher")
      sample(c(-1, 1), n, replace = TRUE) else stats::rnorm(n)
    process <- colSums(score_matrix[, regular, drop = FALSE] * xi) /
      (sqrt(n) * sigma[regular])
    suprema[b] <- max(abs(process))
  }
  critical <- unname(stats::quantile(suprema, conf_level, type = 8))
  zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)
  gaps$pointwise_low <- gaps$gap - zcrit * gaps$se
  gaps$pointwise_high <- gaps$gap + zcrit * gaps$se
  gaps$simultaneous_low <- gaps$gap - critical * gaps$se
  gaps$simultaneous_high <- gaps$gap + critical * gaps$se

  crossing <- function(z, value) {
    z <- z[order(z[[radius_y]]), , drop = FALSE]
    yy <- z[[value]]; xx <- z[[radius_y]]
    ok <- is.finite(yy) & is.finite(xx)
    yy <- yy[ok]; xx <- xx[ok]
    if (!length(yy) || !any(yy <= 0)) return(NA_real_)
    k <- which(yy <= 0)[1L]
    if (k == 1L || !is.finite(yy[k - 1L]) || yy[k] == yy[k - 1L])
      return(xx[k])
    frac <- -yy[k - 1L] / (yy[k] - yy[k - 1L])
    if (!is.finite(frac) || frac < 0 || frac > 1) return(xx[k])
    xx[k - 1L] + frac * (xx[k] - xx[k - 1L])
  }
  group_key <- interaction(gaps$benchmark_id, gaps$endpoint,
                           gaps[[radius_x]], drop = TRUE)
  frontier <- do.call(rbind, lapply(split(gaps, group_key), function(z) {
    vals <- vapply(c("gap", "pointwise_low", "pointwise_high",
      "simultaneous_low", "simultaneous_high"),
      function(nm) crossing(z, nm), numeric(1))
    data.frame(benchmark_id = z$benchmark_id[1L], label = z$label[1L],
      endpoint = z$endpoint[1L], direction = z$direction[1L],
      radius_x_value = z[[radius_x]][1L], estimate = vals[["gap"]],
      pointwise_low = vals[["pointwise_low"]],
      pointwise_high = vals[["pointwise_high"]],
      simultaneous_low = vals[["simultaneous_low"]],
      simultaneous_high = vals[["simultaneous_high"]],
      stringsAsFactors = FALSE)
  }))
  names(frontier)[names(frontier) == "radius_x_value"] <- radius_x
  for (pair in list(c("pointwise_low", "pointwise_high"),
                    c("simultaneous_low", "simultaneous_high"))) {
    lo <- pmin(frontier[[pair[1L]]], frontier[[pair[2L]]], na.rm = TRUE)
    hi <- pmax(frontier[[pair[1L]]], frontier[[pair[2L]]], na.rm = TRUE)
    neither <- !is.finite(frontier[[pair[1L]]]) &
      !is.finite(frontier[[pair[2L]]])
    lo[neither] <- hi[neither] <- NA_real_
    frontier[[pair[1L]]] <- lo; frontier[[pair[2L]]] <- hi
  }
  frontier <- frontier[order(frontier$benchmark_id, frontier$endpoint,
                             frontier[[radius_x]]), , drop = FALSE]
  rownames(frontier) <- NULL
  out <- list(frontier = frontier, gaps = gaps, scores = score_matrix,
    critical_value = critical, suprema = suprema, conf_level = conf_level,
    multiplier = multiplier, endpoint_estimator = endpoint_estimator,
    benchmark_estimator = benchmark_estimator, radius_x = radius_x,
    radius_y = radius_y, B = B, call = match.call())
  class(out) <- "jointbounds_l2_calibration_band"
  out
}

#' Plot calibration-frontier confidence bands
#' @param x Output from `l2_calibration_band()`.
#' @param ... Additional arguments currently ignored.
#' @return A `ggplot2` object.
#' @export
plot.jointbounds_l2_calibration_band <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Plotting calibration bands requires the suggested package ggplot2")
  d <- x$frontier
  if (!nrow(d) || !any(is.finite(d$estimate)))
    stop("No finite calibration-frontier crossings are available to plot")
  rx <- x$radius_x; ry <- x$radius_y
  d$.radius_x <- d[[rx]]
  ggplot2::ggplot(d, ggplot2::aes(x = .radius_x, y = estimate)) +
    ggplot2::geom_ribbon(ggplot2::aes(
      ymin = simultaneous_low, ymax = simultaneous_high),
      fill = "#4C78A8", alpha = .16, na.rm = TRUE) +
    ggplot2::geom_ribbon(ggplot2::aes(
      ymin = pointwise_low, ymax = pointwise_high),
      fill = "#4C78A8", alpha = .28, na.rm = TRUE) +
    ggplot2::geom_line(linewidth = .9, colour = "#173F5F", na.rm = TRUE) +
    ggplot2::facet_wrap(stats::as.formula("~ label")) +
    ggplot2::labs(x = rx, y = paste0("Frontier value of ", ry),
      title = "Estimated calibration frontiers",
      subtitle = paste0(format(100 * x$conf_level, trim = TRUE),
        "% pointwise and simultaneous influence-function bands")) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"))
}

#' Invert a sensitivity surface for many calibration benchmarks
#'
#' @param surface Output from l2_sensitivity_grid or its results data frame.
#' @param benchmarks Output from l2_aipw_benchmarks or its results data frame.
#' @param estimator Which benchmark estimator to invert.
#' @param plot Draw the resulting multi-benchmark frontier.
#' @inheritParams l2_calibration_frontier
#' @return A plot-ready frontier data frame with one curve per benchmark. The
#'   `equal_radius` column is the smallest common budget sufficient to reproduce
#'   the benchmark on the evaluated grid; its summary is also stored as an
#'   attribute and retrieved by `l2_equal_radius()`.
#' @export
l2_benchmark_frontiers <- function(surface, benchmarks,
                                   estimator = c("aipw", "plugin"),
                                   radius_x = NULL, radius_y = NULL,
                                   plot = FALSE) {
  estimator <- match.arg(estimator)
  surface_data <- if (inherits(surface, "jointbounds_l2_grid"))
    surface$results else surface
  if (is.null(radius_x)) radius_x <- if ("delta_M" %in% names(surface_data))
    "delta_M" else "delta_R"
  if (is.null(radius_y)) radius_y <- if ("delta_K" %in% names(surface_data))
    "delta_K" else "delta_A"
  b <- if (inherits(benchmarks, "jointbounds_l2_benchmarks"))
    benchmarks$results else benchmarks
  b <- b[b$estimator == estimator, , drop = FALSE]
  if (!nrow(b)) stop("No benchmarks found for estimator = ", estimator)
  ans <- lapply(seq_len(nrow(b)), function(i) {
    paired_reference <- if ("full_estimate" %in% names(b) &&
      is.finite(b$full_estimate[i])) b$full_estimate[i] else NULL
    z <- l2_calibration_frontier(surface, benchmark = b$delta_bench[i],
      type = "discrepancy", reference = paired_reference,
      radius_x = radius_x, radius_y = radius_y)
    z$benchmark_id <- b$benchmark_id[i]; z$label <- b$label[i]
    z$p_z <- b$p_z[i]; z$estimator <- estimator
    if ("p_x" %in% names(b)) z$p_x <- b$p_x[i]
    if ("n_omitted" %in% names(b)) z$n_omitted <- b$n_omitted[i]
    if ("omitted" %in% names(b)) z$omitted <- b$omitted[i]
    for (nm in c("se", "conf_low", "conf_high", "full_estimate", "full_se",
                 "full_conf_low", "full_conf_high"))
      if (nm %in% names(b)) z[[nm]] <- b[[nm]][i]
    z
  })
  out <- do.call(rbind, ans); rownames(out) <- NULL
  class(out) <- c("jointbounds_l2_benchmark_frontiers", "data.frame")
  attr(out, "radius_x") <- radius_x; attr(out, "radius_y") <- radius_y
  eq <- l2_equal_radius(out)
  key_out <- interaction(out$benchmark_id, out$method, drop = TRUE)
  key_eq <- interaction(eq$benchmark_id, eq$method, drop = TRUE)
  out$equal_radius <- eq$equal_radius[match(key_out, key_eq)]
  attr(out, "equal_radius_summary") <- eq
  if (isTRUE(plot)) print(graphics::plot(out))
  out
}

#' Equal-strength joint calibration value
#'
#' Returns the smallest common radius `x` such that budgets
#' `(delta_M, delta_K) = (x, x)` are sufficient to reproduce each benchmark
#' under the primary model (or the corresponding two radii supplied by the
#' user),
#' using the evaluated sensitivity grid. Equivalently, it minimizes
#' the maximum of the two radii along the compatible frontier.
#'
#' @param x Output from `l2_benchmark_frontiers()`.
#' @return One row per benchmark and bound method, including `equal_radius`.
#' @export
l2_equal_radius <- function(x) {
  if (!inherits(x, "jointbounds_l2_benchmark_frontiers"))
    stop("x must be returned by l2_benchmark_frontiers()")
  rx <- attr(x, "radius_x"); ry <- paste0("minimum_", attr(x, "radius_y"))
  groups <- split(x, interaction(x$benchmark_id, x$method, drop = TRUE))
  ans <- lapply(groups, function(d) {
    ok <- d$compatible_on_grid & is.finite(d[[rx]]) & is.finite(d[[ry]])
    value <- if (any(ok)) min(pmax(d[[rx]][ok], d[[ry]][ok])) else NA_real_
    data.frame(benchmark_id = d$benchmark_id[1L], label = d$label[1L],
      p_z = d$p_z[1L], estimator = d$estimator[1L], method = d$method[1L],
      benchmark = d$benchmark[1L], equal_radius = value,
      p_x = if ("p_x" %in% names(d)) d$p_x[1L] else NA_integer_,
      n_omitted = if ("n_omitted" %in% names(d)) d$n_omitted[1L] else NA_integer_,
      omitted = if ("omitted" %in% names(d)) d$omitted[1L] else NA_character_,
      compatible_on_grid = any(ok), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, ans); rownames(out) <- NULL
  out[order(out$p_z, out$benchmark_id, out$method), , drop = FALSE]
}

#' Overall equal-radius tipping value
#'
#' Finds the smallest common sensitivity budget `x` for which the estimated
#' ATE bound contains `threshold`. The calculation minimizes
#' the maximum of the two chosen radii over tipping grid points and is
#' therefore the grid analogue of moving along their equal-radius diagonal.
#'
#' @param surface Output from `l2_sensitivity_grid()` or its results data frame.
#' @param threshold Null value for the ATE.
#' @param radius_x,radius_y Sensitivity-radius columns.
#' @return One row per bound method.
#' @export
l2_overall_equal_radius <- function(surface, threshold = 0,
                                    radius_x = NULL,
                                    radius_y = NULL) {
  d <- if (inherits(surface, "jointbounds_l2_grid")) surface$results else surface
  if (is.null(radius_x)) radius_x <- if ("delta_M" %in% names(d))
    "delta_M" else "delta_R"
  if (is.null(radius_y)) radius_y <- if ("delta_K" %in% names(d))
    "delta_K" else "delta_A"
  needed <- c("lower", "upper", "method", radius_x, radius_y)
  if (!is.data.frame(d) || !all(needed %in% names(d)))
    stop("surface does not contain the required bound and radius columns")
  groups <- split(d, d$method)
  ans <- lapply(groups, function(z) {
    finite_grid <- is.finite(z[[radius_x]]) & is.finite(z[[radius_y]])
    ok <- is.finite(z$lower) & is.finite(z$upper) & finite_grid &
      z$lower <= threshold & z$upper >= threshold
    value <- if (any(ok)) min(pmax(z[[radius_x]][ok], z[[radius_y]][ok])) else NA_real_
    searched <- if (any(finite_grid))
      max(pmax(z[[radius_x]][finite_grid], z[[radius_y]][finite_grid])) else NA_real_
    data.frame(method = z$method[1L], threshold = threshold,
      overall_equal_radius = value, tipped_on_grid = any(ok),
      maximum_radius_searched = searched,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, ans); rownames(out) <- NULL
  out[order(out$method), , drop = FALSE]
}

l2_grid_equal_radius_value <- function(surface, threshold, method,
                                       radius_x, radius_y, shift = 0) {
  d <- if (inherits(surface, "jointbounds_l2_grid")) surface$results else surface
  d <- d[d$method == method, , drop = FALSE]
  ok <- is.finite(d$lower) & is.finite(d$upper) &
    is.finite(d[[radius_x]]) & is.finite(d[[radius_y]]) &
    d$lower + shift <= threshold & d$upper + shift >= threshold
  if (any(ok)) min(pmax(d[[radius_x]][ok], d[[radius_y]][ok])) else NA_real_
}

l2_benchmark_equal_radius_value <- function(surface, benchmark, method,
                                             radius_x, radius_y,
                                             reference = NULL) {
  z <- l2_calibration_frontier(surface, benchmark = benchmark,
    type = "discrepancy", reference = reference,
    radius_x = radius_x, radius_y = radius_y)
  ry <- paste0("minimum_", radius_y)
  z <- z[z$method == method & z$compatible_on_grid &
    is.finite(z[[radius_x]]) & is.finite(z[[ry]]), , drop = FALSE]
  if (nrow(z)) min(pmax(z[[radius_x]], z[[ry]])) else NA_real_
}

#' Compare calibrated and overall equal-radius values
#'
#' @param surface Output from `l2_sensitivity_grid()`.
#' @param frontiers Output from `l2_benchmark_frontiers()`.
#' @param threshold Null value for the overall ATE.
#' @param uncertainty Add transformed benchmark Wald intervals and an overall
#'   Wald location-shift approximation. Requires `variance = TRUE` when the
#'   benchmarks were estimated.
#' @param conf_level Confidence level used when `uncertainty = TRUE`.
#' @param plot Draw the comparison plot.
#' @return A data frame comparing every calibration benchmark with the overall
#'   tipping value. `strength_ratio` is benchmark divided by overall strength.
#' @export
l2_equal_radius_comparison <- function(surface, frontiers, threshold = 0,
                                       uncertainty = FALSE, conf_level = 0.95,
                                       plot = FALSE) {
  eq <- l2_equal_radius(frontiers)
  overall <- l2_overall_equal_radius(surface, threshold,
    radius_x = attr(frontiers, "radius_x"),
    radius_y = attr(frontiers, "radius_y"))
  out <- merge(eq, overall, by = "method", all.x = TRUE, sort = FALSE)
  out$strength_ratio <- out$equal_radius / out$overall_equal_radius
  out$strength_class <- ifelse(!is.finite(out$strength_ratio), "Unavailable",
    ifelse(out$strength_ratio >= 1, "At or above tipping value",
      "Below tipping value"))
  out$overall_status <- ifelse(out$tipped_on_grid,
    "Tipping value found", "Tipping value beyond evaluated grid")
  out$benchmark_radius_low <- out$benchmark_radius_high <- out$equal_radius
  out$overall_radius_low <- out$overall_radius_high <- out$overall_equal_radius
  if (isTRUE(uncertainty)) {
    if (length(conf_level) != 1L || !is.finite(conf_level) ||
        conf_level <= 0 || conf_level >= 1)
      stop("conf_level must lie strictly between zero and one")
    zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)
    rx <- attr(frontiers, "radius_x"); ry <- attr(frontiers, "radius_y")
    meta <- frontiers[!duplicated(interaction(frontiers$benchmark_id,
      frontiers$method, drop = TRUE)), , drop = FALSE]
    required <- c("benchmark", "se", "full_estimate", "full_se")
    if (!all(required %in% names(meta)) ||
        any(!is.finite(as.matrix(meta[required]))))
      stop("uncertainty = TRUE requires benchmarks estimated with variance = TRUE; recreate bench and frontiers")
    for (i in seq_len(nrow(out))) {
      m <- meta[meta$benchmark_id == out$benchmark_id[i] &
        meta$method == out$method[i], , drop = FALSE][1L, ]
      bench_low <- m$benchmark - zcrit * m$se
      bench_high <- m$benchmark + zcrit * m$se
      full_low <- m$full_estimate - zcrit * m$full_se
      full_high <- m$full_estimate + zcrit * m$full_se
      br <- c(l2_benchmark_equal_radius_value(surface, bench_low,
        out$method[i], rx, ry, reference = m$full_estimate),
        l2_benchmark_equal_radius_value(surface, bench_high,
        out$method[i], rx, ry, reference = m$full_estimate))
      out$benchmark_radius_low[i] <- if (bench_low <= 0 && bench_high >= 0)
        0 else if (any(is.finite(br))) min(br[is.finite(br)]) else NA_real_
      out$benchmark_radius_high[i] <- if (all(!is.finite(br))) NA_real_ else
        max(br, na.rm = TRUE)
      near <- if (m$full_estimate >= threshold) full_low else full_high
      far <- if (m$full_estimate >= threshold) full_high else full_low
      out$overall_radius_low[i] <- l2_grid_equal_radius_value(surface, threshold,
        out$method[i], rx, ry, shift = near - m$full_estimate)
      out$overall_radius_high[i] <- l2_grid_equal_radius_value(surface, threshold,
        out$method[i], rx, ry, shift = far - m$full_estimate)
      oo <- sort(c(out$overall_radius_low[i], out$overall_radius_high[i]))
      if (all(is.finite(oo))) {
        out$overall_radius_low[i] <- oo[1L]
        out$overall_radius_high[i] <- oo[2L]
      }
    }
  }
  attr(out, "uncertainty") <- isTRUE(uncertainty)
  attr(out, "conf_level") <- conf_level
  attr(out, "uncertainty_note") <- if (isTRUE(uncertainty))
    paste("Benchmark intervals transform paired AIPW Wald intervals;",
      "overall intervals are Wald location-shift approximations holding",
      "the estimated sensitivity-width surface fixed.") else "Point estimates"
  class(out) <- c("jointbounds_l2_equal_radius_comparison", "data.frame")
  if (isTRUE(plot)) print(graphics::plot(out))
  out
}

#' @export
plot.jointbounds_l2_equal_radius_comparison <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Plotting equal-radius comparisons requires the suggested package ggplot2")
  d <- x[is.finite(x$equal_radius), , drop = FALSE]
  if (!nrow(d)) stop("No calibrated equal-radius values are available to plot; expand the calibration grid")
  if (all(!is.finite(d$n_omitted)))
    stop("Recreate frontiers with the updated package to retain omitted-covariate counts")
  d$.n_omitted_plot <- d$n_omitted
  refs <- unique(d[is.finite(d$overall_equal_radius),
    c("method", "overall_equal_radius", "overall_radius_low",
      "overall_radius_high"), drop = FALSE])
  reached <- nrow(refs) > 0L
  show_uncertainty <- isTRUE(attr(x, "uncertainty"))
  pos <- ggplot2::position_jitter(width = 0, height = .09, seed = 1)
  p <- ggplot2::ggplot(d, ggplot2::aes(
      x = equal_radius, y = .n_omitted_plot)) +
    ggplot2::geom_rect(data = refs,
      ggplot2::aes(xmin = overall_radius_low,
        xmax = overall_radius_high, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "#D95F5F", alpha = if (show_uncertainty) .16 else 0) +
    ggplot2::geom_vline(data = refs,
      ggplot2::aes(xintercept = overall_equal_radius),
      inherit.aes = FALSE, linewidth = .85, linetype = "dashed",
      colour = "#B33A3A") +
    ggplot2::geom_errorbar(ggplot2::aes(
      xmin = benchmark_radius_low, xmax = benchmark_radius_high),
      orientation = "y", position = pos, width = .08, linewidth = .65,
      alpha = if (show_uncertainty) .8 else 0, colour = "#176B87") +
    ggplot2::geom_point(position = pos, size = 3.4, alpha = .78,
      colour = "#176B87") +
    ggplot2::facet_wrap(stats::as.formula("~ method")) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(.04, .08))) +
    ggplot2::scale_y_reverse(
      breaks = sort(unique(d$n_omitted[is.finite(d$n_omitted)]))) +
    ggplot2::labs(x = "Equal L2 radius", y = "Number of covariates omitted",
      title = "Calibrated strength versus overall tipping strength",
      subtitle = if (reached)
        if (show_uncertainty)
          paste0("Horizontal intervals and red band use ",
            format(100 * attr(x, "conf_level"), trim = TRUE),
            "% Wald uncertainty") else
          "Dashed line: smallest equal radius for which the ATE bounds include the null" else
        "Overall tipping value was not reached on the evaluated sensitivity grid") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(colour = "grey35"),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"))
  p
}

# Axis label for a radius column name, e.g. "delta_M" -> delta[M]
radius_label <- function(name) {
  sub_name <- sub("^delta_", "", name)
  if (identical(sub_name, name)) return(name)
  bquote(delta[.(sub_name)])
}

#' @export
plot.jointbounds_l2_benchmark_frontiers <- function(x, facet_pz = TRUE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Plotting benchmark frontiers requires the suggested package ggplot2")
  rx <- attr(x, "radius_x"); ry <- paste0("minimum_", attr(x, "radius_y"))
  d <- x[x$compatible_on_grid & is.finite(x[[rx]]) & is.finite(x[[ry]]), , drop = FALSE]
  if (!nrow(d)) stop("No compatible frontier points are available to plot")
  d$.radius_R <- d[[rx]]; d$.radius_A <- d[[ry]]
  eq <- attr(x, "equal_radius_summary")
  if (is.null(eq)) eq <- l2_equal_radius(x)
  eq <- eq[is.finite(eq$equal_radius), , drop = FALSE]
  eq$.radius_R <- eq$equal_radius; eq$.radius_A <- eq$equal_radius
  p <- ggplot2::ggplot(d, ggplot2::aes(
      x = .radius_R, y = .radius_A, colour = label,
      linetype = method, group = interaction(label, method))) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey55",
      linewidth = .6, linetype = "dashed") +
    ggplot2::geom_line(linewidth = 1.05, lineend = "round") +
    ggplot2::geom_point(data = eq, ggplot2::aes(
      x = .radius_R, y = .radius_A, colour = label),
      inherit.aes = FALSE, shape = 21, fill = "white", stroke = 1.1,
      size = 3.1) +
    ggplot2::scale_colour_viridis_d(option = "D", end = .85) +
    ggplot2::labs(x = radius_label(attr(x, "radius_x")),
      y = radius_label(attr(x, "radius_y")),
      colour = "Retained adjustment set", linetype = "Bound method",
      title = "Joint calibration frontiers",
      subtitle = "Dots mark the smallest equal MAR and NUC sensitivity radius") +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = .35),
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(colour = "grey35"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom", legend.box = "vertical",
      axis.title = ggplot2::element_text(face = "bold"))
  if (isTRUE(facet_pz))
    p <- p + ggplot2::facet_wrap(stats::as.formula("~ p_z"),
      labeller = ggplot2::label_both)
  p
}
