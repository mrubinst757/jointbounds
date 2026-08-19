#' Multiplier simultaneous bands for L2 sensitivity endpoints
#'
#' Applies the standardized multiplier process described in the companion paper
#' to a collection of cross-fitted endpoint fits evaluated on the same ordered
#' observations.  This function uses the row-level estimated influence scores
#' returned by `l2_cs_crossfit()` and `l2_sharp_crossfit()`.
#'
#' @param fits A cross-fitted L2 fit, or a named list of such fits evaluated on
#'   the same observations.
#' @param estimator Endpoint estimator whose scores should be used.  Defaults to
#'   `"dr"` for sharp fits and `"eif"` for CS fits.
#' @param B Number of multiplier draws.
#' @param conf_level Simultaneous confidence level.
#' @param multiplier Either Rademacher or standard-normal multipliers.
#' @param seed Random seed.
#' @return An object containing long-form pointwise and simultaneous endpoint
#'   intervals, outward identified-set bands, the critical value, and the
#'   multiplier suprema.
#' @export
l2_multiplier_band <- function(fits, estimator = NULL, B = 1000L,
                               conf_level = .95,
                               multiplier = c("rademacher", "normal"),
                               seed = 1L) {
  multiplier <- match.arg(multiplier)
  B <- as.integer(B)
  if (length(B) != 1L || !is.finite(B) || B < 2L)
    stop("B must be an integer of at least two")
  if (length(conf_level) != 1L || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1)
    stop("conf_level must lie strictly between zero and one")
  is_fit <- function(z) is.list(z) && is.data.frame(z$ate) && is.list(z$scores)
  if (is_fit(fits)) fits <- list(fits)
  if (!is.list(fits) || !length(fits) || any(!vapply(fits, is_fit, logical(1))))
    stop("fits must be a cross-fitted L2 fit or a nonempty list of such fits")
  if (is.null(names(fits))) names(fits) <- rep("", length(fits))
  blank <- !nzchar(names(fits))
  names(fits)[blank] <- paste0("theta_", which(blank))
  if (is.null(estimator)) {
    available <- unique(unlist(lapply(fits, function(z) as.character(z$ate$estimator))))
    estimator <- if ("dr" %in% available) "dr" else "eif"
  }
  if (length(estimator) != 1L || !nzchar(estimator))
    stop("estimator must be one nonempty name")

  columns <- list(); rows <- list(); h <- 0L; n <- NULL
  for (j in seq_along(fits)) {
    fit <- fits[[j]]
    for (endpoint in c("lower", "upper")) {
      rr <- fit$ate[fit$ate$estimator == estimator &
                      fit$ate$endpoint == endpoint, , drop = FALSE]
      key <- paste(estimator, endpoint, sep = "_")
      score <- fit$scores[[key]]
      if (nrow(rr) != 1L || is.null(score))
        stop("Fit ", names(fits)[j], " does not contain estimator '",
             estimator, "' and endpoint '", endpoint, "'")
      if (is.null(n)) n <- length(score)
      if (length(score) != n || any(!is.finite(score)))
        stop("All influence-score vectors must be finite and have a common length")
      h <- h + 1L
      score <- score - mean(score)
      columns[[h]] <- score
      rows[[h]] <- data.frame(grid_id = j, label = names(fits)[j],
        endpoint = endpoint, estimator = estimator, estimate = rr$estimate,
        se = stats::sd(score) / sqrt(n), stringsAsFactors = FALSE)
    }
  }
  score_matrix <- do.call(cbind, columns)
  colnames(score_matrix) <- paste0(
    vapply(rows, `[[`, character(1), "label"), "_",
    vapply(rows, `[[`, character(1), "endpoint"))
  sigma <- apply(score_matrix, 2L, stats::sd)
  regular <- is.finite(sigma) & sigma > 0
  if (!any(regular)) stop("All endpoint influence scores have zero variance")
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
  results <- do.call(rbind, rows); rownames(results) <- NULL
  results$pointwise_low <- results$estimate - zcrit * results$se
  results$pointwise_high <- results$estimate + zcrit * results$se
  results$simultaneous_low <- results$estimate - critical * results$se
  results$simultaneous_high <- results$estimate + critical * results$se
  outward <- do.call(rbind, lapply(split(results, results$grid_id), function(z) {
    lo <- z[z$endpoint == "lower", , drop = FALSE]
    hi <- z[z$endpoint == "upper", , drop = FALSE]
    data.frame(grid_id = z$grid_id[1L], label = z$label[1L],
      estimator = estimator, lower = lo$estimate, upper = hi$estimate,
      outward_lower = lo$simultaneous_low,
      outward_upper = hi$simultaneous_high, stringsAsFactors = FALSE)
  }))
  rownames(outward) <- NULL
  out <- list(results = results, outward = outward, critical_value = critical,
    scores = score_matrix, suprema = suprema,
    conf_level = conf_level, multiplier = multiplier,
    estimator = estimator, n = n, B = B, call = match.call())
  class(out) <- "marbounds_l2_multiplier_band"
  out
}

#' Cross-fitted L2 sensitivity surface with simultaneous bands
#'
#' Evaluates either the CS or sharp estimator over a sensitivity grid using a
#' common fold split, then applies `l2_multiplier_band()` jointly over both ATE
#' endpoints and every evaluated grid point.
#'
#' @inheritParams l2_bounds
#' @param grid A nonempty data frame containing the sensitivity-radius columns.
#'   The primary net model uses `delta_M` and `delta_K`; the separated model
#'   uses `delta_R` and `delta_K`. The legacy name `delta_A` is accepted in
#'   place of `delta_K`. Columns may also be arm-specific by ending in `_0`
#'   and `_1`.
#' @param bound_method Use the influence-function CS estimator or the sharp
#'   cross-fitted estimator.
#' @param estimator Estimator passed to `l2_multiplier_band()`; defaults to
#'   `"eif"` for CS and `"dr"` for sharp bounds.
#' @param B Number of multiplier draws.
#' @param conf_level Simultaneous confidence level.
#' @param multiplier Multiplier distribution.
#' @param keep_fits Retain the potentially large list of pointwise fits.
#' @param ... Additional arguments passed to `l2_cs_crossfit()` or
#'   `l2_sharp_crossfit()`.
#' @return A multiplier-band object augmented with the sensitivity grid.
#' @export
l2_sensitivity_band <- function(data, Y, A, C, X, grid,
                                delta = c(1, 1),
                                model = c("net", "separated"),
                                bound_method = c("cs", "sharp"),
                                estimator = NULL, B = 1000L,
                                conf_level = .95,
                                multiplier = c("rademacher", "normal"),
                                seed = 1L, keep_fits = FALSE, ...) {
  model <- match.arg(model); bound_method <- match.arg(bound_method)
  multiplier <- match.arg(multiplier)
  if (!is.data.frame(grid) || !nrow(grid))
    stop("grid must be a nonempty data frame")
  pair <- function(row, nm, default) {
    arm_names <- paste0(nm, c("_0", "_1"))
    if (all(arm_names %in% names(row)))
      return(as.numeric(row[1, arm_names]))
    if (nm %in% names(row)) return(rep(as.numeric(row[[nm]]), 2L))
    default
  }
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots)))))
    stop("All arguments supplied through ... must be named")
  merge_args <- function(base, extra) {
    base[names(extra)] <- NULL
    c(base, extra)
  }
  fits <- vector("list", nrow(grid))
  names(fits) <- paste0("grid_", seq_len(nrow(grid)))
  for (i in seq_len(nrow(grid))) {
    row <- grid[i, , drop = FALSE]
    rr <- pair(row, "delta_R", c(0, 0))
    aa <- pair(row, "delta_K", pair(row, "delta_A", c(0, 0)))
    dd <- pair(row, "delta", delta)
    mm <- pair(row, "delta_M", rr)
    args <- merge_args(list(data = data, Y = Y, A = A, C = C, X = X,
      delta = dd, delta_R = rr, delta_A = aa, delta_M = mm,
      model = model, seed = seed, conf_level = conf_level), dots)
    fits[[i]] <- if (bound_method == "cs")
      do.call(l2_cs_crossfit, args) else do.call(l2_sharp_crossfit, args)
  }
  if (is.null(estimator)) estimator <- if (bound_method == "cs") "eif" else "dr"
  out <- l2_multiplier_band(fits, estimator = estimator, B = B,
    conf_level = conf_level, multiplier = multiplier, seed = seed)
  out$grid <- grid
  out$results <- cbind(grid[out$results$grid_id, , drop = FALSE], out$results)
  out$outward <- cbind(grid[out$outward$grid_id, , drop = FALSE], out$outward)
  out$bound_method <- bound_method; out$model <- model
  if (isTRUE(keep_fits)) out$fits <- fits
  class(out) <- c("marbounds_l2_sensitivity_band", class(out))
  out
}
