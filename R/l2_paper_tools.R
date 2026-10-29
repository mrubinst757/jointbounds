# Internal joint covariance estimates, retaining singular matrices.
l2_joint_covariance <- function(scores) {
  estimators <- unique(sub("_(lower|upper)$", "", names(scores)))
  influence <- lapply(estimators, function(est) {
    z <- cbind(lower = scores[[paste0(est, "_lower")]],
               upper = scores[[paste0(est, "_upper")]])
    stats::cov(z)
  })
  names(influence) <- estimators
  n <- length(scores[[1L]])
  list(influence = influence, estimate = lapply(influence, function(z) z / n))
}

# First point satisfying every oriented inequality (each column <= 0).
l2_frontier_crossing <- function(x, y, interpolate = FALSE) {
  y <- matrix(y, nrow = length(x))
  for (k in seq_along(x)) {
    if (all(is.finite(y[k, ])) && all(y[k, ] <= 0)) return(x[k])
    if (!interpolate || k == length(x) ||
        any(!is.finite(y[c(k, k + 1L), ]))) next
    lo <- 0; hi <- 1
    for (j in seq_len(ncol(y))) {
      a <- y[k, j]; slope <- y[k + 1L, j] - a
      if (slope > 0) hi <- min(hi, -a / slope)
      else if (slope < 0) lo <- max(lo, -a / slope)
      else if (a > 0) hi <- -Inf
    }
    if (lo <= hi) return(x[k] + lo * (x[k + 1L] - x[k]))
  }
  NA_real_
}

#' Confidence bands for a fixed-target tipping frontier
#'
#' Inverts both endpoint contrast bands for a fixed, nonestimated target.
#' Unlike calibration, no benchmark influence score is subtracted.
#' @param surface Output from `l2_sensitivity_band(..., keep_fits = TRUE)`.
#' @param threshold Fixed target, usually zero.
#' @param ... Arguments to `l2_calibration_band()`, such as `B`, `conf_level`,
#'   `radius_x`, `radius_y`, and `interpolate`. Both endpoints are always used.
#' @return A calibration-band object with `target_type = "fixed"`. The
#'   `compatibility` table contains estimated, inner, and outer grid sets.
#' @export
l2_tipping_band <- function(surface, threshold = 0, ...) {
  if (length(threshold) != 1L || !is.finite(threshold))
    stop("threshold must be one finite number")
  if (is.null(surface$fits) || !length(surface$fits))
    stop("surface must retain fits; use keep_fits = TRUE")
  n <- length(surface$fits[[1L]]$scores[[1L]])
  benchmark <- structure(list(results = data.frame(benchmark_id = 1L,
    label = "fixed target", estimator = "aipw", full_estimate = threshold,
    reduced_estimate = threshold), full_scores = list(aipw = numeric(n)),
    target_scores = list(`1_aipw` = numeric(n)), fixed_target = TRUE),
    class = "jointbounds_l2_benchmarks")
  dots <- list(...)
  if (any(names(dots) %in% c("calibration", "alpha", "active_endpoint",
                             "benchmark_estimator", "benchmarks")))
    stop("Tipping bands use a fixed target and both endpoints")
  out <- do.call(l2_calibration_band, c(list(surface = surface,
    benchmarks = benchmark, calibration = "signed", active_endpoint = "both"), dots))
  out$target_type <- "fixed"
  out$threshold <- threshold
  out$call <- match.call()
  out
}
