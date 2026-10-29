#' Select calibrated sensitivity budgets for subsequent bound estimation
#'
#' Selects the first compatible evaluated vertical radius at each requested
#' horizontal radius. Upper-limit selections invert the inner compatibility
#' set using both endpoint inequalities. No interpolation or extrapolation is
#' used, so the selected point itself satisfies the grid-based criterion.
#'
#' @param x Output from `l2_calibration_band()` using both endpoints.
#' @param selection `"estimate"` uses estimated compatibility;
#'   `"pointwise_upper"` uses the upper confidence limit for a prespecified
#'   horizontal radius; `"simultaneous_upper"` uses the simultaneous upper
#'   confidence band, permitting selection across the evaluated horizontal
#'   radii. Selection across benchmarks additionally requires a band computed
#'   with `simultaneous_scope = "all_benchmarks"`.
#' @param radius_x Optional horizontal radius values, already present in the
#'   evaluated grid. NULL returns all horizontal radii.
#' @param benchmark_id Optional benchmark IDs. NULL selects all benchmarks.
#' @return A data frame with benchmark IDs, labels, and the two radius columns,
#'   suitable as `grid` in `l2_sensitivity_band()`. Attributes record the
#'   selection and confidence scope. An unavailable crossing raises an error;
#'   it is never replaced by zero or by the largest evaluated radius.
#' @details Upper-limit calibration enlarges the selected sensitivity region
#'   relative to grid-based point calibration. It does not itself give confidence
#'   coverage for subsequently estimated bounds at data-selected radii. Their
#'   estimation uncertainty and selection require separate inference. Pointwise
#'   limits do not justify inspecting the full curve and choosing a radius.
#' @export
l2_calibrated_grid <- function(x,
    selection = c("estimate", "pointwise_upper", "simultaneous_upper"),
    radius_x = NULL, benchmark_id = NULL) {
  selection <- match.arg(selection)
  if (!inherits(x, "jointbounds_l2_calibration_band"))
    stop("x must be returned by l2_calibration_band()")
  rx <- x$radius_x; ry <- x$radius_y
  z <- x$gaps
  column <- switch(selection, estimate = "gap",
    pointwise_upper = "pointwise_high", simultaneous_upper = "simultaneous_high")
  required <- c("benchmark_id", "label", "endpoint", rx, ry, column)
  if (!is.data.frame(z) || !all(required %in% names(z)))
    stop("x does not retain the required contrast grid")
  if (!is.null(benchmark_id)) {
    if (!length(benchmark_id) || anyNA(benchmark_id) ||
        any(!benchmark_id %in% z$benchmark_id)) stop("Unknown benchmark_id")
    z <- z[z$benchmark_id %in% benchmark_id, , drop = FALSE]
  }
  if (!is.null(radius_x)) {
    if (!is.numeric(radius_x) || !length(radius_x) || any(!is.finite(radius_x)) ||
        any(!radius_x %in% z[[rx]]))
      stop("radius_x must contain values already evaluated in the grid")
    for (id in unique(z$benchmark_id)) {
      if (!all(radius_x %in% z[[rx]][z$benchmark_id == id]))
        stop("Each selected benchmark must include every requested horizontal radius")
    }
    z <- z[z[[rx]] %in% radius_x, , drop = FALSE]
  }
  if (!nrow(z)) stop("No selected grid points")
  groups <- split(z, interaction(z$benchmark_id, z[[rx]], drop = TRUE))
  rows <- lapply(groups, function(g) {
    vertical <- sort(unique(g[[ry]]))
    compatible <- vapply(vertical, function(v) {
      q <- g[g[[ry]] == v, , drop = FALSE]
      if (nrow(q) != 2L || !setequal(q$endpoint, c("lower", "upper")))
        stop("Calibration must retain both endpoints; use active_endpoint = 'both'")
      all(is.finite(q[[column]])) && all(q[[column]] <= 0)
    }, logical(1))
    if (!any(compatible)) stop("No ", selection, " crossing within the grid for benchmark ",
      g$benchmark_id[1L], " at ", rx, " = ", g[[rx]][1L],
      "; extend the sensitivity grid or select another radius")
    row <- g[1L, c("benchmark_id", "label", rx), drop = FALSE]
    row[[ry]] <- vertical[which(compatible)[1L]]
    row
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$benchmark_id, out[[rx]]), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "selection") <- selection
  attr(out, "conf_level") <- x$conf_level
  attr(out, "simultaneous_scope") <- x$simultaneous_scope
  attr(out, "inference_note") <- paste(
    "Calibration uncertainty only; this does not establish confidence coverage",
    "for estimated bounds evaluated at the selected radii.")
  out
}
