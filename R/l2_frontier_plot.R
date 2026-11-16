#' Plot frontier contours in the two-radius plane
#'
#' Compare calibration or tipping curves, with optional uncertainty ribbons.
#' Linear connections are descriptive: finite-grid confidence guarantees do
#' not establish coverage between grid points. Missing crossings are retained
#' as gaps and are never extrapolated or replaced by the largest grid radius.
#'
#' @param x Output from `l2_calibration_band()`, `l2_tipping_band()`,
#'   `l2_benchmark_frontiers()`, or `l2_tipping_frontier()`.
#' @param uncertainty Confidence limits to display: `"simultaneous"`,
#'   `"pointwise"`, `"both"`, or `"none"`. Defaults to simultaneous for band
#'   objects and none for point-only frontier objects.
#' @param interpolate Connect adjacent finite frontier points and interval
#'   endpoints by straight lines. FALSE shows only grid points and error bars.
#'   This controls display along the horizontal radius; use `interpolate` in
#'   `l2_calibration_band()` or `l2_tipping_band()` to interpolate crossings
#'   along the vertical radius before plotting.
#' @param facet Show one panel per benchmark/method instead of overlaying curves.
#' @param points Mark evaluated horizontal grid radii.
#' @return A ggplot object. Its data retain the original grid crossings.
#' @export
l2_plot_frontier <- function(x, uncertainty = NULL, interpolate = TRUE,
                             facet = FALSE, points = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Plotting frontiers requires the suggested package ggplot2")
  band <- inherits(x, "jointbounds_l2_calibration_band")
  if (band) {
    d <- x$frontier; rx <- x$radius_x; ry <- x$radius_y
    d$.curve <- if ("benchmark_id" %in% names(d))
      paste(d$benchmark_id, d$label, sep = ": ") else as.character(d$label)
  } else if (inherits(x, "jointbounds_l2_benchmark_frontiers") ||
             inherits(x, "jointbounds_l2_frontier")) {
    d <- as.data.frame(x); rx <- attr(x, "radius_x"); ry <- attr(x, "radius_y")
    if (is.null(rx)) rx <- if ("delta_M" %in% names(d)) "delta_M" else "delta_R"
    if (is.null(ry)) ry <- if ("delta_K" %in% names(d)) "delta_K" else "delta_A"
    yy <- if (inherits(x, "jointbounds_l2_benchmark_frontiers")) paste0("minimum_", ry) else ry
    d$estimate <- d[[yy]]
    keys <- intersect(c("benchmark_id", "label", "method"), names(d))
    d$.curve <- if (length(keys)) do.call(paste, c(d[keys], sep = ": ")) else "Frontier"
  } else stop("x must be a calibration or tipping frontier object")
  if (is.null(uncertainty)) uncertainty <- if (band) "simultaneous" else "none"
  uncertainty <- match.arg(uncertainty, c("none", "pointwise", "simultaneous", "both"))
  if (!band && uncertainty != "none") stop("This object has no uncertainty limits; use a frontier band object")
  if (!nrow(d) || !any(is.finite(d$estimate))) stop("No finite frontier crossings are available to plot")
  if (isTRUE(interpolate)) warning(
    "Linear interpolation is descriptive; confidence coverage is established only on the evaluated grid, not between points.",
    call. = FALSE)
  d$.radius <- d[[rx]]
  d <- d[order(d$.curve, d$.radius), , drop = FALSE]
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .radius, y = estimate,
    colour = .curve, fill = .curve, group = .curve))
  kinds <- if (uncertainty == "both") c("simultaneous", "pointwise") else
    if (uncertainty == "none") character() else uncertainty
  for (kind in kinds) {
    z <- d; z$.low <- z[[paste0(kind, "_low")]]; z$.high <- z[[paste0(kind, "_high")]]
    # Split at missing limits so ribbons cannot bridge unavailable crossings.
    z$.segment <- ave(as.integer(!is.finite(z$.low) | !is.finite(z$.high)),
                      z$.curve, FUN = cumsum)
    z$.ribbon <- interaction(z$.curve, z$.segment, drop = TRUE)
    z <- z[is.finite(z$.low) & is.finite(z$.high), , drop = FALSE]
    if (isTRUE(interpolate)) p <- p + ggplot2::geom_ribbon(data = z,
      ggplot2::aes(ymin = .low, ymax = .high, group = .ribbon),
      alpha = if (kind == "simultaneous") .15 else .3, colour = NA)
    else p <- p + ggplot2::geom_errorbar(data = z,
      ggplot2::aes(ymin = .low, ymax = .high), width = 0,
      alpha = if (kind == "simultaneous") .5 else 1)
  }
  d$.segment <- ave(as.integer(!is.finite(d$estimate)), d$.curve, FUN = cumsum)
  d$.line <- interaction(d$.curve, d$.segment, drop = TRUE)
  finite <- d[is.finite(d$estimate), , drop = FALSE]
  if (isTRUE(interpolate)) p <- p + ggplot2::geom_line(data = finite,
    ggplot2::aes(group = .line), linewidth = .8)
  if (isTRUE(points) || !isTRUE(interpolate)) p <- p +
    ggplot2::geom_point(data = finite, size = 1.8)
  if (isTRUE(facet)) p <- p + ggplot2::facet_wrap(~ .curve)
  subtitle <- if (uncertainty == "none") "Point-estimate contours" else
    paste0(100 * x$conf_level, "% ", if (uncertainty == "both")
      "pointwise and simultaneous" else uncertainty, " frontier intervals",
      if (uncertainty == "both") " (darker: pointwise)" else "")
  p + ggplot2::labs(x = rx, y = ry, colour = "Frontier", fill = "Frontier",
    subtitle = subtitle, caption = if (isTRUE(interpolate))
      "Lines and ribbons between grid points are descriptive; between-grid coverage is not established." else
      "Evaluated grid points only; missing crossings are beyond the evaluated grid.") +
    ggplot2::theme_minimal(base_size = 12)
}

utils::globalVariables(c(".curve", ".radius", ".low", ".high", ".ribbon", ".line"))
