\name{l2_benchmark_frontiers}
\alias{l2_benchmark_frontiers}
\alias{l2_equal_radius}
\alias{l2_overall_equal_radius}
\alias{l2_equal_radius_comparison}
\alias{plot.marbounds_l2_benchmark_frontiers}
\alias{plot.marbounds_l2_equal_radius_comparison}
\title{Invert sensitivity bounds for multiple calibration benchmarks}
\usage{
l2_benchmark_frontiers(
  surface, benchmarks, estimator = c("aipw", "plugin"),
  radius_x = NULL, radius_y = NULL, plot = FALSE
)
\method{plot}{marbounds_l2_benchmark_frontiers}(x, facet_pz = TRUE, ...)

l2_equal_radius(x)

l2_overall_equal_radius(
  surface, threshold = 0, radius_x = NULL, radius_y = NULL
)

l2_equal_radius_comparison(
  surface, frontiers, threshold = 0, uncertainty = FALSE,
  conf_level = 0.95, plot = FALSE
)

\method{plot}{marbounds_l2_equal_radius_comparison}(x, ...)
}
\description{Computes and plots the joint MAR/NUC radius frontier capable of
producing every supplied leave-covariate-out effect discrepancy.}
\value{A plot-ready data frame with one frontier per benchmark and method.
\code{l2_equal_radius()} returns the smallest common radius sufficient to
reproduce each benchmark on the evaluated grid. The plot method returns a
\code{ggplot} object. \code{l2_overall_equal_radius()} returns the smallest
common radius for which the overall ATE interval includes the null.
\code{l2_equal_radius_comparison()} returns benchmark-to-tipping ratios and a
plot-ready comparison object.}
