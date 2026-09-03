\name{l2_calibration_band}
\alias{l2_calibration_band}
\title{Simultaneous confidence bands for calibration frontiers}
\usage{
l2_calibration_band(
  surface,
  benchmarks,
  benchmark_estimator = "aipw",
  endpoint_estimator = NULL,
  active_endpoint = c("auto", "lower", "upper", "both"),
  B = 1000L,
  conf_level = 0.95,
  multiplier = c("rademacher", "normal"),
  seed = 1L,
  radius_x = NULL,
  radius_y = NULL
)
}
\arguments{
\item{surface}{Output from \code{l2_sensitivity_band(..., keep_fits = TRUE)}.}
\item{benchmarks}{Output from
\code{l2_aipw_benchmarks(..., variance = TRUE)}.}
\item{benchmark_estimator}{Benchmark estimator; combined influence-function
inference currently requires \code{"aipw"}.}
\item{endpoint_estimator}{Endpoint estimator stored in \code{surface}; by
default, \code{"eif"} for CS bounds and \code{"dr"} for sharp bounds.}
\item{active_endpoint}{Use \code{"auto"} to choose the endpoint from the
signed benchmark, or prespecify \code{"lower"}, \code{"upper"}, or
\code{"both"}.}
\item{B}{Number of multiplier draws.}
\item{conf_level}{Pointwise and simultaneous confidence level.}
\item{multiplier}{Rademacher or standard-normal multipliers.}
\item{seed}{Random seed.}
\item{radius_x, radius_y}{Sensitivity-radius columns to place on the two axes.}
}
\value{
\code{l2_calibration_band()} returns the estimated frontier, grid-level gap
estimates and confidence bands, paired influence scores, and the multiplier
critical value.
}
\description{
Combines the sensitivity-endpoint and reduced-adjustment ATE influence scores
observation by observation, constructs a multiplier band for their difference,
and inverts the band to obtain calibration-frontier confidence regions.
}
