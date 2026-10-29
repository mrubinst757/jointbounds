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
  radius_y = NULL,
  calibration = c("toward_zero", "signed"),
  alpha = 1,
  simultaneous_scope = c("all_benchmarks", "per_benchmark"),
  interpolate = FALSE
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
\item{active_endpoint}{\code{"auto"} and \code{"both"} require both
endpoint inequalities for interval inclusion. \code{"lower"} or
\code{"upper"} requests only the corresponding one-sided crossing.}
\item{B}{Number of multiplier draws.}
\item{conf_level}{Pointwise and simultaneous confidence level.}
\item{multiplier}{Rademacher or standard-normal multipliers.}
\item{seed}{Random seed.}
\item{radius_x, radius_y}{Sensitivity-radius columns to place on the two axes.}
\item{calibration}{\code{"toward_zero"} uses
full - sign(full) * alpha * abs(reduced - full). \code{"signed"} uses
full + alpha * (reduced - full).}
\item{interpolate}{FALSE uses evaluated grid radii. TRUE interpolates descriptively and emits a warning; finite-grid coverage does not control interpolation error.}
\item{alpha}{One finite positive benchmark multiplier.}
\item{simultaneous_scope}{\code{"all_benchmarks"} (the existing default)
covers all supplied benchmarks, endpoints, and grid points jointly.
\code{"per_benchmark"} covers the endpoints and full grid separately
for each benchmark, without joint coverage across omitted sets.
Pointwise frontier intervals are unaffected.}
}
\value{
\code{l2_calibration_band()} returns the estimated frontier, grid-level gap
estimates and confidence bands, paired influence scores, and the multiplier
critical value. The \code{compatibility} table records estimated, inner, and outer compatible grid sets for the requested endpoint inequalities.
With \code{simultaneous_scope = "per_benchmark"}, \code{critical_value}
is indexed by benchmark ID; for multiple benchmarks, \code{suprema} has
one column per benchmark. The grid-level \code{gaps} table records the
critical value used for each row.
}
\description{
Combines the sensitivity-endpoint and transformed target influence scores
observation by observation, constructs a multiplier band for their difference,
and inverts the band to obtain calibration-frontier confidence regions.
For toward-zero calibration the target score is the full-X score minus
alpha * sign(full) * sign(reduced - full) times the paired reduced-minus-full
score. Exact zero references or omission discrepancies are rejected because
ordinary target-EIF inference is nonregular there. Near these boundaries,
pointwise normal approximations are not uniform. Both paired full-X and reduced-Z scores are required.
For signed calibration the target score is (1-alpha) times the full-X score
plus alpha times the reduced-Z score; zero references and discrepancies
do not require exclusion for this linear target.
Frontier intervals pointwise in radius_x use a Bonferroni correction across
all searched radius_y values and endpoint inequalities for that benchmark.
Simultaneous bands use a multiplier maximum across the full grid and both
endpoints, including all benchmarks only when \code{simultaneous_scope}
is \code{"all_benchmarks"}. Both options use common observation multipliers.
Crossings use the first compatible evaluated grid radius by default. Optional interpolation has no between-grid coverage guarantee.
An NA edge means that its crossing was not found within the grid;
it is not replaced by the other edge. For coherent calibration, use the
zero-budget surface estimate and its influence scores as the full-X benchmark.
}
