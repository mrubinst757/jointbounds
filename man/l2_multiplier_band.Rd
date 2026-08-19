\name{l2_multiplier_band}
\alias{l2_multiplier_band}
\alias{l2_sensitivity_band}
\title{Multiplier simultaneous bands for L2 sensitivity endpoints}
\usage{
l2_multiplier_band(
  fits, estimator = NULL, B = 1000L, conf_level = 0.95,
  multiplier = c("rademacher", "normal"), seed = 1L
)

l2_sensitivity_band(
  data, Y, A, C, X, grid, delta = c(1, 1),
  model = c("net", "separated"), bound_method = c("cs", "sharp"),
  estimator = NULL, B = 1000L, conf_level = 0.95,
  multiplier = c("rademacher", "normal"), seed = 1L,
  keep_fits = FALSE, ...
)
}
\arguments{
  \item{fits}{One cross-fitted L2 fit or a list evaluated on the same rows.}
  \item{estimator}{The score-bearing estimator; defaults to EIF for CS and DR
  for sharp fits.}
  \item{B}{Number of multiplier draws.}
  \item{conf_level}{Simultaneous confidence level.}
  \item{multiplier}{Rademacher or standard-normal multiplier distribution.}
  \item{seed}{Random seed.}
  \item{data, Y, A, C, X, delta, model}{As in \code{l2_bounds()}.}
  \item{grid}{Sensitivity-parameter grid.}
  \item{bound_method}{CS or sharp cross-fitted endpoints.}
  \item{keep_fits}{Retain all pointwise fits.}
  \item{...}{Arguments passed to the selected cross-fitted estimator.}
}
\value{Pointwise and simultaneous endpoint intervals, outward identified-set
bands, row-level scores, the multiplier critical value, and draw diagnostics.}
\description{Implements the standardized supremum multiplier process over both
ATE endpoints and every evaluated sensitivity-grid point.}
