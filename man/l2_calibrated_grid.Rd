\name{l2_calibrated_grid}
\alias{l2_calibrated_grid}
\title{Select calibrated sensitivity budgets for subsequent bound estimation}
\usage{
l2_calibrated_grid(x,
  selection = c("estimate", "pointwise_upper", "simultaneous_upper"),
  radius_x = NULL, benchmark_id = NULL)
}
\arguments{
\item{x}{Output from \code{l2_calibration_band()} using both endpoints.}
\item{selection}{Point calibration, pointwise upper confidence-limit calibration,
or simultaneous upper confidence-band calibration.}
\item{radius_x}{Optional horizontal radius values already evaluated in the grid.
NULL selects all.}
\item{benchmark_id}{Optional benchmark IDs. NULL selects all.}
}
\description{
Constructs a sensitivity grid for computing bounds at calibrated radii.
Upper-limit selections use the inner compatibility set: both upper limits
on the oriented compatibility contrasts must be nonpositive.
}
\details{
The first compatible evaluated vertical radius is selected for each requested
horizontal radius and benchmark. No interpolation or extrapolation is used.
If the selected crossing is unavailable, the function raises an error instead
of substituting zero or the largest available radius.

Use \code{pointwise_upper} at a prespecified horizontal radius. Use
\code{simultaneous_upper} when selecting among horizontal radii after
inspecting the results. Selection among omitted sets additionally requires
\code{simultaneous_scope = "all_benchmarks"} when constructing the input band.
Both signed and toward-zero targets are supported through the input band.

Upper-limit calibration enlarges the sensitivity region relative to point
calibration on the same grid, but does not itself provide confidence coverage
for subsequently estimated bounds at data-selected radii. Endpoint estimation
uncertainty and selection require separate inference.
}
\value{
A data frame with benchmark IDs, labels, and the two radius columns, suitable
as \code{grid} in \code{l2_sensitivity_band()}. Attributes record the selection,
confidence level, simultaneous scope, and inference limitation.
}
\seealso{\code{\link{l2_calibration_band}}, \code{\link{l2_sensitivity_band}}}
