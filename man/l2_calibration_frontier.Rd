\name{l2_calibration_frontier}
\alias{l2_calibration_frontier}
\title{Invert a sensitivity grid to obtain a joint calibration frontier}
\usage{
l2_calibration_frontier(
  x,
  benchmark,
  type = c("target", "discrepancy"),
  reference = NULL,
  radius_x = NULL,
  radius_y = NULL
)
}
\arguments{
  \item{x}{Output from \code{l2_sensitivity_grid()} or its results data frame.}
  \item{benchmark}{One observed benchmark target or discrepancy.}
  \item{type}{Whether \code{benchmark} is a target or a discrepancy from the reference.}
  \item{reference}{Optional reference value. By default, uses \code{estimate_reference}.}
  \item{radius_x, radius_y}{Names of the two sensitivity-radius columns.}
}
\value{
A data frame containing the minimum compatible second radius at every first
radius and method represented on the grid.
}
\description{
Inverts the joint sensitivity bounds. A radius pair is compatible when its
bound contains the observed benchmark. The frontier therefore describes
jointly sufficient sensitivity budgets and their possible substitution, not an
estimate of the true mechanism-specific radii. Inverting CS bounds gives an
outer approximation; inverting sharp sieve bounds gives the corresponding
empirical sieve frontier.
}
\examples{
grid <- expand.grid(delta_R = c(0, .2), delta_A = c(0, .1, .2))
grid$method <- "sharp_sieve"
grid$estimate_reference <- 1
grid$lower <- 1 - grid$delta_R - grid$delta_A
grid$upper <- 1 + grid$delta_R + grid$delta_A
l2_calibration_frontier(grid, .15, type = "discrepancy", radius_x = "delta_R")
}
