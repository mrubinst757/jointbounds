\name{l2_calibrate}
\alias{l2_calibrate}
\title{Calibrate global L2 radii from benchmark likelihood-ratio tilts}
\usage{
l2_calibrate(R = NULL, K = NULL, M = NULL, weights = NULL,
  benchmark = c("user_supplied", "observed_data"), conditional_on = NULL)
}
\arguments{
  \item{R, K, M}{Benchmark likelihood-ratio vectors. At least one is required.}
  \item{weights}{Optional nonnegative sampling weights.}
  \item{benchmark}{Whether tilts are user supplied or estimated from the observed data.}
  \item{conditional_on}{Use \code{"exchangeability"} with an observed-data \code{R} benchmark, or \code{"MAR"} with an observed-data \code{K} benchmark.}
}
\value{A named vector containing the applicable delta_R, delta_M, and delta_K radii. The legacy delta_A name is also returned for compatibility.}
\description{
Translates estimated or elicited benchmark tilts into the global L2 scales used
by the continuous-outcome sensitivity models. The confounding radius uses M as
the factual-law weight when it is supplied. Reduced-form M calibrates combined
distortion. Separate observed-data R and K benchmarks are not jointly identified.
}
\examples{
R <- c(.8, 1, 1.2)
K <- c(.9, 1, 1.1)
M <- c(.95, 1, 1.05)
l2_calibrate(R, K, M)
}
