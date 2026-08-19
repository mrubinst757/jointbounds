\name{l2_bootstrap}
\alias{l2_bootstrap}
\title{Bootstrap continuous-outcome L2 sensitivity bounds}
\usage{
l2_bootstrap(data, Y, A, C, X, delta = c(1, 1),
  delta_R = c(0, 0), delta_A = c(0, 0), delta_M = delta_R,
  model = c("net", "separated"), method = c("cs", "sharp"),
  folds = 5L, nuisance_method = c("SuperLearner", "glm"),
  sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm",
  sl_lib_outcome = "SL.glm", basis = NULL, seed = 1L, control = list(),
  B = 200L, conf_level = 0.95, delta_K = NULL)
}
\arguments{
  \item{data, Y, A, C, X, delta, delta_R, delta_A, delta_K, delta_M, model,
  method, folds, nuisance_method, sl_lib_prop, sl_lib_miss, sl_lib_outcome,
  basis, seed, control}{As in \code{l2_bounds()}.}
  \item{B}{Number of nonparametric bootstrap samples.}
  \item{conf_level}{Confidence level for percentile intervals.}
}
\value{The original fit, bootstrap endpoint replicates, standard errors, percentile intervals, and the number of successful fits.}
\description{
Refits the full nuisance and bound procedure in nonparametric bootstrap samples.
For nonregular sharp bounds this is a practical diagnostic rather than a claim
of ordinary bootstrap validity at active-set transitions.
}
\examples{
\dontrun{
d <- simulate_l2_data(500, seed = 1)
b <- l2_bootstrap(d, "Y", "A", "C", c("X1", "X2"),
  delta = .5, delta_R = .2, delta_A = .1, B = 100)
b$ci
}
}
