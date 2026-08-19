\name{l2_prevalence_bounds}
\alias{l2_prevalence_bounds}
\title{MAR-only L2 prevalence bounds for binary alternative estimands}
\usage{
l2_prevalence_bounds(
  data, Y, A, C, X, estimand = c("composite", "sde"),
  kappa = c(0, 0), folds = 5L,
  nuisance_method = c("SuperLearner", "glm"),
  sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm",
  sl_lib_outcome = "SL.glm", seed = 1L, conf_level = 0.95,
  keep_nuisance = FALSE
)
}
\arguments{
  \item{data, Y, A, C, X, folds, nuisance_method, sl_lib_prop, sl_lib_miss,
  sl_lib_outcome, seed}{As in \code{l2_bounds()}; Y must be binary.}
  \item{estimand}{The composite effect or separable direct effect.}
  \item{kappa}{One or two L2 radii for the informative-prevalence ratio.}
  \item{conf_level}{Wald confidence level.}
  \item{keep_nuisance}{Retain nuisances and empirical optimal tilts.}
}
\value{Bounds, influence-function standard errors and intervals, row scores,
and capped-support diagnostics.}
\description{Implements the capped L2 support-function formulas under
conditional exchangeability. The influence score includes the envelope term
from the active norm constraint because the norm is defined under the observed
covariate distribution.}
