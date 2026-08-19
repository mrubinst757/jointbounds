\name{l2_linf_bounds}
\alias{l2_linf_bounds}
\title{Continuous-outcome L-infinity sensitivity bounds}
\usage{
l2_linf_bounds(
  data, Y, A, C, X, delta = c(1, 1), delta_R_inf = c(0, 0),
  delta_A_inf = c(0, 0), method = c("both", "outer", "sharp"),
  folds = 5L, nuisance_method = c("SuperLearner", "glm"),
  sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm",
  sl_lib_outcome = "SL.glm", basis = NULL, max_n = 1000L, seed = 1L
)
}
\arguments{
  \item{data, Y, A, C, X, delta, folds, nuisance_method, sl_lib_prop,
  sl_lib_miss, sl_lib_outcome, basis, seed}{As in \code{l2_bounds()}.}
  \item{delta_R_inf, delta_A_inf}{Arm-specific symmetric L-infinity radii.}
  \item{method}{Return the Holder outer bound, sharp empirical-sieve LP, or both.}
  \item{max_n}{Maximum respondent sample size per arm in the sharp LP.}
}
\value{Arm-specific and ATE bounds, parameters, and sharp-LP diagnostics.}
\description{Computes worst-case likelihood-ratio bounds without requiring a
bounded outcome. The sharp LP enforces conditional-normalization sieve moments,
likelihood-ratio boxes, and the nonnegative-mixture lower envelope.}

