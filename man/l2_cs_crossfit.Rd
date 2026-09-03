\name{l2_cs_crossfit}
\alias{l2_cs_crossfit}
\title{Cross-fitted plug-in and EIF estimators of L2 CS outer bounds}
\usage{
l2_cs_crossfit(
  data, Y, A, C, X, delta = c(1, 1), delta_R = c(0.2, 0.2),
  delta_A = c(0.2, 0.2), delta_M = delta_R,
  model = c("net", "separated"), folds = 5L,
  nuisance_method = c("SuperLearner", "glm"),
  sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm",
  sl_lib_outcome = "SL.glm", seed = 1L, conf_level = 0.95,
  nuisance_source = c("estimated", "oracle_error"),
  nuisance_error_rate = 0.25, nuisance_error_scale = 1,
  nuisance_error_sign = c(e = 1, rho = -1, regression = 1),
  delta_K = NULL
)
}
\arguments{
  \item{data, Y, A, C, X, delta, delta_R, delta_A, delta_K, delta_M, model, folds,
  nuisance_method, sl_lib_prop, sl_lib_miss, sl_lib_outcome, seed}{As in
  \code{l2_bounds()}.}
  \item{conf_level}{Wald confidence level.}
  \item{nuisance_source}{Cross-fitted estimated nuisances or controlled
  oracle-error perturbations for simulation. Controlled errors are constructed
  directly and use no sample splitting.}
  \item{nuisance_error_rate, nuisance_error_scale, nuisance_error_sign}{Controls
  for the oracle-error simulation mode; see \code{l2_sharp_crossfit()}.}
}
\description{Estimates the closed-form Cauchy--Schwarz L2 outer endpoints using
cross-fitted plug-in moments and the analytic influence functions for the
reference mean, second moments, and fourth interaction moment.}
\value{Arm and ATE endpoint estimates, standard errors, confidence intervals,
outward identified-set intervals, influence scores, and fold identifiers. For
data returned by \code{simulate_l2_data()}, the output also contains a
\code{true_eif} diagnostic that evaluates the EIF one-step estimator using the
known nuisance functions.}
\details{The EIF estimator is Neyman-orthogonal rather than generically
exactly doubly robust. The manuscript's product-rate conditions, plus moments
bounded away from zero for the square-root delta method, imply root-n
asymptotic linearity. A finite eighth outcome moment is a simple sufficient
condition for finite variance of the fourth-moment endpoint score.
In controlled-error mode, \eqn{e_a}, \eqn{rho_a}, \eqn{mu_a}, and the
conditional residual moments \eqn{nu_{k,a}} receive separate errors of the
specified order; no regressions are fitted.
}
