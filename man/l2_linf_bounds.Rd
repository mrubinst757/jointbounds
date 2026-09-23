\name{l2_linf_bounds}
\alias{l2_linf_bounds}
\title{Continuous-outcome L-infinity sensitivity bounds}
\usage{
l2_linf_bounds(
  data, Y, A, C, X, delta = c(1, 1), delta_R_inf = c(0, 0),
  delta_A_inf = c(0, 0), method = c("both", "outer", "sharp"),
  folds = 5L, nuisance_method = c("SuperLearner", "glm"),
  sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm",
  sl_lib_outcome = "SL.glm", basis = NULL, max_n = 1000L, seed = 1L,
  model = c("net", "separated"), delta_M_inf = delta_R_inf,
  delta_K_inf = NULL
)
}
\arguments{
  \item{data, Y, A, C, X, delta, folds, nuisance_method, sl_lib_prop,
  sl_lib_miss, sl_lib_outcome, basis, seed}{As in \code{l2_bounds()}.
  \code{delta} is used only when \code{model = "separated"}.}
  \item{model}{\code{"net"} (default) imposes the net boxes
  \eqn{|M_a - 1| \le \delta_{M,a}^\infty} and
  \eqn{|K_a - 1| \le \delta_{K,a}^\infty}. \code{"separated"} instead uses
  \eqn{|M_a - 1| \le \delta_a \pi_a(X) \delta_{R,a}^\infty} with the same
  \eqn{K_a} box.}
  \item{delta_M_inf}{Arm-specific L-infinity radii for the net missingness
  tilt \eqn{M_a} (\code{model = "net"}). Defaults to \code{delta_R_inf}.}
  \item{delta_R_inf}{Arm-specific L-infinity radii for the informative-severity
  ratio (\code{model = "separated"}).}
  \item{delta_A_inf, delta_K_inf}{Arm-specific L-infinity radii for the
  confounding tilt \eqn{K_a}; \code{delta_K_inf}, when supplied, overrides the
  legacy name \code{delta_A_inf}.}
  \item{method}{Return the Holder outer bound, sharp empirical-sieve LP, or both.}
  \item{max_n}{Maximum respondent sample size per arm in the sharp LP.}
}
\value{An object of class \code{jointbounds_linf} with arm-specific and ATE
bounds, parameters, and sharp-LP diagnostics.}
\description{Computes worst-case likelihood-ratio bounds without requiring a
bounded outcome. The outer bound uses the sign-aware Holder radius
\eqn{E|Z|\{C^+(X) + C^-(X)\}/2}, whose downward corners respect
\eqn{M_a \ge \rho_a} and \eqn{K_a \ge 0}. The sharp LP enforces
conditional-normalization sieve moments, likelihood-ratio boxes, and the
nonnegative-mixture lower envelope.}
\details{
These are plug-in point estimates of the bounds, without inference. Unlike
\code{l2_sharp_crossfit()} and \code{l2_cs_crossfit()}, no
efficient-influence-function (doubly robust) estimator is implemented for
either L-infinity interval, and no standard errors or confidence intervals are
returned.

\itemize{
  \item The outer interval is centred at the cross-fitted AIPW estimate of the
  arm mean. Its Holder radius is an inverse-probability-weighted plug-in that
  is not orthogonalized, so it is first-order sensitive to errors in the
  treatment, response, and outcome-regression estimates.
  \item The sharp interval is a pure inverse-probability-weighted plug-in: the
  linear program weights respondent outcomes by
  \eqn{1/\{\hat e_a(X)\hat\rho_a(X)\}} and does not use the outcome regression.
  It is consistent only if the treatment and response models are, and it has
  no doubly robust protection.
  \item The sharp program uses at most \code{max_n} respondents per arm, while
  the outer interval uses the full sample. Because the two intervals are
  centred on different estimates, the estimated sharp interval need not lie
  inside the estimated outer interval, even though it does in the population.
  Increase \code{max_n} to reduce this discrepancy.
}

Use these functions as a sensitivity comparison, or, with known nuisance
functions, as a population benchmark. They are not a basis for finite-sample
inference.
}
