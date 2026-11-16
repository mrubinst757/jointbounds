\name{l2_cs_crossfit}
\alias{l2_cs_crossfit}
\title{Cross-fitted plug-in and EIF estimators of L2 CS outer bounds}
\usage{
l2_cs_crossfit (data, Y, A, C, X, delta = c(1, 1), delta_R = c(0.2, 0.2), delta_A = c(0.2, 0.2),
    delta_M = delta_R, model = c("net", "separated"), folds = 5L, nuisance_method = c("SuperLearner",
        "glm"), sl_lib_prop = sl_libs, sl_lib_miss = sl_libs, sl_lib_outcome = sl_libs, seed = 1L,
    conf_level = 0.95, nuisance_source = c("estimated", "oracle_error"), nuisance_error_rate = 0.25,
    nuisance_error_scale = 1, nuisance_error_sign = c(e = 1, rho = -1, regression = 1), delta_K = NULL,
    propensity_clip = clip, moment_strategy = c("three_way", "shared_training"), sl_libs = "SL.glm",
    response_clip = clip, clip = 1e-04, nuisance_fit = NULL, learner_fallback = c("error", "glm"))
}
\arguments{
\item{learner_fallback}{On SuperLearner failure, stop with the original error by default (\code{"error"}). The explicit \code{"glm"} option permits a substitute GLM with a warning. Initial-fit reuse inherits its saved policy. Older caches default to error for new regressions; cached predictions are not recomputed. The legacy benchmark strategy supports error only.}
\item{nuisance_fit}{Reusable output from \code{l2_initial_fit()}. Supplies cached structural predictions, folds, seed, learners, and clipping, overriding downstream settings for these. Requires the same ordered data. Sharp estimation requires the three-way strategy. Dual fitting and pseudo-outcome regressions still run separately.}
\item{clip}{Common two-sided clipping threshold, default 1e-4. Probabilities are restricted to [clip, 1-clip]; no observations are removed.}
\item{propensity_clip, response_clip}{Overrides for treatment and response probabilities. Both must be finite scalars strictly between zero and 0.5. Reused reference fits supply both settings; older fits use response clipping at 1e-4.}

\item{sl_libs}{Common SuperLearner library for all nuisance regressions. Model-specific \code{sl_lib_prop}, \code{sl_lib_miss}, and \code{sl_lib_outcome} override it.}
  \item{moment_strategy}{The default \code{"three_way"} fits structural nuisances,
  residual-moment regressions, and evaluation scores on disjoint folds (at least
  three). \code{"shared_training"} reproduces the earlier implementation,
  reusing structural training observations for residual regressions; this
  requires additional stability assumptions. Ignored in oracle-error mode.}

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
\value{The \code{influence_covariance} and \code{estimate_covariance} lists contain joint lower/upper covariance matrices by estimator; the latter divides by the sample size. Arm and ATE endpoint estimates, standard errors, confidence intervals,
outward identified-set intervals, influence scores, fold identifiers, and
\code{crossfit_roles} describing the training and evaluation folds. For
data returned by \code{simulate_l2_data()}, the output also contains a
\code{true_eif} diagnostic that evaluates the EIF one-step estimator using the
known nuisance functions.}
\details{The \code{moment_diagnostics} table flags nonpositive required estimated moments and truncation to zero. Smooth EIF inference is not justified at those boundary cases; positive estimates do not establish population regularity. The EIF estimator is Neyman-orthogonal rather than generically
exactly doubly robust. The manuscript's product-rate conditions, plus positive population moments for the square-root delta method, imply root-n
asymptotic linearity. A finite eighth outcome moment is a simple sufficient
condition for finite variance of the fourth-moment endpoint score.
In controlled-error mode, \eqn{e_a}, \eqn{rho_a}, \eqn{mu_a}, and the
conditional residual moments \eqn{nu_{k,a}} receive separate errors of the
specified order; no regressions are fitted.
}
