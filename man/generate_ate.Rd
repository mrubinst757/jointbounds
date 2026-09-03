\name{generate_ate}
\alias{generate_ate}
\alias{print.marbounds_ate}
\title{Cross-fitted doubly robust ATE under MAR and no unmeasured confounding}
\usage{
generate_ate(
  data, Y, A, C, X, folds = 5L,
  nuisance_method = c("SuperLearner", "glm"),
  sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm",
  sl_lib_outcome = "SL.glm",
  outcome_family = c("gaussian", "binomial"), seed = 1L,
  trim = 1e-4, conf_level = 0.95, diagnostic_control = list(),
  keep_nuisance = FALSE
)

\method{print}{marbounds_ate}(x, ...)
}
\arguments{
  \item{data}{A data frame.}
  \item{Y, A, C}{Column names. \code{A} and \code{C} must be binary and
  \code{C = 1} denotes a missing outcome.}
  \item{X}{Character vector of baseline covariates. Use \code{character()}
  for an intercept-only analysis.}
  \item{folds}{Number of cross-fitting folds.}
  \item{nuisance_method}{Either \code{"SuperLearner"} or \code{"glm"}.}
  \item{sl_lib_prop, sl_lib_miss, sl_lib_outcome}{SuperLearner libraries for
  the treatment, observation, and outcome regressions.}
  \item{outcome_family}{Either \code{"gaussian"} or \code{"binomial"}.}
  \item{seed}{Random seed used to construct folds.}
  \item{trim}{Lower and upper truncation level for estimated treatment and
  response probabilities.}
  \item{conf_level}{Wald confidence level.}
  \item{diagnostic_control}{Optional list controlling the consolidated warning
  and thresholds for large weights and low weighted effective sample sizes.}
  \item{keep_nuisance}{Return row-level cross-fitted nuisance estimates.}
  \item{x}{An object returned by \code{generate_ate()}.}
  \item{...}{Unused.}
}
\value{An object of class \code{marbounds_ate} containing the ATE and
arm-specific AIPW estimates, standard errors, Wald intervals, centered
influence scores, and positivity diagnostics.}
\description{
Estimates the point-identified average treatment effect under consistency,
positivity, MAR, and no unmeasured confounding. Treatment, response, and outcome
nuisance functions are estimated out of fold using SuperLearner or GLMs.
}
\details{
For arm \eqn{a}, the cross-fitted score is
\deqn{\mu_a(X)+\frac{I(A=a,C=0)}{e_a(X)\rho_a(X)}\{Y-\mu_a(X)\}.}
The estimator is consistent if the outcome regressions are correct or if the
arm-specific joint treatment-response scores \eqn{e_a(X)\rho_a(X)} are
correct. Wald standard errors use the empirical variance of the centered,
cross-fitted influence score.
}
\examples{
d <- simulate_l2_data(500, violation = "none", seed = 1)
fit <- generate_ate(d, "Y", "A", "C", c("X1", "X2"),
  folds = 2, nuisance_method = "glm")
fit$ate
}
