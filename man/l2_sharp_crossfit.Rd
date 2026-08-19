\name{l2_sharp_crossfit}
\alias{l2_sharp_crossfit}
\title{Cross-fitted plug-in and doubly robust sharp L2 endpoints}
\usage{
l2_sharp_crossfit(
  data, Y, A, C, X, delta = c(1, 1), delta_R = c(0.2, 0.2),
  delta_A = c(0.2, 0.2), delta_M = delta_R,
  model = c("net", "separated"), folds = 5L, basis = NULL,
  nuisance_method = c("SuperLearner", "glm"), sl_lib_prop = "SL.glm",
  sl_lib_miss = "SL.glm", sl_lib_outcome = "SL.glm", seed = 1L,
  control = list(), conf_level = 0.95, boundary_tol = 1e-7,
  nuisance_source = c("estimated", "oracle_error"),
  nuisance_error_rate = 0.25, nuisance_error_scale = 1,
  nuisance_error_sign = c(e = 1, rho = -1, regression = 1),
  delta_K = NULL
)
}
\arguments{
  \item{data, Y, A, C, X, delta, delta_R, delta_A, delta_K, delta_M, model, folds, basis,
  nuisance_method, sl_lib_prop, sl_lib_miss, sl_lib_outcome, seed, control}{As
  in \code{l2_bounds()}. Both the separated and reduced-form net-M models are
  supported.}
  \item{conf_level}{Wald confidence level.}
  \item{boundary_tol}{Tolerance for declaring the mixture lower envelope active.}
  \item{nuisance_source}{Use fitted nuisance models, or controlled perturbations
  of oracle nuisance functions. The latter requires simulated DGP columns.}
  \item{nuisance_error_rate}{Exponent alpha in the error magnitude
  \eqn{n^{-\alpha}}. Set this to zero for fixed misspecification.}
  \item{nuisance_error_scale}{Scalar or named vector with entries \code{e},
  \code{rho}, and \code{regression}; the last controls errors in the
  conditional loss and derivative regressions.}
  \item{nuisance_error_sign}{Independently configurable signs, each equal to
  -1 or 1, for the \code{e}, \code{rho}, and \code{regression} errors.}
}
\value{A list containing arm and ATE endpoint estimates, standard errors,
confidence intervals, outward identified-set intervals, influence scores,
fold identifiers, and dual and active-set diagnostics. For simulated DGP data,
\code{true_eif} uses exact
sampling scores and quadrature-based conditional EIF regressions.}
\description{
Fits the empirical sieve dual on training folds and evaluates it on held-out
folds. The plug-in estimator uses inverse respondent weights. The augmented
estimator estimates the conditional dual-loss regression and the conditional
derivative regressions for the structural propensity and response probability.
The response derivative includes the KKT multiplier when the moving mixture
lower envelope binds.
}
\details{
Double robustness applies to evaluation of a fixed correctly specified dual
loss: either its conditional respondent-law regression or the respondent
sampling score must be consistently estimated. It does not protect against an
inconsistent dual optimizer or inconsistent structural nuisance values inside
the dual loss. Wald intervals also require a unique regular optimizer and a
stable active set.
For root-n inference the manuscript requires the second-order products
\eqn{r_m r_\gamma + r_\gamma(r_e+r_\rho) + r_e r_{d_e} +
r_\rho r_{d_\rho} + r_e^2 + r_\rho^2 + r_e r_\rho = o_p(n^{-1/2})},
together with dual excess risk and numerical primal--dual gap of
\eqn{o_p(n^{-1/2})}. Uniform \eqn{o_p(n^{-1/4})} nuisance rates are sufficient
but not necessary.
}
\examples{
d <- simulate_l2_data(400, seed = 2)
fit <- l2_sharp_crossfit(d, "Y", "A", "C", c("X1", "X2"),
  delta = .5, delta_R = .2, delta_A = .1, folds = 2,
  nuisance_method = "glm", control = list(maxit = 200))
fit$ate
}
