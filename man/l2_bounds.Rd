\name{l2_bounds}
\alias{l2_bounds}
\alias{print.jointbounds_l2}
\title{Continuous-outcome global L2 sensitivity bounds}
\usage{
l2_bounds (data, Y, A, C, X, delta = c(1, 1), delta_R = c(0, 0), delta_A = c(0, 0), delta_M = delta_R,
    model = c("net", "separated"), method = c("both", "cs", "sharp"), folds = 5L, nuisance_method = c("SuperLearner",
        "glm"), sl_lib_prop = sl_libs, sl_lib_miss = sl_libs, sl_lib_outcome = sl_libs, basis = NULL,
    seed = 1L, control = list(), diagnostic_control = list(), keep_nuisance = FALSE, return_true_bounds = FALSE,
    delta_K = NULL, sl_libs = "SL.glm", propensity_clip = clip, response_clip = clip, clip = 1e-04)

\method{print}{jointbounds_l2}(x, ...)
}
\arguments{
\item{clip}{Common two-sided clipping threshold, default 1e-4. Probabilities are restricted to [clip, 1-clip]; no observations are removed.}
\item{propensity_clip, response_clip}{Overrides for treatment and response probabilities. Both must be finite scalars strictly between zero and 0.5. Reused reference fits supply both settings; older fits use response clipping at 1e-4.}

\item{sl_libs}{Common SuperLearner library for all nuisance regressions. Model-specific \code{sl_lib_prop}, \code{sl_lib_miss}, and \code{sl_lib_outcome} override it.}
  \item{data}{A data frame.}
  \item{Y, A, C}{Column names. A and C are binary and C = 1 denotes missing.}
  \item{X}{Character vector of baseline covariate names.}
  \item{delta}{One value or arm-specific values bounding the informative fraction of missingness.}
  \item{delta_R, delta_A}{Global L2 radii for informative-missingness severity and unmeasured confounding.}
  \item{delta_K}{Preferred alias for the confounding density-ratio radius.}
  \item{delta_M}{Global L2 radius for the net missingness distortion in the reduced-form model.}
  \item{model}{The separated prevalence/severity model or the reduced-form net-M model.}
  \item{method}{Closed-form Cauchy--Schwarz outer bounds, empirical sharp sieve-dual bounds, or both.}
  \item{folds}{Number of folds for cross-fitted nuisance estimation.}
  \item{nuisance_method}{Cross-fitted SuperLearner or GLM nuisance estimation.}
  \item{sl_lib_prop, sl_lib_miss, sl_lib_outcome}{SuperLearner libraries for treatment, missingness, and continuous-outcome regressions.}
  \item{basis}{One-sided formula for the sharp approximation. NULL balances an intercept and additive covariates in X (means for numeric covariates; indicators for factors). For example, \code{~ X1 + X2 + I(X1^2) + X1:X2} or \code{~ splines::ns(X1, df = 3) + X2}. Only X covariates are allowed; an intercept is required; dot expands to X. Cross-fitting learns transformations and scaling on training folds.}
  \item{seed}{Random seed for fold construction.}
  \item{control}{Control list passed to \code{optim()} for sharp bounds.}
  \item{diagnostic_control}{Optional list controlling sharp-solver
  diagnostics. Set \code{warn = FALSE} to suppress the consolidated warning
  or \code{compute_hessian = FALSE} to skip local-curvature diagnostics.
  Numerical tolerances for normalization, budgets, the primal--dual gap,
  multipliers, curvature, conditioning, and boundary saturation may also be
  changed.}
  \item{keep_nuisance}{Whether row-level nuisance estimates are returned.}
  \item{return_true_bounds}{Return the requested oracle Monte Carlo CS and/or
  sharp-sieve bounds; requires known-DGP columns from
  \code{simulate_l2_data()}.}
  \item{x}{An object returned by \code{l2_bounds()}.}
  \item{...}{Unused.}
}
\value{An object containing arm-specific bounds, ATE bounds, raw diagnostics,
an endpoint-level \code{diagnostic_table}, a one-row
\code{diagnostic_summary}, sensitivity parameters, and an empirical
fourth-moment check.}
\description{
Computes continuous-outcome sensitivity bounds under global L2 restrictions.
Nuisance probabilities and outcome regressions are estimated as functions of X
using cross-fitting. The reference mean is augmented inverse-probability
weighted; nonlinear CS moments are cross-fitted plug-in inverse-weighted moments.
The CS method returns the paper's closed-form outer radius. The sharp method
solves the empirical sieve-dual problem and reports optimization, divergence,
conditional-normalization, budget-activity, multiplier, boundary, and local-
curvature diagnostics. A single warning is printed when any diagnostic is
flagged; the detailed reason remains available in the returned tables.
}
\details{
The sharp result is exact for the empirical conditional-normalization sieve
specified by \code{basis}; increasing the basis enriches the approximation.
A finite empirical fourth moment does not verify the population moment
assumption and is returned only as a diagnostic.
}
\examples{
d <- simulate_l2_data(500, seed = 1)
fit <- l2_bounds(d, "Y", "A", "C", c("X1", "X2"),
  delta = 0.5, delta_R = 0.25, delta_A = 0.15,
  method = "cs", folds = 3)
fit$ate
}
