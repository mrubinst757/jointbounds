\name{l2_compare_dr_plugin}
\alias{l2_compare_dr_plugin}
\title{Compare plug-in and DR sharp-bound estimators in simulation}
\usage{
l2_compare_dr_plugin(
  B = 200L, n = 1000L,
  violation = c("both", "mar", "nuc", "none"),
  density_ratio_link = c("bounded", "linear"),
  confounding_strength = 0.6, truth = NULL,
  truth_n = 20000L,
  seed = 1L, conf_level = 0.95,
  nuisance_simulation = c("oracle_error", "estimated"),
  parameter_evaluation = c("auto", "oracle", "specified"),
  linf_method = c("outer", "sharp", "both"),
  l2_method = c("sharp", "cs", "both"), linf_max_n = 1000L,
  dgp_args = list(), ...
)
}
\arguments{
  \item{B}{Number of Monte Carlo replications.}
  \item{n}{Sample size per replication.}
  \item{violation}{Which identifying assumptions the simulated DGP violates.}
  \item{density_ratio_link}{Bounded or linear latent tilts, passed to the data
  generator.}
  \item{confounding_strength}{Strength of treatment selection on the shared
  potential-outcome disturbance; passed to \code{simulate_l2_data()}.}
  \item{truth}{Optional lower and upper reference endpoints.}
  \item{truth_n}{Reference-sample size when \code{truth} is not supplied.}
  \item{seed}{Random seed.}
  \item{conf_level}{Wald confidence level.}
  \item{nuisance_simulation}{\code{"oracle_error"} directly perturbs known DGP
  nuisance functions without cross-fitting; \code{"estimated"} cross-fits the
  requested learners.}
  \item{parameter_evaluation}{\code{"auto"} uses oracle parameters if no
  sensitivity parameters are supplied; \code{"oracle"} always uses them;
  \code{"specified"} evaluates the supplied or function-default values.}
  \item{linf_method}{Compute the Holder outer bound, the sharp empirical-sieve
  linear program, or both. The sharp options require \code{lpSolve}.}
  \item{l2_method}{Compute sharp-sieve endpoints, Cauchy--Schwarz outer
  endpoints, or both. CS results include plug-in and EIF estimators.}
  \item{linf_max_n}{Maximum calibration-sample size used in the sharp linear
  program.}
  \item{dgp_args}{Named list of additional arguments passed to
  \code{simulate_l2_data()}.}
  \item{...}{Arguments passed to \code{l2_sharp_crossfit()}.}
}
\value{A list containing replication-level results, endpoint bias, RMSE,
empirical standard deviations, mean standard errors, Wald coverage, true-ATE
point-bound and outward-Wald coverage, failure rates, the endpoint reference,
the minimal oracle sensitivity
parameters, and the parameters actually evaluated. It also returns oracle
L-infinity radii, an L-infinity/L1 Holder outer bound,
an optional sharp L-infinity empirical-sieve interval, and a width comparison
with the L2 sharp-sieve interval. The L-infinity intervals are population
benchmarks computed from the known DGP. The sharp one is a
\code{linf_max_n}-draw inverse-weighted LP discretization, so it can fall
slightly outside the outer interval. No DR estimator or inference is provided
for either L-infinity interval; see \code{l2_linf_bounds()}.}
\description{
Runs the continuous-outcome simulator and compares inverse-respondent-weighted
plug-in evaluation, error-contaminated DR evaluation, and the true EIF without
nuisance error. If analytic sharp endpoints
are unavailable, a large Monte Carlo empirical-sieve calculation supplies the
reference and is labeled accordingly.
}
\examples{
\dontrun{
out <- l2_compare_dr_plugin(B = 100, n = 1000, truth_n = 30000,
  delta = .5, delta_R = .2, delta_A = .1, folds = 2,
  nuisance_method = "glm")
out$summary
out$bound_coverage
}
}
