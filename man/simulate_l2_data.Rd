\name{simulate_l2_data}
\alias{simulate_l2_data}
\title{Simulate continuous outcomes with mixed informative missingness}
\usage{
simulate_l2_data(
  n = 1000L, seed = 1L,
  violation = c("both", "mar", "nuc", "none"),
  density_ratio_link = c("bounded", "linear"),
  confounding_strength = 0.6,
  informative_strength = 0.45,
  informative_intercept = -2.4,
  noninformative_intercept = -1.8,
  treatment_scale = 1,
  covariate_complexity = c("linear", "nonlinear"),
  informative_concentration = c("diffuse", "high_variance"),
  outcome_distribution = c("normal", "t"),
  outcome_df = 12,
  outcome_scale = 1
)
}
\arguments{
  \item{n}{Sample size.}
  \item{seed}{Random seed.}
  \item{violation}{Which assumptions are violated: \code{"mar"},
  \code{"nuc"}, \code{"both"}, or \code{"none"}.}
  \item{density_ratio_link}{\code{"bounded"} uses tanh latent effects and has
  finite oracle L-infinity density-ratio radii; \code{"linear"} has infinite
  L-infinity radii when the corresponding mechanism is active.}
  \item{confounding_strength}{Strength of treatment selection on the shared
  latent potential-outcome disturbance. Set to zero for conditional
  exchangeability.}
  \item{informative_strength}{Outcome-dependent informative-missingness tilt.}
  \item{informative_intercept}{Baseline informative-missingness log odds.}
  \item{noninformative_intercept}{Baseline noninformative-missingness log odds.}
  \item{treatment_scale}{Multiplier controlling treatment overlap.}
  \item{covariate_complexity}{Linear or nonlinear nuisance functions.}
  \item{informative_concentration}{Diffuse or covariate-concentrated outcome tilt.}
  \item{outcome_distribution}{Normal or standardized Student t errors.}
  \item{outcome_df}{Degrees of freedom, greater than eight, for t errors.}
  \item{outcome_scale}{Marginal standard deviation of the outcome error.}
}
\value{A data frame containing observed data, both potential outcomes,
potential informative-type and missingness indicators, and true conditional
DGP quantities for outcomes, treatment, informative type, and missingness.}
\description{
Generates Gaussian or finite-eighth-moment Student t potential outcomes.
Informatively missing units are
always missing; noninformative missingness is independent of the outcome
conditional on treatment and covariates.
Treatment depends on a disturbance shared by both potential outcomes, producing
nonzero residual confounding by default while retaining an analytically
tractable likelihood ratio. The \code{violation} argument can turn informative
missingness and residual confounding on independently.
For continuous \code{Y(a)}, the output reports its conditional mean, standard
deviation, and density at the realized potential outcome rather than calling a
density a probability. Missingness probabilities are returned both conditional
on \code{Y(a),X} and marginalized over \code{Y(a)} given X.
}
\examples{
d <- simulate_l2_data(200, seed = 4)
mean(d$C)
mean(d$Y1 - d$Y0)
}
