\name{l2_calibration_simulation}
\alias{l2_calibration_simulation}
\title{Evaluate covariate-omission calibration in simulation}
\usage{
l2_calibration_simulation(
  B = 500L, n = 1000L, truth_n = 100000L, p_z = c(0L, 1L),
  max_per_pz = 2L, simultaneous = TRUE, multiplier_B = 500L,
  seed = 202609L, violation = "both",
  density_ratio_link = "bounded", confounding_strength = 0.6,
  dgp_args = list(), folds = 2L,
  nuisance_method = c("SuperLearner", "glm"), conf_level = 0.95
)
}
\arguments{
  \item{B,n,truth_n}{Replications, sample size, and reference-population size.}
  \item{p_z,max_per_pz}{Retained-set sizes and cap.}
  \item{simultaneous}{Compute simultaneous multiplier coverage.}
  \item{multiplier_B}{Multiplier draws per replication.}
  \item{seed}{Master seed.}
  \item{violation,density_ratio_link,confounding_strength,dgp_args}{DGP settings.}
  \item{folds,nuisance_method}{Nuisance-estimation settings.}
  \item{conf_level}{Confidence level.}
}
\value{A list of reference benchmarks, replication results, pointwise
operating characteristics, simultaneous coverage, and failures.}
\description{Evaluates plug-in and paired cross-fitted AIPW estimators of the
leave-covariate-out calibration discrepancies against an independent
large-population reference.}
\examples{
\dontrun{
l2_calibration_simulation(B = 50, n = 500, truth_n = 10000,
  nuisance_method = "glm")
}
}
