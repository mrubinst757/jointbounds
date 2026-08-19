\name{l2_oracle_sensitivity_parameters}
\alias{l2_oracle_sensitivity_parameters}
\title{Minimal oracle sensitivity parameters for the continuous DGP}
\usage{
l2_oracle_sensitivity_parameters(data, nodes = 200L)
}
\arguments{
  \item{data}{Data returned by \code{simulate_l2_data()}.}
  \item{nodes}{Number of deterministic normal quadrature nodes.}
}
\value{A list containing the arm-specific minimal \code{delta},
\code{delta_M}, \code{delta_R}, \code{delta_K} (also returned under the legacy
name \code{delta_A}), their L-infinity analogues, the covariate-level informative-prevalence
ratios, the calibration-sample empirical envelope, and a numerical-integration
note.}
\description{
Calculates the smallest prevalence envelope and global informative-missingness
L2 radius under the simulator's known law. Because its Gaussian covariate has
unbounded support, the population prevalence envelope is one; the empirical
maximum is also returned. Treatment depends on the shared potential-outcome
disturbance, so the confounding radius is nonzero by default. Setting
\code{confounding_strength = 0} recovers a zero radius.
}
\examples{
d <- simulate_l2_data(500, seed = 1)
l2_oracle_sensitivity_parameters(d)
}
