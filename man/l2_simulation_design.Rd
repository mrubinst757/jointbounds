\name{l2_simulation_design}
\alias{l2_simulation_design}
\title{Reproducible designs for the continuous-outcome simulation study}
\usage{
l2_simulation_design(profile = c("main", "fast"))
}
\arguments{
  \item{profile}{\code{"main"} returns the prespecified paper design;
  \code{"fast"} returns a four-scenario smoke-test design.}
}
\value{A data frame with one row per data-generating scenario.}
\description{Constructs the scenario registry used by the paper simulation,
including identifying-assumption violations, tilt direction, missingness
prevalence and strength, concentration, overlap, nuisance complexity, and
outcome tails.}
\examples{
l2_simulation_design("fast")
}
