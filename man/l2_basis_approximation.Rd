\name{l2_basis_approximation}
\alias{l2_basis_approximation}
\title{Compare sharp-sieve bounds across normalization bases}
\usage{
l2_basis_approximation(data, parameters = NULL, bases = NULL, control = list())
}
\arguments{
  \item{data}{Data returned by \code{simulate_l2_data()}.}
  \item{parameters}{Optional oracle sensitivity parameters.}
  \item{bases}{Named list of normalization-basis formulas.}
  \item{control}{Optimizer control list.}
}
\value{A data frame containing endpoints, width, and excess width relative to
the richest supplied basis.}
\description{Quantifies the outer-relaxation cost of enforcing conditional
normalization only against a finite basis. The default sequence contains an
intercept-only, linear, and cubic/interacted basis.}
\examples{
\dontrun{
d <- simulate_l2_data(10000)
l2_basis_approximation(d)
}
}
