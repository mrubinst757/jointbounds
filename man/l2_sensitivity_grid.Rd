\name{l2_sensitivity_grid}
\alias{l2_sensitivity_grid}
\alias{l2_tipping_frontier}
\alias{plot.marbounds_l2_frontier}
\title{Continuous-outcome L2 sensitivity surfaces and tipping frontiers}
\usage{
l2_sensitivity_grid(data, Y, A, C, X, grid, delta = c(1, 1),
  model = c("net", "separated"), method = c("both", "cs", "sharp"),
  folds = 5L, nuisance_method = c("SuperLearner", "glm"),
  sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm",
  sl_lib_outcome = "SL.glm", basis = NULL, seed = 1L,
  control = list())

l2_tipping_frontier(x, threshold = 0, endpoint = c("lower", "upper"), plot = FALSE)

\method{plot}{marbounds_l2_frontier}(x, ...)
}
\arguments{
  \item{data, Y, A, C, X, delta, model, method, folds, nuisance_method, sl_lib_prop, sl_lib_miss, sl_lib_outcome, basis, seed, control}{As in \code{l2_bounds()}.}
  \item{grid}{Data frame containing delta_M and delta_K for the primary net model, or delta_R and delta_K for the separated model. Legacy delta_A is accepted. Arm-specific columns may use suffixes _0 and _1.}
  \item{x}{A grid object or its results data frame.}
  \item{threshold}{Null threshold defining the tipping frontier.}
  \item{endpoint}{Endpoint used to define tipping.}
  \item{plot}{Draw the frontier before returning it.}
  \item{...}{Additional base-graphics arguments.}
}
\value{\code{l2_sensitivity_grid()} returns a grid object; \code{l2_tipping_frontier()} returns a plot-ready frontier data frame.}
\description{Evaluates arm and ATE bounds over a two-dimensional sensitivity
grid and extracts the first threshold crossing, using linear interpolation
between adjacent grid points. Groups that do not cross are explicitly marked.}
\examples{
d <- simulate_l2_data(300, seed = 2)
g <- expand.grid(delta_R = c(0, .2), delta_A = c(0, .1))
s <- l2_sensitivity_grid(d, "Y", "A", "C", c("X1", "X2"),
  grid = g, delta = .5, method = "cs", folds = 2)
f <- l2_tipping_frontier(s)
if (any(f$tipped_on_grid)) plot(f)
}
