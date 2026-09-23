\name{l2_simulation_study}
\alias{l2_simulation_study}
\alias{print.jointbounds_l2_simulation}
\alias{plot.jointbounds_l2_simulation}
\title{Run and summarize the continuous-outcome simulation study}
\usage{
l2_simulation_study(
  design = l2_simulation_design("main"), B = 1000L,
  truth_n = 100000L, seed = 202608L, conf_level = 0.95,
  l2_method = c("both", "sharp", "cs"),
  linf_method = c("outer", "both", "sharp"), folds = 3L,
  progress = interactive(), continue_on_error = TRUE,
  checkpoint_file = NULL, ...
)

\method{print}{jointbounds_l2_simulation}(x, ...)

\method{plot}{jointbounds_l2_simulation}(x,
  metric = c("coverage", "bias", "width"), ...)
}
\arguments{
  \item{design}{A design returned by \code{l2_simulation_design()}.}
  \item{B}{Monte Carlo replications per scenario.}
  \item{truth_n}{Independent Monte Carlo population size for oracle targets.}
  \item{seed}{Master seed.}
  \item{conf_level}{Wald confidence level.}
  \item{l2_method, linf_method}{Bound methods to compute.}
  \item{folds}{Number of cross-fitting folds. At least three are required when
  sharp one-step estimators are included.}
  \item{progress}{Whether to print scenario progress.}
  \item{continue_on_error}{Whether to record a failure and continue.}
  \item{checkpoint_file}{Optional RDS checkpoint path.}
  \item{x}{A fitted simulation object.}
  \item{metric}{Operating characteristic to plot.}
  \item{...}{Additional arguments passed to \code{l2_compare_dr_plugin()}.}
}
\value{\code{l2_simulation_study()} returns a
\code{jointbounds_l2_simulation} object with replication results, endpoint and
set-coverage summaries, width and mechanism comparisons, and failure
diagnostics. The plot method returns a ggplot object.}
\description{Runs the prespecified paper simulation using minimal oracle
sensitivity radii, independent oracle endpoint references, cross-fitted sharp
and Cauchy--Schwarz estimators, and L-infinity comparisons.}
\examples{
\dontrun{
s <- l2_simulation_study(l2_simulation_design("fast"), B = 100,
  truth_n = 20000, folds = 3)
plot(s, "coverage")
}
}
