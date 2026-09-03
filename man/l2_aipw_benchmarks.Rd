\name{l2_aipw_benchmarks}
\alias{l2_aipw_benchmarks}
\title{Estimate simultaneous leave-covariate-out AIPW benchmarks}
\usage{
l2_aipw_benchmarks(
  data, Y, A, C, X, p_z = NULL, max_per_pz = 10L, subsets = NULL,
  subset_selection = c("random", "first"),
  estimator = c("both", "aipw", "plugin"), variance = FALSE,
  conf_level = 0.95, folds = 5L,
  nuisance_method = c("SuperLearner", "glm"),
  sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm",
  sl_lib_outcome = "SL.glm", seed = 1L
)
}
\arguments{
  \item{data}{A data frame.}
  \item{Y, A, C}{Outcome, treatment, and missingness-indicator column names.}
  \item{X}{Character vector containing the full covariate set.}
  \item{p_z}{Integer vector of retained-set sizes; zero is allowed.}
  \item{max_per_pz}{Maximum number of retained sets at each size.}
  \item{subsets}{Optional list of explicit retained covariate sets.}
  \item{subset_selection}{How to select sets when their number exceeds the cap.}
  \item{estimator}{Return AIPW, plug-in, or both estimators.}
  \item{variance}{Whether to return paired influence-score standard errors and intervals.}
  \item{conf_level}{Confidence level for pointwise intervals.}
  \item{folds}{Number of cross-fitting folds.}
  \item{nuisance_method}{Nuisance-estimation method.}
  \item{sl_lib_prop, sl_lib_miss, sl_lib_outcome}{SuperLearner libraries.}
  \item{seed}{Random seed.}
}
\description{Compares the full-X and reduced-Z ATEs for many retained
covariate sets. The paired reduced-minus-full discrepancy is a combined
MAR/NUC calibration benchmark.}
\value{Benchmark estimates, selected subsets, optional paired discrepancy,
reduced-target, and full-X influence scores, the full-X estimate, and
metadata. The target scores are used by \code{l2_calibration_band()} to retain
covariance with the sensitivity-bound estimators.}
