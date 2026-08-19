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
\description{Compares the full-X and reduced-Z ATEs for many retained
covariate sets. The paired reduced-minus-full discrepancy is a combined
MAR/NUC calibration benchmark.}
\value{Benchmark estimates, selected subsets, optional paired influence
scores, the full-X estimate, and metadata.}
