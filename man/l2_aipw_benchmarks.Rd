\name{l2_aipw_benchmarks}
\alias{l2_aipw_benchmarks}
\title{Estimate simultaneous leave-covariate-out AIPW benchmarks}
\usage{
l2_aipw_benchmarks (data, Y, A, C, X, p_z = NULL, max_per_pz = 10L, subsets = NULL, subset_selection = c("random",
    "first"), estimator = c("both", "aipw", "plugin"), variance = FALSE, conf_level = 0.95,
    folds = 5L, nuisance_method = c("SuperLearner", "glm"), sl_lib_prop = sl_libs, sl_lib_miss = sl_libs,
    sl_lib_outcome = sl_libs, seed = 1L, reference_fit = NULL, reference_strategy = c("three_way",
        "shared_training", "legacy"), propensity_clip = clip, sl_libs = "SL.glm", response_clip = clip,
    clip = 1e-04, learner_fallback = c("error", "glm"))
}
\arguments{
\item{learner_fallback}{On SuperLearner failure, stop with the original error by default (\code{"error"}). The explicit \code{"glm"} option permits a substitute GLM with a warning. Initial-fit reuse inherits its saved policy. Older caches default to error for new regressions; cached predictions are not recomputed. The legacy benchmark strategy supports error only.}
\item{clip}{Common two-sided clipping threshold, default 1e-4. Probabilities are restricted to [clip, 1-clip]; no observations are removed.}
\item{propensity_clip, response_clip}{Overrides for treatment and response probabilities. Both must be finite scalars strictly between zero and 0.5. Reused reference fits supply both settings; older fits use response clipping at 1e-4.}

\item{sl_libs}{Common SuperLearner library for all nuisance regressions. Model-specific \code{sl_lib_prop}, \code{sl_lib_miss}, and \code{sl_lib_outcome} override it.}
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
  \item{reference_fit}{An \code{l2_initial_fit()} object may be supplied directly. Optional zero-radius sharp or CS fit, or a sensitivity band with that fit retained. Reuses the exact full-X AIPW estimate and scores. Reduced fits inherit its folds, training strategy, learner settings, fitting seed and clipping. The analysis columns and observation order must match. Refit objects saved before version 0.4.2.}
  \item{reference_strategy}{Without a supplied fit, use the same rotating three-stage training as the default bound estimators, shared training for two-stage CS, or legacy benchmark fitting.}

  \item{seed}{Random seed.}
}
\description{Compares the full-X and reduced-Z ATEs for many retained
covariate sets. The paired reduced-minus-full discrepancy is a combined
MAR/NUC calibration benchmark.}
\value{Benchmark estimates, selected subsets, optional paired discrepancy,
reduced-target, and full-X influence scores, the full-X estimate, and
metadata. Raw signed discrepancies and reduced-target scores are retained.
\code{l2_calibration_band()} combines the full and reduced scores for its
toward-zero target by default, retaining covariance with bound estimators.}

\section{Predictor diagnostics}{
The \code{predictor_diagnostics} table records training-set design checks,
including context (fold, model and arm where applicable), training size,
design rank, constant predictors, duplicate columns, and a representative
set of linearly dependent design columns. \code{status} is \code{ok},
\code{redundant}, or \code{unchecked}. Factor predictors are contrast-coded;
single-observed-level factors are separately flagged. Numerical rank is
checked after column scaling with QR tolerance 1e-7.

Checks run before each nuisance fit and on the actual screened training subsets
used by SL.glm inside SuperLearner. Sharp and CS fits include pseudo-outcome
regressions; sensitivity-band tables include the grid index. Reused caches
carry previously computed structural diagnostics and must be rebuilt if they
predate this feature. Checks issue warnings but do not remove predictors or
alter estimates. They do not establish absence of separation or reliability
of every prediction. One possible redundant column set is reported; the choice
can depend on column order. Older plug-in fitting paths report warnings but do
not expose this combined table.
}
