\name{l2_initial_fit}
\alias{l2_initial_fit}
\title{Reusable cross-fitted AIPW reference and structural nuisance fits}
\usage{
l2_initial_fit(data, Y, A, C, X, folds = 5L,
  nuisance_method = c("SuperLearner", "glm"), sl_libs = "SL.glm",
  sl_lib_prop = sl_libs, sl_lib_miss = sl_libs, sl_lib_outcome = sl_libs,
  clip = 1e-4, propensity_clip = clip, response_clip = clip,
  seed = 1L, conf_level = .95, strategy = c("three_way", "shared_training"), learner_fallback = c("error", "glm"))
}
\arguments{
\item{learner_fallback}{On SuperLearner failure, stop with the original error by default (\code{"error"}). The explicit \code{"glm"} option permits a substitute GLM with a warning. Initial-fit reuse inherits its saved policy. Older caches default to error for new regressions; cached predictions are not recomputed. The legacy benchmark strategy supports error only.}
\item{data}{Analysis data frame, preserving observation order for downstream reuse.}
\item{Y,A,C}{Column names: outcome, binary treatment, and missingness (C=1 means missing).}
\item{X}{Character vector of adjustment covariate names.}
\item{folds}{Number of folds; at least three for three-way training and two for shared training.}
\item{nuisance_method}{SuperLearner or GLM structural regressions.}
\item{sl_libs}{Common SuperLearner library.}
\item{sl_lib_prop,sl_lib_miss,sl_lib_outcome}{Model-specific library overrides.}
\item{clip}{Common two-sided probability clipping threshold.}
\item{propensity_clip,response_clip}{Treatment and response clipping overrides; strictly between zero and 0.5.}
\item{seed}{Random seed for splitting and structural fitting.}
\item{conf_level}{Wald confidence level.}
\item{strategy}{Three disjoint training/regression/evaluation roles by default. Shared training is supported by CS only.}
}
\description{Fits structural nuisance models once per rotation and computes the AIPW estimate under MAR and no unmeasured confounding, without sensitivity optimization.}
\value{An object with \code{ate} (estimate, SE, CI), \code{scores} (raw AIPW and plug-in scores), \code{influence_function$ate} (centered AIPW scores), \code{nuisance} (held-out predictions), \code{fold_id}, \code{crossfit_roles}, \code{parameters}, and cached \code{splits}.}
\details{Pass as \code{nuisance_fit} to \code{l2_sharp_crossfit}, \code{l2_cs_crossfit}, or \code{l2_sensitivity_band}, and as \code{reference_fit} to \code{l2_aipw_benchmarks}. Cached settings supersede downstream fold, seed, learner, and clipping settings. The same data, columns and row order are required. Radius-specific optimization and pseudo-outcome regression still run. Reduced-covariate benchmarks fit their own models using the same fold scheme. Ordinary \code{generate_ate()} fits cannot substitute for this object.}

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
