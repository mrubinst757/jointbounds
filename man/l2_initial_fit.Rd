\name{l2_initial_fit}
\alias{l2_initial_fit}
\title{Reusable cross-fitted AIPW reference and structural nuisance fits}
\usage{
l2_initial_fit(data, Y, A, C, X, folds = 5L,
  nuisance_method = c("SuperLearner", "glm"), sl_libs = "SL.glm",
  sl_lib_prop = sl_libs, sl_lib_miss = sl_libs, sl_lib_outcome = sl_libs,
  clip = 1e-4, propensity_clip = clip, response_clip = clip,
  seed = 1L, conf_level = .95, strategy = c("three_way", "shared_training"))
}
\arguments{
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
