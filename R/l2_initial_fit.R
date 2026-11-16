#' Fit a reusable cross-fitted AIPW reference and structural nuisances
#'
#' Fits treatment, response, and outcome regressions once for each rotation.
#' The default separates structural training, pseudo-outcome regression, and
#' evaluation folds, matching sharp and CS inference. No sensitivity budgets
#' or dual optimization are needed for this initial fit.
#'
#' @inheritParams l2_sharp_crossfit
#' @param strategy Three-stage training (default), or shared training for CS only.
#' @return An object with `ate` (estimate, SE and Wald CI), `scores` (raw AIPW
#'   scores), `influence_function$ate` (centered scores), `nuisance` (held-out
#'   structural predictions), `fold_id`, `crossfit_roles`, and cached `splits`.
#' @details Pass this object as `nuisance_fit` to `l2_sharp_crossfit()`,
#'   `l2_cs_crossfit()`, or `l2_sensitivity_band()`, and as `reference_fit` to
#'   `l2_aipw_benchmarks()`. The same ordered data and covariates are required.
#'   The cache supplies folds, seed, learner settings and clipping; downstream
#'   settings for these are superseded. Radius-specific dual fitting and
#'   pseudo-outcome regressions still run. Reduced-covariate benchmarks fit
#'   their own structural models using the cached fold scheme and settings.
#'   Ordinary `generate_ate()` objects cannot be substituted for this cache.
#' @export
l2_initial_fit <- function(data, Y, A, C, X, folds = 5L,
    nuisance_method = c("SuperLearner", "glm"), sl_libs = "SL.glm",
    sl_lib_prop = sl_libs, sl_lib_miss = sl_libs, sl_lib_outcome = sl_libs,
    clip = 1e-4, propensity_clip = clip, response_clip = clip,
    seed = 1L, conf_level = .95, strategy = c("three_way", "shared_training")) {
  nuisance_method <- match.arg(nuisance_method); strategy <- match.arg(strategy)
  l2_validate_propensity_clip(propensity_clip)
  l2_validate_propensity_clip(response_clip, "response_clip")
  l2_validate_inputs(data,Y,A,C,X,c(1,1),c(0,0),c(0,0),c(0,0))
  n <- nrow(data); minimum <- if (strategy == "three_way") 3L else 2L
  if (length(folds) != 1L || !is.finite(folds) || folds < 1 || n < minimum)
    stop("A valid fold count and at least ",minimum," observations are required")
  if (length(conf_level) != 1L || !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("conf_level must lie strictly between zero and one")
  folds <- max(minimum,min(as.integer(folds),n))
  set.seed(seed); id <- sample(rep(seq_len(folds),length.out=n))
  splits <- l2_reference_splits(data,Y,A,C,X,id,strategy,nuisance_method,
    sl_lib_prop,sl_lib_miss,sl_lib_outcome,propensity_clip,seed,response_clip)
  ref <- c(l2_reference_scores(data,Y,A,C,id,splits),
    list(signature=l2_reference_signature(data,Y,A,C,X),fold_id=id,
      strategy=strategy,method=nuisance_method,seed=seed,clip=propensity_clip,
      response_clip=response_clip,lib_e=sl_lib_prop,lib_r=sl_lib_miss,lib_y=sl_lib_outcome))
  if (any(!is.finite(ref$scores$aipw))) stop("Initial fitting returned nonfinite AIPW scores")
  ref$centered <- lapply(ref$scores,function(z) z-mean(z))
  se <- stats::sd(ref$scores$aipw)/sqrt(n); est <- ref$estimates[["aipw"]]
  crit <- stats::qnorm(1-(1-conf_level)/2)
  nuisance <- lapply(splits[[1]]$evaluation,function(z) numeric(n))
  for(v in seq_len(folds)) for(nm in names(nuisance))
    nuisance[[nm]][id==v] <- splits[[v]]$evaluation[[nm]]
  roles <- data.frame(evaluation_fold=seq_len(folds),
    regression_fold=if(strategy=="three_way") c(folds,seq_len(folds-1L)) else NA_integer_,
    fit_folds=vapply(seq_len(folds),function(v) {
      excluded <- if(strategy=="three_way") c(v,if(v==1L) folds else v-1L) else v
      paste(setdiff(seq_len(folds),excluded),collapse=",")
    },character(1)))
  structure(list(ate=data.frame(estimator="aipw",estimate=est,se=se,
      conf_low=est-crit*se,conf_high=est+crit*se),
    scores=ref$scores, influence_function=list(ate=ref$centered$aipw),
    nuisance=nuisance,fold_id=id,crossfit_roles=roles,splits=splits,
    benchmark_reference=ref,conf_level=conf_level,
    parameters=list(folds=folds,strategy=strategy,nuisance_method=nuisance_method,
      propensity_clip=propensity_clip,response_clip=response_clip,seed=seed,
      sl_lib_prop=sl_lib_prop,sl_lib_miss=sl_lib_miss,sl_lib_outcome=sl_lib_outcome),
    call=match.call()),
    class="jointbounds_l2_initial_fit")
}

l2_cached_reference <- function(fit,data,Y,A,C,X) {
  if (!inherits(fit,"jointbounds_l2_initial_fit"))
    stop("nuisance_fit must be produced by l2_initial_fit(); ordinary AIPW fits are not reusable here")
  ref <- l2_get_reference(fit,l2_reference_signature(data,Y,A,C,X))
  if (is.null(fit$splits) || length(fit$splits)!=max(ref$fold_id))
    stop("The initial fit is missing its cached fold predictions")
  ref$splits <- fit$splits
  ref
}
