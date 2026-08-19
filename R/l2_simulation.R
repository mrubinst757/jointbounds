#' Reproducible designs for the continuous-outcome simulation study
#'
#' @param profile `"main"` returns the prespecified paper design and `"fast"`
#'   returns a three-scenario smoke-test design.
#' @return A data frame with one row per data-generating scenario.
#' @export
l2_simulation_design <- function(profile = c("main", "fast")) {
  profile <- match.arg(profile)
  base <- data.frame(
    scenario_id = c("null", "mar_only", "nuc_only", "joint_primary",
      "rare_severe", "common_mild", "concentrated",
      "poor_overlap", "heavy_tail", "unbounded_ratio"),
    scenario = c("No violations", "Informative missingness only",
      "Unmeasured confounding only", "Joint violation, n = 1000",
      "Rare, severe informative missingness",
      "Common, mild informative missingness", "Covariate-concentrated missingness tilt",
      "Poor treatment overlap", "Heavy-tailed outcome", "Unbounded density ratios"),
    family = c("identification", "identification", "identification", "joint",
      "missingness", "missingness", "concentration", "overlap",
      "tails", "linfinity"),
    n = rep(1000L, 10L),
    violation = c("none", "mar", "nuc", rep("both", 7L)),
    density_ratio_link = c(rep("bounded", 9L), "linear"),
    confounding_strength = c(0, 0, .6, rep(.6, 7L)),
    informative_strength = c(0, .45, 0, .45, .90, .20, .90,
      .45, .45, .45),
    informative_intercept = c(rep(-2.4, 4L), -3.4, -1.3, rep(-2.4, 4L)),
    noninformative_intercept = rep(-1.8, 10L),
    treatment_scale = c(rep(1, 7L), 1.8, rep(1, 2L)),
    covariate_complexity = rep("linear", 10L),
    informative_concentration = c(rep("diffuse", 6L), "high_variance",
      rep("diffuse", 3L)),
    outcome_distribution = c(rep("normal", 8L), "t", "normal"),
    outcome_df = rep(12, 10L),
    nuisance_simulation = rep("oracle_error", 10L),
    nuisance_error_rate = rep(.25, 10L),
    model = rep("net", 10L),
    stringsAsFactors = FALSE
  )
  size_rows <- base[base$scenario_id == "joint_primary", , drop = FALSE]
  size_rows <- size_rows[rep(1L, 2L), , drop = FALSE]
  size_rows$scenario_id <- c("joint_n500", "joint_n2000")
  size_rows$scenario <- c("Joint, n = 500", "Joint, n = 2000")
  size_rows$family <- "sample_size"
  size_rows$n <- c(500L, 2000L)
  slow <- base[base$scenario_id == "joint_primary", , drop = FALSE]
  slow$scenario_id <- "slow_nuisance"
  slow$scenario <- "Joint, n^(-0.15) nuisance error"
  slow$family <- "nuisance"
  slow$nuisance_error_rate <- .15
  base <- rbind(base, size_rows, slow)
  if (profile == "fast")
    base <- base[base$scenario_id %in% c("null", "joint_primary",
      "heavy_tail"), , drop = FALSE]
  rownames(base) <- NULL
  base
}

l2_simulation_validate_design <- function(design) {
  required <- c("scenario_id", "scenario", "n", "violation",
    "density_ratio_link", "confounding_strength", "informative_strength",
    "informative_intercept", "noninformative_intercept", "treatment_scale",
    "covariate_complexity", "informative_concentration",
    "outcome_distribution", "outcome_df", "nuisance_simulation", "model",
    "nuisance_error_rate")
  missing <- setdiff(required, names(design))
  if (length(missing)) stop("design is missing: ", paste(missing, collapse = ", "))
  if (!nrow(design) || anyDuplicated(design$scenario_id))
    stop("design must have at least one row and unique scenario_id values")
  if (any(!is.finite(design$n)) || any(design$n < 50))
    stop("Every design n must be at least 50")
  design
}

l2_simulation_bind <- function(results, component, scenario_columns) {
  out <- lapply(seq_along(results), function(i) {
    z <- results[[i]][[component]]
    if (is.null(z) || !NROW(z)) return(NULL)
    z <- as.data.frame(z)
    meta <- scenario_columns[i, , drop = FALSE]
    ans <- cbind(meta[rep(1L, nrow(z)), , drop = FALSE], z)
    rownames(ans) <- NULL
    ans
  })
  out <- Filter(Negate(is.null), out)
  if (length(out)) do.call(rbind, out) else data.frame()
}

l2_simulation_mechanism_comparison <- function(data, parameters, basis = NULL,
                                               control = list(),
                                               model = c("net", "separated")) {
  model <- match.arg(model)
  zero <- c(0, 0)
  ref <- l2_oracle_sharp_reference(data, delta = parameters$delta,
    delta_R = zero, delta_A = zero, delta_M = zero, model = model,
    basis = basis, control = control)
  miss <- l2_oracle_sharp_reference(data, delta = parameters$delta,
    delta_R = parameters$delta_R, delta_A = zero, basis = basis,
    delta_M = parameters$delta_M, model = model, control = control)
  nuc <- l2_oracle_sharp_reference(data, delta = parameters$delta,
    delta_R = zero, delta_A = parameters$delta_A, basis = basis,
    delta_M = zero, model = model, control = control)
  joint <- l2_oracle_sharp_reference(data, delta = parameters$delta,
    delta_R = parameters$delta_R, delta_A = parameters$delta_A,
    delta_M = parameters$delta_M, model = model,
    basis = basis, control = control)
  additive <- c(lower = miss[["lower"]] + nuc[["lower"]] - ref[["lower"]],
    upper = miss[["upper"]] + nuc[["upper"]] - ref[["upper"]])
  complete_case_value <- with(data,
    mean(Y[A == 1 & C == 0], na.rm = TRUE) -
      mean(Y[A == 0 & C == 0], na.rm = TRUE))
  complete_case <- c(lower = complete_case_value, upper = complete_case_value)
  vals <- list(complete_case = complete_case, reference = ref,
    missingness_only = miss, confounding_only = nuc,
    additive_separate = additive, joint_sharp = joint)
  do.call(rbind, lapply(names(vals), function(nm) data.frame(
    analysis = nm, lower = vals[[nm]][["lower"]],
    upper = vals[[nm]][["upper"]],
    width = vals[[nm]][["upper"]] - vals[[nm]][["lower"]],
    stringsAsFactors = FALSE)))
}

#' Compare sharp-sieve bounds across normalization bases
#'
#' @param data Data returned by `simulate_l2_data()`.
#' @param parameters Optional output from `l2_oracle_sensitivity_parameters()`.
#' @param bases Named list of formulas. By default, it compares intercept-only,
#'   linear, and cubic/interacted normalization restrictions.
#' @param control Optimizer control list.
#' @return A data frame of oracle sharp endpoints and excess width relative to
#'   the richest supplied basis.
#' @export
l2_basis_approximation <- function(data, parameters = NULL,
                                   bases = NULL, control = list(),
                                   model = c("net", "separated")) {
  model <- match.arg(model)
  if (is.null(parameters)) parameters <- l2_oracle_sensitivity_parameters(data)
  if (is.null(bases)) bases <- list(intercept = ~ 1,
    linear = ~ X1 + X2,
    rich = ~ X1 + X2 + I(X1^2) + I(X2^2) + I(X1^3) + I(X2^3) + I(X1 * X2))
  if (!is.list(bases) || !length(bases) || is.null(names(bases)) ||
      any(!nzchar(names(bases)))) stop("bases must be a nonempty named list")
  ans <- do.call(rbind, lapply(seq_along(bases), function(j) {
    ep <- l2_oracle_sharp_reference(data, delta = parameters$delta,
      delta_R = parameters$delta_R, delta_A = parameters$delta_A,
      delta_M = parameters$delta_M, model = model,
      basis = bases[[j]], control = control)
    data.frame(basis = names(bases)[j], order = j,
      lower = ep[["lower"]], upper = ep[["upper"]],
      width = ep[["upper"]] - ep[["lower"]], stringsAsFactors = FALSE)
  }))
  ans$excess_width_vs_richest <- ans$width - ans$width[nrow(ans)]
  ans
}

#' Run the continuous-outcome simulation study
#'
#' This high-level runner evaluates endpoint bias, empirical standard deviation,
#' standard-error calibration, pointwise endpoint coverage, true-ATE containment,
#' outward-Wald coverage, optimizer failure, and width inflation relative to the
#' oracle sharp sieve. Sensitivity radii are set to the minimal oracle radii of
#' each DGP unless overridden through `...`.
#'
#' @param design A data frame returned by `l2_simulation_design()`, possibly
#'   subsetted or modified.
#' @param B Monte Carlo replications per scenario.
#' @param truth_n Size of the independent Monte Carlo population used for oracle
#'   radii and endpoint references.
#' @param seed Master seed. Scenario `j` uses `seed + 10000*j`.
#' @param conf_level Confidence level for endpoint Wald intervals.
#' @param l2_method Passed to `l2_compare_dr_plugin()`.
#' @param linf_method Passed to `l2_compare_dr_plugin()`. The main finite-sample
#'   study defaults to the inexpensive outer comparison.
#' @param folds Number of cross-fitting folds.
#' @param progress Print scenario progress.
#' @param continue_on_error If true, record a failed scenario and continue.
#' @param checkpoint_file Optional RDS path updated after every scenario.
#' @param ... Additional named arguments passed to `l2_compare_dr_plugin()`.
#' @return An object of class `marbounds_l2_simulation` containing the design,
#'   scenario-specific fits, and tidy summaries.
#' @export
l2_simulation_study <- function(design = l2_simulation_design("main"),
                                B = 1000L, truth_n = 100000L, seed = 202608L,
                                conf_level = .95,
                                l2_method = c("both", "sharp", "cs"),
                                linf_method = c("outer", "both", "sharp"),
                                folds = 2L, progress = interactive(),
                                continue_on_error = TRUE,
                                checkpoint_file = NULL, ...) {
  design <- l2_simulation_validate_design(as.data.frame(design))
  l2_method <- match.arg(l2_method); linf_method <- match.arg(linf_method)
  B <- as.integer(B); truth_n <- as.integer(truth_n); folds <- as.integer(folds)
  if (B < 2L || truth_n < 1000L || folds < 2L)
    stop("B must be at least 2, truth_n at least 1000, and folds at least 2")
  extra <- list(...)
  if (length(extra) && (is.null(names(extra)) || any(!nzchar(names(extra)))))
    stop("All arguments supplied through ... must be named")
  fits <- vector("list", nrow(design)); mechanisms <- vector("list", nrow(design))
  failures <- vector("list", nrow(design))
  for (j in seq_len(nrow(design))) {
    d <- design[j, , drop = FALSE]
    if (isTRUE(progress)) message("[", j, "/", nrow(design), "] ", d$scenario)
    dgp_args <- list(informative_strength = d$informative_strength,
      informative_intercept = d$informative_intercept,
      noninformative_intercept = d$noninformative_intercept,
      treatment_scale = d$treatment_scale,
      covariate_complexity = d$covariate_complexity,
      informative_concentration = d$informative_concentration,
      outcome_distribution = d$outcome_distribution,
      outcome_df = d$outcome_df)
    call_args <- list(B = B, n = d$n, violation = d$violation,
      density_ratio_link = d$density_ratio_link,
      confounding_strength = d$confounding_strength,
      truth_n = truth_n, seed = seed + 10000L * j, conf_level = conf_level,
      nuisance_simulation = d$nuisance_simulation,
      parameter_evaluation = "oracle", linf_method = linf_method,
      l2_method = l2_method, dgp_args = dgp_args, folds = folds,
      model = d$model)
    if (length(extra)) {
      call_args[names(extra)] <- NULL
      call_args <- c(call_args, extra)
    }
    if (is.finite(d$nuisance_error_rate))
      call_args$nuisance_error_rate <- d$nuisance_error_rate
    fit <- try(do.call(l2_compare_dr_plugin, call_args), silent = TRUE)
    if (inherits(fit, "try-error")) {
      failures[[j]] <- data.frame(scenario_id = d$scenario_id,
        error = as.character(fit), stringsAsFactors = FALSE)
      if (!continue_on_error) stop(as.character(fit))
    } else {
      fits[[j]] <- fit
      calibration <- do.call(simulate_l2_data, c(list(n = truth_n,
        seed = seed + 10000L * j + 100000L, violation = d$violation,
        density_ratio_link = d$density_ratio_link,
        confounding_strength = d$confounding_strength), dgp_args))
      mech <- try(l2_simulation_mechanism_comparison(calibration,
        fit$oracle_parameters, basis = extra$basis,
        control = if (is.null(extra$control)) list() else extra$control,
        model = d$model),
        silent = TRUE)
      if (inherits(mech, "try-error")) {
        failures[[j]] <- data.frame(scenario_id = d$scenario_id,
          error = paste("Mechanism comparison:", as.character(mech)),
          stringsAsFactors = FALSE)
      } else mechanisms[[j]] <- mech
    }
    if (!is.null(checkpoint_file)) saveRDS(list(design = design,
      completed = j, fits = fits, mechanisms = mechanisms,
      failures = failures), checkpoint_file)
  }
  ok <- !vapply(fits, is.null, logical(1))
  scenario_columns <- design[, intersect(c("scenario_id", "scenario", "family",
    "n", "violation", "nuisance_simulation"), names(design)), drop = FALSE]
  reps <- l2_simulation_bind(fits, "replicates", scenario_columns)
  summary <- l2_simulation_bind(fits, "summary", scenario_columns)
  coverage <- l2_simulation_bind(fits, "bound_coverage", scenario_columns)
  widths <- l2_simulation_bind(fits, "bound_comparison", scenario_columns)
  mechanism <- lapply(seq_along(mechanisms), function(i) {
    if (is.null(mechanisms[[i]])) return(NULL)
    ans <- cbind(scenario_columns[i, , drop = FALSE][rep(1L,
      nrow(mechanisms[[i]])), , drop = FALSE], mechanisms[[i]])
    rownames(ans) <- NULL
    ans
  })
  mechanism <- Filter(Negate(is.null), mechanism)
  mechanism <- if (length(mechanism)) do.call(rbind, mechanism) else data.frame()
  diagnostics <- do.call(rbind, lapply(which(ok), function(i) data.frame(
    scenario_id = design$scenario_id[i], scenario = design$scenario[i],
    n = design$n[i], n_requested = fits[[i]]$n_requested,
    n_success = fits[[i]]$n_success, failure_rate = fits[[i]]$failure_rate)))
  if (nrow(reps)) {
    sharp_diag <- unique(reps[reps$bound_method == "sharp",
      c("scenario_id", "replication", "optimizer_converged",
        "maximum_primal_dual_gap", "maximum_normalization_error",
        "maximum_budget_violation", "maximum_active_fraction"), drop = FALSE])
    if (nrow(sharp_diag)) {
      opt <- do.call(rbind, lapply(split(sharp_diag, sharp_diag$scenario_id),
        function(z) data.frame(scenario_id = z$scenario_id[1L],
          optimizer_convergence_rate = mean(z$optimizer_converged),
          max_primal_dual_gap = max(z$maximum_primal_dual_gap),
          max_normalization_error = max(z$maximum_normalization_error),
          max_budget_violation = max(z$maximum_budget_violation),
          mean_active_fraction = mean(z$maximum_active_fraction))))
      diagnostics <- merge(diagnostics, opt, by = "scenario_id", all.x = TRUE,
        sort = FALSE)
    }
  }
  failed <- Filter(Negate(is.null), failures)
  failed <- if (length(failed)) do.call(rbind, failed) else data.frame()
  out <- list(design = design, fits = fits, replicates = reps,
    endpoint_summary = summary, set_coverage = coverage,
    width_comparison = widths, mechanism_comparison = mechanism,
    diagnostics = diagnostics, failed_scenarios = failed,
    settings = list(B = B, truth_n = truth_n, seed = seed,
      conf_level = conf_level, folds = folds, l2_method = l2_method,
      linf_method = linf_method), call = match.call())
  class(out) <- c("marbounds_l2_simulation", "list")
  out
}

#' @export
print.marbounds_l2_simulation <- function(x, ...) {
  cat("Continuous-outcome L2 simulation study\n")
  cat("Scenarios:", nrow(x$design), "  completed:", nrow(x$diagnostics),
      "  failed:", nrow(x$failed_scenarios), "\n")
  if (nrow(x$endpoint_summary)) print(x$endpoint_summary, row.names = FALSE)
  invisible(x)
}

#' Plot operating characteristics from an L2 simulation study
#'
#' @param x An object returned by `l2_simulation_study()`.
#' @param metric `"coverage"`, `"bias"`, or `"width"`.
#' @param ... Unused.
#' @export
plot.marbounds_l2_simulation <- function(x,
                                         metric = c("coverage", "bias", "width"),
                                         ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Plotting simulation results requires ggplot2")
  metric <- match.arg(metric)
  if (metric == "coverage") {
    d <- x$endpoint_summary
    p <- ggplot2::ggplot(d, ggplot2::aes(x = scenario, y = coverage,
      colour = estimator, shape = bound_method)) +
      ggplot2::geom_hline(yintercept = x$settings$conf_level,
        linetype = 2, colour = "grey45") + ggplot2::geom_point(size = 2.4) +
      ggplot2::facet_wrap(~ endpoint) +
      ggplot2::labs(x = NULL, y = "Pointwise endpoint coverage",
        colour = "Estimator", shape = "Bound")
  } else if (metric == "bias") {
    d <- x$endpoint_summary
    p <- ggplot2::ggplot(d, ggplot2::aes(x = scenario, y = bias,
      colour = estimator, shape = bound_method)) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey45") +
      ggplot2::geom_point(size = 2.4) + ggplot2::facet_wrap(~ endpoint) +
      ggplot2::labs(x = NULL, y = "Endpoint bias", colour = "Estimator",
        shape = "Bound")
  } else {
    d <- x$width_comparison
    p <- ggplot2::ggplot(d, ggplot2::aes(x = scenario, y = width,
      colour = model)) + ggplot2::geom_point(size = 2.4) +
      ggplot2::labs(x = NULL, y = "Oracle interval width", colour = "Method")
  }
  p + ggplot2::coord_flip() + ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom")
}

#' Monte Carlo evaluation of covariate-omission calibration benchmarks
#'
#' @param B,n,truth_n Replications, sample size, and reference-population size.
#' @param p_z,max_per_pz Retained-set sizes and cap.
#' @param simultaneous Compute a multiplier simultaneous band per replication.
#' @param multiplier_B Number of multiplier draws.
#' @param seed Master seed.
#' @param violation,density_ratio_link,confounding_strength DGP arguments.
#' @param dgp_args Additional named arguments for `simulate_l2_data()`.
#' @param folds,nuisance_method Nuisance-estimation settings.
#' @param conf_level Confidence level.
#' @return Reference values, replication results, operating characteristics,
#'   and simultaneous coverage.
#' @export
l2_calibration_simulation <- function(B = 500L, n = 1000L,
    truth_n = 100000L, p_z = c(0L, 1L), max_per_pz = 2L,
    simultaneous = TRUE, multiplier_B = 500L, seed = 202609L,
    violation = "both", density_ratio_link = "bounded",
    confounding_strength = .6, dgp_args = list(), folds = 2L,
    nuisance_method = c("SuperLearner", "glm"), conf_level = .95) {
  nuisance_method <- match.arg(nuisance_method)
  B <- as.integer(B); if (B < 2L) stop("B must be at least two")
  make_data <- function(nn, ss) do.call(simulate_l2_data, c(list(n = nn,
    seed = ss, violation = violation, density_ratio_link = density_ratio_link,
    confounding_strength = confounding_strength), dgp_args))
  pop <- make_data(truth_n, seed + 100000L)
  ref_fit <- l2_aipw_benchmarks(pop, "Y", "A", "C", c("X1", "X2"),
    p_z = p_z, max_per_pz = max_per_pz, estimator = "aipw",
    variance = FALSE, folds = folds, nuisance_method = nuisance_method,
    seed = seed + 100000L)
  reference <- ref_fit$results[, c("benchmark_id", "label", "p_z",
    "n_omitted", "delta_bench")]
  names(reference)[names(reference) == "delta_bench"] <- "truth"
  rows <- list(); simultaneous_rows <- list(); failures <- character()
  for (b in seq_len(B)) {
    fit <- try(l2_aipw_benchmarks(make_data(n, seed + b), "Y", "A", "C",
      c("X1", "X2"), subsets = ref_fit$subsets, estimator = "both",
      variance = TRUE, conf_level = conf_level, folds = folds,
      nuisance_method = nuisance_method, seed = seed + b), silent = TRUE)
    if (inherits(fit, "try-error")) {
      failures <- c(failures, as.character(fit)); next
    }
    z <- merge(fit$results, reference[, c("benchmark_id", "truth")],
      by = "benchmark_id", all.x = TRUE, sort = FALSE)
    z$replication <- b
    z$covered <- z$conf_low <= z$truth & z$truth <= z$conf_high
    rows[[length(rows) + 1L]] <- z
    if (isTRUE(simultaneous)) {
      band <- try(l2_benchmark_band(fit, B = multiplier_B,
        conf_level = conf_level, seed = seed + 10000L + b), silent = TRUE)
      if (!inherits(band, "try-error")) {
        q <- merge(band$results, reference[, c("benchmark_id", "truth")],
          by = "benchmark_id", all.x = TRUE, sort = FALSE)
        simultaneous_rows[[length(simultaneous_rows) + 1L]] <- data.frame(
          replication = b, covered = all(q$simultaneous_low <= q$truth &
            q$truth <= q$simultaneous_high))
      }
    }
  }
  reps <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (!nrow(reps)) stop("All calibration simulation fits failed")
  groups <- split(reps, interaction(reps$benchmark_id, reps$estimator,
    drop = TRUE))
  summary <- do.call(rbind, lapply(groups, function(z) data.frame(
    benchmark_id = z$benchmark_id[1L], label = z$label[1L],
    p_z = z$p_z[1L], n_omitted = z$n_omitted[1L],
    estimator = z$estimator[1L], n_success = nrow(z),
    bias = mean(z$delta_bench - z$truth),
    rmse = sqrt(mean((z$delta_bench - z$truth)^2)),
    empirical_sd = stats::sd(z$delta_bench), mean_se = mean(z$se),
    coverage = mean(z$covered))))
  sim <- if (length(simultaneous_rows)) do.call(rbind, simultaneous_rows) else
    data.frame()
  list(reference = reference, replicates = reps, summary = summary,
    simultaneous_coverage = if (nrow(sim)) mean(sim$covered) else NA_real_,
    n_requested = B, n_success = length(unique(reps$replication)),
    failure_rate = 1 - length(unique(reps$replication)) / B,
    fit_errors = unique(failures), call = match.call())
}
