# jointbounds

`jointbounds` is the development repository for the `marbounds` R package. It
implements bounds on causal effects under joint sensitivity to informative
missingness and unmeasured confounding, while retaining the original
mixed-missingness methods (see *Bounding causal effects with an unknown mixture
of informative and non-informative missingness*, https://arxiv.org/pdf/2411.16902).

## Key features

- **Generic data interface**: Provide a `data.frame` with columns for outcome `Y`, treatment `A`, missingness indicator `C`, and covariates `X`.
- **Nuisance estimation**: SuperLearner with user-specified libraries and V-fold cross-fitting for propensity score, missingness probabilities, and outcome regressions.
- **Point-identified ATE**: Cross-fitted AIPW estimation and Wald inference under MAR and no unmeasured confounding via `generate_ate()`.
- **User-specified estimand**: Average treatment effect (ATE), composite ATE (Ψ₁), or separable direct effect (Ψ₂).
- **User-specified assumptions**: General bounds, bounded proportion of informative missingness (δ), monotonicity (positive/negative), bounded outcome risk (τ), or point identification under known sensitivity parameters.
- **Multiplier bootstrap**: Simultaneous inference over a grid of sensitivity parameters via multiplier bootstrap.
- **Continuous outcomes**: Global L2 Cauchy--Schwarz bounds, empirical sharp
  sieve-dual bounds, reduced-form net-missingness models, calibration, bootstrap
  uncertainty, sensitivity surfaces, and tipping frontiers.

The continuous module defaults to cross-fitted SuperLearner nuisance models.
Its `sharp` method is an empirical conditional-normalization sieve-dual
optimization; it is not an alias for the CS method.
See [THEORY_AUDIT.md](THEORY_AUDIT.md) for a formula-by-formula map from the
companion manuscript to the exported functions and the remaining nonregular or
theoretical-only cases.

## Installation

Install the current development version from GitHub:

```r
remotes::install_github("mrubinst757/jointbounds")
```

The earlier release remains available from R-universe:

```r
# Install from r-universe
install.packages('marbounds', repos = 'https://mrubinst757.r-universe.dev')
```

Requires: `SuperLearner` (and its dependencies).

## Quick example

Install and load the package, then load SuperLearner (required for nuisance estimation). Simulated data: `Y` outcome, `A` treatment, `C` missingness (1 = missing), `X` covariates.

```r
library(marbounds)
suppressPackageStartupMessages(library(SuperLearner))

n <- 500
set.seed(1)
X <- runif(n, -2, 2)
A <- rbinom(n, 1, plogis(X))
C <- rbinom(n, 1, 0.2 + 0.1 * A)
Y <- ifelse(C == 1, NA, rbinom(n, 1, plogis(0.5 * A + 0.3 * X)))
dat <- data.frame(Y = Y, A = A, C = C, X = X)

# General bounds on ATE (no extra assumptions)
fit <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                  estimand = "ate", assumption = "general",
                  sl_lib_prop = "SL.glm", sl_lib_miss = "SL.glm", sl_lib_outcome = "SL.glm")

fit$result

fit2 <- mar_bounds(dat, Y = "Y", A = "A", C = "C", X = "X",
                   assumption = "bounded_risk",
                   param_grid = list(
                     delta_0u = seq(0.5, 1, 0.1),
                     delta_1u = seq(0.5, 1, 0.1),
                     tau_1 = seq(1.5, 5, 0.5),
                     tau_0 = seq(1.5, 5, 0.5)
                   ),
                   sl_lib = "SL.glm")
# Check the dataframe output
fit2$result
```

When MAR and no unmeasured confounding are maintained, obtain the
point-identified ATE with a cross-fitted doubly robust estimator:

```r
ate_fit <- generate_ate(
  dat, Y = "Y", A = "A", C = "C", X = "X",
  folds = 5, nuisance_method = "SuperLearner"
)
ate_fit$ate
ate_fit$diagnostic_summary
```

## Continuous-outcome L2 example

`simulate_l2_data()` uses Gaussian potential-outcome errors by default, with
configurable finite-eighth-moment t tails, overlap, nonlinear nuisance
functions, missingness prevalence/severity, and concentrated tilts.

```r
dat_cont <- simulate_l2_data(1000, seed = 2026)

# The default DGP has residual confounding. Set confounding_strength = 0 for
# treatment exchangeability given X.
l2_oracle_sensitivity_parameters(dat_cont)

# Select the violated assumptions explicitly:
dat_mar  <- simulate_l2_data(1000, violation = "mar")
dat_nuc  <- simulate_l2_data(1000, violation = "nuc")
dat_both <- simulate_l2_data(1000, violation = "both")

fit_l2 <- l2_bounds(
  dat_cont,
  Y = "Y", A = "A", C = "C", X = c("X1", "X2"),
  delta_M = 0.20,
  delta_K = 0.15,
  method = "both",
  nuisance_method = "SuperLearner",
  return_true_bounds = TRUE
)

fit_l2$ate
fit_l2$diagnostics
fit_l2$true_bounds$cs$ate    # simulation-only oracle CS comparison
fit_l2$true_bounds$sharp$ate # simulation-only oracle sharp-sieve comparison

# Cross-fitted sharp-bound inference. Cyclic three-stage splitting keeps
# structural nuisance/dual fitting, pseudo-outcome regression, and final
# influence-score evaluation on disjoint observations.
fit_cf <- l2_sharp_crossfit(
  dat_cont, "Y", "A", "C", c("X1", "X2"),
  delta_M = 0.20, delta_K = 0.15,
  folds = 5, nuisance_method = "SuperLearner"
)
fit_cf$ate

# The joint framework nests both single-mechanism analyses.
fit_nuc_only <- l2_bounds(
  dat_cont, "Y", "A", "C", c("X1", "X2"),
  delta_M = 0, delta_K = .15, method = "both"
)
fit_mar_only <- l2_bounds(
  dat_cont, "Y", "A", "C", c("X1", "X2"),
  delta_M = .20, delta_K = 0, method = "both"
)

# Monte Carlo comparison of endpoint bias, RMSE, SE calibration, and coverage.
# Supply truth = c(lower = ..., upper = ...) when analytic/reference endpoints
# are available; otherwise a large simulated empirical-sieve reference is used.
sim_cf <- l2_compare_dr_plugin(
  B = 200, n = 1000, truth_n = 30000,
  delta_M = 0.20, delta_K = 0.15,
  nuisance_simulation = "oracle_error",
  nuisance_error_rate = 0.25,
  nuisance_error_scale = c(e = 1, rho = 1, mu = 1, regression = 1),
  nuisance_error_sign = c(e = 1, rho = -1, mu = 1, regression = -1)
)
sim_cf$summary
sim_cf$bound_coverage
sim_cf$oracle_parameters
sim_cf$evaluated_parameters
sim_cf$linf_bounds
sim_cf$linf_sharp_bounds
sim_cf$bound_comparison

# Prespecified paper study (use profile = "fast" before a full run).
design <- l2_simulation_design("fast")
# study <- l2_simulation_study(design, B = 10, truth_n = 2000, folds = 3)
# plot(study, "coverage")

# Separate Monte Carlo assessment of leave-covariate-out AIPW calibration.
# cal_mc <- l2_calibration_simulation(B = 100, n = 1000, truth_n = 30000)

# Compare sharp L2 and CS outer endpoints. The CS rows contain plug-in and EIF
# estimators and are identified by bound_method = "cs".
# sim_cf <- l2_compare_dr_plugin(..., l2_method = "both")

# The summary contains plugin, dr (error-added EIF), and true_eif. In the
# controlled-error mode no regressions are fitted and cross-fitting is not
# used: e, rho, mu, and the second-stage regressions are directly perturbed.

# Set nuisance_error_rate = 0 for fixed misspecification, or use
# nuisance_simulation = "estimated" to refit the requested learners.

# If sensitivity radii are omitted, the simulation defaults to the
# minimal oracle DGP parameters. Use parameter_evaluation = "oracle" to force
# oracle calibration even when sensitivity arguments were also supplied.

# density_ratio_link = "bounded" (the default) gives finite L-infinity
# density-ratio radii. With "linear", active mechanisms have infinite
# L-infinity radii and the corresponding outer interval is unbounded.

# For a sharp-for-sieve comparison (requires lpSolve):
# sim_cf <- l2_compare_dr_plugin(..., linf_method = "both", linf_max_n = 1000)

grid <- expand.grid(
  delta_M = seq(0, 0.4, length.out = 5),
  delta_K = seq(0, 0.3, length.out = 5)
)
surface <- l2_sensitivity_grid(
  dat_cont, "Y", "A", "C", c("X1", "X2"),
  grid = grid, model = "net", method = "cs"
)
frontier <- l2_tipping_frontier(surface)
if (any(frontier$tipped_on_grid)) plot(frontier)

# Simultaneous multiplier bands over both endpoints and the full grid.
# Use keep_fits = TRUE only if the individual pointwise fits are also needed.
surface_band <- l2_sensitivity_band(
  dat_cont, "Y", "A", "C", c("X1", "X2"), grid,
  model = "net", bound_method = "cs", estimator = "eif", B = 1000,
  keep_fits = TRUE
)
surface_band$outward

# Invert the bounds: least delta_K jointly sufficient with each delta_M
# to reproduce an observed benchmark discrepancy.
l2_calibration_frontier(surface, benchmark = 0.1,
                        type = "discrepancy")

# Jointly calibrate total hidden distortion from an observed benchmark tilt.
l2_calibrate(M = c(0.9, 1, 1.1), benchmark = "observed_data")
```

Many leave-covariate-out benchmarks can also be estimated together. The cap is
applied separately within each retained-set size, so a ten-covariate analysis
does not require enumerating every possible subset.

```r
bench <- l2_aipw_benchmarks(
  dat_cont, "Y", "A", "C", c("X1", "X2"),
  p_z = 0:1, max_per_pz = 5,
  estimator = "both", variance = TRUE
)
benchmark_band <- l2_benchmark_band(bench, B = 1000)
calibration_band <- l2_calibration_band(
  surface_band, bench, B = 1000
)
plot(calibration_band)
frontiers <- l2_benchmark_frontiers(
  surface, bench, estimator = "aipw", plot = TRUE
)
l2_equal_radius(frontiers)

comparison <- l2_equal_radius_comparison(
  surface, frontiers, threshold = 0, uncertainty = TRUE,
  conf_level = 0.90, plot = TRUE
)
```

The optional AIPW variance uses the paired reduced-minus-full influence score.
The curves give joint \eqn{(\delta_M,\delta_K)} budgets, with each radius equal
across treatment arms, sufficient to reproduce each signed benchmark. They do
not identify a MAR-versus-NUC decomposition.
The equal-radius summary is the smallest common value \eqn{x} for which
\eqn{(\delta_M,\delta_K)=(x,x)} can reproduce a benchmark on the evaluated
grid; it is an E-value-like robustness summary on the model's L2 scale.
The comparison reports each calibrated radius relative to the smallest common
radius at which the overall ATE bounds include the null.
Its dot plot places the number of omitted covariates on the vertical axis and
uses slight jitter when multiple benchmark sets have the same size.
With `uncertainty = TRUE`, horizontal intervals transform the paired AIPW Wald
intervals for the benchmark discrepancies. The red band is a Wald
location-shift approximation for the overall tipping value that holds the
estimated sensitivity-width surface fixed; joint bootstrap surfaces are needed
for fully simultaneous inference.

Separate observed-data calibrations of informative missingness and unmeasured
confounding are not jointly identified. `l2_calibrate()` requires
`conditional_on = "exchangeability"` for an observed `R` benchmark or
`conditional_on = "MAR"` for an observed `K` benchmark. A reduced-form `M`
benchmark measures their combined distortion without attributing it.

The sharp result is exact for the empirical conditional-normalization sieve
specified by `basis`. Always inspect optimizer convergence, achieved divergence
budgets, and normalization residuals in `fit_l2$diagnostics` before reporting a
sharp endpoint.

## Reference

Methods and notation follow the paper on bounding causal effects under mixed informative and non-informative missingness, with influence-function-based estimation and cross-fitting. Note: by default, the bounded_risk option uses the method that assumes $\mu_a^\star/\mu_a \le \min(1/\mu_a, \tau_a)$, outlined in detail the Appendix. The option in the main paper may be recovered using the ``bounded_risk_unbounded_tau'' option.
