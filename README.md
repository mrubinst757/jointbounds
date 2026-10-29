# jointbounds

`jointbounds` bounds causal effects under joint sensitivity to informative
missingness and unmeasured confounding. It implements the methods in *Joint
Distributional Sensitivity Analysis for Unmeasured Confounding and Mixed
Informative Missingness* and retains the original mixed-missingness methods of
*Bounding causal effects with an unknown mixture of informative and
non-informative missingness* (https://arxiv.org/pdf/2411.16902).

> **`jointbounds` supersedes `marbounds`.** All `marbounds` functionality is
> included here, and `marbounds` will not receive further development.

### Migrating from marbounds

- Replace `library(marbounds)` with `library(jointbounds)`, and
  `marbounds::f()` with `jointbounds::f()`.
- Exported function names and arguments are unchanged. For example,
  `mar_bounds()`, `generate_ate()`, `multiplier_bootstrap_grid()` and the
  `l2_*()` functions all work as before.
- S3 classes are renamed from `marbounds_*` to `jointbounds_*`; for example,
  `marbounds_l2` becomes `jointbounds_l2`. Update any code that calls
  `inherits()` on the old class names.
- Version 0.4.0 fixes several bugs from `marbounds` 0.3.0 that can change
  numerical results. Among them, the default binary analyses previously used
  marginal-mean nuisances. See [NEWS.md](NEWS.md).

## Key features

- **Generic data interface**: Provide a `data.frame` with columns for outcome `Y`, treatment `A`, missingness indicator `C`, and covariates `X`.
- **Nuisance estimation**: SuperLearner with user-specified libraries and V-fold cross-fitting for propensity score, missingness probabilities, and outcome regressions.
- **Point-identified ATE**: Cross-fitted AIPW estimation and Wald inference under MAR and no unmeasured confounding via `generate_ate()`.
- **User-specified estimand**: Average treatment effect (ATE), composite ATE (Ψ₁), or separable direct effect (Ψ₂).
- **User-specified assumptions**: General bounds, bounded proportion of informative missingness (δ), monotonicity (positive/negative), bounded outcome risk (τ), or point identification under known sensitivity parameters.
- **Multiplier bootstrap**: Simultaneous inference over a grid of sensitivity parameters via multiplier bootstrap.
- **Continuous outcomes**:
  - **Bounds**: global L2 Cauchy--Schwarz bounds, empirical sharp sieve-dual
    bounds with cross-fitted EIF inference, and reduced-form net-missingness
    (`model = "net"`) and prevalence--severity (`model = "separated"`)
    models.
  - **L-infinity comparisons** (`l2_linf_bounds()`): sign-aware Hölder outer
    bounds and sharp-sieve LP bounds. These are plug-in estimates without
    inference: the sharp LP is inverse-probability weighted, and neither
    interval has a doubly robust estimator.
  - **Uncertainty**: multiplier simultaneous bands and bootstrap uncertainty.
  - **Sensitivity and calibration**: sensitivity surfaces, tipping frontiers,
    and covariate-omission calibration frontiers.

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

The legacy `marbounds` release is still on R-universe
(`install.packages("marbounds", repos = "https://mrubinst757.r-universe.dev")`),
but new work should use `jointbounds`. The two packages can be installed side
by side, but they export the same function names, so attach only one of them.

Requires: `SuperLearner` (and its dependencies).

## Quick example

Install and load the package, then load SuperLearner (required for nuisance estimation). Simulated data: `Y` outcome, `A` treatment, `C` missingness (1 = missing), `X` covariates.

```r
library(jointbounds)
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

# Data-facing L-infinity comparison. The default net model bounds
# |M_a - 1| <= delta_M_inf and |K_a - 1| <= delta_K_inf; the outer interval
# uses the sign-aware Holder radius. Plug-in estimates only (no DR estimator,
# no standard errors); the sharp LP is IPW and ignores the outcome regression.
linf <- l2_linf_bounds(
  dat_cont, "Y", "A", "C", c("X1", "X2"),
  delta_M_inf = 0.2, delta_K_inf = 0.3, method = "outer"
)

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
  surface_band, bench, calibration = "toward_zero", alpha = 1, B = 1000
)
plot(calibration_band)
signed_band <- l2_calibration_band(
  surface_band, bench, calibration = "signed", alpha = 2, B = 1000
)
frontiers <- l2_benchmark_frontiers(
  surface, bench, estimator = "aipw", calibration = "toward_zero",
  alpha = 1, plot = TRUE
)
signed_frontiers <- l2_benchmark_frontiers(
  surface, bench, estimator = "aipw", calibration = "signed", alpha = 2
)
l2_equal_radius(frontiers)

comparison <- l2_equal_radius_comparison(
  surface, frontiers, threshold = 0, uncertainty = TRUE,
  conf_level = 0.90, plot = TRUE
)
```

The optional AIPW variance uses the paired reduced-minus-full influence score.
The curves give joint \eqn{(\delta_M,\delta_K)} budgets, with each radius equal
across treatment arms, at which the chosen interval contains the target. They do
not identify a MAR-versus-NUC decomposition.

Both `l2_calibration_frontier()` and the benchmark/frontier-band wrappers
support the same two targets, writing `full = psi_X`, `reduced = psi_Z`:

- `calibration = "toward_zero"` (default):
  `full - sign(full) * alpha * abs(reduced - full)`.
- `calibration = "signed"`: `full + alpha * (reduced - full)`.

`alpha` must be positive. At `alpha = 1`, the signed target is exactly the
reduced-covariate estimate. Both modes transform the paired influence scores
using the same target formula. Only toward-zero inference excludes exact zero
references or discrepancies. Explicit `type = "target"` values are unchanged.
For coherent calibration, use the surface's zero-budget estimate and score
as the full-X benchmark. The default band inversion checks both endpoints.

The `pointwise_*` frontier limits are pointwise in the horizontal radius:
they use Bonferroni-adjusted Wald bands over both endpoints and every searched
vertical radius. `simultaneous_*` limits use a common multiplier critical
value across all grid points and supplied benchmarks by default. Set
`simultaneous_scope = "all_benchmarks"` for this joint coverage across omitted
sets, or `simultaneous_scope = "per_benchmark"` for separate coverage of each
benchmark's full frontier. For example:

```r
joint <- l2_calibration_band(surface_band, bench,
                            simultaneous_scope = "all_benchmarks")
separate <- l2_calibration_band(surface_band, bench,
                               simultaneous_scope = "per_benchmark")
```

The latter does not ensure that all benchmark frontiers are covered together.
Neither option changes the point estimates or pointwise frontier intervals.
An `NA` edge means no
crossing was found within the grid. Linear interpolation is a plotting
approximation; grid coverage alone does not guarantee off-grid coverage.
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

To compute bounds using conservative joint calibration, select upper
confidence limits for the frontier radii:

```r
calibration <- l2_calibration_band(surface_band, bench,
  calibration = "toward_zero", alpha = 1,
  simultaneous_scope = "all_benchmarks")
calibrated_grid <- l2_calibrated_grid(calibration,
  selection = "simultaneous_upper", radius_x = 0, benchmark_id = 1)
calibrated_bounds <- l2_sensitivity_band(
  dat_cont, "Y", "A", "C", c("X1", "X2"),
  grid = calibrated_grid, model = "net", bound_method = "cs")
```

Use the same model and bound method as in the calibration surface. Replace
`radius_x = 0` with an evaluated missingness radius, or omit it for all radii.
`selection = "pointwise_upper"` is appropriate for a prespecified horizontal
radius; `"estimate"` selects the point-calibrated grid. Simultaneous limits
permit selection across evaluated horizontal radii; selection across omitted
sets additionally requires `simultaneous_scope = "all_benchmarks"`.
The helper selects evaluated grid points rather than interpolated radii and
stops if a requested crossing is unavailable. This accounts for calibration
uncertainty only: the returned bound estimator's ordinary fixed-grid intervals
do not automatically give valid inference after choosing radii from the data.

Separate observed-data calibrations of informative missingness and unmeasured
confounding are not jointly identified. `l2_calibrate()` requires
`conditional_on = "exchangeability"` for an observed `R` benchmark or
`conditional_on = "MAR"` for an observed `K` benchmark. A reduced-form `M`
benchmark measures their combined distortion without attributing it.

The sharp result is exact for the empirical conditional-normalization sieve
specified by `basis`. Always inspect optimizer convergence, achieved divergence
budgets, and normalization residuals in `fit_l2$diagnostics` before reporting a
sharp endpoint.

## References

- Rubinstein, M., Agniel, D., Horvitz-Lennon, M., & Normand, S.-L. *Joint
  Distributional Sensitivity Analysis for Unmeasured Confounding and Mixed
  Informative Missingness.* Working paper. The `l2_*()` functions and the
  continuous-outcome methods follow this paper.
- Rubinstein, M., Agniel, D., Han, L., Horvitz-Lennon, M., & Normand, S.-L.
  (2026). *Bounding causal effects with an unknown mixture of informative and
  non-informative missingness.* Journal of Causal Inference (accepted).
  https://arxiv.org/pdf/2411.16902. `mar_bounds()` and the binary-outcome
  methods follow this paper; they were formerly distributed as `marbounds`.

Run `citation("jointbounds")` for citation details.

Note: by default, the bounded_risk option uses the method that assumes $\mu_a^\star/\mu_a \le \min(1/\mu_a, \tau_a)$, outlined in detail the Appendix. The option in the main paper may be recovered using the ``bounded_risk_unbounded_tau'' option.

### Frontier inference and direct optimization

`l2_calibration_band()` defaults to finite-grid inversion and returns the
`compatibility` table of estimated, inner, and outer sets. Use
`simultaneous_scope = "all_benchmarks"` for coverage across retained covariate
sets. `interpolate = TRUE` is an optional descriptive interpolation, not a
between-grid confidence guarantee. Both signed and toward-zero alpha-calibration
remain available through the `calibration` and `alpha` arguments.

For a fixed target rather than an estimated benchmark, use
`l2_tipping_band(surface, threshold = 0)` with a surface produced by
`l2_sensitivity_band(..., keep_fits = TRUE)`. Descriptive
`l2_tipping_frontier()` now uses both endpoints by default, including for
negative estimated effects.

For the paper's optional grid-free convex computation, use
`l2_frontier_direct(data, Y, A, C, X, delta_M = .1, target = 0,
 nuisance_method = "glm")`. This requires CVXR and returns the smallest common
confounding radius, arm tilts, solver status, and feasibility residuals. It is a
plug-in empirical sieve calculation without a confidence interval.

Cross-fitted CS and sharp results return `influence_covariance` and
`estimate_covariance`, each indexed by estimator (for example, `$eif`). CS
`moment_diagnostics` flags nonpositive required estimated moments and truncation;
ordinary smooth inference is not justified at these boundaries. See
[THEORY_AUDIT.md](THEORY_AUDIT.md) for the complete functionality map and limits.
Set `interpolate = TRUE` in `l2_calibration_band()` or `l2_tipping_band()` to interpolate crossings. This emits a warning that coverage between evaluated grid points is not guaranteed; the default remains `FALSE`.

## CS sample splitting

`l2_cs_crossfit()` now uses rotating three-stage splitting by default: fit
structural nuisances on the training folds, regress residual powers on a
separate fold, and evaluate scores on the held-out fold. This follows the
paper's recommended construction. The nuisance rates and moment assumptions
remain necessary. Use `moment_strategy = "shared_training"` to reproduce
learned-nuisance results from version 0.4.0; that approach needs additional
stability assumptions because the residuals are constructed in sample.
Controlled `nuisance_source = "oracle_error"` simulations are unchanged.
