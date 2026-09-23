# jointbounds 0.4.0

## Renamed from marbounds

- The package is now `jointbounds` and supersedes `marbounds`. Function names
  and arguments are unchanged; replace `library(marbounds)` with
  `library(jointbounds)`. S3 classes were renamed from `marbounds_*` to
  `jointbounds_*` (for example `jointbounds_l2_sensitivity_band`), so code that tests
  `inherits(x, "marbounds_...")` must be updated. Entries below 0.4.0 refer to
  the package under its former name.

## Continuous-outcome L-infinity models

- `l2_linf_bounds()` gains `model = c("net", "separated")`. The new default
  `"net"` implements the manuscript's net L-infinity model,
  `|M_a - 1| <= delta_M_inf` and `|K_a - 1| <= delta_K_inf`; `"separated"`
  keeps the previous prevalence-capped `delta * pi * delta_R_inf` box.
- The L-infinity outer radius is now sign-aware,
  `E|Z| {C^+(X) + C^-(X)} / 2`, where the downward corners respect
  `M_a >= rho_a` and `K_a >= 0`. It is never wider than the previous symmetric
  radius.
- `l2_oracle_sensitivity_parameters()` also returns the minimal net radii
  `delta_M_inf` and `delta_K_inf`; the oracle L-infinity outer and sharp
  bounds and `l2_compare_dr_plugin()` follow the requested `model`.

## Bug fixes

- `estimate_nuisance(V = 1)` previously had an empty training set and silently
  fell back to marginal means. Because `mar_bounds()` uses `V = 1` whenever all
  libraries are `"SL.glm"`, the default binary analyses were affected. With
  `V = 1` nuisances are now fitted and predicted on the full sample.
- `estimate_nuisance(stratify_mu = FALSE)` now includes treatment in the pooled
  outcome model; previously it forced `mu0 == mu1`.
- The smooth-approximation (`smooth_approximation = TRUE`) bounded-risk ATE
  scores used an uncentered and sign-incorrect chain-rule term, which made the
  bounds diverge as `epsilon` shrank. They now converge to the indicator bounds.
- `mar_bounds()` grids now honour the `delta` and `tau` shorthands for every
  ATE and Psi_1 assumption. The bounded-risk grid defaults `tau_1` to `tau_0`,
  matching the scalar path.
- `multiplier_bootstrap_grid(assumption = "bounded_risk")` now uses the same
  masked `min{1 - mu, (tau - 1) mu}` score as `mar_bounds()`, and it no longer
  resets the global RNG seed.

# marbounds 0.3.0

- Controlled nuisance-error simulations no longer claim or perform
  cross-fitting: they directly perturb DGP nuisance functions on the full
  simulated sample. Outcome regressions and CS conditional residual moments
  now receive explicit errors. Sharp simulations recompute the dual solution
  and dual loss H from the contaminated nuisances, while `true_eif` remains
  the exact-nuisance diagnostic.

- Sharp endpoint estimation now uses the manuscript's rotating three-stage
  cross-fitting: structural nuisances/dual fitting, pseudo-outcome regression,
  and final score evaluation are performed on disjoint folds.
- Added `l2_calibration_band()` for combined influence-function pointwise and
  simultaneous calibration-frontier bands obtained by inverting an evaluated
  sensitivity grid.
- Added explicit single-mechanism tests showing that the sharp sieve and
  Cauchy--Schwarz formulas agree when only one L2 mechanism is active and the
  normalization basis contains the relevant optimizer.
- Corrected the legacy bounded-`delta` Psi2 estimator to apply its sensitivity
  radius, use the proper positive-part influence scores, and preserve coherent
  finite-sample endpoints by projecting nonnegative support estimates to zero.

- Added l2_simulation_design() and l2_simulation_study() to reproduce the
  paper's finite-sample and sharpness experiments, with checkpointing, tidy
  operating-characteristic summaries, and publication-ready ggplot methods.
- The continuous DGP now exposes informative-missingness prevalence and
  strength, treatment overlap, nonlinear nuisance functions, concentrated
  tilts, and finite-eighth-moment Student t outcomes. Oracle quadrature and
  sensitivity radii use the same configurable DGP.
- Simulation output distinguishes point-bound containment from outward-Wald
  confidence-set coverage of the true ATE and reports scenario failure rates.

- Added `l2_sensitivity_band()` and `l2_multiplier_band()` implementing the
  standardized multiplier simultaneous bands in the companion paper, including
  outward bands for the entire identified interval.
- Added `l2_linf_bounds()` as a data-facing Holder and sharp-sieve comparison.
  The sharp LP now always enforces the mixture lower envelope, including when a
  symmetric likelihood-ratio radius is greater than one.
- Added cross-fitted CS and sharp inference for the reduced-form net-M model.
- Added `l2_prevalence_bounds()` for the capped L2 prevalence support-function
  bounds on binary composite and separable direct effects.
- Zero missingness budgets are now imposed exactly as M=1 instead of being
  approximated by a very large dual multiplier. Sharp diagnostics now include
  the empirical primal value, primal-dual gap, normalization error, and budget
  violation.
- Added `l2_sharp_crossfit()` for cross-fitted plug-in and augmented doubly
  robust evaluation of empirical sharp sieve endpoints.
- Implemented the conditional dual-loss regression, structural propensity and
  response derivative corrections, and the moving-boundary KKT multiplier.
- Added `l2_compare_dr_plugin()` to compare endpoint bias, standard-error
  calibration, Wald coverage, and true-ATE bound coverage in simulation.
- Added a fast `oracle_error` simulation mode with smooth nuisance errors of
  user-controlled size `n^(-alpha)`. Fixed perturbation directions expose
  first-order plug-in bias instead of averaging it away across replications.
- Simulation summaries now always include `true_eif`, based on exact DGP
  sampling scores and quadrature-based EIF regressions, alongside the plug-in
  and error-contaminated EIF. Nuisance-error signs are independently settable.
- Added oracle DGP calibration for the minimal prevalence, informative-
  missingness, and confounding sensitivity parameters. Simulations default to
  these parameters when the user does not specify radii.
- The continuous DGP now includes configurable residual confounding through a
  treatment mechanism depending on a shared potential-outcome disturbance;
  `confounding_strength = 0` recovers the earlier exchangeable design.
- Added `violation = "mar"`, `"nuc"`, `"both"`, or `"none"` to switch the two
  identifying-assumption violations on independently in data generation and
  Monte Carlo comparisons.
- Added bounded-link DGPs, oracle L-infinity density-ratio radii, and an
  L-infinity/L1 Holder outer-bound comparison with the L2 sharp-sieve interval.
  Linear-link DGPs correctly report infinite L-infinity radii.
- Added an optional sharp L-infinity empirical-sieve linear program enforcing
  pointwise density-ratio envelopes and basis normalization of both M and L.
- Added cross-fitted plug-in and analytic-EIF estimators of the closed-form L2
  Cauchy--Schwarz outer endpoints, available through `l2_method = "cs"` or
  `"both"` in the Monte Carlo comparison.
- Added simultaneous leave-covariate-out combined MAR/NUC benchmarks with
  cross-fitted plug-in and AIPW estimates, optional paired variance estimates,
  capped subset generation by retained-set size, and multi-benchmark frontier
  inversion and plotting.
- Added a styled ggplot display for joint benchmark frontiers and an
  equal-radius summary giving the smallest common L2 budget sufficient to
  reproduce each benchmark.
- Added an overall equal-radius tipping value for the ATE and a ggplot
  dot-plot comparison of calibrated benchmark strength against that tipping
  value. The vertical axis records the number of omitted covariates, with
  slight jitter to reveal overlapping benchmarks.
- Added optional transformed Wald intervals for calibrated equal-radius values
  and a location-shift uncertainty band for the overall tipping radius.
- Corrected the nuisance-error experiment to hold the oracle sieve-dual
  solution fixed across plug-in, contaminated-EIF, and true-EIF evaluations.
  This matches the envelope-theorem comparison and avoids changing the target
  by re-optimizing under deliberately contaminated probabilities.

# marbounds 0.2.0

- Added true conditional DGP quantities and potential missingness variables to
  `simulate_l2_data()`.
- Added cross-fitted SuperLearner nuisance estimation and simulation-only oracle
  Monte Carlo CS comparisons to `l2_bounds()`.
- Added plot-ready tipping frontiers with an S3 plot method.
- Clarified that `method = "sharp"` is the empirical sieve-dual optimizer.

* Added continuous-outcome global L2 sensitivity analysis via `l2_bounds()`.
* Added closed-form Cauchy--Schwarz outer bounds and empirical sharp sieve-dual
  bounds for jointly nonignorable missingness and unmeasured confounding.
* Added the reduced-form net-missingness model based on an L2 restriction on
  `M - 1`.
* Added sensitivity surfaces, tipping frontiers, benchmark calibration, and a
  full nonparametric bootstrap.
* Added `simulate_l2_data()` and a continuous-outcome vignette using Gaussian
  errors that satisfy the paper's polynomial moment restrictions.
* Preserved the original binary-outcome API.

# marbounds 0.1.0

Initial release of the marbounds package for bounding causal effects under mixed informative and non-informative missingness.

## Features

* `mar_bounds()` - Main function for computing bounds on causal effects
* Support for multiple estimands: ATE, composite ATE (Ψ₁), separable direct effect (Ψ₂)
* Multiple assumption types:
  - General bounds (no assumptions)
  - Bounded proportion of informative missingness (δ)
  - Monotonicity (positive/negative)
  - Bounded outcome risk (τ)
  - Point identification under known sensitivity parameters
* SuperLearner integration for nuisance parameter estimation with V-fold cross-fitting
* Multiplier bootstrap for simultaneous inference over parameter grids
* Influence function-based estimators with asymptotic standard errors

## Reference

Rubinstein, M., Agniel, D., Han, L., Horvitz-Lennon, M., & Normand, S.-L. (2026).
Bounding causal effects with an unknown mixture of informative and non-informative missingness.
Journal of Causal Inference (Accepted). https://arxiv.org/pdf/2411.16902
