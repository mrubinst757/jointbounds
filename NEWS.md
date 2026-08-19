# marbounds 0.3.0

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
