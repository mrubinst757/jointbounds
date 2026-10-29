# Paper-to-package audit (2026-10-29)

Reviewed against bounds/main.tex (Overleaf revision 9c3b2bc). In the software,
C = 1 denotes missingness; the paper uses R = 1 for observation, so R = 1 - C.
The current package and manuscript use the name jointbounds. This audit maps
implemented procedures to the theory; it does not establish the population
regularity assumptions for a particular dataset.

| Paper functionality | R interface |
|---|---|
| Net and separated L2 models; empirical sharp sieve and CS outer bounds | `l2_bounds(model = "net" / "separated", method = "sharp" / "cs")` |
| Zero-budget submodels | Zero `delta_M`/`delta_R` or `delta_K` (alias `delta_A`) |
| Sharp cross-fitted one-step inference | `l2_sharp_crossfit()`; rotating three-stage sample splitting |
| CS cross-fitted one-step inference | `l2_cs_crossfit()`; rotating three-stage residual-moment estimation by default |
| Joint lower/upper asymptotic covariance | Both cross-fit results: `influence_covariance` and `estimate_covariance`, keyed by estimator |
| Finite-grid simultaneous endpoint inference | `l2_sensitivity_band()`, `l2_multiplier_band()` |
| Sensitivity surfaces and descriptive tipping | `l2_sensitivity_grid()`, `l2_tipping_frontier()` |
| Fixed-target tipping confidence bands | `l2_tipping_band()` on retained grid fits |
| Covariate-omission benchmarks with paired scores | `l2_aipw_benchmarks()` |
| Signed and toward-zero alpha-calibration | `l2_benchmark_frontiers()`, `l2_calibration_frontier()`, `l2_calibration_band()` |
| Inner/outer compatible grid sets and frontier confidence limits | `l2_calibration_band()`: `compatibility`, `frontier`; grid inversion by default |
| Simultaneity across retained covariate sets | `simultaneous_scope = "all_benchmarks"` (or `"per_benchmark"`) |
| Bounds at point or upper-confidence calibration radii | `l2_calibrated_grid()` |
| Direct convex empirical net frontier optimization | `l2_frontier_direct()`; optional CVXR; fixed missingness budgets and target, common confounding radius |
| L-infinity comparisons | `l2_linf_bounds()` |
| Simulation and oracle comparisons | `simulate_l2_data()`, `l2_simulation_design()`, `l2_oracle_sensitivity_parameters()`, `l2_compare_dr_plugin()` |

## Interpretation and limitations

- A sharp empirical sieve result imposes finite basis moments. Approximation of
  unrestricted sharp bounds, negligible sieve bias, nuisance rates, uniqueness,
  and population regularity remain assumptions, not software guarantees.
- Numerical diagnostics report balance, budgets, multiplier behavior, and the
  objective gap. The gap alone does not certify feasibility or solver accuracy.
- CS `moment_diagnostics` identifies estimated nonpositive required moments and
  clipping. Ordinary smooth EIF inference is not justified at those cases;
  positive estimates alone cannot verify the assumptions at the true law.
- Both calibration directions use paired endpoint and benchmark influence
  scores. Toward-zero calibration is nonregular at zero full-adjustment effect
  or zero benchmark discrepancy; ordinary Wald/multiplier theory does not cover
  these boundaries.
- Confidence-frontier inversion defaults to actual evaluated grid points.
  `interpolate = TRUE` is descriptive and does not supply between-grid coverage.
  Missing upper crossings stay unavailable, never zero or the last grid radius.
- `l2_tipping_frontier()` checks both interval endpoints by default, including
  negative effects. Its interpolation is descriptive; use `l2_tipping_band()`
  for finite-grid confidence limits.
- `l2_calibrated_grid()` is a radius-selection utility. Its use of an upper
  calibration confidence limit is not a theorem giving confidence coverage for
  bounds at a data-selected radius.
- Direct convex frontier optimization is a plug-in finite-basis calculation;
  no direct-frontier EIF or confidence interval is asserted.
- Legacy binary/prevalence functions remain available but are not replacements
  for the paper's continuous-outcome L2 estimators.

## Validation added in this review

Regression tests cover grid versus interpolated inversion, inner/outer set
inclusion, fixed-target paired scores, negative-effect tipping, joint covariance
and singular covariance cases, and direct optimization against analytic feasible
and infeasible examples. Existing tests cover both calibration directions,
cross-benchmark simultaneity, upper-radius selection, solver diagnostics, net and
separated models, oracle estimators, and L-infinity comparisons.
## October 29 alignment

The sharp optimized-loss derivatives, CS primitive moment and radius scores,
ATE endpoint signs, joint covariance scaling, paired calibration scores,
Bonferroni adjustment over searched radii, multiplier maxima, and inner/outer
grid inversion were checked against the current manuscript. Zero-budget
submodels and finite-basis versus unrestricted sharp targets remain distinct.

The substantive change is CS sample splitting. Residual-power regressions now
use observations separate from structural-nuisance training and score
evaluation, matching the practical nuisance-estimation remark following the
CS EIF proposition. `moment_strategy = "shared_training"` retains the earlier
implementation for reproducibility and requires additional stability
conditions. Existing simulation output is not retroactively changed. The
controlled oracle-error mode does not fit regressions and is unchanged.

A regression test records the actual rows passed to structural fitting,
residual regression, and score evaluation and checks that all three sets are
disjoint in every rotation. Population moment/rate assumptions, regularity of
active budget multipliers, and negligible sieve bias cannot be certified by
these tests or by numerical optimizer diagnostics.

The complete test suite passed 638 assertions with no failures or errors;
one slow legacy coverage test was skipped. Twelve warnings concerned small
legacy test samples, dependencies built under a newer R patch version, and
CVXR/Matrix class registration. Direct-frontier analytic tests passed despite
those dependency warnings. A separate integration check verified forwarding
of `moment_strategy = "shared_training"` through the joint simulation helper.
R CMD build and R CMD check completed with zero errors after using a
Windows-supported locale. Installation, namespace/code checks, documentation
consistency, and examples passed; tests were run separately, and PDF manual
and vignette rebuilding were excluded. The sole check warning identified
obsolete duplicate help aliases. Those aliases were removed, all 51 aliases
were verified unique, and all 45 Rd files then passed tools::checkRd.
The initial check attempt failed because Windows R rejected the inherited
C.UTF-8 locale; the repeated check used LC_ALL=C.
