# Theory-to-software audit

This audit maps the current package to the sensitivity models and estimators in
the companion manuscript `bounds/main.tex`. The package uses `C = 1` for a
missing outcome, whereas the manuscript uses `R = 1` for an observed outcome;
thus every implementation uses `R = 1 - C`.

## Implemented manuscript functionality

| Manuscript component | Package implementation | Theory check |
|---|---|---|
| Separated likelihood-ratio model | `l2_bounds(model = "separated")` | Uses `bar_pi = delta * (1-rho)`, `M = 1 + bar_pi(R-1)`, `L = MK`, the two conditional-normalization constraints, the mixture lower envelope, and the two stated divergences. |
| CS outer radius | `l2_bounds(method = "cs")` | Implements `delta_R sqrt(V_R) + delta_A sqrt(V_A + delta_R sqrt(W_RA))` arm by arm and combines arm endpoints in the correct worst-case direction for the ATE. |
| Reduced-form net-M model | `l2_bounds(model = "net")` | Enforces `M >= rho`, uses the unweighted `M-1` divergence, and uses `V_M` and `W_MA` in the outer radius. |
| Sharp empirical sieve program | `l2_bounds(method = "sharp")` | The dual KKT updates match the manuscript. Basis moments enforce normalization of both M and L. Diagnostics include divergence use, balance residuals, primal value, dual value, primal-dual gap, and maximum constraint violations. |
| MAR-only and NUC-only special cases | Set `delta_A = 0` or the relevant missingness radius to zero | Zero missingness budgets now impose `M = 1` exactly rather than approximating the equality with a diverging multiplier. |
| CS EIF estimator | `l2_cs_crossfit()` | Implements the weighted-central-moment EIF, including the `-4 nu_3 Z` correction for the fourth moment and all derivatives through `e` and `rho`. Both separated and net-M models are supported. |
| Sharp augmented estimator | `l2_sharp_crossfit()` | Implements the respondent-law regression of H, inverse respondent augmentation, `d_e = E[sigma Y(M-L)|X]`, and the interior plus active-boundary terms in `d_rho`. Structural nuisance/dual fitting, pseudo-outcome regression, and final score evaluation use disjoint folds in a rotating three-stage construction. Both separated and net-M models are supported. |
| Pointwise and outward Wald inference | `l2_cs_crossfit()` and `l2_sharp_crossfit()` | Returns endpoint intervals, row-level influence scores, and an outward interval for the whole identified set. Plug-in standard errors are explicitly labeled as conditional on fitted nuisances. |
| Simultaneous sensitivity bands | `l2_sensitivity_band()` and `l2_multiplier_band()` | Uses the standardized supremum multiplier process jointly over both endpoints and every evaluated grid point. |
| L-infinity comparison | `l2_linf_bounds()` | Provides centered Holder outer bounds and a sharp empirical-sieve LP. The LP enforces both nonnegative likelihood ratios and `M >= 1-bar_pi`, including radii greater than one. |
| Sensitivity surfaces and tipping | `l2_sensitivity_grid()` and `l2_tipping_frontier()` | The tipping frontier uses the first actual threshold crossing and linear interpolation; it no longer labels the merely closest non-crossing point as a tip. |
| Joint covariate-omission calibration | `l2_aipw_benchmarks()`, `l2_benchmark_frontiers()`, and equal-radius helpers | The benchmark is reduced-adjustment minus full-adjustment. Its AIPW score is paired on the same observations. Frontier inversion uses the signed benchmark on the estimand scale. |
| Calibration uncertainty | `l2_benchmark_band()`, `l2_calibration_band()`, and `l2_equal_radius_comparison()` | `l2_calibration_band()` combines endpoint and reduced-adjustment AIPW influence scores observation by observation and jointly inverts pointwise and multiplier bands over the sensitivity grid. `l2_equal_radius_comparison()` remains a simpler transformed-radius approximation. |
| Binary composite and SDE prevalence bounds | `l2_prevalence_bounds()` | Implements the capped L2 support optimizer and the manuscript's fixed-q EIFs. It also includes the active norm-constraint multiplier term required because the norm is under the estimated `P_X`. |
| Simulation and oracle checks | `simulate_l2_data()`, `l2_oracle_sensitivity_parameters()`, and `l2_compare_dr_plugin()` | Supports MAR, NUC, both, or neither; bounded or unbounded selection links; oracle sensitivity radii; sharp, CS, and L-infinity comparisons; and true-EIF nuisance-error experiments. |

## Important interpretation and inference boundaries

- A `sharp_sieve` result is sharp for the stated empirical basis restrictions,
  not automatically for unrestricted conditional normalization. Always inspect
  the returned numerical diagnostics.
- Ordinary Wald and multiplier results presume a unique regular optimizer,
  stable active set, nonzero smooth moments, and the rate conditions in the
  manuscript. The package does not relabel a nonregular result as regular.
- `l2_equal_radius_comparison(uncertainty = TRUE)` transforms pointwise Wald
  intervals while holding the estimated grid surface fixed. It is intentionally
  not described as the exact implicit-frontier EIF calculation.
- The manuscript's simultaneous NUC/MAR composite extension and joint-state SDE
  formulation are theoretical sharp programs rather than claimed turnkey
  software routines. The latter requires profiling or a globally certified
  nonconvex solver, as the manuscript itself notes.
- The legacy `mar_bounds()` interface implements the earlier binary-outcome
  pointwise models. It should not be interpreted as an alternative interface to
  the new likelihood-ratio L2 model.

## Corrections made during this audit

1. Enforced the nonnegative-R mixture envelope in the sharp L-infinity LP.
2. Imposed zero missingness budgets exactly.
3. Added primal-dual and constraint diagnostics.
4. Added net-M cross-fitted EIF/DR estimation.
5. Added continuous multiplier bands and outward identified-set intervals.
6. Corrected tipping-frontier crossing logic.
7. Made `return_true_bounds` return the requested oracle method rather than
   silently returning CS bounds for a sharp request.
8. Added data-facing L-infinity bounds and capped-L2 alternative-estimand bounds.
9. Corrected the manuscript's support-function EIF by including the active
   L2-norm multiplier contribution.
10. Made frontier inversion add each benchmark discrepancy to its paired
    full-X AIPW estimate, so the calibration target is exactly the reduced-Z
    estimate rather than a separately fitted surface reference.
11. Added rotating three-stage cross-fitting for the sharp one-step estimator.
12. Added exact combined-score calibration-frontier bands and corresponding
    grid inversion.
13. Added tests for the single-mechanism sharp/CS equality when the sieve basis
    represents the Cauchy--Schwarz optimizer.
14. Corrected the legacy bounded-prevalence Psi2 influence-score implementation
    and imposed a finite-sample coherence projection on its nonnegative support
    terms.
