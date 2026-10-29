test_that("calibration targets and derivatives cover every sign combination", {
  for (reference in c(-2, 2)) for (delta in c(-.5, .5)) {
    alpha <- 1.7
    z <- l2_calibration_target(reference, delta, alpha)
    expect_equal(z$target, reference - sign(reference) * alpha * abs(delta))
    reduced <- reference + delta
    h <- 1e-6
    f <- function(x, y) l2_calibration_target(x, y - x, alpha)$target
    expect_equal((f(reference + h, reduced) - f(reference - h, reduced)) /
      (2 * h), 1 - z$coefficient, tolerance = 1e-7)
    expect_equal((f(reference, reduced + h) - f(reference, reduced - h)) /
      (2 * h), z$coefficient, tolerance = 1e-7)
  }
  expect_equal(l2_calibration_target(.2, .4)$target, -.2)
  expect_equal(l2_calibration_target(0, .4)$target, 0)
  expect_equal(l2_calibration_target(2, 0)$target, 2)
  expect_equal(l2_calibration_target(2, .5, 2, "signed")$target, 3)
  for (bad in list(0, -1, Inf, NA_real_, c(1, 2)))
    expect_error(l2_calibration_target(1, .2, bad), "alpha")
})

test_that("asymmetric frontier inversion uses the requested direction", {
  grid <- data.frame(delta_M = 0, delta_K = c(0, .1, .2, .4),
    method = "sharp_sieve", estimate_reference = 1,
    lower = 1 - c(0, .1, .2, .4), upper = 1 + 2 * c(0, .1, .2, .4))
  toward <- l2_calibration_frontier(grid, .2, type = "discrepancy")
  signed <- l2_calibration_frontier(grid, .2, type = "discrepancy",
    calibration = "signed")
  expect_equal(toward$target, .8)
  expect_equal(toward$minimum_delta_K, .2)
  expect_equal(signed$target, 1.2)
  expect_equal(signed$minimum_delta_K, .1)
  explicit <- l2_calibration_frontier(grid, 1.2, type = "target", alpha = 2)
  expect_equal(explicit$target, 1.2)
  b <- data.frame(estimator = "aipw", benchmark_id = 1, label = "omit",
    p_z = 1, full_estimate = 1, delta_bench = .2)
  expect_equal(l2_benchmark_frontiers(grid, b, alpha = 2)$target, .6)
})

test_that("frontier bands transform paired scores and reject nonregular ties", {
  n <- 12L
  full <- rep(c(-.1, .1), 6)
  red <- rep(c(-.2, 0, .2), 4)
  lower <- rep(c(-.3, .3, .1), 4)
  upper <- -lower
  lower <- lower - mean(lower); upper <- upper - mean(upper)
  grid <- expand.grid(delta_M = 0, delta_K = c(0, .5))
  for (ref in c(-1, 1)) for (delta in c(-.2, .2)) {
    surface <- structure(list(grid = grid, bound_method = "cs",
      fits = lapply(grid$delta_K, function(w) list(
        ate = data.frame(endpoint = c("lower", "upper"), estimator = "eif",
          estimate = c(ref - w, ref + w)),
        scores = list(eif_lower = lower, eif_upper = upper)))),
      class = "jointbounds_l2_sensitivity_band")
    b <- structure(list(results = data.frame(benchmark_id = 1, label = "omit",
      estimator = "aipw", full_estimate = ref, reduced_estimate = ref + delta),
      target_scores = list(`1_aipw` = red), full_scores = list(aipw = full)),
      class = "jointbounds_l2_benchmarks")
    z <- l2_calibration_band(surface, b, alpha = 1.5, B = 19)
    target_score <- full - 1.5 * sign(ref) * sign(delta) * (red - full)
    expected <- lower - target_score
    expect_equal(z$scores[, 1], expected - mean(expected))
    expect_equal(z$gaps$target, rep(ref - sign(ref) * 1.5 * abs(delta), 4))
    expect_equal(unique(z$gaps$endpoint), c("lower", "upper"))
    legacy <- l2_calibration_band(surface, b, calibration = "signed", B = 19)
    expect_equal(legacy$gaps$target, rep(ref + delta, 4))
    b$results$reduced_estimate <- ref
    expect_error(l2_calibration_band(surface, b, B = 19), "nonzero")
    b$results$full_estimate <- 0
    expect_error(l2_calibration_band(surface, b, B = 19), "nonzero")
  }
})

test_that("joint inversion rejects false zero inclusion and preserves missing edges", {
  score <- rep(c(-.001, .001), 10)
  surface <- structure(list(grid = data.frame(delta_M = 0, delta_K = c(0, 1)),
    bound_method = "cs", fits = lapply(c(0, 1), function(w) list(
      ate = data.frame(endpoint = c("lower", "upper"), estimator = "eif",
        estimate = c(1 - w, 1 + w)),
      scores = list(eif_lower = score, eif_upper = score)))),
    class = "jointbounds_l2_sensitivity_band")
  b <- structure(list(results = data.frame(benchmark_id = 1, label = "omit",
    estimator = "aipw", full_estimate = 2, reduced_estimate = 1.5),
    target_scores = list(`1_aipw` = -score), full_scores = list(aipw = score)),
    class = "jointbounds_l2_benchmarks")
  expect_warning(z <- l2_calibration_band(surface, b, B = 29, interpolate = TRUE), "coverage between evaluated grid points")
  expect_equal(z$frontier$estimate, .5)
  expect_gt(z$frontier$pointwise_low, 0)
  expect_gt(z$frontier$pointwise_high, z$frontier$pointwise_low)
  b$target_scores$`1_aipw` <- score * 10000
  expect_warning(z <- l2_calibration_band(surface, b, B = 29, interpolate = TRUE), "coverage between evaluated grid points")
  expect_equal(z$frontier$pointwise_low, 0)
  expect_true(is.na(z$frontier$pointwise_high))
  expect_true(is.na(z$frontier$simultaneous_high))
})

test_that("zero-budget sharp endpoints share the reference AIPW score", {
  d <- simulate_l2_data(120, seed = 802, violation = 'both', density_ratio_link = 'bounded')
  for (source in c('estimated', 'oracle_error')) {
    sharp <- l2_sharp_crossfit(d, Y='Y', A='A', C='C', X=c('X1','X2'),
      delta_M=0, delta_K=0, model='net', nuisance_method='glm',
      nuisance_source=source, folds=3, seed=12)
    z <- sharp$ate[sharp$ate$estimator == 'dr', ]
    expect_equal(z$estimate[1], z$estimate[2], tolerance=1e-12)
    expect_equal(sharp$scores$dr_lower, sharp$scores$dr_upper, tolerance=1e-12)
    if (source == 'oracle_error') {
      cs <- l2_cs_crossfit(d,Y='Y',A='A',C='C',X=c('X1','X2'),
        delta_M=0,delta_K=0,model='net',nuisance_method='glm',
        nuisance_source=source,folds=3,seed=12)
      expect_equal(z$estimate, cs$ate$estimate[cs$ate$estimator=='eif'], tolerance=1e-12)
      expect_equal(sharp$scores$dr_lower, cs$scores$eif_lower, tolerance=1e-12)
    }
  }
})

test_that("signed calibration permits zero discrepancies and scales paired scores", {
  n <- 20
  full <- rep(c(-.1,.1),10); red <- rep(c(-.2,0,.2,0),5)
  grid <- data.frame(delta_M=0,delta_K=c(0,1))
  surface <- structure(list(grid=grid,bound_method='cs',fits=lapply(c(0,1),function(w)
    list(ate=data.frame(endpoint=c('lower','upper'),estimator='eif',estimate=c(-w,w)),
      scores=list(eif_lower=full,eif_upper=-full)))),class='jointbounds_l2_sensitivity_band')
  b <- structure(list(results=data.frame(benchmark_id=1,label='omit',estimator='aipw',
    full_estimate=0,reduced_estimate=.2),full_scores=list(aipw=full),
    target_scores=list(`1_aipw`=red)),class='jointbounds_l2_benchmarks')
  expect_warning(z <- l2_calibration_band(surface,b,calibration='signed',alpha=2,B=19,interpolate=TRUE), 'coverage between evaluated grid points')
  expect_equal(unique(z$gaps$target),.4)
  expect_equal(z$scores[,1], 2*full-2*red)
  expect_equal(z$frontier$estimate,.4)
  expect_equal((z$gaps$pointwise_high-z$gaps$gap)/z$gaps$se,
    rep(qnorm(1-.05/8),4))
  b$results$reduced_estimate <- 0
  expect_warning(z <- l2_calibration_band(surface,b,calibration='signed',alpha=2,B=19,interpolate=TRUE),
    'equals its full-X')
  expect_equal(z$frontier$estimate,0)
})
