test_that("radius-scale benchmarks return finite arm and combined radii", {
  set.seed(41)
  n <- 800
  x1 <- rnorm(n); x2 <- rnorm(n)
  A <- rbinom(n, 1, plogis(.3 * x1 - .2 * x2))
  R <- rbinom(n, 1, plogis(.7 - .2 * A + .4 * x1 + .3 * x2))
  y_full <- .5 * A + x1 + .7 * x2 + rnorm(n)
  dat <- data.frame(Y = ifelse(R == 1, y_full, NA_real_), A = A,
                    C = 1 - R, X1 = x1, X2 = x2)
  ans <- l2_radius_benchmarks(dat, "Y", "A", "C", c("X1", "X2"),
    subsets = list(`omit X2` = "X1"), folds = 3, seed = 9)
  expect_s3_class(ans, "jointbounds_l2_radius_benchmarks")
  expect_equal(nrow(ans$arm_results), 2L)
  expect_equal(nrow(ans$results), 1L)
  expect_true(all(is.finite(ans$arm_results$delta_M)))
  expect_true(all(is.finite(ans$arm_results$delta_K)))
  expect_true(all(ans$arm_results$delta_M >= 0))
  expect_true(all(ans$arm_results$delta_K >= 0))
  expect_equal(ans$results$delta_M, max(ans$arm_results$delta_M))
  expect_equal(ans$results$delta_K, max(ans$arm_results$delta_K))
  expect_identical(ans$metadata$uncertainty, "not estimated")
})

test_that("default radius benchmarks omit each covariate and all covariates", {
  set.seed(12)
  n <- 500
  dat <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  dat$A <- rbinom(n, 1, .5)
  R <- rbinom(n, 1, .8)
  dat$C <- 1 - R
  y <- dat$A + dat$X1 + dat$X2 + rnorm(n)
  dat$Y <- ifelse(R == 1, y, NA_real_)
  ans <- l2_radius_benchmarks(dat, "Y", "A", "C", c("X1", "X2"),
                              folds = 2, seed = 3)
  expect_equal(nrow(ans$results), 3L)
  expect_setequal(ans$results$label, c("omit X1", "omit X2", "omit all"))
})

test_that("radius benchmark input validation is informative", {
  dat <- data.frame(Y = 1:8, A = rep(0:1, 4), C = 0, X = rnorm(8))
  expect_error(l2_radius_benchmarks(dat, "Y", "A", "C", "X",
    subsets = list("X")), "must omit at least one")
  expect_error(l2_radius_benchmarks(dat, "Y", "A", "C", "X",
    probability_clip = .6), "strictly between")
})
