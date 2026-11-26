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
  expect_gt(max(ans$arm_results$delta_M, ans$arm_results$delta_K), 0)
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

test_that("radius benchmarks sample retained sets within p_z", {
  set.seed(18)
  n <- 600
  dat <- as.data.frame(matrix(rnorm(n * 4), ncol = 4))
  names(dat) <- paste0("X", 1:4)
  dat$A <- rbinom(n, 1, plogis(.2 * dat$X1 - .1 * dat$X2))
  R <- rbinom(n, 1, plogis(1 + .1 * dat$X3))
  dat$C <- 1 - R
  yf <- dat$A + dat$X1 + dat$X2 + rnorm(n)
  dat$Y <- ifelse(R == 1, yf, NA_real_)
  args <- list(data = dat, Y = "Y", A = "A", C = "C",
    X = paste0("X", 1:4), p_z = c(1, 2), max_per_pz = 2,
    subset_selection = "random", folds = 2, seed = 27)
  a <- do.call(l2_radius_benchmarks, args)
  b <- do.call(l2_radius_benchmarks, args)
  expect_equal(nrow(a$results), 4L)
  expect_equal(as.integer(table(a$results$p_z)), c(2L, 2L))
  expect_equal(a$subsets, b$subsets)
  expect_true(all(a$results$n_omitted == 4L - a$results$p_z))
})

test_that("radius benchmark input validation is informative", {
  dat <- data.frame(Y = 1:8, A = rep(0:1, 4), C = 0, X = rnorm(8))
  expect_error(l2_radius_benchmarks(dat, "Y", "A", "C", "X",
    subsets = list("X")), "must omit at least one")
  expect_error(l2_radius_benchmarks(dat, "Y", "A", "C", "X",
    probability_clip = .6), "strictly between")
  expect_error(l2_radius_benchmarks(dat, "Y", "A", "C", "X",
    p_z = 1), "between 0")
})

test_that("radius benchmark plot combines sharp and CS frontier bands", {
  skip_if_not_installed("ggplot2")
  bench <- structure(list(results = data.frame(
    label = c("Omit X1", "Omit X2"), omitted=c("X1", "X2"), n_omitted=c(1L,1L),
    delta_M = c(.05, .12), delta_K = c(.18, .08))),
    class = "jointbounds_l2_radius_benchmarks")
  make_band <- function(shift, estimator) structure(list(
    frontier = data.frame(delta_M = c(0, .1, .2),
      estimate = c(.4, .3, .2) + shift,
      pointwise_low = c(.35, .25, .15) + shift,
      pointwise_high = c(.45, .35, .25) + shift,
      simultaneous_low = c(.32, .22, .12) + shift,
      simultaneous_high = c(.48, .38, .28) + shift),
    radius_x = "delta_M", radius_y = "delta_K",
    endpoint_estimator = estimator, conf_level = .95),
    class = "jointbounds_l2_calibration_band")
  sharp <- make_band(.05, "dr"); cs <- make_band(0, "eif")
  p <- plot(bench, frontier = list(`Sharp sieve` = sharp, `CS outer` = cs),
            uncertainty = "simultaneous")
  expect_s3_class(p, "ggplot")
  q <- l2_plot_radius_benchmarks(list(cs=cs,sharp=sharp),bench,
    frontier_type="both",uncertainty="both",color_by="variables_omitted",
    xlim=c(0,.15),ylim=c(0,.5))
  expect_equal(q$coordinates$limits$x,c(0,.15))
  expect_equal(q$coordinates$limits$y,c(0,.5))
  expect_setequal(q$data$.bound,c("CS outer frontier","Sharp sieve approximation frontier"))
  expect_equal(q$layers[[1]]$aes_params$linetype,"dashed")
  expect_length(unique(ggplot2::ggplot_build(q)$data[[2]]$group),2)
  expect_error(l2_plot_radius_benchmarks(cs,bench,frontier_type="both"),"requires both")
  bench$metadata <- list(selection_source="sizes",X=c("X1","X2","X3"),p_z=1,max_per_pz=2,subset_selection="random")
  bench$results$omitted <- c("X1,X2","X2,X3")
  q <- l2_plot_radius_benchmarks(cs,bench,color_by="variables_omitted")
  expect_match(q$labels$caption,"up to 2 subsets.*at random")
  expect_match(q$labels$caption,"remaining 3 - p_z")
  expect_equal(tail(q$layers,1)[[1]]$data$label,c("X1, X2","X2, X3"))
  expect_gte(length(p$layers), 4L)
  expect_error(l2_plot_radius_benchmarks(
    structure(data.frame(delta_M = 0, delta_K = 0),
      class = c("jointbounds_l2_frontier", "data.frame")),
    bench, uncertainty = "simultaneous"), "requires every object")
})
