test_that("prepare_data returns correct structure and validates inputs", {
  dat <- make_test_data(n = 100)
  prep <- jointbounds:::prepare_data(dat, "Y", "A", "C", "X")
  expect_type(prep, "list")
  expect_equal(prep$n, 100)
  expect_equal(length(prep$Y_vec), 100)
  expect_equal(length(prep$A_vec), 100)
  expect_equal(length(prep$C_vec), 100)
  expect_equal(nrow(prep$X_mat), 100)
  expect_equal(ncol(prep$X_mat), 1)
  expect_equal(prep$X_names, "X")
})

test_that("prepare_data errors on invalid inputs", {
  dat <- make_test_data(n = 100)
  expect_error(jointbounds:::prepare_data(as.list(dat), "Y", "A", "C", "X"), "data.frame")
  expect_error(jointbounds:::prepare_data(dat, "Y", "A", "C", "NotACol"), "exist|undefined columns")
  dat_bad <- dat
  dat_bad$A <- dat_bad$A + 1L
  expect_error(jointbounds:::prepare_data(dat_bad, "Y", "A", "C", "X"), "0/1")
})

test_that("clip_probs keeps values in (0, 1)", {
  p <- c(0, 0.5, 1, -0.1, 1.2)
  out <- jointbounds:::clip_probs(p)
  expect_true(all(out > 0 & out < 1))
  expect_equal(out[2], 0.5)
})

test_that("expit is inverse of logit", {
  x <- seq(-2, 2, by = 0.5)
  expect_equal(jointbounds:::expit(x), 1 / (1 + exp(-x)))
  expect_equal(jointbounds:::expit(c(-Inf, Inf)), c(0, 1), ignore_attr = TRUE)
})

test_that("%||% returns second arg when first is NULL", {
  expect_equal(NULL %||% 5, 5)
  expect_equal(0 %||% 5, 0)
  expect_equal(NA %||% 5, NA)
})
