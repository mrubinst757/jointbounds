test_that("generate_ate returns the cross-fitted AIPW score", {
  d <- simulate_l2_data(600, violation = "none", seed = 501)
  fit <- generate_ate(d, "Y", "A", "C", c("X1", "X2"),
    folds = 2, nuisance_method = "glm", seed = 7,
    diagnostic_control = list(warn = FALSE), keep_nuisance = TRUE)
  expect_s3_class(fit, "jointbounds_ate")
  expect_true(all(is.finite(unlist(fit$ate[c("estimate", "se",
    "conf_low", "conf_high")]))))
  expect_equal(fit$ate$estimate,
    fit$arm$estimate[fit$arm$arm == 1] - fit$arm$estimate[fit$arm$arm == 0])
  expect_equal(mean(fit$influence_function$ate), 0, tolerance = 1e-12)

  y <- d$Y
  y[!is.finite(y)] <- 0
  score <- vector("list", 2L)
  for (a in 0:1) {
    e <- fit$nuisance[[paste0("e", a)]]
    rho <- fit$nuisance[[paste0("rho", a)]]
    mu <- fit$nuisance[[paste0("mu", a)]]
    S <- as.numeric(d$A == a & d$C == 0)
    score[[a + 1L]] <- mu + S / (e * rho) * (y - mu)
  }
  expect_equal(fit$ate$estimate, mean(score[[2L]] - score[[1L]]))
  expect_lt(abs(fit$ate$estimate - .7), .3)
})

test_that("generate_ate supports intercept-only and SuperLearner fits", {
  d <- simulate_l2_data(180, violation = "none", seed = 502)
  intercept <- generate_ate(d, "Y", "A", "C", character(), folds = 2,
    nuisance_method = "glm", diagnostic_control = list(warn = FALSE))
  expect_true(is.finite(intercept$ate$estimate))

  skip_if_not_installed("SuperLearner")
  sl <- generate_ate(d, "Y", "A", "C", c("X1", "X2"), folds = 2,
    nuisance_method = "SuperLearner", sl_lib_prop = "SL.glm",
    sl_lib_miss = "SL.glm", sl_lib_outcome = "SL.glm",
    diagnostic_control = list(warn = FALSE))
  expect_true(is.finite(sl$ate$estimate))
  expect_equal(length(sl$influence_function$ate), nrow(d))
})

test_that("generate_ate consolidates positivity warnings", {
  d <- simulate_l2_data(160, violation = "none", seed = 503)
  expect_warning(
    fit <- generate_ate(d, "Y", "A", "C", c("X1", "X2"), folds = 2,
      nuisance_method = "glm",
      diagnostic_control = list(max_weight = .01)),
    "ATE diagnostics flagged"
  )
  expect_false(fit$diagnostic_summary$diagnostic_ok)
  expect_true(any(grepl("large_inverse_probability_weight",
                        fit$diagnostic_table$diagnostic_codes)))
})
