test_that("failed requested learners stop by default and preserve root cause", {
  local_mocked_bindings(l2_superlearner_fit=function(...) stop("ranger dependency unavailable"),
    .package="jointbounds")
  xx <- data.frame(x=seq(-1,1,length.out=30)); yy <- rep(0:1,15)
  expect_error(l2_learner_predict(yy,xx,xx,"binomial","SuperLearner","SL.ranger",
    context="fold 2 response arm 1"),"fold 2 response arm 1.*ranger dependency unavailable.*No substitute")
  expect_warning(z <- l2_learner_predict(yy,xx,xx,"binomial","SuperLearner","SL.ranger",
    learner_fallback="glm"),"explicitly requested GLM fallback")
  expect_equal(z,as.numeric(predict(glm(yy~x,data=xx,family=binomial()),type="response")))
})

test_that("invalid or failed library predictions cannot silently pass", {
  xx <- data.frame(x=1:10); yy <- seq_len(10)
  local_mocked_bindings(l2_superlearner_fit=function(...) list(SL.predict=rep(NA_real_,10)),.package="jointbounds")
  expect_error(l2_learner_predict(yy,xx,xx,"gaussian","SuperLearner","SL.mean"),"nonfinite")
  local_mocked_bindings(l2_superlearner_fit=function(...) list(SL.predict=rep(0,10),errorsInCVLibrary=TRUE),.package="jointbounds")
  expect_error(l2_learner_predict(yy,xx,xx,"gaussian","SuperLearner","SL.mean"),"All requested")
})

test_that("public APIs supply fold and model context and preserve the policy", {
  d <- simulate_l2_data(100,seed=829)[,c("Y","A","C","X1","X2")]
  local_mocked_bindings(l2_superlearner_fit=function(...) stop("deliberate fit failure"),.package="jointbounds")
  expect_error(l2_initial_fit(d,"Y","A","C",c("X1","X2"),sl_libs="SL.ranger"),
    "fold 1 treatment propensity.*deliberate fit failure")
  expect_error(l2_sensitivity_band(d,"Y","A","C",c("X1","X2"),
    grid=data.frame(delta_M=.1,delta_K=.1),sl_libs="SL.ranger"),"fold 1 treatment propensity")
  z <- suppressWarnings(l2_initial_fit(d,"Y","A","C",c("X1","X2"),sl_libs="SL.ranger",learner_fallback="glm"))
  expect_equal(z$parameters$learner_fallback,"glm")
  expect_equal(z$benchmark_reference$learner_fallback,"glm")
})

test_that("SuperLearner wrapper lookup works without attaching the package", {
  xx <- data.frame(x=seq(-1,1,length.out=40)); yy <- rep(0:1,20)
  expect_no_error(z <- l2_learner_predict(yy,xx,xx,"binomial","SuperLearner","SL.mean"))
  expect_equal(z,rep(.5,40))
})

test_that("SL.ranger can run through the cross-fitted helper", {
  skip_if_not_installed("ranger")
  xx <- data.frame(x=seq(-1,1,length.out=40)); yy <- rep(0:1,20)
  z <- l2_learner_predict(yy,xx,xx,"binomial","SuperLearner","SL.ranger")
  expect_length(z,40)
  expect_true(all(is.finite(z) & z>=0 & z<=1))
})


test_that("legacy nuisance fits report errors instead of returning marginal means", {
  d <- simulate_l2_data(100,seed=828)
  local_mocked_bindings(l2_quiet_superlearner=function(...) stop("original SL failure"),.package="jointbounds")
  expect_error(estimate_nuisance(as.matrix(d[c("X1","X2")]),d$A,d$C,d$Y),
    "treatment propensity; fold 1.*original SL failure")
})
