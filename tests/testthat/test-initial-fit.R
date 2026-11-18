test_that("initial reference reuses structural fits throughout the pipeline", {
  d <- simulate_l2_data(120,seed=218)[,c("Y","A","C","X1","X2")]
  x <- c("X1","X2")
  init <- l2_initial_fit(d,"Y","A","C",x,nuisance_method="glm",folds=3,
    seed=19,clip=.03,response_clip=.05)
  expect_equal(init$ate$estimate,mean(init$scores$aipw))
  expect_equal(init$ate$se,sd(init$scores$aipw)/sqrt(nrow(d)))
  expect_equal(init$influence_function$ate,init$scores$aipw-init$ate$estimate)
  expect_equal(length(init$nuisance$e1),nrow(d))
  for(v in 1:3) {
    i <- init$fold_id==v
    expect_equal(init$nuisance$rho0[i],init$splits[[v]]$evaluation$rho0)
  }
  original <- l2_split_nuisance
  calls <- 0L
  local_mocked_bindings(l2_split_nuisance=function(...) {
    calls <<- calls+1L;original(...)
  },.package="jointbounds")
  for(method in c("sharp","cs")) {
    surf <- suppressWarnings(l2_sensitivity_band(d,"Y","A","C",x,
      grid=data.frame(delta_M=c(0,.05),delta_K=c(0,.05)),bound_method=method,
      nuisance_fit=init,keep_fits=TRUE,B=12))
    expect_equal(calls,0L)
    expect_true(all(vapply(surf$fits,`[[`,logical(1),"reused_nuisance_fit")))
    z <- surf$fits[[1L]]
    est <- if(method=="sharp") "dr" else "eif"
    expect_equal(subset(z$ate,estimator==est)$estimate,rep(init$ate$estimate,2))
    expect_equal(z$scores[[paste0(est,"_lower")]],init$influence_function$ate)
    expect_identical(z$fold_id,init$fold_id)
    expect_equal(z$parameters$response_clip,.05)
    b <- l2_aipw_benchmarks(d,"Y","A","C",x,subsets=list(full=x),
      variance=TRUE,reference_fit=surf)
    expect_equal(calls,0L)
    expect_identical(b$full_scores$aipw,init$influence_function$ate)
    expect_identical(b$full[["aipw"]],init$ate$estimate)
  }
  b <- l2_aipw_benchmarks(d,"Y","A","C",x,subsets=list(reduced="X2"),
    variance=TRUE,reference_fit=init)
  expect_equal(calls,3L) # reduced set refits once per rotation; full set never refits
  expect_identical(b$full_scores$aipw,init$influence_function$ate)
})

test_that("cache validation protects data order and cross-fitting roles", {
  d <- simulate_l2_data(80,seed=217)[,c("Y","A","C","X1","X2")]
  init <- l2_initial_fit(d,"Y","A","C",c("X1","X2"),nuisance_method="glm",folds=2,
    strategy="shared_training")
  expect_error(l2_sharp_crossfit(d,"Y","A","C",c("X1","X2"),nuisance_fit=init),"three_way")
  z <- l2_cs_crossfit(d,"Y","A","C",c("X1","X2"),nuisance_fit=init,delta_M=0,delta_K=0)
  expect_equal(z$parameters$moment_strategy,"shared_training")
  expect_error(l2_cs_crossfit(d,"Y","A","C",c("X1","X2"),nuisance_fit=init,
    moment_strategy="three_way"),"match")
  dd <- d[nrow(d):1,];rownames(dd)<-NULL
  expect_error(l2_cs_crossfit(dd,"Y","A","C",c("X1","X2"),nuisance_fit=init),"order")
  dd <- d;dd$X1[1]<-dd$X1[1]+1
  expect_error(l2_cs_crossfit(dd,"Y","A","C",c("X1","X2"),nuisance_fit=init),"match")
  expect_error(l2_cs_crossfit(d,"Y","A","C",c("X1","X2"),nuisance_fit=init,
    nuisance_source="oracle_error"),"oracle_error")
  expect_error(l2_cs_crossfit(d,"Y","A","C",c("X1","X2"),nuisance_fit=list()),"l2_initial_fit")
})

test_that("cached structural fits do not override second-stage learners", {
  d <- simulate_l2_data(120, seed = 221)[, c("Y", "A", "C", "X1", "X2")]
  x <- c("X1", "X2")
  init <- suppressWarnings(l2_initial_fit(d, "Y", "A", "C", x,
    nuisance_method = "glm", folds = 3, seed = 17))
  seen_method <- character()
  seen_library <- list()
  local_mocked_bindings(l2_pseudo_regression = function(train, X, pseudo,
      use, test, method, library, ...) {
    seen_method <<- c(seen_method, method)
    seen_library[[length(seen_library) + 1L]] <<- library
    rep(mean(pseudo), nrow(test))
  }, .package = "jointbounds")
  lib <- c("SL.mean", "SL.ranger")
  sharp <- suppressWarnings(l2_sharp_crossfit(d, "Y", "A", "C", x,
    delta_M = .04, delta_K = .03, nuisance_fit = init,
    second_stage_method = "SuperLearner", sl_libs = lib))
  expect_true(length(seen_method) > 0L)
  expect_true(all(seen_method == "SuperLearner"))
  expect_true(all(vapply(seen_library, identical, logical(1), lib)))
  expect_identical(sharp$parameters$second_stage_method, "SuperLearner")
  expect_identical(sharp$parameters$sl_lib_second_stage, lib)

  seen_method <- character(); seen_library <- list()
  cs <- suppressWarnings(l2_cs_crossfit(d, "Y", "A", "C", x,
    delta_M = .04, delta_K = .03, nuisance_fit = init,
    second_stage_method = "SuperLearner", sl_libs = lib))
  expect_true(length(seen_method) > 0L)
  expect_true(all(seen_method == "SuperLearner"))
  expect_true(all(vapply(seen_library, identical, logical(1), lib)))
  expect_identical(cs$parameters$second_stage_method, "SuperLearner")
  expect_identical(cs$parameters$sl_lib_second_stage, lib)
})

test_that("pseudo regressions diagnose nonfinite responses", {
  d <- data.frame(X = 1:2)
  expect_error(jointbounds:::l2_pseudo_regression(d, "X", c(1, Inf),
    c(TRUE, TRUE), d, "glm", "SL.glm"), "pseudo-outcomes are nonfinite")
})
