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
