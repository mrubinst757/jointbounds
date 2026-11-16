test_that("shared clipping and overrides affect only the intended nuisance", {
  d <- simulate_l2_data(90,seed=123)
  all <- rep(TRUE,nrow(d))
  args <- list(data=d,X=c("X1","X2"),tr=all,te=all,rate=.1,
    scale=c(e=4,rho=4,mu=.5,regression=.5),sign=c(e=-1,rho=-1,mu=1,regression=1))
  old <- do.call(l2_split_nuisance_oracle_error,args)$test
  z <- do.call(l2_split_nuisance_oracle_error,c(args,list(clip=.15,propensity_clip=.05)))$test
  expect_equal(z$e1,clip_probs(old$e1,.05))
  expect_equal(z$e0,1-z$e1)
  for(a in 0:1) expect_equal(z[[paste0("rho",a)]],clip_probs(old[[paste0("rho",a)]],.15))
  expect_equal(z$mu0,old$mu0)
  expect_true(any(z$rho0!=old$rho0))
  for(fun in list(l2_cs_crossfit,l2_sharp_crossfit)) {
    expect_error(fun(NULL,"Y","A","C","X",response_clip=0),"response_clip")
    expect_error(fun(NULL,"Y","A","C","X",clip=.5),"propensity_clip")
    base <- suppressWarnings(fun(d,"Y","A","C",c("X1","X2"),nuisance_source="oracle_error"))
    alt <- suppressWarnings(fun(d,"Y","A","C",c("X1","X2"),nuisance_source="oracle_error",clip=.1))
    expect_equal(alt$parameters$response_clip,.1)
    expect_equal(subset(base$ate,estimator=="true_eif"),subset(alt$ate,estimator=="true_eif"))
  }
})

test_that("reference reuse inherits both clipping thresholds", {
  d <- simulate_l2_data(100,seed=9)[,c("Y","A","C","X1","X2")]
  for(fun in list(l2_cs_crossfit,l2_sharp_crossfit)) {
    z <- fun(d,"Y","A","C",c("X1","X2"),delta_M=0,delta_K=0,
      nuisance_method="glm",clip=.08,propensity_clip=.03)
    b <- l2_aipw_benchmarks(d,"Y","A","C",c("X1","X2"),subsets=list("X2"),variance=TRUE,reference_fit=z)
    expect_equal(b$metadata$response_clip,.08)
    expect_equal(b$metadata$propensity_clip,.03)
    expect_identical(b$full_scores$aipw,z$benchmark_reference$centered$aipw)
    standalone <- l2_aipw_benchmarks(d,"Y","A","C",c("X1","X2"),subsets=list("X2"),variance=TRUE,nuisance_method="glm",clip=.08,propensity_clip=.03)
    expect_equal(b$results,standalone$results)
    z$benchmark_reference$response_clip <- NULL
    legacy <- l2_aipw_benchmarks(d,"Y","A","C",c("X1","X2"),subsets=list("X2"),reference_fit=z)
    expect_equal(legacy$metadata$response_clip,1e-4)
  }
})

test_that("GLM and SL fitting paths apply both clipping settings", {
  d <- simulate_l2_data(100,seed=20)
  for(method in c("glm","SuperLearner")) {
    z <- suppressWarnings(generate_ate(d,"Y","A","C",c("X1","X2"),
      nuisance_method=method,folds=2,clip=.2,propensity_clip=.1,keep_nuisance=TRUE))
    expect_true(all(z$nuisance$e1>=.1 & z$nuisance$e1<=.9))
    expect_true(all(z$nuisance$rho0>=.2 & z$nuisance$rho0<=.8))
    expect_true(all(z$nuisance$rho1>=.2 & z$nuisance$rho1<=.8))
  }
  args <- list(data=d,Y="Y",A="A",C="C",X=c("X1","X2"),nuisance_method="glm",method="cs",seed=3)
  grid <- do.call(l2_sensitivity_grid,c(args,list(grid=data.frame(delta_M=.1,delta_K=.1),clip=.2)))
  single <- do.call(l2_bounds,c(args,list(delta_M=.1,delta_K=.1,clip=.2)))
  expect_equal(grid$results$lower,single$ate$lower)
  expect_equal(grid$results$upper,single$ate$upper)
  legacy <- l2_aipw_benchmarks(d,"Y","A","C",c("X1","X2"),subsets=list("X2"),
    nuisance_method="glm",reference_strategy="legacy",clip=.2)
  expect_equal(legacy$metadata$response_clip,.2)
})
