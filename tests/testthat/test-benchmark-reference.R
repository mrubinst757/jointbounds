test_that("standalone and reused benchmarks align to zero-radius one-step scores", {
  x <- c("X1","X2")
  d <- simulate_l2_data(150, seed=114)[,c("Y","A","C",x)]
  common <- list(data=d,Y="Y",A="A",C="C",X=x,folds=5,seed=23,
                 propensity_clip=.02,nuisance_method="glm")
  args <- c(common,list(subsets=list("X2",character()),variance=TRUE))
  b <- do.call(l2_aipw_benchmarks,args)
  for (method in c("sharp","cs")) {
    surf <- do.call(l2_sensitivity_band,c(common,list(grid=data.frame(delta_M=0,delta_K=0),
      bound_method=method,keep_fits=TRUE,B=12)))
    z <- surf$fits[[1]]
    key <- if(method=="sharp") "dr_lower" else "eif_lower"
    expect_equal(b$full[["aipw"]],surf$results$estimate[1])
    expect_equal(b$full_scores$aipw,z$scores[[key]])
    alt <- args; alt$seed <- 50; alt$folds <- 2; alt$reference_fit <- surf
    reused <- do.call(l2_aipw_benchmarks,alt)
    expect_identical(reused$full_scores$aipw,z$scores[[key]])
    expect_equal(reused$results,b$results,tolerance=1e-12)
    expect_identical(reused$metadata$fold_id,z$fold_id)
    sc <- reused$target_scores[["1_aipw"]]-z$scores[[key]]
    expect_equal(reused$scores[["1_aipw"]],sc)
    expect_equal(subset(reused$results,benchmark_id==1 & estimator=="aipw")$se,
                 sd(sc)/sqrt(nrow(d)))
    alt$data <- d[rev(seq_len(nrow(d))),]; rownames(alt$data)<-NULL
    expect_error(do.call(l2_aipw_benchmarks,alt),"order")
  }
})

test_that("stochastic structural predictions are shared at zero", {
  original <- l2_learner_predict
  local_mocked_bindings(l2_learner_predict=function(...) {
    original(...) + stats::runif(1,-.0001,.0001)
  }, .package="jointbounds")
  x <- c("X1","X2");d <- simulate_l2_data(90,seed=121)[,c("Y","A","C",x)]
  z <- l2_sharp_crossfit(d,"Y","A","C",x,folds=3,delta_M=0,delta_K=0,
                         nuisance_method="glm",seed=5)
  b <- l2_aipw_benchmarks(d,"Y","A","C",x,subsets=list("X2"),variance=TRUE,
                          nuisance_method="glm",seed=5,folds=3)
  expect_equal(z$scores$dr_lower,z$scores$dr_upper)
  expect_equal(b$full_scores$aipw,z$scores$dr_lower)
  z <- l2_cs_crossfit(d,"Y","A","C",x,folds=2,delta_M=0,delta_K=0,
                      nuisance_method="glm",seed=5,moment_strategy="shared_training")
  b <- l2_aipw_benchmarks(d,"Y","A","C",x,subsets=list(x),reference_fit=z,variance=TRUE)
  expect_identical(b$full_scores$aipw,z$scores$eif_lower)
  expect_equal(b$results$se,c(0,0))
})
