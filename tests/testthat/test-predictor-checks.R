test_that("diagnostics distinguish constants, duplicates, and linear dependence", {
  x <- data.frame(a=1:12,b=(1:12)^2,copy=1:12,sum=(1:12)+(1:12)^2,
                  constant=2,group=factor(rep("only",12)))
  original <- x
  expect_warning(z <- l2_predictor_check(x,"fold 2 response arm 1"),"fold 2 response arm 1.*redundancy")
  expect_identical(x,original)
  expect_equal(z$rank,3)
  expect_match(z$constant_predictors,"constant")
  expect_match(z$constant_predictors,"group")
  expect_match(z$duplicate_columns,"copy")
  expect_match(z$dependent_columns,"sum")
  expect_equal(z$status,"redundant")
  expect_no_warning(good <- l2_predictor_check(x[c("a","b")],"good"))
  expect_equal(good$rank,good$design_columns)
  expect_equal(good$status,"ok")
})

test_that("rank checks use each model's training subset", {
  set.seed(81)
  n <- 150; d <- data.frame(X1=rnorm(n),X2=rnorm(n),A=rep(0:1,75),C=rep(c(0,0,1),50))
  d$Y <- d$X1+rnorm(n); d$Y[d$C==1] <- NA
  d$arm_only <- d$A*d$X1
  original <- d
  z <- suppressWarnings(l2_initial_fit(d,"Y","A","C",c("X1","X2","arm_only"),
    nuisance_method="glm",folds=3))
  expect_identical(d,original)
  tab <- z$predictor_diagnostics
  expect_equal(nrow(tab),15)
  expect_true(all(tab$status[grepl("treatment",tab$context)]=="ok"))
  expect_true(all(tab$status[grepl("arm 0",tab$context)]=="redundant"))
  expect_true(all(grepl("arm_only",tab$constant_predictors[grepl("arm 0",tab$context)])))
  for(v in 1:3) expect_true(any(grepl(paste("fold",v),tab$context)))
  surf <- suppressWarnings(l2_sensitivity_band(d,"Y","A","C",c("X1","X2","arm_only"),
    grid=data.frame(delta_M=.1,delta_K=.1),bound_method="cs",nuisance_fit=z,B=10))
  expect_true(any(grepl("CS residual",surf$predictor_diagnostics$context)))
  expect_true(all(surf$predictor_diagnostics$grid_id==1))
})

test_that("SL.glm internal training subsets are checked and retained", {
  xx <- data.frame(x=seq(-1,1,length.out=40),duplicate=seq(-1,1,length.out=40))
  sink <- l2_predictor_sink()
  z <- suppressWarnings(l2_learner_predict(rep(0:1,20),xx,xx,"binomial","SuperLearner",c("SL.glm","SL.mean"),
    context="fold 1 treatment",diagnostic_sink=sink))
  tab <- l2_predictor_table(sink)
  expect_true(any(grepl("SL.glm internal",tab$context)))
  expect_true(any(tab$n_train<40))
  expect_true(all(tab$status=="redundant"))
  expect_length(z,40)
})
