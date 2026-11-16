test_that("common libraries propagate and individual overrides win", {
  calls <- list()
  local_mocked_bindings(l2_learner_predict = function(y, xtrain, xnew, family, method, library, ...) {
    calls[[length(calls)+1L]] <<- library
    rep(mean(y), nrow(xnew))
  }, .package = "jointbounds")
  d <- simulate_l2_data(120, seed=8)
  l2_sharp_crossfit(d, "Y", "A", "C", c("X1","X2"), delta_M=0, delta_K=0,
    folds=3, sl_libs=c("SL.glm","SL.mean"), sl_lib_prop="SL.mean")
  expect_true(all(vapply(calls, function(x) identical(x,"SL.mean") ||
    identical(x,c("SL.glm","SL.mean")), logical(1))))
  expect_true(any(vapply(calls, identical, logical(1), "SL.mean")))
  expect_true(any(vapply(calls, identical, logical(1), c("SL.glm","SL.mean"))))
  funcs <- c("generate_ate","estimate_nuisance","mar_bounds","l2_bounds",
    "l2_cs_crossfit","l2_sharp_crossfit","l2_aipw_benchmarks","l2_sensitivity_grid",
    "l2_bootstrap","l2_frontier_direct","l2_linf_bounds","l2_prevalence_bounds")
  expect_true(all(vapply(funcs, function(f) "sl_libs" %in% names(formals(get(f))), logical(1))))
})

test_that("basis defaults and formula prediction preserve training transformations", {
  d <- data.frame(X1=seq(-2,2,length.out=40), X2=rep(0:1,20), Y=rnorm(40))
  expect_equal(l2_basis(d,c("X1","X2"),NULL), l2_basis(d,c("X1","X2"),~ X1+X2))
  expect_equal(l2_basis(d,c("X1","X2"),~ .), l2_basis(d,c("X1","X2"),NULL))
  tr <- 1:30; te <- 31:40
  ff <- ~ splines::ns(X1,df=3)+X2+X1:X2
  got <- l2_basis_split(d,c("X1","X2"),ff,tr,te)
  mf <- model.frame(ff,d[tr,]); tt <- attr(mf,"terms")
  raw <- model.matrix(tt,mf); pred <- model.matrix(tt,model.frame(tt,d[te,]))
  for(j in 2:ncol(raw)) pred[,j] <- (pred[,j]-mean(raw[,j]))/sd(raw[,j])
  expect_equal(got$test,pred,ignore_attr=TRUE)
  expect_error(l2_basis(d,c("X1","X2"),~ Y),"covariates")
  expect_error(l2_basis(d,c("X1","X2"),~ X1-1),"intercept")
  chars <- data.frame(X1=1:6, X2=c("a","b","a","b","b","b"))
  cb <- l2_basis_split(chars,c("X1","X2"),NULL,1:4,5:6)
  expect_identical(colnames(cb$train),colnames(cb$test))
  expect_equal(nrow(cb$test),2)
  d$X1[1] <- NA
  expect_error(l2_basis(d,c("X1","X2"),NULL),"missing")
})

test_that("frontier plots switch intervals and do not bridge missing crossings", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(benchmark_id=1, label="Z",delta_M=0:4,estimate=c(4,3,NA,1,0),
    simultaneous_low=c(3,2,NA,0,0),simultaneous_high=c(5,4,NA,2,1),
    pointwise_low=c(3.5,2.5,NA,.5,0),pointwise_high=c(4.5,3.5,NA,1.5,.5))
  x <- structure(list(frontier=d,radius_x="delta_M",radius_y="delta_K",conf_level=.95),
    class="jointbounds_l2_calibration_band")
  expect_warning(p <- l2_plot_frontier(x),"between points")
  b <- ggplot2::ggplot_build(p)
  expect_equal(length(unique(b$data[[1]]$group)),2)
  expect_equal(length(unique(b$data[[2]]$group)),2)
  expect_equal(nrow(p$data),5)
  expect_no_warning(q <- plot(x,uncertainty="none",interpolate=FALSE))
  expect_length(q$layers,1)
  expect_no_warning(q <- plot(x,uncertainty="both",interpolate=FALSE))
  expect_length(q$layers,3)
  expect_warning(q <- l2_plot_frontier(x,uncertainty="none",facet=TRUE),"between points")
  expect_length(q$layers,2)
})
