test_that('simultaneous scope controls coverage across benchmarks only', {
  n <- 40
  a <- rep(c(-1,1),20); b <- rep(c(-1,-1,1,1),10)
  grid <- expand.grid(delta_M=c(0,.1),delta_K=c(0,.5,1))
  surface <- structure(list(grid=grid,bound_method='cs',fits=lapply(seq_len(nrow(grid)),function(i) {
    w <- sum(grid[i,])
    list(ate=data.frame(endpoint=c('lower','upper'),estimator='eif',estimate=c(1-w,1+w)),
      scores=list(eif_lower=.1*a+w*.1*b,eif_upper=.1*a-w*.1*b))
  })),class='jointbounds_l2_sensitivity_band')
  bench <- structure(list(results=data.frame(benchmark_id=c(7,12),label=c('a','b'),
    estimator='aipw',full_estimate=1,reduced_estimate=c(.7,1.2)),
    full_scores=list(aipw=.05*a), target_scores=list(`7_aipw`=.2*a,`12_aipw`=.2*b)),
    class='jointbounds_l2_benchmarks')
  joint <- l2_calibration_band(surface,bench,B=199)
  separate <- l2_calibration_band(surface,bench,B=199,simultaneous_scope='per_benchmark')
  expect_equal(joint$simultaneous_scope,'all_benchmarks')
  expect_equal(joint$suprema,apply(separate$suprema,1,max))
  expect_true(all(joint$critical_value >= separate$critical_value))
  expect_equal(names(separate$critical_value),c('7','12'))
  expect_equal(joint$frontier$estimate,separate$frontier$estimate)
  expect_equal(joint$frontier$pointwise_low,separate$frontier$pointwise_low)
  expect_equal(joint$frontier$pointwise_high,separate$frontier$pointwise_high)
  for (id in c(7,12)) {
    one <- bench
    one$results <- bench$results[bench$results$benchmark_id==id,,drop=FALSE]
    single <- l2_calibration_band(surface,one,B=199)
    expect_equal(unname(separate$critical_value[as.character(id)]),single$critical_value)
    expect_equal(separate$gaps$simultaneous_critical_value[separate$gaps$benchmark_id==id],
      rep(single$critical_value,nrow(single$gaps)))
  }
  expect_error(l2_calibration_band(surface,bench,simultaneous_scope='invalid'),'arg')
})
