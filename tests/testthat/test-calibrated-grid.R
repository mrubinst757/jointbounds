test_that('calibrated grids select confidence limits and preserve scope', {
  z <- expand.grid(endpoint=c('lower','upper'),delta_K=c(0,.2,.4,.6),
    delta_M=c(0,.1),benchmark_id=c(2,7), stringsAsFactors=FALSE)
  z$label <- paste('omit',z$benchmark_id)
  z$gap <- ifelse(z$endpoint=='lower',.2-z$delta_K,-.2-z$delta_K)
  z$pointwise_high <- z$gap+.1
  z$simultaneous_high <- z$gap+.3
  band <- structure(list(gaps=z,radius_x='delta_M',radius_y='delta_K',
    conf_level=.95, simultaneous_scope='all_benchmarks'),class='jointbounds_l2_calibration_band')
  p <- l2_calibrated_grid(band,selection='estimate')
  u <- l2_calibrated_grid(band,selection='pointwise_upper')
  s <- l2_calibrated_grid(band,selection='simultaneous_upper')
  expect_equal(p$delta_K,rep(.2,4))
  expect_equal(u$delta_K,rep(.4,4))
  expect_equal(s$delta_K,rep(.6,4))
  expect_equal(attr(s,'simultaneous_scope'),'all_benchmarks')
  one <- l2_calibrated_grid(band,'simultaneous_upper',radius_x=.1,benchmark_id=7)
  expect_equal(nrow(one),1L)
  expect_equal(one$delta_M,.1)
  expect_equal(one$benchmark_id,7)
  expect_error(l2_calibrated_grid(band,radius_x=.15),'already evaluated')
  expect_error(l2_calibrated_grid(band,benchmark_id=8),'Unknown')
  band$gaps$simultaneous_high <- band$gaps$simultaneous_high+1
  expect_error(l2_calibrated_grid(band,'simultaneous_upper'),'extend the sensitivity grid')
  band$gaps <- subset(band$gaps,endpoint=='lower')
  expect_error(l2_calibrated_grid(band),'both endpoints')
})

test_that('upper calibration uses both inequalities rather than target direction', {
  z <- expand.grid(endpoint=c('lower','upper'),delta_K=c(0,.25,.5))
  z$benchmark_id <- 1; z$label <- 'omit'; z$delta_M <- 0
  z$gap <- ifelse(z$endpoint=='lower',-.2-z$delta_K,.2-z$delta_K)
  z$pointwise_high <- z$gap+.1
  z$simultaneous_high <- z$gap+.1
  x <- structure(list(gaps=z,radius_x='delta_M',radius_y='delta_K'),class='jointbounds_l2_calibration_band')
  expect_equal(l2_calibrated_grid(x)$delta_K,.25)
  expect_equal(l2_calibrated_grid(x,'pointwise_upper')$delta_K,.5)
})

test_that('selected calibration grid can be used to compute bounds', {
  d <- simulate_l2_data(150,seed=913,violation='both',density_ratio_link='bounded')
  surface <- l2_sensitivity_band(d,'Y','A','C',c('X1','X2'),
    grid=data.frame(delta_M=0,delta_K=c(0,.5,1)),model='net',bound_method='cs',
    nuisance_method='glm',nuisance_source='oracle_error',keep_fits=TRUE,B=19,seed=4)
  f <- surface$fits[[1]]
  ref <- f$ate$estimate[f$ate$estimator=='eif' & f$ate$endpoint=='lower']
  score <- f$scores$eif_lower
  bench <- structure(list(results=data.frame(benchmark_id=1,label='omit',
    estimator='aipw',full_estimate=ref,reduced_estimate=ref+.1),
    full_scores=list(aipw=score),target_scores=list(`1_aipw`=score*.9)),
    class='jointbounds_l2_benchmarks')
  band <- l2_calibration_band(surface,bench,calibration='signed',B=19)
  grid <- l2_calibrated_grid(band,selection='simultaneous_upper',radius_x=0)
  expect_equal(nrow(grid),1L)
  expect_gte(grid$delta_K,l2_calibrated_grid(band)$delta_K)
  fit <- l2_sensitivity_band(d,'Y','A','C',c('X1','X2'),grid=grid,
    model='net',bound_method='cs',nuisance_method='glm',nuisance_source='oracle_error',B=19,seed=4)
  expect_true(all(is.finite(fit$results$estimate)))
  expect_equal(fit$grid$delta_K,grid$delta_K)
})
