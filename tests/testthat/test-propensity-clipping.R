test_that('treatment clipping preserves complements and other nuisance estimates', {
  d <- simulate_l2_data(120, seed=914, violation='both')
  all <- rep(TRUE,nrow(d))
  args <- list(d,c('X1','X2'),all,all,.1,c(e=3.5,rho=.5,mu=.5,regression=.5),c(e=-1,rho=1,mu=1,regression=-1))
  old <- do.call(jointbounds:::l2_split_nuisance_oracle_error,args)
  new <- do.call(jointbounds:::l2_split_nuisance_oracle_error,c(args,list(propensity_clip=.01)))
  for(part in c('train','test')) {
    expect_equal(new[[part]]$e1,pmin(pmax(old[[part]]$e1,.01),.99))
    expect_equal(new[[part]]$e0,1-new[[part]]$e1)
    keep <- setdiff(names(old[[part]]),c('e0','e1'))
    expect_equal(new[[part]][keep],old[[part]][keep])
  }
  expect_true(any(abs(new$test$e1-old$test$e1)>0))
})
test_that('public estimators validate clipping and leave oracle scores unchanged', {
  for(bad in list(NA,0,.5,Inf,c(.01,.02),'0.01')) {
    expect_error(l2_sharp_crossfit(NULL,'Y','A','C','X',propensity_clip=bad),'propensity_clip')
    expect_error(l2_cs_crossfit(NULL,'Y','A','C','X',propensity_clip=bad),'propensity_clip')
  }
  d<-simulate_l2_data(100,seed=916,violation='both')
  args<-list(data=d,Y='Y',A='A',C='C',X=c('X1','X2'),nuisance_source='oracle_error',nuisance_error_rate=.1,nuisance_error_scale=c(e=3.5,rho=.5,mu=.5,regression=.5))
  for(fun in list(l2_sharp_crossfit,l2_cs_crossfit)) {
    old<-do.call(fun,args); new<-do.call(fun,c(args,list(propensity_clip=.01)))
    expect_equal(new$parameters$propensity_clip,.01)
    expect_equal(subset(old$ate,estimator=='true_eif'),subset(new$ate,estimator=='true_eif'))
    expect_true(all(is.finite(new$ate$estimate)))
  }
})
