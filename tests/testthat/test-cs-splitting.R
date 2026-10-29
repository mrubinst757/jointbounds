test_that("CS uses disjoint data for nuisance fitting, residual regression, and evaluation", {
  d <- simulate_l2_data(180, seed = 126)
  d <- d[, c('Y', 'A', 'C', 'X1', 'X2')]
  rownames(d) <- as.character(seq_len(nrow(d)))
  structural <- list(); moments <- list()
  original_split <- jointbounds:::l2_split_nuisance
  original_regression <- jointbounds:::l2_pseudo_regression
  local_mocked_bindings(
    l2_split_nuisance = function(data, Y, A, C, X, train, test, ...) {
      structural[[length(structural) + 1L]] <<- which(train)
      original_split(data, Y, A, C, X, train, test, ...)
    },
    l2_pseudo_regression = function(data, X, response, obs, newdata, ...) {
      moments[[length(moments) + 1L]] <<- list(
        reg = as.integer(rownames(data)), eval = as.integer(rownames(newdata)))
      original_regression(data, X, response, obs, newdata, ...)
    }, .package = 'jointbounds')
  z <- l2_cs_crossfit(d, 'Y', 'A', 'C', c('X1','X2'),
    folds = 3, nuisance_method = 'glm')
  expect_length(structural, 6)
  expect_length(moments, 18)
  for (j in seq_along(moments)) {
    a <- structural[[ceiling(j/3)]]; b <- moments[[j]]$reg; c <- moments[[j]]$eval
    expect_length(intersect(a,b), 0)
    expect_length(intersect(a,c), 0)
    expect_length(intersect(b,c), 0)
    expect_equal(sort(c(a,b,c)), seq_len(nrow(d)))
  }
  expect_match(z$parameters$sample_splitting, 'three-stage')
  expect_true(all(is.finite(z$ate$estimate)))
})

test_that("CS shared-training mode remains explicit and oracle errors ignore splitting", {
  d <- simulate_l2_data(120, seed=127)
  args <- list(data=d,Y='Y',A='A',C='C',X=c('X1','X2'),
    nuisance_method='glm',folds=3)
  old <- do.call(l2_cs_crossfit,c(args,list(moment_strategy='shared_training')))
  expect_match(old$parameters$sample_splitting,'shared-training')
  args$nuisance_source <- 'oracle_error'
  a <- do.call(l2_cs_crossfit,args)
  b <- do.call(l2_cs_crossfit,c(args,list(moment_strategy='shared_training')))
  expect_equal(a$ate,b$ate)
  expect_equal(a$scores,b$scores)
})
