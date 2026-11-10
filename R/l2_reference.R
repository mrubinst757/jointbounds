# Structural predictions for a zero-budget AIPW reference. Each fold is
# fitted once so both signed endpoints and both arms use the same predictions.
l2_reference_splits <- function(data, Y, A, C, X, id, strategy, method,
                                lib_e, lib_r, lib_y, clip, seed) {
  set.seed(seed)
  V <- max(id)
  lapply(seq_len(V), function(v) {
    tr <- id != v
    if (strategy == "three_way") tr <- tr & id != (if (v == 1L) V else v - 1L)
    z <- l2_split_nuisance(data, Y, A, C, X, tr, !tr, method,
                           lib_e, lib_r, lib_y, propensity_clip = clip)
    z$evaluation <- lapply(z$test, `[`, id[!tr] == v)
    z
  })
}

l2_reference_signature <- function(data, Y, A, C, X) {
  digest::digest(list(Y, A, C, X, rownames(data),
    lapply(c(Y,A,C,X), function(nm) data[[nm]])), algo = "xxhash64")
}

l2_reference_scores <- function(data, Y, A, C, id, splits) {
  n <- nrow(data)
  plugin <- aipw <- numeric(n)
  for (v in seq_along(splits)) {
    i <- which(id == v); z <- splits[[v]]$evaluation
    for (a in 0:1) {
      mu <- z[[paste0("mu",a)]]
      S <- data[[A]][i] == a & data[[C]][i] == 0
      aug <- numeric(length(i))
      aug[S] <- (data[[Y]][i[S]] - mu[S]) /
        (z[[paste0("e",a)]][S] * z[[paste0("rho",a)]][S])
      sg <- 2*a-1
      plugin[i] <- plugin[i] + sg*mu
      aipw[i] <- aipw[i] + sg*(mu+aug)
    }
  }
  scores <- list(plugin=plugin,aipw=aipw)
  list(scores=scores, estimates=vapply(scores,mean,numeric(1)))
}

l2_get_reference <- function(fit, signature) {
  if (inherits(fit,"jointbounds_l2_sensitivity_band")) {
    ii <- which(vapply(fit$fits, function(z) !is.null(z$benchmark_reference), logical(1)))
    if (length(ii) != 1L) stop("Expected one retained zero-radius fit; use keep_fits=TRUE and include (0,0)")
    fit <- fit$fits[[ii]]
  }
  ref <- fit$benchmark_reference
  if (is.null(ref)) stop("A zero-radius estimated-nuisance fit made with jointbounds >= 0.4.2 is required")
  if (!identical(ref$signature, signature)) stop("Reference data, column names or observation order do not match")
  ref
}
