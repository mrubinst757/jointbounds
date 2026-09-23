l2_quiet_superlearner <- function(...) {
  # SuperLearner defaults to parent.frame() when resolving `All` and learner
  # wrappers.  A namespace-qualified call from jointbounds would otherwise look
  # in the jointbounds namespace and fail with "object 'All' not found".
  suppressWarnings(SuperLearner::SuperLearner(
    ..., verbose = FALSE, env = asNamespace("SuperLearner")
  ))
}

#' Estimate nuisance functions using SuperLearner with cross-fitting
#'
#' Estimates e(X)=P(A=1|X), pi_0(X)=P(C=1|X,A=0), pi_1(X)=P(C=1|X,A=1),
#' mu_0(X)=E(Y|X,A=0,C=0), mu_1(X)=E(Y|X,A=1,C=0). Uses V-fold cross-fitting:
#' for each fold, fits nuisances on the training part and predicts on the held-out part.
#'
#' @param X Matrix of covariates.
#' @param A Treatment indicator (0/1).
#' @param C Missingness indicator (1=missing, 0=observed).
#' @param Y Outcome (NA or value when C=0); only used when C=0 for outcome models.
#' @param V Number of cross-fitting folds (default 2; the calling function may use
#'   \code{V=1} when a simple single-library SuperLearner such as \code{SL.glm} is used,
#'   in which case nuisances are fitted and predicted on the full sample).
#' @param sl_lib_prop Character vector of SuperLearner library for propensity e (default "SL.glm").
#' @param sl_lib_miss Character vector of SuperLearner library for missingness pi_0, pi_1 (default "SL.glm").
#' @param sl_lib_outcome Character vector of SuperLearner library for outcome mu_0, mu_1 (default "SL.glm").
#' @param stratify_mu Logical; if TRUE (default), estimate separate outcome models mu_0 and mu_1 stratified by treatment A. If FALSE, estimate a single pooled model E(Y|X,A,C=0) with A as a covariate and predict at A=0 and A=1.
#' @param family_Y Character; family for outcome model ("gaussian" or "binomial"). Default "gaussian".
#' @param seed Optional seed for fold splits.
#' @return List with components: e, pi0, pi1, mu0, mu1 (each length n), and fold_id (fold index per row).
#' @export
estimate_nuisance <- function(X, A, C, Y,
                              V = 2,
                              sl_lib_prop = "SL.glm",
                              sl_lib_miss = "SL.glm",
                              sl_lib_outcome = "SL.glm",
                              stratify_mu = TRUE,
                              family_Y = "gaussian",
                              seed = NULL) {
  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    stop("Package 'SuperLearner' is required. Install with install.packages('SuperLearner').")
  }
  n <- length(A)
  if (!is.matrix(X)) X <- as.matrix(X)
  family_Y_object <- switch(family_Y,
    gaussian = stats::gaussian(),
    binomial = stats::binomial(),
    stop("family_Y must be 'gaussian' or 'binomial'")
  )
  V <- min(V, max(1L, n))
  if (is.null(seed)) seed <- 1L
  set.seed(seed)
  fold_id <- sample(rep(seq_len(V), length.out = n))

  e <- numeric(n)
  pi0 <- numeric(n)
  pi1 <- numeric(n)
  mu0 <- numeric(n)
  mu1 <- numeric(n)
  successful_sl_fits <- 0L
  sl_errors <- character()

  for (v in seq_len(V)) {
    # V = 1 means no sample splitting: fit and predict on the full sample.
    train <- if (V == 1L) rep(TRUE, n) else fold_id != v
    eval <- fold_id == v
    n_train <- sum(train)
    Xtrain <- X[train, , drop = FALSE]
    Xeval <- X[eval, , drop = FALSE]
    # Propensity e(X) = P(A=1|X) — skip SuperLearner if training set empty or too small
    if (n_train < 2L || length(unique(A[train])) < 2L) {
      e[eval] <- clip_probs(mean(A))
    } else {
      fit_e <- tryCatch(
        {
          l2_quiet_superlearner(
              Y = A[train],
              X = as.data.frame(Xtrain),
              newX = as.data.frame(Xeval),
              family = stats::binomial(),
              SL.library = sl_lib_prop
          )
        },
        error = function(e) {
          sl_errors <<- c(sl_errors, conditionMessage(e)); NULL
        }
      )
      if (!is.null(fit_e)) successful_sl_fits <- successful_sl_fits + 1L
      if (is.null(fit_e)) {
        e[eval] <- clip_probs(mean(A[train]))
      } else {
        e[eval] <- clip_probs(fit_e$SL.predict)
      }
    }

    # Missingness: fit P(C=1|X,A) by fitting separately in A=0 and A=1
    for (a in 0:1) {
      ia <- train & (A == a)
      n_ia <- sum(ia)
      if (n_ia < 2L || length(unique(C[ia])) < 2L) {
        if (a == 0) pi0[eval] <- clip_probs(mean(C[A == 0])) else pi1[eval] <- clip_probs(mean(C[A == 1]))
        next
      }
      fit_pi <- tryCatch(
        {
          l2_quiet_superlearner(
              Y = C[ia],
              X = as.data.frame(X[ia, , drop = FALSE]),
              newX = as.data.frame(Xeval),
              family = stats::binomial(),
              SL.library = sl_lib_miss
          )
        },
        error = function(e) {
          sl_errors <<- c(sl_errors, conditionMessage(e)); NULL
        }
      )
      if (!is.null(fit_pi)) successful_sl_fits <- successful_sl_fits + 1L
      if (is.null(fit_pi)) {
        if (a == 0) pi0[eval] <- clip_probs(mean(C[ia])) else pi1[eval] <- clip_probs(mean(C[ia]))
      } else {
        if (a == 0) pi0[eval] <- clip_probs(fit_pi$SL.predict) else pi1[eval] <- clip_probs(fit_pi$SL.predict)
      }
    }

    # Outcome: E(Y|X,A=a,C=0) for a=0,1
    if (stratify_mu) {
      # Stratified: separate models for A=0 and A=1
      for (a in 0:1) {
        ia <- train & (A == a) & (C == 0)
        n_ia <- sum(ia)
        Y_ia <- Y[ia]
        # Check for sufficient data and variation
        has_variation <- (family_Y == "gaussian") || (length(unique(Y_ia[!is.na(Y_ia)])) >= 2L)
        if (n_ia < 2L || !has_variation) {
          if (a == 0) mu0[eval] <- mean(Y[A == 0 & C == 0], na.rm = TRUE) else mu1[eval] <- mean(Y[A == 1 & C == 0], na.rm = TRUE)
          next
        }
        fit_mu <- tryCatch(
          {
            l2_quiet_superlearner(
                Y = Y_ia,
                X = as.data.frame(X[ia, , drop = FALSE]),
                newX = as.data.frame(Xeval),
                family = family_Y_object,
                SL.library = sl_lib_outcome
            )
          },
          error = function(e) {
            sl_errors <<- c(sl_errors, conditionMessage(e)); NULL
          }
        )
        if (!is.null(fit_mu)) successful_sl_fits <- successful_sl_fits + 1L
        if (is.null(fit_mu)) {
          if (a == 0) mu0[eval] <- mean(Y_ia, na.rm = TRUE) else mu1[eval] <- mean(Y_ia, na.rm = TRUE)
        } else {
          pred <- fit_mu$SL.predict
          if (a == 0) mu0[eval] <- pred else mu1[eval] <- pred
        }
      }
    } else {
      # Pooled: single model E(Y|X, C=0)
      io <- train & (C == 0)
      n_io <- sum(io)
      Y_io <- Y[io]
      has_variation <- (family_Y == "gaussian") || (length(unique(Y_io[!is.na(Y_io)])) >= 2L)
      if (n_io < 2L || !has_variation) {
        mu_pooled <- mean(Y[C == 0], na.rm = TRUE)
        mu0[eval] <- mu_pooled
        mu1[eval] <- mu_pooled
      } else {
        # A enters the pooled design so that mu0 and mu1 can differ; predict
        # both arms by setting A = 0 and A = 1 on the evaluation fold.
        XA_train <- data.frame(X[io, , drop = FALSE], .A = A[io])
        n_eval <- sum(eval)
        XA_eval <- rbind(
          data.frame(Xeval, .A = rep(0, n_eval)),
          data.frame(Xeval, .A = rep(1, n_eval))
        )
        fit_mu_pooled <- tryCatch(
          {
            l2_quiet_superlearner(
                Y = Y_io,
                X = XA_train,
                newX = XA_eval,
                family = family_Y_object,
                SL.library = sl_lib_outcome
            )
          },
          error = function(e) {
            sl_errors <<- c(sl_errors, conditionMessage(e)); NULL
          }
        )
        if (!is.null(fit_mu_pooled)) successful_sl_fits <- successful_sl_fits + 1L
        if (is.null(fit_mu_pooled)) {
          mu_all <- mean(Y_io, na.rm = TRUE)
          mu0[eval] <- if (any(A[io] == 0)) mean(Y_io[A[io] == 0], na.rm = TRUE) else mu_all
          mu1[eval] <- if (any(A[io] == 1)) mean(Y_io[A[io] == 1], na.rm = TRUE) else mu_all
        } else {
          pred <- as.numeric(fit_mu_pooled$SL.predict)
          mu0[eval] <- pred[seq_len(n_eval)]
          mu1[eval] <- pred[n_eval + seq_len(n_eval)]
        }
      }
    }
  }

  if (successful_sl_fits == 0L && ncol(X) > 0L) {
    detail <- if (length(sl_errors)) paste(unique(sl_errors), collapse = "; ") else
      "no error message was returned"
    warning("All SuperLearner fits failed; nuisance estimates used marginal fallbacks. ",
            "Underlying error: ", detail)
  }
  list(e = e, pi0 = pi0, pi1 = pi1, mu0 = mu0, mu1 = mu1, fold_id = fold_id)
}
