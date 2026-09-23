## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)


## -----------------------------------------------------------------------------
library(jointbounds)
suppressPackageStartupMessages(library(SuperLearner))


## ---- eval = FALSE------------------------------------------------------------
## # See all available algorithms
## listWrappers()


## -----------------------------------------------------------------------------
# Example data
set.seed(20260313)
n <- 500
X <- runif(n, -2, 2)
A <- rbinom(n, 1, plogis(X))
C <- rbinom(n, 1, 0.2 + 0.1 * A)
Y <- ifelse(C == 1, NA, rbinom(n, 1, plogis(0.5 * A + 0.3 * X)))
dat <- data.frame(Y = Y, A = A, C = C, X = X)

# Simple ensemble
fit_ensemble <- mar_bounds(
  dat,
  Y = "Y",
  A = "A",
  C = "C",
  X = "X",
  estimand = "ate",
  assumption = "general",
  sl_lib = c("SL.glm", "SL.gam", "SL.mean"),
  V = 2,
  seed = 1
)

# Check SuperLearner weights (stored in nuisance object)
fit_ensemble$nuisance$SL_weights_prop  # Propensity score weights


## ---- eval = FALSE------------------------------------------------------------
## fit_custom <- mar_bounds(
##   dat,
##   Y = "Y",
##   A = "A",
##   C = "C",
##   X = "X",
##   estimand = "ate",
##   assumption = "general",
##   sl_lib_prop = c("SL.glm", "SL.gam"),           # Propensity score
##   sl_lib_miss = c("SL.glm", "SL.ranger"),        # Missingness mechanism
##   sl_lib_outcome = c("SL.glm", "SL.xgboost"),    # Outcome regression
##   V = 5,  # More folds for stable SuperLearner
##   seed = 1
## )


## -----------------------------------------------------------------------------
# V = 2 (default): Split data in half
fit_v2 <- mar_bounds(
  dat, Y = "Y", A = "A", C = "C", X = "X",
  estimand = "ate", assumption = "general",
  V = 2, sl_lib = "SL.glm", seed = 1
)

# V = 5: More folds (better for larger samples)
fit_v5 <- mar_bounds(
  dat, Y = "Y", A = "A", C = "C", X = "X",
  estimand = "ate", assumption = "general",
  V = 5, sl_lib = "SL.glm", seed = 1
)

# V = 1: No cross-fitting (not recommended for inference)
fit_v1 <- mar_bounds(
  dat, Y = "Y", A = "A", C = "C", X = "X",
  estimand = "ate", assumption = "general",
  V = 1, sl_lib = "SL.glm", seed = 1
)


## -----------------------------------------------------------------------------
# Check fold assignments
table(fit_v2$nuisance$fold_id)

# Nuisance estimates are out-of-sample within each fold
head(fit_v2$nuisance$fold_id)


## -----------------------------------------------------------------------------
phi <- fit_v2$phi
dim(phi)  # n x 6 matrix
colnames(phi)


## -----------------------------------------------------------------------------
# Standard error calculation
se_custom <- function(phi_col) {
  sqrt(var(phi_col) / length(phi_col))
}

# Custom confidence interval
estimate <- fit_v2$lower
se <- se_custom(phi[, 2] + phi[, 3])  # Combine relevant components
ci_lower <- estimate - 1.96 * se
ci_upper <- estimate + 1.96 * se

cat("Custom 95% CI for lower bound: [", ci_lower, ",", ci_upper, "]\n")


## ---- eval = FALSE------------------------------------------------------------
## # Extract fitted object components
## phi <- fit_v2$phi
## nuisance <- fit_v2$nuisance
## 
## # Compute bounds at different δ values using stored influence functions
## # (This is what multiplier_bootstrap_grid does internally)
## new_bounds <- jointbounds:::compute_bounds(
##   phi,
##   estimand = "ate",
##   assumption = "general",
##   delta_0 = 0.6,
##   delta_1 = 0.6
## )


## -----------------------------------------------------------------------------
# Stratified (default): E[Y|A=a,X] fit separately for a=0 and a=1
fit_stratified <- mar_bounds(
  dat, Y = "Y", A = "A", C = "C", X = "X",
  estimand = "ate", assumption = "general",
  stratify_mu = TRUE,
  V = 2, sl_lib = "SL.glm", seed = 1
)

# Pooled: E[Y|A,X] fit once with A as covariate
fit_pooled <- mar_bounds(
  dat, Y = "Y", A = "A", C = "C", X = "X",
  estimand = "ate", assumption = "general",
  stratify_mu = FALSE,
  V = 2, sl_lib = "SL.glm", seed = 1
)


## -----------------------------------------------------------------------------
fit_general <- mar_bounds(
  dat, Y = "Y", A = "A", C = "C", X = "X",
  estimand = "ate",
  assumption = "general",
  V = 2, sl_lib = "SL.glm", seed = 1
)


## -----------------------------------------------------------------------------
fit_bounded_delta <- mar_bounds(
  dat, Y = "Y", A = "A", C = "C", X = "X",
  estimand = "psi1",
  assumption = "bounded_delta",
  delta_0u = 0.8,  # Upper bound on proportion
  delta_1u = 0.8,
  V = 2, sl_lib = "SL.glm", seed = 1
)


## -----------------------------------------------------------------------------
# Positive monotonicity: P(C=1|Y=1,A,X) >= P(C=1|Y=0,A,X)
fit_mono_pos <- mar_bounds(
  dat, Y = "Y", A = "A", C = "C", X = "X",
  estimand = "ate",
  assumption = "monotonicity_pos",
  delta_0u = 0.8,
  delta_1u = 0.8,
  V = 2, sl_lib = "SL.glm", seed = 1
)

# Negative monotonicity: reverse inequality
fit_mono_neg <- mar_bounds(
  dat, Y = "Y", A = "A", C = "C", X = "X",
  estimand = "ate",
  assumption = "monotonicity_neg",
  delta_0u = 0.8,
  delta_1u = 0.8,
  V = 2, sl_lib = "SL.glm", seed = 1
)


## -----------------------------------------------------------------------------
fit_bounded_risk <- mar_bounds(
  dat, Y = "Y", A = "A", C = "C", X = "X",
  estimand = "ate",
  assumption = "bounded_risk",
  delta_0u = 0.8,
  delta_1u = 0.8,
  tau_0 = 2,  # Risk ratio bound in control
  tau_1 = 2,  # Risk ratio bound in treatment
  V = 2, sl_lib = "SL.glm", seed = 1
)


## -----------------------------------------------------------------------------
fit_point <- mar_bounds(
  dat, Y = "Y", A = "A", C = "C", X = "X",
  estimand = "ate",
  assumption = "point_ate",
  delta_0 = 0.5,  # Known value
  delta_1 = 0.5,
  tau = 2,
  V = 2, sl_lib = "SL.glm", seed = 1
)


## ---- eval = FALSE------------------------------------------------------------
## # Rademacher (default): Random ±1
## fit_rademacher <- mar_bounds(
##   dat, Y = "Y", A = "A", C = "C", X = "X",
##   param_grid = list(delta_0 = c(0.5, 0.7, 0.9),
##                     delta_1 = c(0.5, 0.7, 0.9)),
##   multiplier = "rademacher",
##   B = 500, seed = 1
## )
## 
## # Gaussian: Random N(0,1)
## fit_gaussian <- mar_bounds(
##   dat, Y = "Y", A = "A", C = "C", X = "X",
##   param_grid = list(delta_0 = c(0.5, 0.7, 0.9),
##                     delta_1 = c(0.5, 0.7, 0.9)),
##   multiplier = "gaussian",
##   B = 500, seed = 1
## )


## ---- eval = FALSE------------------------------------------------------------
## # Quick exploration
## fit_quick <- mar_bounds(..., B = 100)
## 
## # Standard analysis
## fit_standard <- mar_bounds(..., B = 500)
## 
## # Publication-quality
## fit_publication <- mar_bounds(..., B = 1000)


## ---- eval = FALSE------------------------------------------------------------
## fit_small <- mar_bounds(
##   dat,
##   Y = "Y", A = "A", C = "C", X = "X",
##   sl_lib = "SL.glm",  # Parametric only
##   V = 2,              # Minimal cross-fitting
##   B = 200,            # Fewer bootstrap replicates
##   seed = 1
## )


## ---- eval = FALSE------------------------------------------------------------
## fit_large <- mar_bounds(
##   dat,
##   Y = "Y", A = "A", C = "C", X = "X",
##   sl_lib = c("SL.glm", "SL.gam", "SL.ranger"),  # Can use complex models
##   V = 5,              # More folds for stability
##   B = 1000,           # More bootstrap replicates
##   seed = 1
## )


## ---- eval = FALSE------------------------------------------------------------
## # Only compute one bound
## fit_lower <- mar_bounds(..., grid_bounds = "lower")
## fit_upper <- mar_bounds(..., grid_bounds = "upper")
## 
## # Or split into smaller grids
## grid1 <- list(delta_0 = seq(0.3, 0.6, 0.1), delta_1 = seq(0.3, 0.6, 0.1))
## grid2 <- list(delta_0 = seq(0.7, 1.0, 0.1), delta_1 = seq(0.7, 1.0, 0.1))


## ---- eval = FALSE------------------------------------------------------------
## # Solution 1: Use simpler library
## mar_bounds(..., sl_lib = "SL.glm")
## 
## # Solution 2: Increase V to avoid small training samples
## mar_bounds(..., V = 2)  # instead of V = 5
## 
## # Solution 3: Add fallback algorithms
## mar_bounds(..., sl_lib = c("SL.glm", "SL.mean"))


## ---- eval = FALSE------------------------------------------------------------
## # Correct coding
## dat$Y[dat$C == 1] <- NA
## 
## # Check consistency
## stopifnot(all(is.na(dat$Y) == (dat$C == 1)))


## ---- eval = FALSE------------------------------------------------------------
## fit <- mar_bounds(..., seed = 20260313)  # YYYYMMDD format recommended

