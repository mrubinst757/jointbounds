## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)


## ---- eval = FALSE------------------------------------------------------------
## # Install from source
## install.packages("path/to/jointbounds", repos = NULL, type = "source")


## -----------------------------------------------------------------------------
library(jointbounds)
suppressPackageStartupMessages(library(SuperLearner))


## -----------------------------------------------------------------------------
set.seed(20260313)
n <- 500

# Covariates
X <- runif(n, -2, 2)

# Treatment (depends on X)
A <- rbinom(n, 1, plogis(X))

# Missingness (depends on treatment)
C <- rbinom(n, 1, 0.2 + 0.1 * A)

# Outcome (only observed when C = 0)
Y <- ifelse(C == 1, NA, rbinom(n, 1, plogis(0.5 * A + 0.3 * X)))

# Create data frame
dat <- data.frame(Y = Y, A = A, C = C, X = X)

# Summary
summary(dat)
cat("\nMissing outcomes:", sum(is.na(dat$Y)), "/", n, "\n")


## -----------------------------------------------------------------------------
fit_general <- mar_bounds(
  dat,
  Y = "Y",
  A = "A",
  C = "C",
  X = "X",
  estimand = "ate",
  assumption = "general",
  sl_lib = "SL.glm",
  V = 2,
  seed = 1
)

# View results
fit_general$result


## -----------------------------------------------------------------------------
fit_bounded_delta <- mar_bounds(
  dat,
  Y = "Y",
  A = "A",
  C = "C",
  X = "X",
  estimand = "ate",
  assumption = "general",
  delta_0 = 0.8,  # At most 80% informative in control group
  delta_1 = 0.8,  # At most 80% informative in treatment group
  sl_lib = "SL.glm",
  V = 2,
  seed = 1
)

fit_bounded_delta$result


## -----------------------------------------------------------------------------
fit_point <- mar_bounds(
  dat,
  Y = "Y",
  A = "A",
  C = "C",
  X = "X",
  estimand = "ate",
  assumption = "point_ate",
  delta_0 = 0.5,
  delta_1 = 0.5,
  tau = 2,
  sl_lib = "SL.glm",
  V = 2,
  seed = 1
)

fit_point$result


## -----------------------------------------------------------------------------
results <- data.frame(
  assumption = c("General", "Delta=0.8", "Delta=0.5", "Point ID"),
  lower = c(
    fit_general$result$estimate[1],
    fit_bounded_delta$result$estimate[1],
    NA,  # Will compute below
    NA   # Point estimate doesn't have bounds
  ),
  upper = c(
    fit_general$result$estimate[2],
    fit_bounded_delta$result$estimate[2],
    NA,
    NA
  ),
  point = c(NA, NA, NA, fit_point$result$estimate)
)

# Add delta=0.5 bounds
fit_delta_05 <- mar_bounds(
  dat, Y = "Y", A = "A", C = "C", X = "X",
  estimand = "ate", assumption = "general",
  delta_0 = 0.5, delta_1 = 0.5,
  sl_lib = "SL.glm", V = 2, seed = 1
)
results$lower[3] <- fit_delta_05$result$estimate[1]
results$upper[3] <- fit_delta_05$result$estimate[2]

print(results)


## -----------------------------------------------------------------------------
fit_psi1 <- mar_bounds(
  dat,
  Y = "Y",
  A = "A",
  C = "C",
  X = "X",
  estimand = "psi1",
  assumption = "bounded_delta",
  delta_0u = 0.8,
  delta_1u = 0.8,
  sl_lib = "SL.glm",
  V = 2,
  seed = 1
)

fit_psi1$result


## -----------------------------------------------------------------------------
fit_psi2 <- mar_bounds(
  dat,
  Y = "Y",
  A = "A",
  C = "C",
  X = "X",
  estimand = "psi2",
  assumption = "point_psi2",
  delta_0 = 0.5,
  sl_lib = "SL.glm",
  V = 2,
  seed = 1
)

fit_psi2$result


## ---- eval = FALSE------------------------------------------------------------
## fit_custom <- mar_bounds(
##   dat,
##   Y = "Y",
##   A = "A",
##   C = "C",
##   X = "X",
##   estimand = "ate",
##   assumption = "general",
##   sl_lib_prop = c("SL.glm", "SL.mean"),      # Propensity score
##   sl_lib_miss = c("SL.glm", "SL.gam"),       # Missingness probabilities
##   sl_lib_outcome = c("SL.glm", "SL.ranger"), # Outcome regression
##   V = 2,
##   seed = 1
## )

