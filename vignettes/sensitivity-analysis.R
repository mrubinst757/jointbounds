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


## -----------------------------------------------------------------------------
set.seed(20260313)
n <- 1000

X <- runif(n, -2, 2)
A <- rbinom(n, 1, plogis(X))
C <- rbinom(n, 1, 0.2 + 0.1 * A)
Y <- ifelse(C == 1, NA, rbinom(n, 1, plogis(0.5 * A + 0.3 * X)))

dat <- data.frame(Y = Y, A = A, C = C, X = X)


## -----------------------------------------------------------------------------
fit_delta_grid <- mar_bounds(
  dat,
  Y = "Y",
  A = "A",
  C = "C",
  X = "X",
  estimand = "ate",
  assumption = "general",
  param_grid = list(
    delta_0 = seq(0.3, 1.0, by = 0.1),
    delta_1 = seq(0.3, 1.0, by = 0.1)
  ),
  sl_lib = "SL.glm",
  V = 2,
  alpha = 0.05,
  B = 500,  # Number of bootstrap replicates
  seed = 1
)


## -----------------------------------------------------------------------------
# Lower bound grid
head(fit_delta_grid$lower_grid$grid)

# Upper bound grid
head(fit_delta_grid$upper_grid$grid)


## -----------------------------------------------------------------------------
cat("Lower bound critical value:", fit_delta_grid$lower_grid$critical_value, "\n")
cat("Upper bound critical value:", fit_delta_grid$upper_grid$critical_value, "\n")
cat("Standard normal quantile:", qnorm(0.975), "\n")


## -----------------------------------------------------------------------------
library(ggplot2)

grid <- fit_delta_grid$lower_grid$grid

# When both deltas are equal
diag_subset <- grid[grid$delta_0 == grid$delta_1, ]

ggplot(diag_subset, aes(x = delta_0)) +
  geom_line(aes(y = estimate), linewidth = 1) +
  geom_ribbon(aes(ymin = ci_lower_pointwise, ymax = ci_upper_pointwise),
              alpha = 0.2, fill = "blue") +
  geom_ribbon(aes(ymin = ci_lower_uniform, ymax = ci_upper_uniform),
              alpha = 0.2, fill = "red") +
  labs(
    title = "Lower Bound Sensitivity Analysis",
    subtitle = "Blue = pointwise CIs, Red = simultaneous CIs",
    x = "Proportion of informative missingness (δ)",
    y = "Lower bound on ATE"
  ) +
  theme_minimal()


## -----------------------------------------------------------------------------
ggplot(grid, aes(x = delta_0, y = delta_1, fill = estimate)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Lower\nBound") +
  labs(
    title = "Lower Bound Heatmap",
    x = "δ₀ (control group)",
    y = "δ₁ (treatment group)"
  ) +
  theme_minimal() +
  coord_fixed()


## -----------------------------------------------------------------------------
upper_grid <- fit_delta_grid$upper_grid$grid
upper_diag <- upper_grid[upper_grid$delta_0 == upper_grid$delta_1, ]

ggplot(upper_diag, aes(x = delta_0)) +
  geom_line(aes(y = estimate), linewidth = 1) +
  geom_ribbon(aes(ymin = ci_lower_pointwise, ymax = ci_upper_pointwise),
              alpha = 0.2, fill = "blue") +
  geom_ribbon(aes(ymin = ci_lower_uniform, ymax = ci_upper_uniform),
              alpha = 0.2, fill = "red") +
  labs(
    title = "Upper Bound Sensitivity Analysis",
    subtitle = "Blue = pointwise CIs, Red = simultaneous CIs",
    x = "Proportion of informative missingness (δ)",
    y = "Upper bound on ATE"
  ) +
  theme_minimal()


## ---- eval = FALSE------------------------------------------------------------
## fit_risk_grid <- mar_bounds(
##   dat,
##   Y = "Y",
##   A = "A",
##   C = "C",
##   X = "X",
##   estimand = "ate",
##   assumption = "bounded_risk",
##   param_grid = list(
##     delta_0u = seq(0.5, 1.0, by = 0.1),
##     delta_1u = seq(0.5, 1.0, by = 0.1),
##     tau_0 = c(1.5, 2.0, 2.5),
##     tau_1 = c(1.5, 2.0, 2.5)
##   ),
##   sl_lib = "SL.glm",
##   V = 2,
##   B = 500,
##   seed = 1
## )
## 
## # This creates a 6 x 6 x 3 x 3 = 324 point grid
## nrow(fit_risk_grid$lower_grid$grid)


## -----------------------------------------------------------------------------
# Only lower bound
fit_lower_only <- mar_bounds(
  dat,
  Y = "Y",
  A = "A",
  C = "C",
  X = "X",
  estimand = "ate",
  assumption = "general",
  param_grid = list(
    delta_0 = seq(0.5, 1.0, by = 0.1),
    delta_1 = seq(0.5, 1.0, by = 0.1)
  ),
  grid_bounds = "lower",  # Only compute lower bound
  sl_lib = "SL.glm",
  V = 2,
  B = 500,
  seed = 1
)

names(fit_lower_only)


## -----------------------------------------------------------------------------
fit_psi2_grid <- mar_bounds(
  dat,
  Y = "Y",
  A = "A",
  C = "C",
  X = "X",
  estimand = "psi2",
  assumption = "point_psi2",
  param_grid = list(
    delta_0 = seq(0.3, 0.9, by = 0.1)
  ),
  sl_lib = "SL.glm",
  V = 2,
  B = 500,
  seed = 1
)

# View point estimate grid
print(fit_psi2_grid$estimate_grid$grid)


## -----------------------------------------------------------------------------
est_grid <- fit_psi2_grid$estimate_grid$grid

ggplot(est_grid, aes(x = delta_0)) +
  geom_line(aes(y = estimate), linewidth = 1) +
  geom_ribbon(aes(ymin = ci_lower_pointwise, ymax = ci_upper_pointwise),
              alpha = 0.2, fill = "blue") +
  geom_ribbon(aes(ymin = ci_lower_uniform, ymax = ci_upper_uniform),
              alpha = 0.2, fill = "red") +
  labs(
    title = "Ψ₂ Sensitivity Analysis",
    subtitle = "Blue = pointwise CIs, Red = simultaneous CIs",
    x = "δ₀",
    y = "Separable direct effect estimate"
  ) +
  theme_minimal()


## -----------------------------------------------------------------------------
grid <- fit_delta_grid$lower_grid$grid
diag <- grid[grid$delta_0 == grid$delta_1, ]

cat("Average pointwise CI width:",
    mean(diag$ci_upper_pointwise - diag$ci_lower_pointwise), "\n")
cat("Average uniform CI width:",
    mean(diag$ci_upper_uniform - diag$ci_lower_uniform), "\n")
cat("Ratio:",
    mean(diag$ci_upper_uniform - diag$ci_lower_uniform) /
    mean(diag$ci_upper_pointwise - diag$ci_lower_pointwise), "\n")


## -----------------------------------------------------------------------------
grid <- fit_delta_grid$lower_grid$grid

# Using pointwise CIs
significant_pointwise <- grid[grid$ci_lower_pointwise > 0, ]
cat("Pointwise: ", nrow(significant_pointwise), "out of", nrow(grid),
    "grid points exclude zero\n")

# Using uniform CIs
significant_uniform <- grid[grid$ci_lower_uniform > 0, ]
cat("Uniform: ", nrow(significant_uniform), "out of", nrow(grid),
    "grid points exclude zero\n")


## ---- eval = FALSE------------------------------------------------------------
## # Quick exploration: B = 100-200
## fit_quick <- mar_bounds(..., B = 100)
## 
## # Final analysis: B = 500-1000
## fit_final <- mar_bounds(..., B = 1000)


## ---- eval = FALSE------------------------------------------------------------
## # Coarse grid first
## fit_coarse <- mar_bounds(
##   ...,
##   param_grid = list(
##     delta_0 = seq(0.3, 1.0, by = 0.2),
##     delta_1 = seq(0.3, 1.0, by = 0.2)
##   )
## )
## 
## # Fine grid in region of interest
## fit_fine <- mar_bounds(
##   ...,
##   param_grid = list(
##     delta_0 = seq(0.5, 0.8, by = 0.05),
##     delta_1 = seq(0.5, 0.8, by = 0.05)
##   )
## )

