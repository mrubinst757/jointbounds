#' Coefficient vectors b such that bound = mean(b' phi)
#'
#' theta = (E(mu0), E(mu1), E(pi0), E(pi1), E(mu0*pi0), E(mu1*pi1)) = (m1_0, m1_1, m2_0, m2_1, m3_0, m3_1)
#' Order of phi columns: phi_1_0, phi_1_1, phi_2_0, phi_2_1, phi_3_0, phi_3_1
#'
#' All bounds/estimands are linear: value = b %*% theta_hat = mean(b %*% t(phi)) = mean(phi %*% b).

coef_general_lower_ate <- function() {
  c(-1, 1, -1, 0, 1, -1)  # ell_0
}

coef_general_upper_ate <- function() {
  c(-1, 1, 0, 1, 1, -1)   # u_0
}

coef_delta_lower_ate <- function(delta_0u, delta_1u) {
  # ell_0,delta = Psi_tilde - delta_1u*m3_1 - delta_0u*(m2_0 - m3_0)? No.
  # ell_0,delta = Psi_tilde - [ delta_1u*pi_1*mu_1 + delta_0u*pi_0*(1-mu_0) ] = m1_1 - m1_0 - delta_1u*m3_1 - delta_0u*E(pi_0*(1-mu_0))
  # E(pi_0*(1-mu_0)) = E(pi_0) - E(pi_0*mu_0) = m2_0 - m3_0. So b = (-1, 1, -delta_0u, 0, delta_0u, -delta_1u).
  c(-1, 1, -delta_0u, 0, delta_0u, -delta_1u)
}

coef_delta_upper_ate <- function(delta_0u, delta_1u) {
  # u_0,delta = Psi_tilde + delta_1u*pi_1*(1-mu_1) + delta_0u*pi_0*mu_0 = m1_1 - m1_0 + delta_1u*(m2_1 - m3_1) + delta_0u*m3_0
  c(-1, 1, 0, delta_1u, delta_0u, -delta_1u)
}

coef_mono_pos_lower_ate <- function(delta_0u) {
  c(-1, 1, -delta_0u, 0, delta_0u, 0)
}

coef_mono_pos_upper_ate <- function(delta_1u) {
  c(-1, 1, 0, delta_1u, 0, -delta_1u)
}

coef_mono_neg_lower_ate <- function(delta_1u) {
  c(-1, 1, 0, 0, 0, -delta_1u)
}

coef_mono_neg_upper_ate <- function(delta_0u) {
  c(-1, 1, 0, 0, delta_0u, 0)
}

coef_bounded_risk_lower_ate <- function(delta_0u, delta_1u, tau_0, tau_1) {
  # tilde_ell = Psi_tilde + delta_1u*(1/tau_1 - 1)*m3_1 - delta_0u*(tau_0 - 1)*m3_0
  c(-1, 1, 0, 0, -delta_0u * (tau_0 - 1), delta_1u * (1/tau_1 - 1))
}

coef_bounded_risk_upper_ate <- function(delta_0u, delta_1u, tau_0, tau_1) {
  # tilde_u = Psi_tilde + delta_1u*(tau_1 - 1)*m3_1 - delta_0u*(1/tau_0 - 1)*m3_0
  c(-1, 1, 0, 0, -delta_0u * (1/tau_0 - 1), delta_1u * (tau_1 - 1))
}

coef_point_ate <- function(delta_0, delta_1, tau) {
  # Psi_0 = Psi_tilde + (tau-1)[ delta_1*m3_1 - delta_0*m3_0 ]
  c(-1, 1, 0, 0, -((tau - 1) * delta_0), (tau - 1) * delta_1)
}

# Psi_1 bounds (Proposition 8)
coef_lower_psi1 <- function(delta_1l, delta_0u) {
  # ell_1 = Psi_tilde + delta_1l*(m2_1 - m3_1) - delta_0u*(m2_0 - m3_0)
  c(-1, 1, -delta_0u, delta_1l, delta_0u, -delta_1l)
}

coef_upper_psi1 <- function(delta_1u, delta_0l) {
  c(-1, 1, -delta_0l, delta_1u, delta_0l, -delta_1u)
}

# Psi_1 point (Proposition 9)
coef_point_psi1 <- function(delta_0, delta_1) {
  c(-1, 1, -delta_0, delta_1, delta_0, -delta_1)
}

# Psi_2 point (Proposition 9): Psi_2 = Psi_tilde - delta_0*(E(mu_1*pi_0) - m3_0). E(mu_1*pi_0) uses phi_3,01.
# So estimate = mean(phi_1_1 - phi_1_0 - delta_0*(phi_3_01 - phi_3_0)) where phi_3_01 = phi_1_1*pi0 + phi_2_0*mu1 - mu1*pi0.
# This is not a linear b'*phi with the 6-vector; we need the extra term. So we compute it explicitly in the caller.

#' Compute all requested bound/point estimates from influence matrix and parameters
#'
#' @param phi n x 6 influence matrix.
#' @param estimand "ate", "psi1" (composite ATE), "psi2" (SDE).
#' @param assumption "general", "bounded_delta", "monotonicity_pos", "monotonicity_neg", "bounded_risk", "point_ate", "point_psi1", "point_psi2".
#' @param delta_0u,delta_1u Upper bounds on proportion informative missingness (for bounded/mono bounds).
#' @param delta_0l,delta_1l Lower bounds (for Psi_1 bounds).
#' @param delta_0,delta_1 Point values (for point identification).
#' @param tau_0,tau_1 Bounded risk ratio params (for bounded_risk).
#' @param tau Single sensitivity parameter for point_ate (risk ratio in both arms).
#' @param mu0,mu1 Outcome nuisance functions (n-vectors).
#' @param pi0,pi1 Missingness probability nuisance functions (n-vectors).
#' @param smooth_approximation Logical; for psi2 bounded_delta and bounded_risk, use smooth approximation (TRUE) or indicator method (FALSE).
#' @param epsilon Smoothing parameter for smooth approximation (default 0.001).
#' @return List with elements estimate, lower, upper (as appropriate), and optional se_* for asymptotic SEs.
compute_bounds <- function(phi,
                          estimand = c("ate", "psi1", "psi2"),
                          assumption = c("general", "bounded_delta", "monotonicity_pos", "monotonicity_neg", "bounded_risk", "bounded_risk_unbounded_tau", "point_ate", "point_psi1", "point_psi2"),
                          delta_0u = 1, delta_1u = 1,
                          delta_0l = 0, delta_1l = 0,
                          delta_0 = NULL, delta_1 = NULL,
                          tau_0 = NULL, tau_1 = NULL,
                          tau = NULL,
                          mu0 = NULL, mu1 = NULL,
                          pi0 = NULL, pi1 = NULL,
                          smooth_approximation = FALSE,
                          epsilon = 0.001) {
  estimand <- match.arg(estimand)
  assumption <- match.arg(assumption)
  n <- nrow(phi)
  th <- theta_hat(phi)
  V_phi <- phi_var(phi)

  get_est <- function(b) {
    est <- mean(phi %*% b)
    se <- sqrt(drop(t(b) %*% V_phi %*% b) / n)
    list(est = est, se = se)
  }

  if (estimand == "ate") {
    if (assumption == "general") {
      # Corollary 1: when delta parameters are provided, general bounds shrink
      if (!is.null(delta_0u) && !is.null(delta_1u) && (delta_0u != 1 || delta_1u != 1)) {
        bl <- coef_delta_lower_ate(delta_0u, delta_1u)
        bu <- coef_delta_upper_ate(delta_0u, delta_1u)
      } else {
        bl <- coef_general_lower_ate()
        bu <- coef_general_upper_ate()
      }
      L <- get_est(bl)
      U <- get_est(bu)
      return(list(naive = th[2] - th[1], lower = L$est, upper = U$est, se_lower = L$se, se_upper = U$se))
    }
    if (assumption == "bounded_delta") {
      bl <- coef_delta_lower_ate(delta_0u, delta_1u)
      bu <- coef_delta_upper_ate(delta_0u, delta_1u)
      L <- get_est(bl)
      U <- get_est(bu)
      return(list(naive = th[2] - th[1], lower = L$est, upper = U$est, se_lower = L$se, se_upper = U$se))
    }
    if (assumption == "monotonicity_pos") {
      bl <- coef_mono_pos_lower_ate(delta_0u)
      bu <- coef_mono_pos_upper_ate(delta_1u)
      L <- get_est(bl)
      U <- get_est(bu)
      return(list(naive = th[2] - th[1], lower = L$est, upper = U$est, se_lower = L$se, se_upper = U$se))
    }
    if (assumption == "monotonicity_neg") {
      bl <- coef_mono_neg_lower_ate(delta_1u)
      bu <- coef_mono_neg_upper_ate(delta_0u)
      L <- get_est(bl)
      U <- get_est(bu)
      return(list(naive = th[2] - th[1], lower = L$est, upper = U$est, se_lower = L$se, se_upper = U$se))
    }
    if (assumption == "bounded_risk_unbounded_tau") {
      bl <- coef_bounded_risk_lower_ate(delta_0u, delta_1u, tau_0, tau_1)
      bu <- coef_bounded_risk_upper_ate(delta_0u, delta_1u, tau_0, tau_1)
      L <- get_est(bl)
      U <- get_est(bu)
      return(list(naive = th[2] - th[1], lower = L$est, upper = U$est, se_lower = L$se, se_upper = U$se))
    }
    if (assumption == "bounded_risk") {
      if (is.null(mu0) || is.null(mu1) || is.null(pi0) || is.null(pi1)) {
        stop("For bounded_risk provide mu0, mu1, pi0, and pi1.")
      }

      phi1_0 <- phi[, 1]
      phi1_1 <- phi[, 2]
      phi2_0 <- phi[, 3]
      phi2_1 <- phi[, 4]
      phi3_0 <- phi[, 5]
      phi3_1 <- phi[, 6]

      if (smooth_approximation) {
        # Smooth approximation using Φ_ε
        # Φ_ε(x) = Φ(x/ε) where Φ is standard normal CDF
        phi_eps_tau1_mu1 <- stats::pnorm((tau_1 * mu1 - 1) / epsilon)  # Φ_ε(τ_1*μ_1 - 1)
        phi_eps_tau0_mu0 <- stats::pnorm((tau_0 * mu0 - 1) / epsilon)  # Φ_ε(τ_0*μ_0 - 1)
        phi_eps_1_tau1_mu1 <- stats::pnorm((1 - tau_1 * mu1) / epsilon)  # Φ_ε(1 - τ_1*μ_1)
        phi_eps_1_tau0_mu0 <- stats::pnorm((1 - tau_0 * mu0) / epsilon)  # Φ_ε(1 - τ_0*μ_0)

        # φ_ε(x) = (1/ε)φ(x/ε) where φ is standard normal PDF
        phi_density_tau1_mu1 <- stats::dnorm((tau_1 * mu1 - 1) / epsilon) / epsilon
        phi_density_1_tau1_mu1 <- stats::dnorm((1 - tau_1 * mu1) / epsilon) / epsilon
        phi_density_tau0_mu0 <- stats::dnorm((tau_0 * mu0 - 1) / epsilon) / epsilon
        phi_density_1_tau0_mu0 <- stats::dnorm((1 - tau_0 * mu0) / epsilon) / epsilon

        # Upper bound: uses min{1-mu_1, mu_1*(tau_1-1)} for a=1
        # Base smooth formula
        ub_vec <- (phi1_1 - phi1_0) +
          delta_1u * ((phi2_1 - phi3_1) * phi_eps_tau1_mu1 + (tau_1 - 1) * phi3_1 * (1 - phi_eps_tau1_mu1)) -
          delta_0u * ((tau_0^{-1} - 1) * phi3_0)

        # Additional correction term for a=1
        correction_upper <- delta_1u * tau_1 * phi1_1 * (
          pi1 * (1 - mu1) * phi_density_tau1_mu1 +
          (tau_1 - 1) * pi1 * mu1 * phi_density_1_tau1_mu1
        )
        ub_vec <- ub_vec + correction_upper

        # Lower bound: uses min{1-mu_0, mu_0*(tau_0-1)} for a=0
        # Base smooth formula
        lb_vec <- (phi1_1 - phi1_0) +
          delta_1u * ((tau_1^{-1} - 1) * phi3_1) -
          delta_0u * ((phi2_0 - phi3_0) * phi_eps_tau0_mu0 + (tau_0 - 1) * phi3_0 * (1 - phi_eps_tau0_mu0))

        # Additional correction term for a=0
        correction_lower <- delta_0u * tau_0 * phi1_0 * (
          pi0 * (1 - mu0) * phi_density_tau0_mu0 +
          (tau_0 - 1) * pi0 * mu0 * phi_density_1_tau0_mu0
        )
        lb_vec <- lb_vec - correction_lower  # Note: subtract for lower bound

        lower <- mean(lb_vec)
        upper <- mean(ub_vec)
        se_lower <- sqrt(stats::var(lb_vec) / n)
        se_upper <- sqrt(stats::var(ub_vec) / n)
        return(list(naive = th[2] - th[1], lower = lower, upper = upper, se_lower = se_lower, se_upper = se_upper))
      } else {
        # Indicator function method (default)
        mask1 <- (tau_1 * mu1 > 1)
        mask0 <- (tau_0 * mu0 > 1)

        ub_vec <- (phi1_1 - phi1_0) +
          delta_1u * ( (phi2_1 - phi3_1) * mask1 + (tau_1 - 1) * phi3_1 * (!mask1) ) -
          delta_0u * ( (tau_0^{-1} - 1) * phi3_0 )

        lb_vec <- (phi1_1 - phi1_0) +
          delta_1u * ( (tau_1^{-1} - 1) * phi3_1 ) -
          delta_0u * ( (phi2_0 - phi3_0) * mask0 + (tau_0 - 1) * phi3_0 * (!mask0) )

        lower <- mean(lb_vec)
        upper <- mean(ub_vec)
        se_lower <- sqrt(stats::var(lb_vec) / n)
        se_upper <- sqrt(stats::var(ub_vec) / n)
        return(list(naive = th[2] - th[1], lower = lower, upper = upper, se_lower = se_lower, se_upper = se_upper))
      }
    }
    if (assumption == "point_ate") {
      b <- coef_point_ate(delta_0, delta_1, tau)
      P <- get_est(b)
      return(list(estimate = P$est, se = P$se, naive = th[2] - th[1]))
    }
  }

  if (estimand == "psi1") {
    if (assumption == "general") {
      # same as bounded with delta_au=1, delta_al=0
      bl <- coef_lower_psi1(0, 1)
      bu <- coef_upper_psi1(1, 0)
      L <- get_est(bl)
      U <- get_est(bu)
      return(list(naive = th[2] - th[1], lower = L$est, upper = U$est, se_lower = L$se, se_upper = U$se))
    }
    if (assumption == "bounded_delta") {
      bl <- coef_lower_psi1(delta_1l, delta_0u)
      bu <- coef_upper_psi1(delta_1u, delta_0l)
      L <- get_est(bl)
      U <- get_est(bu)
      return(list(naive = th[2] - th[1], lower = L$est, upper = U$est, se_lower = L$se, se_upper = U$se))
    }
    if (assumption == "point_psi1") {
      b <- coef_point_psi1(delta_0, delta_1)
      P <- get_est(b)
      return(list(estimate = P$est, se = P$se, naive = th[2] - th[1]))
    }
  }

  if (estimand == "psi2") {
    if (assumption == "point_psi2") {
      if (is.null(pi0) || is.null(mu1)) stop("For point_psi2 provide nuisance pi0 and mu1.")
      phi_3_01 <- phi[, 2] * pi0 + phi[, 3] * mu1 - mu1 * pi0
      psi2_vec <- phi[, 2] - phi[, 1] - delta_0 * (phi_3_01 - phi[, 5])
      est <- mean(psi2_vec)
      se <- sqrt(stats::var(psi2_vec) / n)
      return(list(estimate = est, se = se, naive = th[2] - th[1]))
    }
    if (assumption == "bounded_delta") {
      if (is.null(mu0) || is.null(mu1) || is.null(pi0)) {
        stop("For psi2 bounded_delta provide mu0, mu1, and pi0.")
      }

      if (smooth_approximation) {
        # Smooth x_+ by q_epsilon(x) = {x + sqrt(x^2 + epsilon^2)}/2.
        # Unlike x Phi(x/epsilon), this approximation remains nonnegative.
        phi1_0 <- phi[, 1]
        phi1_1 <- phi[, 2]
        phi2_0 <- phi[, 3]
        contrast <- mu1 - mu0
        reference_score <- phi1_1 - phi1_0
        contrast_score <- reference_score - contrast
        pi_score <- phi2_0 - pi0

        smooth_support_score <- function(x, x_score) {
          root <- sqrt(x^2 + epsilon^2)
          q <- (x + root) / 2
          q_prime <- (1 + x / root) / 2
          pi0 * q + pi0 * q_prime * x_score + q * pi_score
        }
        support_pos <- smooth_support_score(contrast, contrast_score)
        support_neg <- smooth_support_score(-contrast, -contrast_score)
        # A one-step estimate of a nonnegative support functional can be
        # negative in a small sample.  Projecting it to zero preserves the
        # defining order of the endpoints and is asymptotically inactive away
        # from the zero-support boundary.
        support_pos <- if (mean(support_pos) > 0) support_pos else rep(0, n)
        support_neg <- if (mean(support_neg) > 0) support_neg else rep(0, n)
        lb_vec <- reference_score - delta_0u * support_pos
        ub_vec <- reference_score + delta_0u * support_neg

        lower <- mean(lb_vec)
        upper <- mean(ub_vec)
        se_lower <- sqrt(stats::var(lb_vec) / n)
        se_upper <- sqrt(stats::var(ub_vec) / n)

        return(list(naive = th[2] - th[1], lower = lower, upper = upper, se_lower = se_lower, se_upper = se_upper))
      } else {
        # Exact positive-part support functions.  The lower endpoint subtracts
        # delta_0u E[pi0 (mu1-mu0)_+]; the upper endpoint adds
        # delta_0u E[pi0 (mu0-mu1)_+].
        phi1_0 <- phi[, 1]
        phi1_1 <- phi[, 2]
        phi2_0 <- phi[, 3]
        contrast <- mu1 - mu0
        reference_score <- phi1_1 - phi1_0
        contrast_score <- reference_score - contrast
        pi_score <- phi2_0 - pi0

        support_pos <- pi0 * pmax(contrast, 0) +
          pi0 * (contrast > 0) * contrast_score +
          pmax(contrast, 0) * pi_score
        support_neg <- pi0 * pmax(-contrast, 0) -
          pi0 * (contrast < 0) * contrast_score +
          pmax(-contrast, 0) * pi_score
        support_pos <- if (mean(support_pos) > 0) support_pos else rep(0, n)
        support_neg <- if (mean(support_neg) > 0) support_neg else rep(0, n)
        lb_vec <- reference_score - delta_0u * support_pos
        ub_vec <- reference_score + delta_0u * support_neg

        lower <- mean(lb_vec)
        upper <- mean(ub_vec)
        se_lower <- sqrt(stats::var(lb_vec) / n)
        se_upper <- sqrt(stats::var(ub_vec) / n)

        return(list(naive = th[2] - th[1], lower = lower, upper = upper, se_lower = se_lower, se_upper = se_upper))
      }
    }
    return(list(naive = th[2] - th[1], message = "Psi_2 bounds only supported for bounded_delta assumption"))
  }

  stop("Unsupported estimand/assumption combination.")
}
