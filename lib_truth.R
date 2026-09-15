###########################################################
## lib_truth.R — Monte Carlo truth computation
## Faithful to Chunk 3 of original Rmd.
###########################################################

compute_truth <- function(scenario) {
  set.seed(SIM_CONFIG$seed + 9999)
  N <- SIM_CONFIG$N_truth; L <- scenario$L; P <- SIM_CONFIG$P_conf
  coefs <- get_coefs(L)
  C <- mvrnorm(N, rep(0, P), make_Sigma_C(P))
  Z <- C %*% matrix(DGP_COEFS$theta_CZ, P, L) + mvrnorm(N, rep(0, L), make_Sigma_Z(L))
  astar <- apply(Z, 2, quantile, 0.25)
  a     <- apply(Z, 2, quantile, 0.75)

  cMm <- function(z, Cm) {
    Zf <- matrix(z, nrow(Cm), L, byrow = TRUE)
    if (scenario$linearity == "linear") {
      DGP_COEFS$alpha_0 + Zf %*% coefs$alpha_Z + Cm %*% DGP_COEFS$alpha_C
    } else {
      DGP_COEFS$alpha_0 + f_M_nonlinear(Zf, coefs$alpha_Z) + Cm %*% DGP_COEFS$alpha_C
    }
  }

  cYm <- function(z, Mv, Cm) {
    Zf <- matrix(z, nrow(Cm), L, byrow = TRUE)
    if (scenario$linearity == "linear") {
      eta <- DGP_COEFS$beta_0 + Zf %*% coefs$beta_Z + DGP_COEFS$beta_M * Mv + Cm %*% DGP_COEFS$beta_C
    } else {
      eta <- DGP_COEFS$beta_0 + f_Y_nonlinear(Zf, coefs$beta_Z) + DGP_COEFS$beta_M * Mv + Cm %*% DGP_COEFS$beta_C
    }
    if (scenario$interaction == "present") eta <- eta + g_ZM_interaction(Zf, Mv, coefs$beta_ZM)
    eta
  }

  N_mc <- min(N, 100000); idx <- sample(N, N_mc); Cm <- C[idx, , drop = FALSE]
  Ma <- as.numeric(cMm(astar, Cm) + rnorm(N_mc, 0, DGP_COEFS$sigma_M))
  Mb <- as.numeric(cMm(a,     Cm) + rnorm(N_mc, 0, DGP_COEFS$sigma_M))

  if (scenario$outcome_type == "continuous") {
    Ya <- as.numeric(cYm(a,     Ma, Cm))
    Yb <- as.numeric(cYm(astar, Ma, Cm))
    Yc <- as.numeric(cYm(a,     Mb, Cm))
  } else {
    Ya <- plogis(as.numeric(cYm(a,     Ma, Cm)))
    Yb <- plogis(as.numeric(cYm(astar, Ma, Cm)))
    Yc <- plogis(as.numeric(cYm(a,     Mb, Cm)))
  }

  mix_idx <- rowSums(Z[idx, ]); dm <- median(mix_idx)

  list(
    TE  = mean(Yc) - mean(Yb),
    NDE = mean(Ya) - mean(Yb),
    NIE = mean(Yc) - mean(Ya),
    a = a, astar = astar,
    a_deepmed     = apply(Z[idx, ][mix_idx >  dm, , drop = FALSE], 2, mean),
    astar_deepmed = apply(Z[idx, ][mix_idx <= dm, , drop = FALSE], 2, mean),
    scenario_id = scenario$scenario_id
  )
}
