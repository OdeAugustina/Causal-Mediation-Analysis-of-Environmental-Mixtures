###########################################################
## lib_dgp.R — Data Generating Process functions
## Faithful to Chunk 2 of original Rmd.
###########################################################

get_coefs <- function(L) {
  list(alpha_Z  = DGP_COEFS$alpha_Z_base[1:L],
       beta_Z   = DGP_COEFS$beta_Z_base[1:L],
       beta_ZM  = DGP_COEFS$beta_ZM_base[1:L])
}

make_Sigma_Z <- function(L, rho = SIM_CONFIG$rho_Z) {
  S <- matrix(rho, L, L); diag(S) <- 1; S
}

make_Sigma_C <- function(P = SIM_CONFIG$P_conf, rho = SIM_CONFIG$rho_C) {
  S <- matrix(rho, P, P); diag(S) <- 1; S
}

f_M_nonlinear <- function(Z, alpha_Z) {
  L <- ncol(Z); r <- rep(0, nrow(Z))
  r <- r + alpha_Z[1] * Z[,1] + 0.15 * Z[,1]^2
  if (L >= 2) r <- r + alpha_Z[2] * sqrt(abs(Z[,2])) * sign(Z[,2])
  if (L >= 3) r <- r + Z[, 3:L, drop = FALSE] %*% alpha_Z[3:L]
  r
}

f_Y_nonlinear <- function(Z, beta_Z) {
  L <- ncol(Z); r <- rep(0, nrow(Z))
  r <- r + beta_Z[1] * Z[,1] + 0.1 * Z[,1]^2
  if (L >= 2) r <- r + beta_Z[2] * sqrt(abs(Z[,2])) * sign(Z[,2])
  if (L >= 2) r <- r + 0.15 * Z[,1] * Z[,2]
  if (L >= 3) r <- r + Z[, 3:L, drop = FALSE] %*% beta_Z[3:L]
  r
}

g_ZM_interaction <- function(Z, M, beta_ZM) {
  rowSums(sweep(Z, 2, beta_ZM, "*")) * M
}

generate_data <- function(scenario, n = NULL) {
  if (is.null(n)) n <- scenario$n
  L <- scenario$L; P <- SIM_CONFIG$P_conf; coefs <- get_coefs(L)
  C <- mvrnorm(n, rep(0, P), make_Sigma_C(P))
  Z <- C %*% matrix(DGP_COEFS$theta_CZ, P, L) + mvrnorm(n, rep(0, L), make_Sigma_Z(L))
  if (scenario$linearity == "linear") {
    M_mean <- DGP_COEFS$alpha_0 + Z %*% coefs$alpha_Z + C %*% DGP_COEFS$alpha_C
  } else {
    M_mean <- DGP_COEFS$alpha_0 + f_M_nonlinear(Z, coefs$alpha_Z) + C %*% DGP_COEFS$alpha_C
  }
  M <- as.numeric(M_mean + rnorm(n, 0, DGP_COEFS$sigma_M))
  if (scenario$linearity == "linear") {
    eta <- DGP_COEFS$beta_0 + Z %*% coefs$beta_Z + DGP_COEFS$beta_M * M + C %*% DGP_COEFS$beta_C
  } else {
    eta <- DGP_COEFS$beta_0 + f_Y_nonlinear(Z, coefs$beta_Z) + DGP_COEFS$beta_M * M + C %*% DGP_COEFS$beta_C
  }
  if (scenario$interaction == "present") eta <- eta + g_ZM_interaction(Z, M, coefs$beta_ZM)
  if (scenario$outcome_type == "continuous") {
    Y <- as.numeric(eta + rnorm(n, 0, DGP_COEFS$sigma_Y))
  } else {
    Y <- rbinom(n, 1, plogis(as.numeric(eta)))
  }
  list(Z = Z, M = M, Y = Y, C = C, scenario = scenario, coefs = coefs, n = n, L = L)
}
