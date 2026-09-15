#!/usr/bin/env Rscript
###########################################################
## lib_dgp_plasmode.R (v2) -- NHANES-anchored plasmode cells
##
## CELLS
##   13  ("null")       fitted coefficients; true NIE ~ 0.
##                      Tests null recovery under real support geometry.
##   14  ("detectable") alpha_Z amplified, beta_M solved so the true
##                      Q25->Q75 NIE = 0.10 x SD(observed Y).
##
## BASIS
##   Franklin et al. (2014, Comput Stat Data Anal): resample observed
##   covariates/exposures unmodified to preserve their associations,
##   then simulate outcomes under an investigator-chosen true effect.
##   Schreck et al. (2024, Stat Med) review this as standard practice.
##   See also Shaw et al. (2025, arXiv:2504.11740) for causal-inference
##   specific cautions on the Franklin "Sample Treatment" framework.
##
## KEY INVARIANTS
##   1. pm_mu_M() / pm_mu_Y() are the ONLY structural equations.  The
##      generator and the truth both call them; they cannot drift.
##   2. The mediator is CENTRED in the outcome equation.  Without this
##      the interaction contributes ((a-astar).beta_ZM)*E[M] to the NDE,
##      and E[DII] = 1.604, which inflates the NDE by orders of
##      magnitude.  The synthetic DGP has E[M] ~ 0 so this never arose.
##   3. Contrast is FIXED (empirical Q25/Q75 of the full base exposure
##      distribution).  Truth is REPLICATE-SPECIFIC.
##   4. Var(M) is held at the observed DII variance when alpha_Z is
##      amplified, by absorbing the change into sigma_M.
##
## NOTE ON NIE TRUTH VARIABILITY
##   With Y linear in M and a fixed contrast, NIE truth is analytically
##   invariant to the resample (verified numerically: sd = 0 across 200
##   resamples).  TE and NDE do vary.  Truth is still computed per
##   replicate -- correct by construction, and non-degenerate for TE/NDE.
##
## USAGE
##   source("R/lib_setup.R"); source("R/lib_dgp_plasmode.R")
##   pm_build_spec(13); pm_build_spec(14)
##   dat <- generate_data_plasmode(rep_id = 1, sid = 13)
##   tr  <- pm_truth_all_reps(sid = 13, R = 100)
###########################################################

PM_CFG <- list(
  data_path = "/work/users/o/d/odeaug/BART_CMA/completed_liver_pfas_ALL.xlsx",
  sheet     = "Imputation_1",          # 20 sheets: Imputation_1 .. _20

  exposures = c("log_Lead","log_Cadmium","log_Mercury","log_PFOA","log_PFOS"),
  mediator  = "DII",
  outcome   = "LiverFibrosis",
  log_outcome = TRUE,                  # raw in file; applied pipeline logs it
  log_exposures = FALSE,               # ALREADY logged in the file
  scale_Z   = TRUE,

  confounders = c("Age","BMI","Gender_Male","Race_NHWhite",
                  "Education","Income","Smoker","Drinker"),

  n_plasmode = 1000L,                  # resample WITHOUT replacement (base 1819)
  lambda_int = 0.50,                   # beta_ZM = lambda * beta_M * v
  n_int      = 2000L,                  # MC draws for the truth check
  spec_dir   = "results",

  ## per-cell settings
  cells = list(
    "13" = list(label = "plasmode_null",       amp_alpha = 1.0, target_nie_sd = NA),
    "14" = list(label = "plasmode_detectable", amp_alpha = 3.0, target_nie_sd = 0.10)
  )
)

pm_spec_path <- function(sid, cfg = PM_CFG)
  file.path(cfg$spec_dir, sprintf("plasmode_spec_%02d.rds", sid))

## ==================== BASE DATA ====================
pm_load_base <- function(cfg = PM_CFG) {
  if (!requireNamespace("readxl", quietly = TRUE)) stop("need package 'readxl'")
  if (!file.exists(cfg$data_path)) stop("base data not found: ", cfg$data_path)
  df <- as.data.frame(readxl::read_excel(cfg$data_path, sheet = cfg$sheet))

  need <- c(cfg$exposures, cfg$mediator, cfg$outcome, cfg$confounders)
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("columns missing: ", paste(miss, collapse = ", "),
                         "\n available: ", paste(names(df), collapse = ", "))
  df <- df[stats::complete.cases(df[, need, drop = FALSE]), , drop = FALSE]

  Z <- as.matrix(df[, cfg$exposures, drop = FALSE])
  if (cfg$log_exposures) Z <- log(Z)
  if (cfg$scale_Z) Z <- scale(Z)
  Z <- matrix(as.numeric(Z), nrow(Z), ncol(Z),
              dimnames = list(NULL, cfg$exposures))
  C <- as.matrix(df[, cfg$confounders, drop = FALSE])
  storage.mode(C) <- "double"

  Yo <- as.numeric(df[[cfg$outcome]])
  if (cfg$log_outcome) {
    if (any(Yo <= 0)) stop("non-positive outcome values; cannot log")
    Yo <- log(Yo)
  }
  list(Z = Z, C = C, M_obs = as.numeric(df[[cfg$mediator]]), Y_obs = Yo,
       n_base = nrow(Z), L = ncol(Z), P = ncol(C))
}

## ==================== STRUCTURAL EQUATIONS ====================
pm_mu_M <- function(Z, C, sp)
  as.numeric(sp$alpha_0 + Z %*% sp$alpha_Z + C %*% sp$alpha_C)

## Mediator CENTRED at sp$m_center -- see invariant 2.
pm_mu_Y <- function(Z, M, C, sp) {
  as.numeric(sp$beta_0 + Z %*% sp$beta_Z + C %*% sp$beta_C +
               (sp$beta_M + as.numeric(Z %*% sp$beta_ZM)) * (M - sp$m_center))
}

## ==================== TRUTH ON A SUPPORT ====================
pm_truth_on_support <- function(Z, C, sp, n_int = PM_CFG$n_int, seed = NULL) {
  n <- nrow(Z)
  Za  <- matrix(sp$a,     n, ncol(Z), byrow = TRUE); colnames(Za)  <- colnames(Z)
  Zas <- matrix(sp$astar, n, ncol(Z), byrow = TRUE); colnames(Zas) <- colnames(Z)
  mu_hi <- pm_mu_M(Za, C, sp); mu_lo <- pm_mu_M(Zas, C, sp)

  A_hh <- mean(pm_mu_Y(Za,  mu_hi, C, sp))
  A_ll <- mean(pm_mu_Y(Zas, mu_lo, C, sp))
  A_hl <- mean(pm_mu_Y(Za,  mu_lo, C, sp))

  out <- list(TE = A_hh - A_ll, NDE = A_hl - A_ll, NIE = A_hh - A_hl)
  if (!is.null(n_int) && n_int > 1L) {
    if (!is.null(seed)) set.seed(as.integer(seed))
    s_hh <- s_ll <- s_hl <- 0
    for (k in seq_len(n_int)) {
      m_hi <- mu_hi + stats::rnorm(n, 0, sp$sigma_M)
      m_lo <- mu_lo + stats::rnorm(n, 0, sp$sigma_M)
      s_hh <- s_hh + mean(pm_mu_Y(Za,  m_hi, C, sp))
      s_ll <- s_ll + mean(pm_mu_Y(Zas, m_lo, C, sp))
      s_hl <- s_hl + mean(pm_mu_Y(Za,  m_lo, C, sp))
    }
    out$TE_mc  <- s_hh/n_int - s_ll/n_int
    out$NDE_mc <- s_hl/n_int - s_ll/n_int
    out$NIE_mc <- s_hh/n_int - s_hl/n_int
    out$max_abs_diff <- max(abs(c(out$TE - out$TE_mc, out$NDE - out$NDE_mc,
                                  out$NIE - out$NIE_mc)))
  } else out$max_abs_diff <- NA_real_
  out
}

## ==================== SPEC BUILD ====================
pm_build_spec <- function(sid, cfg = PM_CFG, verbose = TRUE) {
  cell <- cfg$cells[[as.character(sid)]]
  if (is.null(cell)) stop("no cell config for scenario ", sid)
  base <- pm_load_base(cfg); Z <- base$Z; C <- base$C

  fm <- stats::lm(base$M_obs ~ Z + C); cm <- stats::coef(fm)
  alpha_0 <- unname(cm[1]); alpha_Z <- unname(cm[2:(1+ncol(Z))])
  alpha_C <- unname(cm[(2+ncol(Z)):(1+ncol(Z)+ncol(C))])
  sigma_M <- stats::sd(stats::residuals(fm))

  fy <- stats::lm(base$Y_obs ~ Z + base$M_obs + C); cy <- stats::coef(fy)
  beta_0 <- unname(cy[1]); beta_Z <- unname(cy[2:(1+ncol(Z))])
  beta_Mf <- unname(cy[2+ncol(Z)])
  beta_C <- unname(cy[(3+ncol(Z)):(2+ncol(Z)+ncol(C))])
  sigma_Y <- stats::sd(stats::residuals(fy))
  alpha_C[is.na(alpha_C)] <- 0; beta_C[is.na(beta_C)] <- 0
  if (is.na(beta_Mf)) stop("fitted beta_M is NA")

  a     <- apply(Z, 2, stats::quantile, probs = 0.75)
  astar <- apply(Z, 2, stats::quantile, probs = 0.25)
  v     <- beta_Z / max(abs(beta_Z))
  varM_obs <- stats::var(base$M_obs)

  ## amplify the Z->M leg, holding Var(M) fixed by absorbing into sigma_M
  cA <- cell$amp_alpha
  alpha_Z2 <- alpha_Z * cA
  varZ <- as.numeric(t(alpha_Z2) %*% stats::cov(Z) %*% alpha_Z2)
  varC <- stats::var(as.numeric(C %*% alpha_C))
  s2   <- varM_obs - varZ - varC
  if (s2 <= 0)
    stop(sprintf("amp_alpha = %.2f exhausts the mediator variance budget ",
                 cA), sprintf("(Z share %.1f%%). Lower it.", 100*varZ/varM_obs))
  sigma_M2 <- sqrt(s2)
  alpha_0b <- mean(base$M_obs) - mean(Z %*% alpha_Z2 + C %*% alpha_C)

  d <- as.numeric((a - astar) %*% alpha_Z2)

  if (is.na(cell$target_nie_sd)) {          # cell 13: fitted magnitudes
    beta_M <- beta_Mf
  } else {                                   # cell 14: solve for the target
    tgt <- cell$target_nie_sd * stats::sd(base$Y_obs) * sign(beta_Mf * d)
    beta_M <- tgt / (d * (1 + cfg$lambda_int * as.numeric(a %*% v)))
  }
  beta_ZM <- cfg$lambda_int * beta_M * v

  sp <- list(alpha_0 = alpha_0b, alpha_Z = alpha_Z2, alpha_C = alpha_C,
             sigma_M = sigma_M2, beta_0 = beta_0, beta_Z = beta_Z,
             beta_M = beta_M, beta_ZM = beta_ZM, beta_C = beta_C,
             sigma_Y = sigma_Y, m_center = mean(base$M_obs),
             a = a, astar = astar, scenario_id = sid, cell = cell,
             sd_Y_obs = stats::sd(base$Y_obs), cfg = cfg,
             base = list(Z = Z, C = C, n_base = base$n_base),
             fitted = list(alpha_Z = alpha_Z, beta_M = beta_Mf,
                           varZ_share = varZ/varM_obs))

  full <- pm_truth_on_support(Z, C, sp, n_int = 200L, seed = 20260807L)
  mix  <- rowSums(Z)
  sp$diagnostics <- list(corr_mix_M_obs = stats::cor(mix, base$M_obs),
                         truth_full = full, varZ_share = varZ/varM_obs)

  if (verbose) {
    cat(sprintf("\n===== PLASMODE SPEC: scenario %d (%s) =====\n", sid, cell$label))
    cat(sprintf("base %d rows, sheet %s | L %d | P %d\n",
                base$n_base, cfg$sheet, base$L, base$P))
    cat(sprintf("SD(Y obs) %.4f | Corr(mixture, DII) %+.4f\n",
                sp$sd_Y_obs, sp$diagnostics$corr_mix_M_obs))
    cat(sprintf("alpha_Z amp %.2fx -> Z explains %.1f%% of Var(M) (fitted %.1f%%)\n",
                cA, 100*varZ/varM_obs, 100*sp$fitted$varZ_share))
    cat(sprintf("Q25->Q75 mediator shift %.4f (%.3f SD_M)\n", d, d/sqrt(varM_obs)))
    cat(sprintf("beta_M %.5f (fitted %.5f) | max|beta_ZM| %.5f\n",
                beta_M, beta_Mf, max(abs(beta_ZM))))
    cat(sprintf("TRUTH  TE %+.5f  NDE %+.5f  NIE %+.5f\n",
                full$TE, full$NDE, full$NIE))
    cat(sprintf("       NIE/SD_Y %+.4f   NDE/SD_Y %+.4f\n",
                full$NIE/sp$sd_Y_obs, full$NDE/sp$sd_Y_obs))
    cat(sprintf("analytic vs MC diff %.2e (want < 1e-2)\n", full$max_abs_diff))
    cat("==========================================\n\n")
  }
  dir.create(cfg$spec_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(sp, pm_spec_path(sid, cfg))
  cat("spec written: ", pm_spec_path(sid, cfg), "\n", sep = "")
  invisible(sp)
}

pm_load_spec <- function(sid, cfg = PM_CFG) {
  p <- pm_spec_path(sid, cfg)
  if (!file.exists(p)) stop("spec missing: ", p, " -- run pm_build_spec(", sid, ")")
  readRDS(p)
}

## ==================== RESAMPLE / GENERATE / TRUTH ====================
pm_sample_support <- function(rep_id, sp, cfg = PM_CFG) {
  set.seed(as.integer(SIM_CONFIG$seed + sp$scenario_id * 1000L + rep_id))
  nb <- nrow(sp$base$Z)
  if (cfg$n_plasmode > nb) stop("n_plasmode exceeds base rows")
  i <- sample.int(nb, cfg$n_plasmode, replace = FALSE)
  list(Z = sp$base$Z[i, , drop = FALSE], C = sp$base$C[i, , drop = FALSE], idx = i)
}

generate_data_plasmode <- function(rep_id, sid, sp = NULL, cfg = PM_CFG) {
  if (is.null(sp)) sp <- pm_load_spec(sid, cfg)
  s <- pm_sample_support(rep_id, sp, cfg); Z <- s$Z; C <- s$C; n <- nrow(Z)
  M <- pm_mu_M(Z, C, sp) + stats::rnorm(n, 0, sp$sigma_M)
  Y <- pm_mu_Y(Z, M, C, sp) + stats::rnorm(n, 0, sp$sigma_Y)
  list(Y = Y, M = M, Z = Z, C = C, n = n, L = ncol(Z), rep_id = rep_id,
       scenario = list(scenario_id = sp$scenario_id, n = n, L = ncol(Z),
                       linearity = "plasmode", interaction = "present",
                       outcome_type = "continuous", group = "ext_plasmode"))
}

pm_truth_all_reps <- function(sid, R = SIM_CONFIG$R, sp = NULL,
                              cfg = PM_CFG, verbose = TRUE) {
  if (is.null(sp)) sp <- pm_load_spec(sid, cfg)
  TE <- NDE <- NIE <- numeric(R)
  for (r in seq_len(R)) {
    s  <- pm_sample_support(r, sp, cfg)
    tt <- pm_truth_on_support(s$Z, s$C, sp, n_int = NULL)
    TE[r] <- tt$TE; NDE[r] <- tt$NDE; NIE[r] <- tt$NIE
  }
  out <- list(TE = mean(TE), NDE = mean(NDE), NIE = mean(NIE),
              TE_rep = TE, NDE_rep = NDE, NIE_rep = NIE,
              a = sp$a, astar = sp$astar,
              replicate_specific = TRUE, R = R, scenario_id = sid,
              sd_across_reps = c(TE = stats::sd(TE), NDE = stats::sd(NDE),
                                 NIE = stats::sd(NIE)))
  if (verbose)
    cat(sprintf("\nscenario %d truth over %d reps:\n  TE %+.5f (sd %.5f)\n  NDE %+.5f (sd %.5f)\n  NIE %+.5f (sd %.5f)\n\n",
                sid, R, out$TE, out$sd_across_reps["TE"], out$NDE,
                out$sd_across_reps["NDE"], out$NIE, out$sd_across_reps["NIE"]))
  out
}

pm_support_diagnostic <- function(rep_id = 1L, sid, sp = NULL, cfg = PM_CFG) {
  d <- generate_data_plasmode(rep_id, sid, sp, cfg)
  mix <- rowSums(scale(d$Z))
  hi <- mix > stats::quantile(mix, 0.75); lo <- d$M < stats::quantile(d$M, 0.25)
  c(corr_mix_M = stats::cor(mix, d$M),
    cross_world_pct = 100 * mean(hi & lo), n = d$n)
}

cat("[plasmode] lib_dgp_plasmode.R v2 loaded (cells 13 null / 14 detectable)\n")
