###########################################################
## lib_methods.R — BKMR-CMA + BART-CMA + rep runner + metrics
## Faithful to Chunks 4, 5, 7, 8 of original Rmd.
###########################################################

## ---- BKMR-CMA (manual mediation via SamplePred) ----
run_bkmr_cma <- function(dat, truth) {
  tryCatch({
    suppressPackageStartupMessages({
      library(bkmr); library(causalbkmr)
    })
    Z <- dat$Z; M <- dat$M; Y <- dat$Y; C <- dat$C; L <- dat$L; n <- dat$n
    scenario <- dat$scenario

    fit.m <- kmbayes(y = M, Z = Z, X = C,
                     iter = METHOD_CONFIG$bkmr_iter,
                     verbose = FALSE, varsel = TRUE)
    fit.y <- kmbayes(y = Y, Z = cbind(Z, M), X = C,
                     iter = METHOD_CONFIG$bkmr_iter,
                     verbose = FALSE, varsel = TRUE,
                     family = if (scenario$outcome_type == "binary") "binomial" else "gaussian")
    fit.y.TE <- kmbayes(y = Y, Z = Z, X = C,
                        iter = METHOD_CONFIG$bkmr_iter,
                        verbose = FALSE, varsel = TRUE,
                        family = if (scenario$outcome_type == "binary") "binomial" else "gaussian")

    sel <- seq(METHOD_CONFIG$bkmr_burnin + 1, METHOD_CONFIG$bkmr_iter, by = 1)
    K   <- METHOD_CONFIG$bkmr_K_med
    a <- truth$a; astar <- truth$astar
    X.predict <- matrix(colMeans(C), nrow = 1)

    set.seed(123)
    TE_samp <- NDE_samp <- NIE_samp <- numeric(length(sel))

    for (j in seq_along(sel)) {
      s <- sel[j]

      Znew_a  <- matrix(a,     nrow = 1)
      Znew_as <- matrix(astar, nrow = 1)
      Y_a  <- SamplePred(fit.y.TE, Znew = Znew_a,  Xnew = X.predict, sel = s)
      Y_as <- SamplePred(fit.y.TE, Znew = Znew_as, Xnew = X.predict, sel = s)
      TE_samp[j] <- Y_a - Y_as

      M_astar_samp <- numeric(K)
      for (k in 1:K) {
        M_pred <- SamplePred(fit.m, Znew = Znew_as, Xnew = X.predict, sel = s)
        M_astar_samp[k] <- rnorm(1, M_pred, sqrt(fit.m$sigsq.eps[s]))
      }

      Y_a_Mas <- Y_as_Mas <- numeric(K)
      for (k in 1:K) {
        Znew_a_m  <- matrix(c(a,     M_astar_samp[k]), nrow = 1)
        Znew_as_m <- matrix(c(astar, M_astar_samp[k]), nrow = 1)
        Y_a_Mas[k]  <- SamplePred(fit.y, Znew = Znew_a_m,  Xnew = X.predict, sel = s)
        Y_as_Mas[k] <- SamplePred(fit.y, Znew = Znew_as_m, Xnew = X.predict, sel = s)
      }
      NDE_samp[j] <- mean(Y_a_Mas) - mean(Y_as_Mas)
      NIE_samp[j] <- TE_samp[j] - NDE_samp[j]
    }

    pips <- ExtractPIPs(fit.m)
    list(
      method = "BKMR-CMA",
      NIE_est = mean(NIE_samp), NIE_lo = quantile(NIE_samp, 0.025), NIE_hi = quantile(NIE_samp, 0.975),
      NDE_est = mean(NDE_samp), NDE_lo = quantile(NDE_samp, 0.025), NDE_hi = quantile(NDE_samp, 0.975),
      TE_est  = mean(TE_samp),  TE_lo  = quantile(TE_samp,  0.025), TE_hi  = quantile(TE_samp,  0.975),
      PIPs = pips, success = TRUE
    )
  }, error = function(e) list(method = "BKMR-CMA", success = FALSE, error = conditionMessage(e)))
}

## ---- BART-CMA ----
run_bart_cma <- function(dat, truth) {
  tryCatch({
    suppressPackageStartupMessages(library(dbarts))
    Z <- dat$Z; M <- dat$M; Y <- dat$Y; C <- dat$C; L <- dat$L; n <- dat$n
    scenario <- dat$scenario; n_post <- METHOD_CONFIG$bart_ndpost

    fit_m <- bart(x.train = cbind(Z, C), y.train = M,
                  ntree = METHOD_CONFIG$bart_ntree,
                  ndpost = n_post, nskip = METHOD_CONFIG$bart_nskip,
                  keeptrees = TRUE, verbose = FALSE)
    fit_y <- bart(x.train = cbind(Z, M, C), y.train = Y,
                  ntree = METHOD_CONFIG$bart_ntree,
                  ndpost = n_post, nskip = METHOD_CONFIG$bart_nskip,
                  keeptrees = TRUE, verbose = FALSE)

    a <- truth$a; astar <- truth$astar
    Za  <- matrix(a,     n, L, byrow = TRUE)
    Zas <- matrix(astar, n, L, byrow = TRUE)
    Mas_draws <- predict(fit_m, newdata = cbind(Zas, C))
    Ma_draws  <- predict(fit_m, newdata = cbind(Za,  C))

    TE_d <- NDE_d <- NIE_d <- numeric(n_post)
    for (d in 1:n_post) {
      Mas_d <- Mas_draws[d, ] + rnorm(n, 0, fit_m$sigma[d])
      Ma_d  <- Ma_draws[d, ]  + rnorm(n, 0, fit_m$sigma[d])
      if (scenario$outcome_type == "binary") {
        Ya_Mas  <- mean(pnorm(predict(fit_y, newdata = cbind(Za,  Mas_d, C))[d, ]))
        Yas_Mas <- mean(pnorm(predict(fit_y, newdata = cbind(Zas, Mas_d, C))[d, ]))
        Ya_Ma   <- mean(pnorm(predict(fit_y, newdata = cbind(Za,  Ma_d,  C))[d, ]))
      } else {
        Ya_Mas  <- mean(predict(fit_y, newdata = cbind(Za,  Mas_d, C))[d, ])
        Yas_Mas <- mean(predict(fit_y, newdata = cbind(Zas, Mas_d, C))[d, ])
        Ya_Ma   <- mean(predict(fit_y, newdata = cbind(Za,  Ma_d,  C))[d, ])
      }
      TE_d[d]  <- Ya_Ma - Yas_Mas
      NDE_d[d] <- Ya_Mas - Yas_Mas
      NIE_d[d] <- Ya_Ma - Ya_Mas
    }

    list(
      method = "BART-CMA",
      NIE_est = mean(NIE_d), NIE_lo = quantile(NIE_d, 0.025), NIE_hi = quantile(NIE_d, 0.975),
      NDE_est = mean(NDE_d), NDE_lo = quantile(NDE_d, 0.025), NDE_hi = quantile(NDE_d, 0.975),
      TE_est  = mean(TE_d),  TE_lo  = quantile(TE_d,  0.025), TE_hi  = quantile(TE_d,  0.975),
      var_importance = colMeans(fit_m$varcount) / sum(colMeans(fit_m$varcount)),
      success = TRUE
    )
  }, error = function(e) list(method = "BART-CMA", success = FALSE, error = conditionMessage(e)))
}

## ---- Single-rep runner (BART + BKMR only; DeepMed runs separately) ----
run_one_rep <- function(rep_id, scenario, truth,
                        methods_to_run = c("BKMR-CMA", "BART-CMA")) {
set.seed(as.integer(SIM_CONFIG$seed + scenario$scenario_id * 1000L + rep_id))
  dat <- generate_data(scenario)
  results <- list()
  if ("BKMR-CMA" %in% methods_to_run) {
    t0 <- proc.time(); r <- run_bkmr_cma(dat, truth)
    r$runtime <- (proc.time() - t0)[3]; results[["BKMR-CMA"]] <- r
  }
  if ("BART-CMA" %in% methods_to_run) {
    t0 <- proc.time(); r <- run_bart_cma(dat, truth)
    r$runtime <- (proc.time() - t0)[3]; results[["BART-CMA"]] <- r
  }
  list(rep_id = rep_id, results = results, data = dat)
}

## ---- Performance metrics ----
compute_metrics <- function(rep_results_list, truth, method_name, estimand = "NIE") {
  reps <- lapply(rep_results_list, function(r) r$results[[method_name]])
  reps <- Filter(function(x) !is.null(x) && isTRUE(x$success), reps)
  ns <- length(reps)
  if (ns < 10) return(list(method = method_name, estimand = estimand,
                           n_success = ns, valid = FALSE))
  tv <- truth[[estimand]]
  fe <- paste0(estimand, "_est"); fl <- paste0(estimand, "_lo"); fh <- paste0(estimand, "_hi")
  est <- sapply(reps, function(x) x[[fe]])
  lo  <- sapply(reps, function(x) x[[fl]])
  hi  <- sapply(reps, function(x) x[[fh]])
  rt  <- sapply(reps, function(x) x$runtime)
  bias <- mean(est) - tv
  rb   <- if (abs(tv) > 1e-6) abs(bias / tv) * 100 else NA
  mse  <- mean((est - tv)^2); rmse <- sqrt(mse); ese <- sd(est)
  vc   <- !is.na(lo) & !is.na(hi)
  cov  <- if (sum(vc) > 0) mean((lo[vc] <= tv) & (hi[vc] >= tv)) else NA
  ciw  <- if (sum(vc) > 0) mean(hi[vc] - lo[vc], na.rm = TRUE) else NA
  list(method = method_name, estimand = estimand, true_val = tv, n_success = ns,
       bias = bias, rel_bias = rb, mse = mse, rmse = rmse, emp_se = ese,
       coverage = cov, ci_width = ciw, mean_runtime = mean(rt, na.rm = TRUE),
       valid = TRUE)
}
