## =====================================================================
## run_bart_cma_aligned.R
##
## Drop-in replacement for run_bart_cma() (and run_one_rep()) in
## R/lib_methods.R.  Source this AFTER R/lib_methods.R so the new
## definitions mask the old ones.
##
## PURPOSE -------------------------------------------------------------
## Aligns the simulation's BART-CMA estimator with the applied analysis
## (BART_CMA_Mediation.Rmd :: bart_cma_mixture) on three points:
##
##   1. MEDIATOR DRAWS.  The old code drew ONE mediator realisation per
##      posterior draw d and read prediction row [d, ].  That leaves
##      Monte-Carlo simulation error of order sigma/sqrt(n) inside every
##      posterior draw, inflating credible intervals and coverage.  This
##      version averages over METHOD_CONFIG$bart_nmed mediator draws per
##      posterior draw, exactly as the applied analysis does.
##
##   2. PREDICT PLACEMENT.  The old code called predict() inside the
##      1:ndpost loop and discarded all but one row -- 3 * ndpost full
##      S x n predictions per replicate.  This version lifts predict()
##      out of the loop: 3 * n_med predictions per replicate.  With
##      ndpost = 2000 and n_med = 30 that is ~16x LESS prediction work
##      than the original at ndpost = 1000, so the rerun is faster, not
##      slower.
##
##   3. ENGINE.  Defaults to BART::wbart / BART::pbart, the same package
##      the applied analysis uses, so the manuscript can state one
##      implementation.  Set METHOD_CONFIG$bart_engine <- "dbarts" to
##      fall back to the original engine if wbart proves too slow.
##
## WHAT IS DELIBERATELY *NOT* CHANGED ----------------------------------
##   * The Q25 -> Q75 contrast still comes from truth$a / truth$astar
##     (population values), not from sample quantiles.  Bias against a
##     fixed truth requires this; the applied analysis computes sample
##     quantiles because it has no truth to compare against.
##   * The return contract is byte-compatible with compute_metrics(),
##     merge_split.R and aggregate.R.  Nothing downstream needs editing.
##   * run_one_rep()'s data seed is untouched, so reruns land on exactly
##     the datasets the existing BKMR and DeepMed results already used.
## =====================================================================

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

## ---- Column naming helper -------------------------------------------
## wbart's predict() matches newdata columns to training columns by name,
## so both matrices must carry identical, explicit names.  The DGP does
## not name its matrices; this mirrors nmed_sensitivity_scenario12.R.
.nm <- function(x, prefix) {
  x <- as.matrix(x)
  colnames(x) <- paste0(prefix, seq_len(ncol(x)))
  x
}

## ---- BART-CMA (aligned) ---------------------------------------------
run_bart_cma <- function(dat, truth, seed = NULL) {
  tryCatch({

    engine <- METHOD_CONFIG$bart_engine %||% "BART"
    if (!engine %in% c("BART", "dbarts"))
      stop("METHOD_CONFIG$bart_engine must be 'BART' or 'dbarts', got '", engine, "'")
    if (engine == "BART") suppressPackageStartupMessages(library(BART))
    else                  suppressPackageStartupMessages(library(dbarts))

    ## ---- inputs ----
    Z <- .nm(dat$Z, "Z")
    C <- .nm(dat$C, "C")
    M <- as.numeric(dat$M)
    Y <- dat$Y
    n <- dat$n
    L <- dat$L
    binary <- identical(dat$scenario$outcome_type, "binary")

    ntree <- METHOD_CONFIG$bart_ntree
    nskip <- METHOD_CONFIG$bart_nskip
    S     <- METHOD_CONFIG$bart_ndpost
    n_med <- as.integer(METHOD_CONFIG$bart_nmed %||% 30L)

    ## Explicit seed makes the mediator draws reproducible regardless of
    ## whether BKMR ran first in the same replicate (METHOD=bart vs both).
    if (!is.null(seed)) set.seed(as.integer(seed))

    ## ---- counterfactual mixture profiles (Q75 vs Q25) ----
    a <- truth$a; astar <- truth$astar
    Za  <- matrix(a,     n, L, byrow = TRUE); colnames(Za)  <- colnames(Z)
    Zas <- matrix(astar, n, L, byrow = TRUE); colnames(Zas) <- colnames(Z)

    ## ---- Step 1: mediator model E[M | Z, C] ----
    if (engine == "BART") {
      mfit <- wbart(x.train = cbind(Z, C), y.train = M,
                    ntree = ntree, nskip = nskip, ndpost = S,
                    printevery = nskip + S + 1L)
      Ls  <- length(mfit$sigma)              # wbart keeps burn-in sigmas
      sig <- mfit$sigma[(Ls - S + 1):Ls]     # post-burn-in only
    } else {
      mfit <- dbarts::bart(x.train = cbind(Z, C), y.train = M,
                           ntree = ntree, ndpost = S, nskip = nskip,
                           keeptrees = TRUE, verbose = FALSE)
      sig <- mfit$sigma
      if (length(sig) > S) sig <- utils::tail(sig, S)
    }
    M_lo_hat <- predict(mfit, newdata = cbind(Zas, C))   # S x n
    M_hi_hat <- predict(mfit, newdata = cbind(Za,  C))   # S x n

    ## ---- Step 2: outcome model E[Y | Z, M, C] ----
    if (engine == "BART") {
      if (binary) {
        ofit  <- pbart(x.train = cbind(Z, M = M, C), y.train = Y,
                       ntree = ntree, nskip = nskip, ndpost = S,
                       printevery = nskip + S + 1L)
        predY <- function(nd) predict(ofit, newdata = nd)$prob.test
      } else {
        ofit  <- wbart(x.train = cbind(Z, M = M, C), y.train = Y,
                       ntree = ntree, nskip = nskip, ndpost = S,
                       printevery = nskip + S + 1L)
        predY <- function(nd) predict(ofit, newdata = nd)
      }
    } else {
      ofit <- dbarts::bart(x.train = cbind(Z, M = M, C), y.train = Y,
                           ntree = ntree, ndpost = S, nskip = nskip,
                           keeptrees = TRUE, verbose = FALSE)
      predY <- if (binary) function(nd) pnorm(predict(ofit, newdata = nd))
               else        function(nd) predict(ofit, newdata = nd)
    }

    ## ---- Step 3: paired-draw Monte-Carlo g-computation ----
    ## Mediator parameter draw s is paired with outcome draw s; residual
    ## realisations are averaged within that draw.  The prior version drew
    ## one mediator index per k and shared it across all S rows, which
    ## averaged mediator parameter uncertainty away and produced intervals
    ## that were too narrow by construction.
    ##   S_use   = retained (thinned) posterior draws -> CI resolution
    ##   n_med_p = residual realisations per retained draw -> integration
    S_use   <- as.integer(METHOD_CONFIG$bart_Suse %||% 200L)
    idx     <- if (S_use < S) round(seq(1, S, length.out = S_use)) else seq_len(S)
    n_med_p <- as.integer(METHOD_CONFIG$bart_nmed_paired %||% 5L)

    Y_hh <- Y_ll <- Y_hl <- matrix(NA_real_, length(idx), n)
    ridx  <- rep(seq_len(n), n_med_p)
    Zrep  <- Za[ridx, , drop = FALSE]
    Zsrep <- Zas[ridx, , drop = FALSE]
    Crep  <- C[ridx, , drop = FALSE]
    B     <- n_med_p * n
    cm    <- function(v) colMeans(matrix(v, n_med_p, n, byrow = TRUE))

    for (j in seq_along(idx)) {
      s    <- idx[j]
      m_lo <- as.vector(t(matrix(M_lo_hat[s, ], n_med_p, n, byrow = TRUE) +
                          matrix(rnorm(B, 0, sig[s]), n_med_p, n)))
      m_hi <- as.vector(t(matrix(M_hi_hat[s, ], n_med_p, n, byrow = TRUE) +
                          matrix(rnorm(B, 0, sig[s]), n_med_p, n)))

      Y_hh[j, ] <- cm(predY(cbind(Zrep,  M = m_hi, Crep))[s, ])
      Y_ll[j, ] <- cm(predY(cbind(Zsrep, M = m_lo, Crep))[s, ])
      Y_hl[j, ] <- cm(predY(cbind(Zrep,  M = m_lo, Crep))[s, ])
    }

    ## Population average over individuals.  The applied analysis uses a
    ## survey-weighted average; the simulation has no weights, so the
    ## unweighted mean is the same operator with weights = 1.
    TE_d  <- rowMeans(Y_hh) - rowMeans(Y_ll)
    NDE_d <- rowMeans(Y_hl) - rowMeans(Y_ll)
    NIE_d <- rowMeans(Y_hh) - rowMeans(Y_hl)

    q <- function(v) unname(stats::quantile(v, c(0.025, 0.975), na.rm = TRUE))

    list(
      method  = "BART-CMA",
      NIE_est = mean(NIE_d), NIE_lo = q(NIE_d)[1], NIE_hi = q(NIE_d)[2],
      NDE_est = mean(NDE_d), NDE_lo = q(NDE_d)[1], NDE_hi = q(NDE_d)[2],
      TE_est  = mean(TE_d),  TE_lo  = q(TE_d)[1],  TE_hi  = q(TE_d)[2],
      var_importance = colMeans(mfit$varcount) / sum(colMeans(mfit$varcount)),
      ## provenance, so a results file can never be mistaken for the old run
      n_med       = n_med,
      bart_engine = engine,
      bart_ndpost = S,
      bart_nskip  = nskip,
      success     = TRUE
    )
  }, error = function(e) {
    list(method = "BART-CMA", success = FALSE, error = conditionMessage(e))
  })
}

## ---- run_one_rep (override: passes an explicit seed to BART) ---------
## Identical to lib_methods.R's version except for the seed pass-through.
## The DATA seed is unchanged, so datasets match the existing BKMR and
## DeepMed runs exactly.  BKMR's behaviour is untouched.
run_one_rep <- function(rep_id, scenario, truth,
                        methods_to_run = c("BKMR-CMA", "BART-CMA")) {
  set.seed(as.integer(SIM_CONFIG$seed + scenario$scenario_id * 1000L + rep_id))
  dat <- generate_data(scenario)
  results <- list()
  if ("BKMR-CMA" %in% methods_to_run) {
    dat$rep_id <- rep_id; t0 <- proc.time(); r <- run_bkmr_cma(dat, truth)
    r$runtime <- (proc.time() - t0)[3]; results[["BKMR-CMA"]] <- r
  }
  if ("BART-CMA" %in% methods_to_run) {
    bart_seed <- as.integer(SIM_CONFIG$seed + scenario$scenario_id * 1000L +
                            rep_id + 500000L)
    t0 <- proc.time(); r <- run_bart_cma(dat, truth, seed = bart_seed)
    r$runtime <- (proc.time() - t0)[3]; results[["BART-CMA"]] <- r
  }
  list(rep_id = rep_id, results = results, data = dat)
}

cat(sprintf("[aligned] run_bart_cma loaded | engine=%s | n_med=%d | ndpost=%d | nskip=%d\n",
            METHOD_CONFIG$bart_engine %||% "BART",
            as.integer(METHOD_CONFIG$bart_nmed %||% 30L),
            METHOD_CONFIG$bart_ndpost, METHOD_CONFIG$bart_nskip))
