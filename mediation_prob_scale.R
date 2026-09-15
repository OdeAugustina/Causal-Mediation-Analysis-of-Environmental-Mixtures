###########################################################
## mediation_prob_scale.R -- probability-scale BKMR-CMA for binary Y
##
## PROBLEM
##   The DGP draws Y ~ Bernoulli(inv_link(eta)), inv_link = pnorm(eta/1.7),
##   and lib_truth.R computes truth as differences in inv_link(cYm(.)) --
##   i.e. PROBABILITY differences.
##
##   kmbayes(family = "binomial") uses probit data augmentation:
##       Y* = h(Z) + X'beta + e ,  Y = 1(Y* > 0)
##   bkmr::SamplePred() returns h + X'beta, the LATENT value, and
##   mediation.bkmr() differences those.  So stock BKMR-CMA reports
##   LATENT LIABILITY differences.
##
##   Measured consequence (scenario 9, 100 reps):
##       TE  +1.476 / truth +0.395 = 3.74
##       NDE +0.749 / truth +0.199 = 3.76
##       NIE +0.727 / truth +0.195 = 3.72
##   A uniform multiplicative factor across all three estimands is the
##   signature of a scale mismatch, not estimation bias.  BART-CMA is
##   unaffected: it applies pnorm() to the outcome predictions BEFORE
##   contrasting (run_bart_cma_aligned.R), giving TE ratios 0.91 / 1.03.
##
## FIX
##   Apply pnorm() to each posterior prediction BEFORE differencing, so
##   BKMR-CMA reports probability differences and aligns with both the
##   truth definition and BART-CMA.
##
##   ORDER MATTERS for the cross-world term.  E[Y(a, M(astar))] is
##   E_M[ P(Y = 1 | a, M) ], so pnorm() is applied to each of the K
##   mediator draws and THEN averaged.  pnorm(mean(.)) != mean(pnorm(.)).
##
## SCOPE
##   Binary cells only (9, 10).  Continuous cells use family="gaussian"
##   with an identity link -- no scale to mismatch -- and are delegated
##   to causalbkmr::mediation.bkmr() unchanged.
##
## ASSUMPTION TO VERIFY
##   That SamplePred() on a binomial fit returns the latent scale.  The
##   uniform 3.7x ratio is strong evidence, and verify_prob_scale.R
##   confirms it directly by checking the ratio moves to ~1.
###########################################################

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

.post_summary <- function(v, alpha = 0.05) {
  c(mean   = mean(v, na.rm = TRUE),
    median = stats::median(v, na.rm = TRUE),
    lower  = unname(stats::quantile(v, alpha / 2,     na.rm = TRUE)),
    upper  = unname(stats::quantile(v, 1 - alpha / 2, na.rm = TRUE)),
    sd     = stats::sd(v, na.rm = TRUE))
}

mediation_bkmr_scaled <- function(a, astar, e.m = NULL, e.y = NULL,
                                  fit.m, fit.y, fit.y.TE,
                                  X.predict.M, X.predict.Y,
                                  alpha = 0.05, sel, seed, K,
                                  binary_y = FALSE) {

  ## ---- continuous: delegate, unchanged (astar patch already applied) ----
  if (!isTRUE(binary_y)) {
    return(causalbkmr::mediation.bkmr(
      a = a, astar = astar, e.m = e.m, e.y = e.y,
      fit.m = fit.m, fit.y = fit.y, fit.y.TE = fit.y.TE,
      X.predict.M = X.predict.M, X.predict.Y = X.predict.Y,
      alpha = alpha, sel = sel, seed = seed, K = K))
  }

  ## ---- binary: probability scale ----
  set.seed(seed)
  z.y  <- c(a, e.y)
  z.ys <- c(astar, e.y)

  ## Ya and Yastar from the total-effect fit; pnorm BEFORE differencing
  Ya_lat  <- as.vector(bkmr::SamplePred(fit.y.TE, Znew = z.y,
                                        Xnew = X.predict.Y, sel = sel))
  Yas_lat <- as.vector(bkmr::SamplePred(fit.y.TE, Znew = z.ys,
                                        Xnew = X.predict.Y, sel = sel))
  Ya  <- stats::pnorm(Ya_lat)
  Yas <- stats::pnorm(Yas_lat)

  ## cross-world: mediator drawn under astar (the astar bug, fixed here too)
  Mastar <- as.vector(bkmr::SamplePred(fit.m, Znew = c(astar, e.m),
                                       Xnew = X.predict.M, sel = sel))
  sig    <- sqrt(fit.m$sigsq.eps[sel])
  Ms     <- Mastar + sig * matrix(stats::rnorm(length(sel) * K),
                                  nrow = length(sel), ncol = K)

  YaMas <- numeric(length(sel))
  for (j in seq_along(sel)) {
    nd <- cbind(matrix(z.y, nrow = K, ncol = length(z.y), byrow = TRUE), Ms[j, ])
    lat <- as.vector(bkmr::SamplePred(fit.y, Znew = nd,
                                      Xnew = X.predict.Y, sel = sel[j]))
    YaMas[j] <- mean(stats::pnorm(lat))     ## pnorm INSIDE, then average
  }

  TE.samp  <- Ya    - Yas
  NDE.samp <- YaMas - Yas
  NIE.samp <- Ya    - YaMas

  est <- rbind(TE  = .post_summary(TE.samp,  alpha),
               NDE = .post_summary(NDE.samp, alpha),
               NIE = .post_summary(NIE.samp, alpha))
  colnames(est) <- c("mean", "median", "lower", "upper", "sd")

  list(est = est, TE.samp = TE.samp, NDE.samp = NDE.samp, NIE.samp = NIE.samp,
       scale = "probability", Ya = mean(Ya), Yastar = mean(Yas),
       YaMastar = mean(YaMas))
}

cat("[scale] mediation_bkmr_scaled() loaded (probability scale for binary Y)\n")
