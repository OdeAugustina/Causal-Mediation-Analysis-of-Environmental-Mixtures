###########################################################
## fix_causalbkmr_astar.R
##
## PACKAGE BUG (causalbkmr 0.1.0, commit 4a33c0b):
##   causalbkmr:::YaMastar.SamplePred() is documented and named to compute
##   the cross-world term E[Y(a, M(astar))].  Its body begins
##
##       z.m <- c(a, e.m)
##
##   using the COMPARATIVE exposure a where the REFERENCE exposure astar
##   belongs.  `astar` is accepted as an argument and never used anywhere
##   in the function body.  The counterfactual mediator is therefore drawn
##   under a, not astar, so the function returns E[Y(a, M(a))] = Ya.
##
##   Downstream in mediation.bkmr():
##       NDE <- YaMastar - Yastar   ==>  Ya - Yastar  ==  TE
##       NIE <- Ya       - YaMastar ==>  Ya - Ya      ==  0
##
##   Observed signature: NDE == TE to 3-4 dp and NIE == 0 in every
##   scenario regardless of sample size, nonlinearity, or interaction.
##   TE and the CDEs are UNAFFECTED (different code paths).
##
## FIX: one line, z.m <- c(astar, e.m).  Everything else is the package's
##      own code, reproduced verbatim so behaviour is otherwise identical.
##
## Source this AFTER library(causalbkmr) and INSIDE every parallel worker.
###########################################################

.YaMastar.SamplePred.fixed <- function(a, astar, e.m = NULL, e.y,
                                       fit.m, fit.y,
                                       X.predict.M, X.predict.Y,
                                       sel, seed, K) {
  set.seed(seed)
  z.m <- c(astar, e.m)                 ## <<< THE FIX (package has c(a, e.m))
  EM.samp <- bkmr::SamplePred(fit.m, Znew = z.m, Xnew = X.predict.M, sel = sel)
  Mastar  <- as.vector(EM.samp)
  sigma.samp  <- sqrt(fit.m$sigsq.eps[sel])
  random.samp <- matrix(stats::rnorm(length(sel) * K),
                        nrow = length(sel), ncol = K)
  Mastar.samp <- Mastar + sigma.samp * random.samp
  YaMastar.samp.mat <- matrix(NA_real_, nrow = length(sel), ncol = K)
  z.y <- c(a, e.y)
  for (j in seq_along(sel)) {
    Mastar.j  <- Mastar.samp[j, ]
    aMastar.j <- cbind(matrix(z.y, nrow = K, ncol = length(z.y), byrow = TRUE),
                       Mastar.j)
    YaMastar.j <- bkmr::SamplePred(fit.y, Znew = aMastar.j,
                                   Xnew = X.predict.Y, sel = sel[j])
    YaMastar.samp.mat[j, ] <- as.vector(YaMastar.j)
  }
  apply(YaMastar.samp.mat, 1, mean)
}

## ---- apply the patch -------------------------------------------------
.apply_causalbkmr_fix <- function(verbose = TRUE) {
  if (!requireNamespace("causalbkmr", quietly = TRUE))
    stop("causalbkmr not installed")
  orig <- utils::getFromNamespace("YaMastar.SamplePred", "causalbkmr")
  src  <- paste(deparse(body(orig)), collapse = " ")
  already <- grepl("c(astar, e.m)", src, fixed = TRUE)
  if (already) {
    if (verbose) cat("[fix] causalbkmr already patched or upstream fixed.\n")
    return(invisible(FALSE))
  }
  if (!grepl("c(a, e.m)", src, fixed = TRUE))
    warning("[fix] expected 'z.m <- c(a, e.m)' not found; package version may ",
            "differ. Patch applied anyway -- VERIFY before trusting results.")
  utils::assignInNamespace("YaMastar.SamplePred",
                           .YaMastar.SamplePred.fixed, ns = "causalbkmr")
  if (verbose) cat("[fix] causalbkmr:::YaMastar.SamplePred patched (astar bug).\n")
  invisible(TRUE)
}

.apply_causalbkmr_fix()
