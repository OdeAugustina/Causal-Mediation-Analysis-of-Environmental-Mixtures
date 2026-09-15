#!/usr/bin/env Rscript
###########################################################
## knots_helper.R (v2) -- Gaussian predictive process knots for BKMR
##
## WHY: at n = 1000 the exact n x n kernel is numerically singular for
## binary outcomes ("H must be positive definite", 100/100 replicates
## in scenario 10, under both varsel settings and both logit and probit
## DGPs).  The predictive process reduces the kernel to n x m, m << n.
##
## DESIGN: n_knots = 100 fixed for EVERY BKMR-CMA fit in the grid, so
## the approximation is not confounded with outcome type, n, or L.
##
## KERNEL SPACES -- knots must live in the same space as the kernel:
##   fit.m      kernel = Z            -> knots from Z
##   fit.y.TE   kernel = Z            -> knots from Z
##   fit.y      kernel = cbind(Z, M)  -> knots from cbind(Z, M)
## Building knots on Z and pinning M at a constant is INVALID for
## fit.y: it places every knot on a single mediator plane and cannot
## represent the joint exposure-mediator surface.
##
## REPRODUCIBILITY: cover.design() is stochastic.  The seed is derived
## from (scenario_id, rep_id, kernel space) so the same replicate always
## yields the same knot set, and the Z-space and ZM-space designs within
## a replicate are independent.
###########################################################

KNOTS_CFG <- list(
  n_knots   = as.integer(Sys.getenv("BKMR_NKNOTS", "100")),
  seed_base = 77777L,           # offset from the data-generation stream
  enabled   = as.logical(Sys.getenv("BKMR_USE_KNOTS", "TRUE"))
)

## space_id: 1 = Z space (fit.m, fit.y.TE), 2 = cbind(Z, M) space (fit.y)
knots_seed <- function(scenario_id, rep_id, space_id) {
  as.integer(KNOTS_CFG$seed_base + scenario_id * 1000L +
               rep_id * 10L + space_id)
}

make_knots <- function(Zk, scenario_id, rep_id, space_id,
                       n_knots = KNOTS_CFG$n_knots) {
  if (!requireNamespace("fields", quietly = TRUE))
    stop("package 'fields' required for cover.design()")
  Zk <- as.matrix(Zk)
  nk <- min(n_knots, nrow(Zk) - 1L)
  sd <- knots_seed(scenario_id, rep_id, space_id)
  set.seed(sd)
  kn <- as.matrix(fields::cover.design(R = Zk, nd = nk)$design)
  colnames(kn) <- colnames(Zk)
  attr(kn, "knots_seed") <- sd
  attr(kn, "n_knots")    <- nk
  attr(kn, "space_dim")  <- ncol(Zk)
  kn
}

knots_info <- function(kn, label) {
  if (is.null(kn)) return(sprintf("%s: exact kernel (no knots)", label))
  sprintf("%s: %d knots x %d dims, seed %d",
          label, nrow(kn), ncol(kn), attr(kn, "knots_seed"))
}

cat(sprintf("[knots] helper v2 | n_knots=%d | enabled=%s\n",
            KNOTS_CFG$n_knots, KNOTS_CFG$enabled))
