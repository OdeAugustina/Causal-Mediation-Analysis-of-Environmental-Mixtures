#!/usr/bin/env python3
"""
patch_applied.py -- fix the mediator-draw pairing in the applied analysis.

THE DEFECT (Step 3 of bart_cma_mixture):
    r <- sample.int(S, 1)
draws ONE mediator index per loop pass and applies it to ALL S posterior
rows. Every row therefore sees the same mediator values, so mediator
uncertainty is averaged out of the interval instead of propagating into it.
Intervals come out too narrow.

THE FIX: pair mediator draw s with outcome draw s, and average the residual
noise WITHIN each draw. Survey weighting via wavg() is preserved unchanged.

Run from the folder holding the .Rmd:  python3 patch_applied.py <file.Rmd>
"""
import re, shutil, sys, os, time

f = sys.argv[1] if len(sys.argv) > 1 else "BART_CMA_Mediation__2_.Rmd"
if not os.path.exists(f):
    print("ABORT: %s not found" % f); sys.exit(1)
s = open(f).read()
if "n_med_paired" in s:
    print("ABORT: already patched"); sys.exit(1)

old_start = s.find("  ## Step 3: vectorised Monte-Carlo g-computation")
old_end   = s.find("  Y_hl <- acc_hl / n_med_draws")
if old_start < 0 or old_end < 0:
    print("ABORT: could not find the Step 3 block"); sys.exit(1)
old_end += len("  Y_hl <- acc_hl / n_med_draws")

block = s[old_start:old_end]
if "sample.int" not in block:
    print("ABORT: block has no sample.int -- wrong block"); sys.exit(1)

new = '''  ## Step 3: paired-draw Monte-Carlo g-computation
  ## Mediator draw s is paired with outcome draw s; residual realisations
  ## are averaged within the draw. The previous version drew one mediator
  ## index per k and shared it across all S rows, which averaged mediator
  ## uncertainty out of the posterior and produced intervals that were too
  ## narrow. Survey weighting (wavg) is unchanged.
  S_use   <- min(S, 200L)
  idx     <- if (S_use < S) round(seq(1, S, length.out = S_use)) else seq_len(S)
  n_med_paired <- 5L

  Y_hh <- Y_ll <- Y_hl <- matrix(NA_real_, length(idx), n)
  ridx  <- rep(seq_len(n), n_med_paired)
  Zh_r  <- Z_at_high[ridx, , drop = FALSE]
  Zl_r  <- Z_at_low[ridx,  , drop = FALSE]
  X_r   <- X[ridx, , drop = FALSE]
  B     <- n_med_paired * n
  cm    <- function(v) colMeans(matrix(v, n_med_paired, n, byrow = TRUE))

  for (j in seq_along(idx)) {
    sdx    <- idx[j]
    m_low  <- as.vector(t(matrix(M_low_hat[sdx, ],  n_med_paired, n, byrow = TRUE) +
                          matrix(rnorm(B, 0, sig[sdx]), n_med_paired, n)))
    m_high <- as.vector(t(matrix(M_high_hat[sdx, ], n_med_paired, n, byrow = TRUE) +
                          matrix(rnorm(B, 0, sig[sdx]), n_med_paired, n)))
    Y_hh[j, ] <- cm(predict_Y(cbind(Zh_r, M = m_high, X_r))[sdx, ])
    Y_ll[j, ] <- cm(predict_Y(cbind(Zl_r, M = m_low,  X_r))[sdx, ])
    Y_hl[j, ] <- cm(predict_Y(cbind(Zh_r, M = m_low,  X_r))[sdx, ])
  }'''

s = s[:old_start] + new + s[old_end:]
shutil.copy2(f, f + ".bak." + time.strftime("%Y%m%d_%H%M%S"))
open(f, "w").write(s)
print("PATCHED", f)
print("  Step 3 replaced. wavg(), TE/NDE/NIE/PM and pooling unchanged.")
