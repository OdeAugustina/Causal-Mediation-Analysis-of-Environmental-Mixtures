#!/usr/bin/env python3
"""
patch_pairing.py -- replace the unpaired Step 3 block in
R/run_bart_cma_aligned.R with the paired-draw version, and add the
S_use / n_med_paired knobs to METHOD_CONFIG in R/lib_setup.R.

Anchored on text, not line numbers. Refuses to run twice.
Run from the project root:  python3 patch_pairing.py
"""
import re, shutil, sys, os, time

BART = "R/run_bart_cma_aligned.R"
SETUP = "R/lib_setup.R"
STAMP = time.strftime("%Y-%m-%d_%H%M%S")

NEW_STEP3 = '''    ## ---- Step 3: paired-draw Monte-Carlo g-computation ----
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
'''

def die(msg):
    print("ABORT: " + msg)
    sys.exit(1)

for f in (BART, SETUP):
    if not os.path.exists(f):
        die("%s not found. Run from the project root." % f)

# ---------------- 1. Step 3 block ----------------
src = open(BART).read()

if "n_med_p" in src or "bart_nmed_paired" in src:
    die("%s already contains the paired-draw code. Nothing to do." % BART)

lines = src.split("\n")
start = end = None
for i, ln in enumerate(lines):
    if start is None and "---- Step 3:" in ln:
        start = i
    if start is not None and re.match(r"\s*Y_hl\s*<-\s*acc_hl\s*/\s*n_med\s*$", ln):
        end = i
        break

if start is None:
    die("could not find the '---- Step 3:' anchor.")
if end is None:
    die("found Step 3 header at line %d but no 'Y_hl <- acc_hl / n_med' after it." % (start + 1))

block = "\n".join(lines[start:end + 1])
if "sample.int" not in block:
    die("block at lines %d-%d has no sample.int -- wrong block, refusing." % (start + 1, end + 1))
if block.count("acc_hh") < 3 or "acc_hl" not in block:
    die("block at lines %d-%d does not look like the accumulator loop; refusing."
        % (start + 1, end + 1))

shutil.copy2(BART, "%s.bak.%s" % (BART, STAMP))
out = lines[:start] + NEW_STEP3.rstrip("\n").split("\n") + lines[end + 1:]
open(BART, "w").write("\n".join(out))
print("PATCHED %s : replaced lines %d-%d (%d lines -> %d lines)"
      % (BART, start + 1, end + 1, end - start + 1, len(NEW_STEP3.rstrip("\n").split("\n"))))

# ---------------- 2. METHOD_CONFIG knobs ----------------
s2 = open(SETUP).read()
if "bart_Suse" in s2:
    print("SKIP %s : bart_Suse already present." % SETUP)
else:
    m = re.search(r"^(\s*)bart_nmed\s*=\s*30\s*,\s*bart_engine\s*=\s*\"BART\"\s*,\s*$",
                  s2, re.M)
    if not m:
        die("could not find the 'bart_nmed = 30, bart_engine = \"BART\",' line in %s."
            % SETUP)
    ind = m.group(1)
    new = (m.group(0) + "\n" + ind +
           "bart_Suse = 200L, bart_nmed_paired = 5L,")
    s2 = s2[:m.start()] + new + s2[m.end():]

    n_sub = 0
    s2, n_sub = re.subn(r"bart_ndpost\s*=\s*2000", "bart_ndpost = 1000", s2)
    if n_sub != 1:
        die("expected exactly one 'bart_ndpost = 2000', found %d." % n_sub)

    shutil.copy2(SETUP, "%s.bak.%s" % (SETUP, STAMP))
    open(SETUP, "w").write(s2)
    print("PATCHED %s : added bart_Suse/bart_nmed_paired, ndpost 2000 -> 1000" % SETUP)

print("\nBackups written with suffix .bak.%s" % STAMP)
