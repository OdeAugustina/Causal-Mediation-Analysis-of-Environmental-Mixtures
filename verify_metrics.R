#!/usr/bin/env Rscript
## verify_metrics.R -- confirm the patch changed nothing for scenarios 1-12
## and works correctly for the plasmode cells 13/14.
##
## Run BEFORE and AFTER patching; the "before" numbers are saved so the
## "after" run can compare automatically.
##
## Usage:  Rscript R/verify_metrics.R

source("R/lib_setup.R"); source("R/lib_methods.R")
SNAP <- "results/metrics_snapshot_prepatch.rds"

grab <- function() {
  out <- list()
  for (s in sprintf("%03d", 1:14)) {
    for (mth in c("bart", "bkmr_pkg")) {
      f <- sprintf("results/raw_scenario_%s_%s.rds", s, mth)
      if (!file.exists(f)) next
      x  <- readRDS(f)
      tf <- sprintf("results/truth/truth_%s.rds", s)
      if (!file.exists(tf)) next
      tr <- readRDS(tf)
      nm <- if (mth == "bart") "BART-CMA" else "BKMR-CMA"
      for (est in c("TE", "NDE", "NIE")) {
        m <- tryCatch(compute_metrics(x, tr, nm, estimand = est),
                      error = function(e) list(valid = FALSE, err = conditionMessage(e)))
        key <- paste(s, mth, est, sep = "|")
        out[[key]] <- if (isTRUE(m$valid))
          c(n = m$n_success, bias = m$bias, rmse = m$rmse, cov = m$coverage,
            width = m$ci_width, mcse = m$mcse %||% NA_real_,
            repspec = as.numeric(isTRUE(m$replicate_specific_truth)))
        else c(n = NA, bias = NA, rmse = NA, cov = NA, width = NA, mcse = NA, repspec = NA)
      }
    }
  }
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a
cur <- grab()

if (!file.exists(SNAP)) {
  saveRDS(cur, SNAP)
  cat("\nPRE-PATCH snapshot saved:", SNAP, "\n")
  cat("cells captured:", length(cur), "\n")
  cat("Now apply the patch, then run this script again.\n\n")
} else {
  old <- readRDS(SNAP)
  cat("\n=== scenarios 1-12: bias / rmse / coverage must be UNCHANGED ===\n")
  bad <- 0
  for (k in names(old)) {
    if (!k %in% names(cur)) { cat("  MISSING now:", k, "\n"); bad <- bad + 1; next }
    s <- as.integer(strsplit(k, "\\|")[[1]][1])
    if (s > 12) next
    o <- old[[k]]; n <- cur[[k]]
    for (f in c("bias", "rmse", "cov", "width")) {
      if (!isTRUE(all.equal(o[[f]], n[[f]], tolerance = 1e-10))) {
        cat(sprintf("  CHANGED %s %s: %.8f -> %.8f\n", k, f, o[[f]], n[[f]]))
        bad <- bad + 1
      }
    }
  }
  cat(if (bad == 0) "  all identical -- backward compatible\n"
      else sprintf("  %d differences -- INVESTIGATE\n", bad))

  cat("\n=== plasmode cells 13-14: should now use replicate-specific truth ===\n")
  for (k in names(cur)) {
    s <- as.integer(strsplit(k, "\\|")[[1]][1])
    if (s < 13) next
    n <- cur[[k]]
    cat(sprintf("  %-22s n=%3.0f bias %+8.5f cov %.2f width %.4f rep-specific=%s\n",
        k, n[["n"]], n[["bias"]], n[["cov"]], n[["width"]],
        if (isTRUE(n[["repspec"]] == 1)) "YES" else "no"))
  }
  cat("\nTE and NDE should show rep-specific=YES; NIE may show either\n")
  cat("(its truth is constant, so both paths give the same answer).\n\n")
}
