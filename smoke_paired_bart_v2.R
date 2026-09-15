#!/usr/bin/env Rscript
###########################################################
## smoke_paired_bart.R  (v2)
##
## Fixes over v1:
##   - results are nested: res[[i]]$results[["BART-CMA"]]$NIE_est
##   - old raw file uses flattened names: results.BART.CMA.NIE_est
##   - saves the raw list IMMEDIATELY after compute, before any
##     summarising, so a formatting bug can never cost the run again
##   - if the raw file already exists, skips compute and re-summarises
##
## Usage:  Rscript R/smoke_paired_bart.R 12 8
###########################################################

suppressPackageStartupMessages(library(parallel))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
SID  <- if (length(args) >= 1) as.integer(args[1]) else 12L
NREP <- if (length(args) >= 2) as.integer(args[2]) else 8L

source("R/lib_setup.R"); source("R/lib_dgp.R"); source("R/lib_truth.R")
source("R/lib_methods.R"); source("R/run_bart_cma_aligned.R")

rawfile <- file.path(RESULTS_DIR, sprintf("smoke_paired_s%02d_raw.rds", SID))

sc <- as.list(design_grid[design_grid$scenario_id == SID, ])
tf <- file.path(TRUTH_DIR, sprintf("truth_%03d.rds", SID))
truth <- readRDS(tf)
methods_to_run <- "BART-CMA"

cat(sprintf("\nscenario %d | reps %d | ndpost %s | S_use %s | n_med_paired %s\n",
            SID, NREP, METHOD_CONFIG$bart_ndpost,
            METHOD_CONFIG$bart_Suse %||% "UNSET",
            METHOD_CONFIG$bart_nmed_paired %||% "UNSET"))
cat(sprintf("truth: TE %+.4f  NDE %+.4f  NIE %+.4f\n\n",
            truth$TE, truth$NDE, truth$NIE))

## ---------------- compute (or reuse) ----------------
if (file.exists(rawfile)) {
  cat("raw results already present, skipping compute: ", rawfile, "\n\n", sep = "")
  res <- readRDS(rawfile)
} else {
  cat(sprintf("running %d reps on %d cores...\n", NREP, NREP))
  t0 <- Sys.time()
  res <- mclapply(seq_len(NREP), function(rep_id) {
    tr  <- Sys.time()
    out <- tryCatch(run_one_rep(rep_id, sc, truth, methods_to_run),
                    error = function(e) list(error = conditionMessage(e)))
    attr(out, "secs") <- as.numeric(difftime(Sys.time(), tr, units = "secs"))
    out
  }, mc.cores = NREP)
  saveRDS(res, rawfile)                       # SAVE FIRST
  cat(sprintf("compute done in %.1f min; raw saved to %s\n\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), rawfile))
}

secs <- vapply(res, function(x) attr(x, "secs") %||% NA_real_, numeric(1))
if (any(is.finite(secs)))
  cat(sprintf("per-rep sec: median %.0f (min %.0f / max %.0f)\n",
              median(secs, na.rm = TRUE), min(secs, na.rm = TRUE),
              max(secs, na.rm = TRUE)))

## ---------------- extraction ----------------
FLD <- c("NIE_est","NIE_lo","NIE_hi","NDE_est","NDE_lo","NDE_hi",
         "TE_est","TE_lo","TE_hi")

pick <- function(x, rep_id) {
  if (!is.null(x$error))
    return(data.frame(rep_id = rep_id, err = x$error, stringsAsFactors = FALSE))
  b <- x$results[["BART-CMA"]] %||% x$results[[1]]
  if (is.null(b))
    return(data.frame(rep_id = rep_id, err = "no BART-CMA element",
                      stringsAsFactors = FALSE))
  d <- as.data.frame(lapply(FLD, function(f) as.numeric(b[[f]] %||% NA_real_)))
  names(d) <- FLD
  cbind(data.frame(rep_id = rep_id), d)
}

new <- do.call(rbind, lapply(seq_along(res), function(i) pick(res[[i]], i)))
if ("err" %in% names(new)) {
  cat("\n!!! ERRORS !!!\n"); print(new); quit(status = 1)
}

## ---------------- old file, same reps ----------------
oldf <- file.path(RESULTS_DIR, sprintf("raw_scenario_%03d_bart.rds", SID))
old <- NULL
if (file.exists(oldf)) {
  o <- readRDS(oldf)
  if (!is.data.frame(o)) o <- do.call(rbind, lapply(o, as.data.frame))
  names(o) <- sub("^results\\.BART[._]CMA\\.", "", names(o))
  keep <- intersect(c("rep_id", FLD), names(o))
  o <- o[, keep, drop = FALSE]
  if ("rep_id" %in% names(o)) o <- o[o$rep_id %in% seq_len(NREP), , drop = FALSE]
  else o <- head(o, NREP)
  for (f in FLD) if (f %in% names(o)) o[[f]] <- as.numeric(o[[f]])
  old <- o
}

summ <- function(d, lab) {
  if (is.null(d) || !nrow(d)) return(NULL)
  data.frame(version = lab, n = nrow(d),
    NIE_est = mean(d$NIE_est, na.rm = TRUE),
    NIE_bias = mean(d$NIE_est, na.rm = TRUE) - truth$NIE,
    NIE_width = mean(d$NIE_hi - d$NIE_lo, na.rm = TRUE),
    NDE_width = mean(d$NDE_hi - d$NDE_lo, na.rm = TRUE),
    TE_width  = mean(d$TE_hi  - d$TE_lo,  na.rm = TRUE),
    NIE_cover = mean(d$NIE_lo <= truth$NIE & truth$NIE <= d$NIE_hi, na.rm = TRUE),
    stringsAsFactors = FALSE)
}

cat("\n=========== COMPARISON (identical datasets) ===========\n")
print(rbind(summ(old, "OLD unpaired"), summ(new, "NEW paired")),
      row.names = FALSE, digits = 4)

if (!is.null(old) && nrow(old)) {
  wo <- mean(old$NIE_hi - old$NIE_lo, na.rm = TRUE)
  wn <- mean(new$NIE_hi - new$NIE_lo, na.rm = TRUE)
  cat(sprintf("\nNIE width ratio new/old : %.3f   (want > 1; expect ~1.2-1.9)\n", wn/wo))
  cat(sprintf("NIE bias shift          : %+.4f  (want ~0)\n",
              mean(new$NIE_est, na.rm = TRUE) - mean(old$NIE_est, na.rm = TRUE)))
}

cat("\n=========== PER-REP (new) ===========\n")
print(new[, c("rep_id","NIE_est","NIE_lo","NIE_hi")], row.names = FALSE, digits = 4)
if (!is.null(old)) {
  cat("\n=========== PER-REP (old) ===========\n")
  print(old[, intersect(c("rep_id","NIE_est","NIE_lo","NIE_hi"), names(old))],
        row.names = FALSE, digits = 4)
}

m <- median(secs, na.rm = TRUE)
if (is.finite(m))
  cat(sprintf("\nPROJECTION 100 reps: %.1f h @40 cores | %.1f h @60 cores\n",
              m*ceiling(100/40)/3600, m*ceiling(100/60)/3600))
cat("\nDONE\n")
