#!/usr/bin/env Rscript
###########################################################
## smoke_paired_bart.R -- validate the paired-draw Step 3 fix.
##
## Runs a few reps of ONE scenario with the patched estimator and
## compares them against the SAME reps in the existing raw results.
## Seeding is seed + scenario_id*1000 + rep_id, so rep k here is the
## same dataset as rep k in the old file -- a paired comparison.
##
## Usage (from project root, on a COMPUTE node):
##   Rscript R/smoke_paired_bart.R 12 8
##      arg1 = scenario id (default 12)
##      arg2 = number of reps (default 8), also the core count
###########################################################

suppressPackageStartupMessages({ library(parallel) })

args   <- commandArgs(trailingOnly = TRUE)
SID    <- if (length(args) >= 1) as.integer(args[1]) else 12L
NREP   <- if (length(args) >= 2) as.integer(args[2]) else 8L

source("R/lib_setup.R")
source("R/lib_dgp.R")
source("R/lib_truth.R")
source("R/lib_methods.R")
source("R/run_bart_cma_aligned.R")

`%||%` <- function(a, b) if (is.null(a)) b else a

cat("\n=========== CONFIG ===========\n")
cat(sprintf("scenario   : %d\n", SID))
cat(sprintf("reps       : %d\n", NREP))
cat(sprintf("ndpost     : %s\n", METHOD_CONFIG$bart_ndpost))
cat(sprintf("nskip      : %s\n", METHOD_CONFIG$bart_nskip))
cat(sprintf("S_use      : %s\n", METHOD_CONFIG$bart_Suse %||% "UNSET"))
cat(sprintf("n_med_pair : %s\n", METHOD_CONFIG$bart_nmed_paired %||% "UNSET"))
if (is.null(METHOD_CONFIG$bart_Suse))
  stop("bart_Suse not in METHOD_CONFIG -- lib_setup.R patch did not apply.")

cat("\n=========== SIGNATURES ===========\n")
cat("run_one_rep : "); print(args(run_one_rep))
cat("run_bart_cma: "); print(args(run_bart_cma))

## ---- scenario + truth, exactly as run_scenario.R does ----
sc <- as.list(design_grid[design_grid$scenario_id == SID, ])
tf <- file.path(TRUTH_DIR, sprintf("truth_%03d.rds", SID))
if (!file.exists(tf)) stop("missing truth file: ", tf)
truth <- readRDS(tf)
methods_to_run <- "BART-CMA"

cat("\n=========== TRUTH ===========\n")
str(truth)

## ---- run the reps ----
cat(sprintf("\n=========== RUNNING %d REPS (%d cores) ===========\n", NREP, NREP))
t0 <- Sys.time()
res <- mclapply(seq_len(NREP), function(rep_id) {
  tr <- Sys.time()
  out <- tryCatch(run_one_rep(rep_id),
                  error = function(e) list(error = conditionMessage(e)))
  attr(out, "secs") <- as.numeric(difftime(Sys.time(), tr, units = "secs"))
  out
}, mc.cores = NREP)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

secs <- vapply(res, function(x) attr(x, "secs") %||% NA_real_, numeric(1))
cat(sprintf("\nwall elapsed        : %.1f min (%d reps in parallel)\n", elapsed, NREP))
cat(sprintf("per-rep sec (median): %.0f  (min %.0f / max %.0f)\n",
            median(secs, na.rm = TRUE), min(secs, na.rm = TRUE), max(secs, na.rm = TRUE)))

cat("\n=========== ONE RAW RESULT ===========\n")
str(res[[1]], max.level = 2)

## ---- flatten to a data frame ----
to_df <- function(x, rep_id) {
  if (is.data.frame(x)) { x$rep_id <- rep_id; return(x) }
  if (!is.null(x$error)) return(data.frame(rep_id = rep_id, error = x$error))
  d <- as.data.frame(x[!vapply(x, function(e) length(e) != 1, logical(1))],
                     stringsAsFactors = FALSE)
  d$rep_id <- rep_id
  d
}
new <- do.call(rbind, lapply(seq_along(res), function(i) to_df(res[[i]], i)))
if ("error" %in% names(new)) {
  cat("\n!!! ERRORS !!!\n"); print(new[!is.na(new$error), ]); quit(status = 1)
}
new <- new[grepl("BART", new$method %||% "BART-CMA"), , drop = FALSE]

## ---- pull the same reps from the OLD file ----
oldf <- file.path(RESULTS_DIR, sprintf("raw_scenario_%03d_bart.rds", SID))
cat("\n=========== OLD FILE ===========\n")
if (file.exists(oldf)) {
  old <- readRDS(oldf)
  if (!is.data.frame(old)) old <- do.call(rbind, lapply(old, as.data.frame))
  cat("columns: ", paste(names(old), collapse = ", "), "\n")
  old <- old[grepl("BART", old$method), , drop = FALSE]
  if ("rep_id" %in% names(old)) old <- old[old$rep_id %in% seq_len(NREP), , drop = FALSE]
  else old <- head(old, NREP)
} else { cat("not found: ", oldf, "\n"); old <- NULL }

## ---- compare ----
summ <- function(d, lab) {
  if (is.null(d) || !nrow(d)) return(NULL)
  w <- d$NIE_hi - d$NIE_lo
  data.frame(
    version   = lab,
    n         = nrow(d),
    NIE_mean  = mean(d$NIE_est, na.rm = TRUE),
    NIE_bias  = mean(d$NIE_est, na.rm = TRUE) - (truth$NIE %||% NA_real_),
    NIE_width = mean(w, na.rm = TRUE),
    NDE_width = mean(d$NDE_hi - d$NDE_lo, na.rm = TRUE),
    TE_width  = mean(d$TE_hi  - d$TE_lo,  na.rm = TRUE),
    stringsAsFactors = FALSE)
}

cat("\n=========== COMPARISON (same datasets) ===========\n")
cmp <- rbind(summ(old, "OLD unpaired"), summ(new, "NEW paired"))
print(cmp, row.names = FALSE, digits = 4)

if (!is.null(old) && nrow(old)) {
  wo <- mean(old$NIE_hi - old$NIE_lo, na.rm = TRUE)
  wn <- mean(new$NIE_hi - new$NIE_lo, na.rm = TRUE)
  cat(sprintf("\nNIE width ratio new/old : %.2f\n", wn / wo))
  cat(sprintf("NIE bias shift          : %+.4f\n",
              mean(new$NIE_est, na.rm = TRUE) - mean(old$NIE_est, na.rm = TRUE)))
}

cat("\n=========== PER-REP DETAIL (new) ===========\n")
print(new[, intersect(c("rep_id","NIE_est","NIE_lo","NIE_hi"), names(new))],
      row.names = FALSE, digits = 4)

cat(sprintf("\nPROJECTION: 100 reps / 40 cores = %.1f h ; / 60 cores = %.1f h\n",
            median(secs, na.rm = TRUE) * ceiling(100/40) / 3600,
            median(secs, na.rm = TRUE) * ceiling(100/60) / 3600))
cat("\nDONE\n")
