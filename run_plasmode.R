#!/usr/bin/env Rscript
###########################################################
## run_plasmode.R -- standalone runner for scenarios 13 / 14
##
## Deliberately does NOT touch design_grid, run_scenario.R, or
## run_one_rep.  Cells 13/14 are produced independently and merged at
## aggregation time.  Output file names and object shapes match the
## existing pipeline exactly:
##
##   results/raw_scenario_0{13,14}_bart.rds
##   results/raw_scenario_0{13,14}_bkmr_pkg.rds
##   results/data_for_deepmed/scenario_0{13,14}.rds
##
##   list of R elements, each list(rep_id = <int>,
##                                 results = list("<METHOD>" = <fields>))
##
## Usage:
##   METHOD=bart     Rscript R/run_plasmode.R 13
##   METHOD=bkmr_pkg Rscript R/run_plasmode.R 13
##   METHOD=data     Rscript R/run_plasmode.R 13     # DeepMed inputs only
##
## Env:
##   BKMR_USE_KNOTS=TRUE|FALSE   (default TRUE)
##   BKMR_NKNOTS=100
##   SIM_CORES from SLURM_CPUS_PER_TASK
###########################################################

suppressPackageStartupMessages(library(parallel))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
SID  <- if (length(args) >= 1) as.integer(args[1]) else
        as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "13"))
METHOD <- tolower(Sys.getenv("METHOD", "bart"))
if (!METHOD %in% c("bart", "bkmr_pkg", "data"))
  stop("METHOD must be bart | bkmr_pkg | data")
if (!SID %in% c(13L, 14L)) stop("run_plasmode.R handles scenarios 13 and 14 only")

source("R/lib_setup.R")
source("R/lib_dgp.R")
source("R/lib_truth.R")
source("R/lib_methods.R")
source("R/run_bart_cma_aligned.R")          # patched paired-draw estimator
source("R/lib_dgp_plasmode.R")
source("R/knots_helper.R")
if (METHOD == "bkmr_pkg") {
  source("R/lib_methods_bkmr_pkg.R")        # includes astar fix + knots
}

R_REPS <- SIM_CONFIG$R
sp     <- pm_load_spec(SID)
truth  <- readRDS(file.path(TRUTH_DIR, sprintf("truth_%03d.rds", SID)))
stopifnot(isTRUE(truth$replicate_specific), length(truth$NIE_rep) >= R_REPS)

method_name <- if (METHOD == "bart") "BART-CMA" else "BKMR-CMA"
outfile <- file.path(RESULTS_DIR, sprintf("raw_scenario_%03d_%s.rds", SID, METHOD))
dmfile  <- file.path(DATA_DIR,    sprintf("scenario_%03d.rds", SID))

cat(sprintf("\n=== plasmode scenario %d | METHOD=%s | reps=%d | cores=%d ===\n",
            SID, METHOD, R_REPS, SIM_CORES))
cat(sprintf("cell: %s | truth TE %+.5f  NDE %+.5f  NIE %+.5f (NIE sd %.2e)\n",
            sp$cell$label, truth$TE, truth$NDE, truth$NIE,
            stats::sd(truth$NIE_rep)))
if (METHOD == "bkmr_pkg")
  cat(sprintf("BKMR knots: %s (n_knots=%d)\n",
              KNOTS_CFG$enabled, KNOTS_CFG$n_knots))

## ---- one replicate -------------------------------------------------
one_rep <- function(rep_id) {
  set.seed(as.integer(SIM_CONFIG$seed + SID * 1000L + rep_id))
  d <- generate_data_plasmode(rep_id, SID, sp)

  ## replicate-specific truth, passed in the shape the methods expect
  tr_r <- list(TE = truth$TE_rep[rep_id], NDE = truth$NDE_rep[rep_id],
               NIE = truth$NIE_rep[rep_id], a = truth$a, astar = truth$astar,
               scenario_id = SID)

  dat <- list(Y = d$Y, M = d$M, Z = d$Z, C = d$C,
              n = d$n, L = d$L, rep_id = rep_id, scenario = d$scenario)

  out <- if (METHOD == "bart") {
    tryCatch(run_bart_cma(dat, tr_r,
               seed = as.integer(SIM_CONFIG$seed + SID*1000L + rep_id)),
             error = function(e) list(method = "BART-CMA", success = FALSE,
                                      error = conditionMessage(e)))
  } else if (METHOD == "bkmr_pkg") {
    tryCatch(run_bkmr_cma(dat, tr_r),
             error = function(e) list(method = "BKMR-CMA", success = FALSE,
                                      error = conditionMessage(e)))
  } else NULL

  list(rep_id = rep_id,
       results = if (is.null(out)) list() else setNames(list(out), method_name),
       data = if (METHOD == "data") dat else NULL)
}

## ---- DeepMed inputs only -------------------------------------------
if (METHOD == "data") {
  dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
  dl <- lapply(seq_len(R_REPS), function(r) {
    set.seed(as.integer(SIM_CONFIG$seed + SID * 1000L + r))
    d <- generate_data_plasmode(r, SID, sp)
    list(Y = d$Y, M = d$M, Z = d$Z, C = d$C, rep_id = r,
         n = d$n, L = d$L, scenario = d$scenario)
  })
  saveRDS(dl, dmfile)
  cat("wrote DeepMed inputs: ", dmfile, "\n", sep = "")
  cat("NOTE: verify this matches the existing format, e.g.\n",
      "  Rscript -e 'str(readRDS(\"results/data_for_deepmed/scenario_001.rds\")[[1]], max.level=1)'\n")
  quit(save = "no", status = 0)
}

## ---- run ------------------------------------------------------------
cl <- makeCluster(SIM_CORES, type = "PSOCK")
on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(MASS); library(dbarts); library(BART)
    library(bkmr); library(causalbkmr); library(fields) })
})
clusterExport(cl, c("SID", "METHOD", "sp", "truth", "method_name",
                    "SIM_CONFIG", "METHOD_CONFIG", "one_rep", "%||%"),
              envir = environment())
clusterEvalQ(cl, {
  source("R/lib_setup.R"); source("R/lib_dgp.R"); source("R/lib_truth.R")
  source("R/lib_methods.R"); source("R/run_bart_cma_aligned.R")
  source("R/lib_dgp_plasmode.R")
  source("R/knots_helper.R")
  if (METHOD == "bkmr_pkg") source("R/lib_methods_bkmr_pkg.R")
  NULL
})

t0 <- Sys.time()
res <- parLapply(cl, seq_len(R_REPS), one_rep)
el <- as.numeric(difftime(Sys.time(), t0, units = "hours"))

nfail <- sum(vapply(res, function(r) {
  x <- r$results[[method_name]]
  is.null(x) || !isTRUE(x$success) }, logical(1)))

saveRDS(res, outfile)
cat(sprintf("\n[scenario %d | %s] done. fails=%d/%d wall=%.2f hours\n",
            SID, METHOD, nfail, R_REPS, el))
cat("wrote: ", outfile, "\n", sep = "")

## ---- quick read ------------------------------------------------------
e <- sapply(res, function(r) {
  x <- r$results[[method_name]]
  if (is.null(x) || !isTRUE(x$success)) NA_real_ else as.numeric(x$NIE_est) })
ok <- is.finite(e)
if (sum(ok) > 1) {
  tv <- truth$NIE_rep[which(ok)]
  cat(sprintf("NIE: mean %+.5f | truth %+.5f | bias %+.5f | MCSE %.5f | n=%d\n",
              mean(e[ok]), mean(tv), mean(e[ok] - tv),
              stats::sd(e[ok] - tv)/sqrt(sum(ok)), sum(ok)))
}
