#!/usr/bin/env Rscript
###########################################################
## run_scenario.R — Run BART-CMA and/or BKMR-CMA for ONE scenario
## Scenario from SLURM_ARRAY_TASK_ID; method from METHOD env
## (bkmr | bart | both; default both). Parallelizes reps
## across SIM_CORES (= SLURM_CPUS_PER_TASK).
###########################################################

source("R/lib_setup.R")
source("R/lib_dgp.R")
source("R/lib_truth.R")
source("R/lib_methods.R")
if (tolower(Sys.getenv("METHOD", "both")) == "bkmr_pkg") source("R/lib_methods_bkmr_pkg.R")
source("R/run_bart_cma_aligned.R")

library(parallel)

`%||%` <- function(a, b) if (is.null(a)) b else a

sid <- get_scenario_id()
sc  <- as.list(design_grid[design_grid$scenario_id == sid, ])
if (length(sc$scenario_id) == 0) stop(sprintf("Unknown scenario_id: %d", sid))

## ---- Which method(s) this job runs ----
method_arg <- tolower(Sys.getenv("METHOD", "both"))
if (!method_arg %in% c("bkmr", "bart", "both", "bkmr_pkg"))
  stop("METHOD must be one of: bkmr, bart, both (got '", method_arg, "')")
methods_to_run <- switch(method_arg,
  bkmr = "BKMR-CMA",
  bkmr_pkg = "BKMR-CMA",
  bart = "BART-CMA",
  both = c("BART-CMA", "BKMR-CMA"))

## ---- Load precomputed truth ----
tf <- file.path(TRUTH_DIR, sprintf("truth_%03d.rds", sid))
if (!file.exists(tf)) {
  stop(sprintf("Truth file missing for scenario %d: %s\n", sid, tf),
       "Run truth job (slurm/01_truth.sh) first.")
}
truth <- readRDS(tf)

## ---- Output paths (method-tagged so split jobs don't clobber) ----
out_tag    <- if (method_arg == "both") "" else paste0("_", method_arg)
outfile    <- file.path(RESULTS_DIR, sprintf("raw_scenario_%03d%s.rds", sid, out_tag))
deepmed_in <- file.path(DATA_DIR,    sprintf("scenario_%03d.rds", sid))
## Only the bart/both job owns the shared per-rep data DeepMed reads.
write_deepmed_in <- method_arg %in% c("bart", "both")

## ---- Skip if already done ----
done <- file.exists(outfile) && (!write_deepmed_in || file.exists(deepmed_in))
if (done) {
  cat(sprintf("[scenario %d | %s] already complete. Skipping.\n", sid, method_arg))
  quit(save = "no", status = 0)
}

cat(sprintf("\n--- Scenario %d [%s]: %s/%s/%s/n=%d/L=%d  reps=%d cores=%d ---\n",
            sid, method_arg, sc$linearity, sc$interaction, sc$outcome_type,
            sc$n, sc$L, SIM_CONFIG$R, SIM_CORES))

## ---- Build PSOCK cluster ----
t0 <- Sys.time()
cl <- makeCluster(SIM_CORES, type = "PSOCK")

clusterExport(cl,
  varlist = c("sc", "truth", "methods_to_run",
              "generate_data", "make_Sigma_Z", "make_Sigma_C", "get_coefs",
              "f_M_nonlinear", "f_Y_nonlinear", "g_ZM_interaction",
              "run_one_rep", "run_bkmr_cma", "run_bart_cma", ".nm", "%||%",
              "DGP_COEFS", "SIM_CONFIG", "METHOD_CONFIG", "inv_link"),
  envir = environment())

clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(MASS); library(dbarts); library(BART)
    library(bkmr); library(causalbkmr)
  })
})

## Data is now seeded deterministically PER REP inside run_one_rep, so
## datasets are identical across methods and core counts. This stream
## seed only governs any residual non-seeded randomness.
if (method_arg == "bkmr_pkg") clusterExport(cl, "BKMR_PKG", envir = environment())
clusterSetRNGStream(cl, iseed = SIM_CONFIG$seed + sid)

## ---- Run all reps in parallel ----
reps <- tryCatch(
  parLapply(cl, 1:SIM_CONFIG$R, function(r) {
    run_one_rep(r, sc, truth, methods_to_run = methods_to_run)
  }),
  error = function(e) {
    cat(sprintf("[scenario %d] FATAL parLapply error: %s\n", sid, conditionMessage(e)))
    NULL
  },
  finally = tryCatch(stopCluster(cl), error = function(e) NULL)
)

if (is.null(reps)) {
  stop(sprintf("Scenario %d failed catastrophically. See log above.", sid))
}

## ---- Per-rep error reporting (only for methods this job ran) ----
err_dir <- file.path(ERRORS_DIR, sprintf("scenario_%03d", sid))
dir.create(err_dir, showWarnings = FALSE, recursive = TRUE)
n_fail_bkmr <- 0; n_fail_bart <- 0
for (r in seq_along(reps)) {
  rep_res <- reps[[r]]$results
  if ("BKMR-CMA" %in% methods_to_run && !isTRUE(rep_res[["BKMR-CMA"]]$success)) {
    n_fail_bkmr <- n_fail_bkmr + 1
    writeLines(
      paste("BKMR-CMA failure for rep", r, ":",
            rep_res[["BKMR-CMA"]]$error %||% "unknown"),
      file.path(err_dir, sprintf("rep_%03d_bkmr.txt", r)))
  }
  if ("BART-CMA" %in% methods_to_run && !isTRUE(rep_res[["BART-CMA"]]$success)) {
    n_fail_bart <- n_fail_bart + 1
    writeLines(
      paste("BART-CMA failure for rep", r, ":",
            rep_res[["BART-CMA"]]$error %||% "unknown"),
      file.path(err_dir, sprintf("rep_%03d_bart.txt", r)))
  }
}

## ---- Save outputs ----
if (write_deepmed_in) saveRDS(reps, deepmed_in)   # full data for DeepMed
reps_light <- reps                                # light copy for metrics
for (r in seq_along(reps_light)) reps_light[[r]]$data <- NULL
saveRDS(reps_light, outfile)

cat(sprintf("[scenario %d | %s] done. BKMR fails=%d BART fails=%d wall=%s\n",
            sid, method_arg, n_fail_bkmr, n_fail_bart, format(Sys.time() - t0)))
