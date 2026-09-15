#!/usr/bin/env Rscript
###########################################################
## run_truth.R — Compute truth for ONE scenario
## Reads SLURM_ARRAY_TASK_ID (or first commandline arg).
## Faithful to Chunk 9 of original Rmd, parallelized per-scenario.
###########################################################

## Working directory must be the project root (SLURM scripts handle this via `cd`).
source("R/lib_setup.R")
source("R/lib_dgp.R")
source("R/lib_truth.R")

sid <- get_scenario_id()
sc  <- as.list(design_grid[design_grid$scenario_id == sid, ])
if (length(sc$scenario_id) == 0) stop(sprintf("Unknown scenario_id: %d", sid))

tf <- file.path(TRUTH_DIR, sprintf("truth_%03d.rds", sid))
if (file.exists(tf)) {
  cat(sprintf("[truth] Scenario %d: cached at %s. Skipping.\n", sid, tf))
  quit(save = "no", status = 0)
}

cat(sprintf("[truth] Scenario %d (%s/%s/%s/n=%d/L=%d) computing...\n",
            sid, sc$linearity, sc$interaction, sc$outcome_type, sc$n, sc$L))
t0 <- Sys.time()
truth <- compute_truth(sc)
saveRDS(truth, tf)
cat(sprintf("[truth] Scenario %d: TE=%.4f NDE=%.4f NIE=%.4f (wall=%s)\n",
            sid, truth$TE, truth$NDE, truth$NIE, format(Sys.time() - t0)))
