#!/usr/bin/env Rscript
## =====================================================================
## run_deepmed.R  —  production DeepMed driver (one scenario per call)
## =====================================================================
## v2: BATCHED CLUSTER TEARDOWN.
##
## v1 held N_OUTER persistent PSOCK workers for all 100 reps. On GPU,
## each DeepMed() call leaks ~12 GB of HOST RAM that .clear_tf() does
## not reclaim (it frees the TF session, not the process's heap). Ten
## workers x rep 2 => 120 GB => SLURM OUT_OF_MEMORY.
##
## Fix: a worker now lives for exactly ONE rep. We run reps in batches
## of N_OUTER, building a fresh cluster per batch and stopping it after.
## Process death is the only reliable free(). Per-rep part files make
## this cost nothing: work is never repeated, and the job resumes.
##
## Data: reads results/data_for_deepmed/scenario_<ID>.rds — the SAME
## datasets BART/BKMR consumed. No regeneration, no reseeding.
## All DeepMed logic lives in lib_deepmed.R::run_deepmed().
##
## Usage:  Rscript R/run_deepmed.R <SCENARIO_ID>
## Output: results/deepmed/deepmed_<ID>.rds
## =====================================================================

suppressMessages(library(parallel))

PROJ  <- Sys.getenv("SIM_PROJ", "/work/users/o/d/odeaug/sim_proj")
R_DIR <- file.path(PROJ, "R")
setwd(PROJ)

source(file.path(R_DIR, "lib_setup.R"))   # design_grid, SIM_CONFIG, METHOD_CONFIG, dirs

DEEPMED_PY <- Sys.getenv("RETICULATE_PYTHON",
                         "/work/users/o/d/odeaug/envs/deepmed_py/bin/python")

## 5 concurrent reps, not 10: peak host RAM ~= N_OUTER x 12 GB.
N_OUTER <- as.integer(Sys.getenv("DM_OUTER", "5"))

sid <- get_scenario_id()
sc  <- as.list(design_grid[design_grid$scenario_id == sid, ])
if (length(sc$scenario_id) == 0) stop(sprintf("Unknown scenario_id: %d", sid))

in_file  <- file.path(DATA_DIR,    sprintf("scenario_%03d.rds", sid))
out_file <- file.path(DEEPMED_DIR, sprintf("deepmed_%03d.rds",  sid))

if (!file.exists(in_file))
  stop(sprintf("Input data missing for scenario %d: %s", sid, in_file))

if (file.exists(out_file)) {
  cat(sprintf("[scenario %d] DeepMed already complete. Skipping.\n", sid))
  quit(save = "no", status = 0)
}

PARTS <- file.path(DEEPMED_DIR, sprintf("parts_s%03d", sid))
dir.create(PARTS, recursive = TRUE, showWarnings = FALSE)

sim_data <- readRDS(in_file)
R_REPS   <- min(length(sim_data), SIM_CONFIG$R)

cat(sprintf("=== DeepMed | scenario %d (n=%d, L=%d, %s) | %d reps | %d concurrent (batched) ===\n",
            sid, sc$n, sc$L, sc$outcome_type, R_REPS, N_OUTER))
cat(sprintf("[data] %s  (%d reps on disk)\n", in_file, length(sim_data)))
cat(sprintf("[py]   %s\n", DEEPMED_PY))

## ------------------------- per-rep worker -----------------------------
## Runs in a FRESH process that is destroyed immediately afterwards.
worker <- function(r) {
  part <- file.path(PARTS, sprintf("rep_%03d.rds", r))
  if (file.exists(part)) return(TRUE)

  Sys.setenv(TF_CPP_MIN_LOG_LEVEL         = "3",
             TF_USE_LEGACY_KERAS          = "1",
             OMP_NUM_THREADS              = "1",
             RETICULATE_PYTHON            = DEEPMED_PY,
             TF_FORCE_GPU_ALLOW_GROWTH    = "true")

  suppressMessages({
    library(reticulate)
    reticulate::use_python(DEEPMED_PY, required = TRUE)
    library(tensorflow); library(keras); library(DeepMed)
    library(foreach)
  })

  ## Hard GPU assert: fail loudly rather than crawl on CPU.
  gpus <- tryCatch(length(tensorflow::tf$config$list_physical_devices("GPU")),
                   error = function(e) 0L)
  if (gpus < 1 && !nzchar(Sys.getenv("DM_ALLOW_CPU")))
    stop(sprintf("rep %d: TensorFlow sees NO GPU. Refusing to run on CPU.", r))

  try(tensorflow::tf$config$threading$set_intra_op_parallelism_threads(1L), silent = TRUE)
  try(tensorflow::tf$config$threading$set_inter_op_parallelism_threads(1L), silent = TRUE)
  foreach::registerDoSEQ()   # CUDA contexts are NOT fork-safe. Never fork here.

  source(file.path(R_DIR, "lib_setup.R"))    # METHOD_CONFIG
  source(file.path(R_DIR, "lib_deepmed.R"))  # run_deepmed()

  dat <- sim_data[[r]]$data   # exactly what BART/BKMR received
  res <- run_deepmed(dat)
  res$rep_id   <- r
  res$scenario <- sid

  saveRDS(res, part)
  TRUE
}

## ------------------------------ run -----------------------------------
t0    <- Sys.time()
todo  <- Filter(function(r) !file.exists(file.path(PARTS, sprintf("rep_%03d.rds", r))),
                seq_len(R_REPS))
done0 <- R_REPS - length(todo)
if (done0 > 0) cat(sprintf("[resume] %d/%d reps already on disk.\n", done0, R_REPS))

## Batch the reps. One worker == one rep == one process lifetime.
batches <- split(todo, ceiling(seq_along(todo) / N_OUTER))

for (b in seq_along(batches)) {
  reps <- batches[[b]]
  cat(sprintf("[batch %d/%d] reps %s ... ", b, length(batches),
              paste(range(reps), collapse = "-")))
  flush.console()

  cl <- makeCluster(length(reps), type = "PSOCK")
  clusterExport(cl,
    c("R_DIR", "PARTS", "DEEPMED_PY", "sim_data", "sid", "worker"),
    envir = environment())

  ok <- tryCatch(clusterApplyLB(cl, reps, worker),
                 error = function(e) { cat("BATCH ERROR:", conditionMessage(e), "\n"); NULL })

  stopCluster(cl)   ## <-- the actual fix: processes die, RAM returns to the OS
  gc(FALSE)

  n_now <- sum(file.exists(file.path(PARTS, sprintf("rep_%03d.rds", seq_len(R_REPS)))))
  cat(sprintf("done (%d/%d total, %.0f min elapsed)\n", n_now, R_REPS,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

## ---------------------------- assemble --------------------------------
parts   <- file.path(PARTS, sprintf("rep_%03d.rds", seq_len(R_REPS)))
results <- lapply(parts, function(p) if (file.exists(p)) readRDS(p) else NULL)
names(results) <- as.character(seq_len(R_REPS))
saveRDS(results, out_file)

n_ok    <- sum(vapply(results, function(x) isTRUE(x$success), logical(1)))
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("Done scenario %d: %d/%d reps ok | %.1f min | -> %s\n",
            sid, n_ok, R_REPS, elapsed, out_file))

errs <- Filter(function(x) !is.null(x) && !isTRUE(x$success), results)
if (length(errs))
  cat(sprintf("[errors] %d reps failed. First: %s\n", length(errs), errs[[1]]$error))
