#!/usr/bin/env Rscript
## =====================================================================
## run_deepmed.R  —  production DeepMed driver (one scenario per call)
## =====================================================================
## REWRITTEN. The previous version reimplemented DeepMed inline and
## REGENERATED its own data with its own seed, which silently broke the
## paired comparison with BART/BKMR. It also referenced SIM_GRID (an
## object that does not exist anywhere in the project).
##
## This version is a THIN DRIVER:
##   - reads the SAME datasets BART/BKMR consumed, from
##     results/data_for_deepmed/scenario_<ID>.rds  (no regeneration,
##     no reseeding => byte-identical data => valid paired comparison)
##   - delegates every DeepMed detail to run_deepmed() in lib_deepmed.R
##     (mixture->D median split, X = cbind(C,Z), METHOD_CONFIG grid,
##      .clear_tf() memory fix, corrected result parsing)
##
## Usage:  Rscript R/run_deepmed.R <SCENARIO_ID>
## Output: results/deepmed/deepmed_<ID>.rds
##         named list, one element per rep, each the list returned by
##         run_deepmed(): NIE_est/NDE_est/TE_est/... /success
##
## Parallel: N_OUTER PSOCK workers over reps (each an isolated R process
## => isolated TensorFlow session). Per-rep part files allow resume.
## =====================================================================

suppressMessages(library(parallel))

PROJ    <- Sys.getenv("SIM_PROJ", "/work/users/o/d/odeaug/sim_proj")
R_DIR   <- file.path(PROJ, "R")
setwd(PROJ)

source(file.path(R_DIR, "lib_setup.R"))   # design_grid, SIM_CONFIG, METHOD_CONFIG, dirs

## ---- the venv that has tensorflow[and-cuda]; NOT the old conda env ----
DEEPMED_PY <- Sys.getenv("RETICULATE_PYTHON",
                         "/work/users/o/d/odeaug/envs/deepmed_py/bin/python")

N_OUTER <- as.integer(Sys.getenv("DM_OUTER", "10"))
N_INNER <- as.integer(Sys.getenv("DM_INNER", "8"))

## --------------------------- scenario id ------------------------------
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

## ---- the SAME data BART and BKMR used -------------------------------
sim_data <- readRDS(in_file)
R_REPS   <- min(length(sim_data), SIM_CONFIG$R)

cat(sprintf("=== DeepMed | scenario %d (n=%d, L=%d, %s) | %d reps | %d outer x %d inner ===\n",
            sid, sc$n, sc$L, sc$outcome_type, R_REPS, N_OUTER, N_INNER))
cat(sprintf("[data] %s  (%d reps on disk)\n", in_file, length(sim_data)))
cat(sprintf("[py]   %s\n", DEEPMED_PY))

## ------------------------- per-rep worker -----------------------------
worker <- function(r) {
  part <- file.path(PARTS, sprintf("rep_%03d.rds", r))
  if (file.exists(part)) return(part)              # resume

  Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3",
             TF_USE_LEGACY_KERAS  = "1",
             OMP_NUM_THREADS      = "1",
             RETICULATE_PYTHON    = DEEPMED_PY)

  suppressMessages({
    library(reticulate)
    reticulate::use_python(DEEPMED_PY, required = TRUE)   # venv, not conda
    library(tensorflow); library(keras); library(DeepMed)
    library(foreach);    library(doParallel)
  })

  ## HARD GPU ASSERT: fail loudly rather than silently crawling on CPU.
  ## Set DM_ALLOW_CPU=1 to override.
  gpus <- tryCatch(length(tensorflow::tf$config$list_physical_devices("GPU")),
                   error = function(e) 0L)
  if (gpus < 1 && !nzchar(Sys.getenv("DM_ALLOW_CPU"))) {
    stop(sprintf("rep %d: TensorFlow sees NO GPU (python=%s). Refusing to run on CPU.",
                 r, reticulate::py_config()$python))
  }

  try(tensorflow::tf$config$threading$set_intra_op_parallelism_threads(1L), silent = TRUE)
  try(tensorflow::tf$config$threading$set_inter_op_parallelism_threads(1L), silent = TRUE)
  registerDoParallel(cores = N_INNER)              # DeepMed's hyper_grid pool

  source(file.path(R_DIR, "lib_setup.R"),   local = FALSE)  # METHOD_CONFIG
  source(file.path(R_DIR, "lib_deepmed.R"), local = FALSE)  # run_deepmed()

  ## the rep's data, exactly as BART/BKMR received it
  dat <- sim_data[[r]]$data

  res <- run_deepmed(dat)                         # ALL DeepMed logic lives here
  res$rep_id   <- r
  res$scenario <- sid

  saveRDS(res, part)
  part
}

## ------------------------------ run -----------------------------------
t0 <- Sys.time()
already <- sum(file.exists(file.path(PARTS, sprintf("rep_%03d.rds", seq_len(R_REPS)))))
if (already > 0) cat(sprintf("[resume] %d/%d reps already on disk.\n", already, R_REPS))

cl <- makeCluster(N_OUTER, type = "PSOCK")
on.exit(stopCluster(cl), add = TRUE)
clusterExport(cl,
  c("R_DIR", "PARTS", "DEEPMED_PY", "N_INNER", "sim_data", "worker"),
  envir = environment())

invisible(clusterApplyLB(cl, seq_len(R_REPS), worker))

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
if (length(errs)) {
  cat(sprintf("[errors] %d reps failed. First: %s\n",
              length(errs), errs[[1]]$error))
}
