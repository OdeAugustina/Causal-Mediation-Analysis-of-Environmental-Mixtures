#!/usr/bin/env Rscript
## =====================================================================
## run_method.R  —  driver for the R-only methods: BKMR-CMA & BART-CMA
## =====================================================================
## DeepMed has its own driver (run_deepmed.R) because TensorFlow is not
## fork-safe. These two methods are pure R/C++, so they parallelize
## cleanly and far faster with mclapply (FORK) — no PSOCK, no conda.
##
## Usage:
##   METHOD=bkmr Rscript R/run_method.R <SCENARIO_ID>
##   METHOD=bart Rscript R/run_method.R <SCENARIO_ID>
##   (or under SLURM array: reads SLURM_ARRAY_TASK_ID)
##
## Submit 01_bkmr.sh AND 02_bart.sh together -> the two methods run as
## independent, concurrent jobs (this is the "in parallel" you want).
##
## Output: results/<method>/<method>_scenario_<ID>.rds
##   tidy df: method | rep | estimand | effect | se | ci_lower | ci_upper
##            | runtime | error
##   estimands (fixed order): TE, NDE_treated, NDE_control,
##                            NIE_treated, NIE_control
##
## PAIRED DESIGN (important):
##   The per-rep seed below is IDENTICAL to run_deepmed.R, so all three
##   methods are scored on the SAME simulated datasets (common random
##   numbers). That makes every method-to-method comparison paired and
##   sharply reduces the Monte-Carlo variance of the differences — the
##   whole point of an aligned benchmark.
## =====================================================================

suppressMessages(library(parallel))

METHOD <- tolower(Sys.getenv("METHOD", ""))
if (!METHOD %in% c("bkmr", "bart"))
  stop("Set METHOD=bkmr or METHOD=bart (DeepMed uses run_deepmed.R).")

## ----------------------------- config --------------------------------
PROJ    <- Sys.getenv("SIM_PROJ", "/work/users/o/d/odeaug/sim_proj")
R_DIR   <- file.path(PROJ, "R")
RES_DIR <- file.path(PROJ, "results", METHOD)
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

N_CORES   <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",
                          Sys.getenv("N_CORES", "8")))
R_REPS    <- as.integer(Sys.getenv("DM_REPS", "100"))
BASE_SEED <- as.integer(Sys.getenv("DM_SEED", "20260101"))

## --------------------------- scenario id ------------------------------
args <- commandArgs(trailingOnly = TRUE)
SID <- if (length(args) >= 1) suppressWarnings(as.integer(args[1])) else
         suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "")))
if (is.na(SID)) stop("No scenario id: pass as an argument or set SLURM_ARRAY_TASK_ID.")

out_file <- file.path(RES_DIR, sprintf("%s_scenario_%02d.rds", METHOD, SID))
if (file.exists(out_file)) {
  message(sprintf("[skip] %s scenario %02d already complete -> %s",
                  METHOD, SID, out_file))
  quit(save = "no", status = 0)
}
PARTS <- file.path(RES_DIR, sprintf("parts_s%02d", SID))
dir.create(PARTS, showWarnings = FALSE)

## --------- load project code once; forked workers inherit it ----------
## NOTE: confirm these names against your lib_setup.R —
##   SIM_GRID       : the 12-row scenario design data.frame
##   METHOD_CONFIG  : list of method hyperparameters (bkmr_iter,
##                    bkmr_burnin, bkmr_K_med, bart_ntree, bart_ndpost,
##                    bart_nskip). If they live in SIM_CONFIG instead,
##                    change METHOD_CONFIG -> SIM_CONFIG below.
suppressMessages({
  source(file.path(R_DIR, "lib_setup.R"))     # SIM_GRID / METHOD_CONFIG
  source(file.path(R_DIR, "lib_dgp.R"))       # generate_data()
  source(file.path(R_DIR, "lib_methods.R"))   # run_bkmr_cma(), run_bart_cma()
})
USE_MOCK <- FALSE   # lib_methods.R defaults this to TRUE — force real methods

method_fun   <- if (METHOD == "bkmr") run_bkmr_cma else run_bart_cma
scen         <- SIM_GRID[SID, ]
outcome_type <- if (!is.null(scen$outcome_type)) as.character(scen$outcome_type) else "continuous"
ESTIMANDS    <- c("TE", "NDE_treated", "NDE_control", "NIE_treated", "NIE_control")

cat(sprintf("=== %s | scenario %02d | %d reps | %d cores | outcome=%s ===\n",
            toupper(METHOD), SID, R_REPS, N_CORES, outcome_type))

## ------------------------- per-rep worker -----------------------------
run_one_rep <- function(r) {
  part <- file.path(PARTS, sprintf("rep_%03d.rds", r))
  if (file.exists(part)) return(part)                 # resume: already done

  set.seed(BASE_SEED + SID * 100000L + r)             # SAME seed as run_deepmed.R
  dat <- generate_data(scenario = scen, n = scen$n)

  res <- tryCatch(
    method_fun(Y = dat$Y, D = dat$D, M = dat$M, Z = dat$Z, X = dat$X,
               outcome_type = outcome_type, config = METHOD_CONFIG),
    error = function(e) structure(list(err = conditionMessage(e)), class = "m_err"))

  if (inherits(res, "m_err")) {
    df <- data.frame(method = toupper(METHOD), rep = r, estimand = ESTIMANDS,
                     effect = NA_real_, se = NA_real_,
                     ci_lower = NA_real_, ci_upper = NA_real_,
                     runtime = NA_real_, error = res$err,
                     stringsAsFactors = FALSE)
  } else {
    df <- data.frame(method = toupper(METHOD), rep = r, estimand = ESTIMANDS,
                     effect   = as.numeric(res$estimates),
                     se       = NA_real_,
                     ci_lower = as.numeric(res$ci_lower),
                     ci_upper = as.numeric(res$ci_upper),
                     runtime  = as.numeric(res$runtime)[1],
                     error    = NA_character_,
                     stringsAsFactors = FALSE)
  }
  saveRDS(df, part)
  part
}

## --------------------------- run in parallel --------------------------
## FORK is safe here (no TensorFlow). mc.preschedule=FALSE load-balances
## the heterogeneous rep times (BKMR reps vary a lot).
t0 <- Sys.time()
already <- sum(file.exists(file.path(PARTS, sprintf("rep_%03d.rds", seq_len(R_REPS)))))
if (already > 0) cat(sprintf("Resuming: %d/%d reps already on disk.\n", already, R_REPS))

RNGkind("L'Ecuyer-CMRG")
invisible(mclapply(seq_len(R_REPS), run_one_rep,
                   mc.cores = N_CORES, mc.preschedule = FALSE))

## ---------------------------- assemble --------------------------------
parts <- file.path(PARTS, sprintf("rep_%03d.rds", seq_len(R_REPS)))
results <- do.call(rbind, lapply(parts, function(p) if (file.exists(p)) readRDS(p)))
saveRDS(results, out_file)

elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
n_ok <- length(unique(results$rep[is.na(results$error)]))
cat(sprintf("Done %s scenario %02d: %d/%d reps ok | %.1f min | -> %s\n",
            toupper(METHOD), SID, n_ok, R_REPS, elapsed, out_file))
