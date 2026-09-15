###########################################################
## lib_setup.R — Configuration sourced by all scripts
## Faithful to original Chunk 1 settings (FULL PRODUCTION).
###########################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(MASS); library(tidyr)
})
filter <- dplyr::filter; select <- dplyr::select

## ---- Paths (override RESULTS_DIR via env var if desired) ----
RESULTS_DIR <- Sys.getenv("SIM_RESULTS_DIR",
                          unset = file.path(getwd(), "results"))
TRUTH_DIR    <- file.path(RESULTS_DIR, "truth")
DATA_DIR     <- file.path(RESULTS_DIR, "data_for_deepmed")
DEEPMED_DIR  <- file.path(RESULTS_DIR, "deepmed")
ERRORS_DIR   <- file.path(RESULTS_DIR, "errors")
FIGURES_DIR  <- file.path(RESULTS_DIR, "figures")

for (d in c(RESULTS_DIR, TRUTH_DIR, DATA_DIR, DEEPMED_DIR, ERRORS_DIR, FIGURES_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

## ---- Cores (SLURM-aware) ----
SIM_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",
                                   unset = max(1, parallel::detectCores() - 1)))

## ---- Simulation config (matches original) ----
SIM_CONFIG <- list(seed = 2026, R = 100, N_truth = 1e6,
                   rho_Z = 0.3, rho_C = 0.2, P_conf = 3)

DGP_COEFS <- list(
  theta_CZ = 0.1, alpha_C = c(0.3,-0.2,0.1), beta_C = c(0.2,-0.1,0.15),
  alpha_0 = 0.5,
  alpha_Z_base = c(0.6,0.4,0.3,0.2,0.1,0.5,0.35,0.25,0.15,0.05),
  sigma_M = 1.0,
  beta_0 = 1.0,
  beta_Z_base = c(0.4,0.3,0.2,0.1,0.05,0.35,0.25,0.15,0.08,0.03),
  beta_M = 0.5,
  beta_ZM_base = c(0.2,0.15,0.1,0.05,0.0,0.18,0.12,0.08,0.03,0.0),
  sigma_Y = 1.0
)

METHOD_CONFIG <- list(
  bkmr_iter = 3000, bkmr_burnin = 1500, bkmr_K_med = 100,
  bart_ntree = 200, bart_ndpost = 1000, bart_nskip = 250,
  dm_method = "DNN",
  dm_epochs = 50, dm_batch_size = 100, dm_trim = 0.05,
  dm_l1 = c(0,0.05), dm_layers = c(2), dm_units = c(20)
)

## ---- 12-scenario factorial design grid ----
design_grid <- tribble(
  ~scenario_id, ~linearity,  ~interaction, ~n,   ~L, ~outcome_type, ~group,
  1,  "linear",    "none",    300,  5,  "continuous", "core",
  2,  "linear",    "none",    1000, 5,  "continuous", "core",
  3,  "linear",    "present", 300,  5,  "continuous", "core",
  4,  "linear",    "present", 1000, 5,  "continuous", "core",
  5,  "nonlinear", "none",    300,  5,  "continuous", "core",
  6,  "nonlinear", "none",    1000, 5,  "continuous", "core",
  7,  "nonlinear", "present", 300,  5,  "continuous", "core",
  8,  "nonlinear", "present", 1000, 5,  "continuous", "core",
  9,  "nonlinear", "present", 300,  5,  "binary",     "ext_outcome",
  10, "nonlinear", "present", 1000, 5,  "binary",     "ext_outcome",
  11, "nonlinear", "present", 300,  3,  "continuous", "ext_mixture",
  12, "nonlinear", "present", 300,  10, "continuous", "ext_mixture"
)

## ---- Helper: read SLURM_ARRAY_TASK_ID or commandline arg ----
get_scenario_id <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  sid_env <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA)
  if (length(args) >= 1 && nzchar(args[1])) {
    as.integer(args[1])
  } else if (!is.na(sid_env) && nzchar(sid_env)) {
    as.integer(sid_env)
  } else {
    stop("No scenario_id provided. Pass as commandline arg or set SLURM_ARRAY_TASK_ID.")
  }
}

cat(sprintf("[setup] RESULTS_DIR = %s\n", RESULTS_DIR))
cat(sprintf("[setup] SIM_CORES   = %d\n", SIM_CORES))
