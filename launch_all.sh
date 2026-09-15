#!/bin/bash
## =====================================================================
## launch_all.sh — pre-flight + production launch for the 3-method study
## Run from the LOGIN node:
##   bash launch_all.sh                          # BKMR + BART; holds DeepMed
##   bash launch_all.sh --with-deepmed           # also submit DeepMed
##   bash launch_all.sh --with-deepmed --chain   # also queue aggregate (dependency)
## =====================================================================
set -uo pipefail

PROJ="/work/users/o/d/odeaug/sim_proj"
cd "$PROJ" || { echo "Cannot cd to $PROJ"; exit 1; }

WITH_DEEPMED=0; CHAIN=0
for a in "$@"; do
  case "$a" in
    --with-deepmed) WITH_DEEPMED=1 ;;
    --chain)        CHAIN=1 ;;
    *) echo "Unknown option: $a"; exit 1 ;;
  esac
done

echo "==> [1/4] Stripping Windows carriage returns"
sed -i 's/\r$//' R/run_method.R R/run_deepmed.R slurm/*.sh 2>/dev/null || true
chmod +x slurm/*.sh 2>/dev/null || true

echo "==> [2/4] Pre-flight wiring check"
module load r/4.5.0 2>/dev/null || true
Rscript -e '
ok <- TRUE
src <- function(f) if (file.exists(f)) source(f) else { message("FAIL: missing ", f); ok <<- FALSE }
src("R/lib_setup.R"); src("R/lib_dgp.R"); src("R/lib_methods.R")

if (!exists("SIM_GRID")) { message("FAIL: SIM_GRID not found in lib_setup.R"); ok <- FALSE
} else if (nrow(SIM_GRID) != 12) message(sprintf("WARN: SIM_GRID has %d rows (expected 12)", nrow(SIM_GRID)))

if (!exists("generate_data")) { message("FAIL: generate_data() not found"); ok <- FALSE }
if (!all(sapply(c("run_bkmr_cma","run_bart_cma"), exists)))
  { message("FAIL: a method wrapper is missing"); ok <- FALSE }
if (!exists("METHOD_CONFIG"))
  message("WARN: METHOD_CONFIG not found — if your hyperparameters live in SIM_CONFIG, edit run_method.R")

if (exists("SIM_GRID") && exists("generate_data")) {
  d <- tryCatch(generate_data(scenario = SIM_GRID[1, ], n = SIM_GRID$n[1]),
                error = function(e) { message("FAIL: generate_data() errored: ", conditionMessage(e)); NULL })
  if (is.null(d)) ok <- FALSE else {
    miss <- setdiff(c("Y","D","M","Z","X"), names(d))
    if (length(miss)) { message("FAIL: generate_data() is missing field(s): ", paste(miss, collapse=", ")); ok <- FALSE
    } else message("OK: generate_data() returns Y, D, M, Z, X")
  }
}
if (!ok) { message("\nPRE-FLIGHT FAILED — fix the items above, then re-run."); quit(status = 1) }
message("Pre-flight passed.")
' || { echo "Aborting: pre-flight failed."; exit 1; }

echo "==> [3/4] Clearing scenario-01 smoke-test remnants"
echo "    (a 3-rep smoke .rds would otherwise make production SKIP scenario 1)"
rm -f  results/bkmr/bkmr_scenario_01.rds \
       results/bart/bart_scenario_01.rds \
       results/deepmed/dm_scenario_01.rds 2>/dev/null || true
rm -rf results/bkmr/parts_s01 results/bart/parts_s01 results/deepmed/parts_s01 2>/dev/null || true

echo "==> [4/4] Submitting jobs"
BKMR_ID=$(sbatch --parsable slurm/01_bkmr.sh); echo "    BKMR  array job: $BKMR_ID"
BART_ID=$(sbatch --parsable slurm/02_bart.sh); echo "    BART  array job: $BART_ID"
DEPS="afterok:${BKMR_ID}:${BART_ID}"

if [ "$WITH_DEEPMED" -eq 1 ]; then
  DM_ID=$(sbatch --parsable slurm/03_deepmed.sh); echo "    DeepMed array job: $DM_ID"
  DEPS="${DEPS}:${DM_ID}"
else
  echo "    DeepMed: HELD until you've confirmed res\$results (runbook Stage 4a)."
  echo "             Release with:  sbatch slurm/03_deepmed.sh"
fi

echo
echo "Monitor:  squeue -u odeaug"
echo "Logs:     tail -f results/bart_logs/bart_1_*.out"
if [ "$CHAIN" -eq 1 ]; then
  AGG_ID=$(sbatch --parsable --dependency="${DEPS}" slurm/04_aggregate.sh)
  echo "Aggregate: queued as $AGG_ID (auto-runs after the method jobs succeed)"
else
  echo "Aggregate (when methods finish):"
  echo "    sbatch --dependency=${DEPS} slurm/04_aggregate.sh"
fi
echo
echo "Done submitting."
