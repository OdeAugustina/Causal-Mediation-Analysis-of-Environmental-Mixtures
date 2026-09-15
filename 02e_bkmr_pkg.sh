#!/bin/bash
## =====================================================================
## 02e_bkmr_pkg.sh — BKMR-CMA via causalbkmr::mediation.bkmr()
##
##   sbatch SLURM/02e_bkmr_pkg.sh              # all 12
##   sbatch --array=1,3,5,7 SLURM/02e_bkmr_pkg.sh   # continuous, n=300
##
## Runs through run_scenario.R so the per-rep seeding
## (2026 + scenario_id*1000 + rep_id) is IDENTICAL to the existing
## BKMR / BART / DeepMed results. Datasets stay paired.
##
## Writes to results/raw_scenario_XXX_bkmr_pkg.rds -- the manual-
## implementation results are NOT overwritten, so you can compare.
##
## NOTE ON WALLTIME: 02a_bkmr.sh used 18h and scenarios 2/4/6/8 (n=1000)
## needed the separate 11-day 02c job. Same applies here -- launch
## 1,3,5,7,9,10,11,12 with this script and mirror 02c for the big-n four.
## =====================================================================
#SBATCH --job-name=bkmrpkg
#SBATCH --array=1-12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=64g
#SBATCH --time=11-00:00:00
#SBATCH --partition=general
#SBATCH --output=results/bkmr_pkg_logs/bkmrpkg_%a_%j.out
#SBATCH --error=results/bkmr_pkg_logs/bkmrpkg_%a_%j.err

set -euo pipefail
cd /work/users/o/d/odeaug/sim_proj
mkdir -p results/bkmr_pkg_logs results/bkmr_pkg_raw

module purge
module load r/4.5.0

export SIM_PROJ=/work/users/o/d/odeaug/sim_proj
export METHOD=bkmr_pkg

## ---- Variable selection: OFF for both arms ----
## varsel=FALSE => BKMR-CMA (base variant, Devick et al. 2022).
## varsel=TRUE  => BKMR-CMA-VS (spike-and-slab on kernel weights).
##
## FALSE is used here because this study's DGP is DENSE: every element
## of alpha_Z_base and beta_Z_base is nonzero, so a sparsity prior is
## misspecified against the truth. It is also the closer match to the
## comparators -- neither dbarts BART nor DeepMed imposes a point mass
## at zero. Devick's finding that VS wins at L=10 came from a SPARSE
## DGP and does not transfer here.
##
## To run the VS arm later: set CONT=TRUE and change the out tag so it
## does not overwrite these results.
export BKMR_VARSEL_CONT=TRUE
export BKMR_VARSEL_BINARY=FALSE

## ---- Runtime levers ----
## mediation.bkmr() cost is O(length(sel) * K).
## Unthinned (by=1) with K=100 is ~150,000 predictions per rep.
## THIN=5 with K=30 is ~9,000 -- roughly 17x cheaper.
## K=30 is supported by the scenario-12 n_med experiment: bias and
## interval width were flat across K in {30,100,300,1000}.
export BKMR_SEL_THIN=5

## Keep rep 1's raw mediation.bkmr() object per scenario for auditing.
export BKMR_SAVE_RAW_REP1=TRUE

echo "BKMR-pkg scenario ${SLURM_ARRAY_TASK_ID} on $(hostname) | ${SLURM_CPUS_PER_TASK} cpus | $(date)"
Rscript R/run_scenario.R "${SLURM_ARRAY_TASK_ID}"
echo "Finished scenario ${SLURM_ARRAY_TASK_ID} | $(date)"
