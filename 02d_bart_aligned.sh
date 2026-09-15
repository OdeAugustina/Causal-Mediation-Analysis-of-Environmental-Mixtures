#!/bin/bash
## =====================================================================
## 02d_bart_aligned.sh — rerun BART-CMA ONLY, with the aligned estimator
##
##   sbatch slurm/02d_bart_aligned.sh
##
## Runs through run_scenario.R (NOT run_method.R) so the per-rep data
## seed stays 2026 + scenario_id*1000 + rep_id.  That is the seed the
## existing raw_scenario_*_bkmr.rds and results/data_for_deepmed/ files
## were built under, so the reruns land on identical datasets and the
## paired design is preserved.
##
## Writes results/raw_scenario_<sid>_bart.rds.  Recombine afterwards:
##   Rscript R/merge_split.R
##   sbatch  slurm/04_aggregate.sh
##
## PREREQUISITE — delete the stale BART files first, or every task will
## skip.  run_scenario.R short-circuits when its outfile already exists:
##   cd /work/users/o/d/odeaug/sim_proj
##   mkdir -p results/_archive_bart_nmed1
##   mv results/raw_scenario_*_bart.rds results/_archive_bart_nmed1/
##   cp results/all_metrics.rds results/_archive_bart_nmed1/
## =====================================================================
#SBATCH --job-name=bart_aligned
#SBATCH --array=1-12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=64g
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --output=results/bart_logs/bart_aligned_%a_%j.out
#SBATCH --error=results/bart_logs/bart_aligned_%a_%j.err

set -euo pipefail
cd /work/users/o/d/odeaug/sim_proj
mkdir -p results/bart_logs

module purge
module load r/4.5.0

export METHOD=bart          # run_scenario.R: BART only, writes *_bart.rds
export SIM_PROJ=/work/users/o/d/odeaug/sim_proj

echo "BART-CMA (aligned) scenario ${SLURM_ARRAY_TASK_ID} on $(hostname)"
echo "  cpus : ${SLURM_CPUS_PER_TASK}"
echo "  start: $(date)"

Rscript R/run_scenario.R "${SLURM_ARRAY_TASK_ID}"

echo "  end  : $(date)"
