#!/bin/bash
#SBATCH --job-name=cma_bart_bkmr
#SBATCH --array=1-12
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=48g
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --output=logs/bart_bkmr_%A_%a.out
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=aoodediran@aggies.ncat.edu

set -euo pipefail

module purge
module load r/4.4.0

cd "${SLURM_SUBMIT_DIR}"

echo "===== BART+BKMR job ====="
echo "Host:           $(hostname)"
echo "Job ID:         ${SLURM_JOB_ID}"
echo "Array Task:     ${SLURM_ARRAY_TASK_ID}"
echo "CPUs per task:  ${SLURM_CPUS_PER_TASK}"
echo "Working Dir:    $(pwd)"
echo "Start:          $(date)"
echo

Rscript R/run_scenario.R "${SLURM_ARRAY_TASK_ID}"

echo
echo "End: $(date)"
