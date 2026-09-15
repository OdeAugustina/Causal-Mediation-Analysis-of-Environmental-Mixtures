#!/bin/bash
#SBATCH --job-name=cma_bkmr
#SBATCH --array=1-12
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=48g
#SBATCH --time=18:00:00
#SBATCH --partition=general
#SBATCH --output=logs/bkmr_%A_%a.out
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=aoodediran@aggies.ncat.edu
set -euo pipefail
module purge
module load r/4.5.0
cd "${SLURM_SUBMIT_DIR}"
export METHOD=bkmr
echo "===== BKMR-CMA job ====="
echo "Host:          $(hostname)"
echo "Array Task:    ${SLURM_ARRAY_TASK_ID}"
echo "CPUs:          ${SLURM_CPUS_PER_TASK}"
echo "Method:        ${METHOD}"
echo "Start:         $(date)"
echo
Rscript R/run_scenario.R "${SLURM_ARRAY_TASK_ID}"
echo
echo "End: $(date)"
