#!/bin/bash
#SBATCH --job-name=cma_truth
#SBATCH --array=1-12
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --time=01:00:00
#SBATCH --partition=general
#SBATCH --output=logs/truth_%A_%a.out
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=aoodediran@aggies.ncat.edu

set -euo pipefail

module purge
module load r/4.4.0

cd "${SLURM_SUBMIT_DIR}"

echo "===== Truth job ====="
echo "Host:        $(hostname)"
echo "Job ID:      ${SLURM_JOB_ID}"
echo "Array Task:  ${SLURM_ARRAY_TASK_ID}"
echo "Working Dir: $(pwd)"
echo "Start:       $(date)"
echo

Rscript R/run_truth.R "${SLURM_ARRAY_TASK_ID}"

echo
echo "End: $(date)"
