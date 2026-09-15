#!/bin/bash
#SBATCH --job-name=cma_aggregate
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16g
#SBATCH --time=01:00:00
#SBATCH --partition=general
#SBATCH --output=logs/aggregate_%j.out
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=aoodediran@aggies.ncat.edu

set -euo pipefail

module purge
module load r/4.4.0

cd "${SLURM_SUBMIT_DIR}"

echo "===== Aggregate job ====="
echo "Host:        $(hostname)"
echo "Job ID:      ${SLURM_JOB_ID}"
echo "Working Dir: $(pwd)"
echo "Start:       $(date)"
echo

Rscript R/aggregate.R

echo
echo "End: $(date)"
