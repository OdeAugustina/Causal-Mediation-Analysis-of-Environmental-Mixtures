#!/bin/bash
#SBATCH --job-name=cma_bkmr_bign
#SBATCH --array=2,4,6,8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=64g
#SBATCH --time=11-00:00:00
#SBATCH --partition=general
#SBATCH --output=logs/bkmr_bign_%A_%a.out
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=aoodediran@aggies.ncat.edu
set -euo pipefail
module purge
module load r/4.5.0
cd "${SLURM_SUBMIT_DIR}"
export METHOD=bkmr
echo "===== BKMR-CMA n=1000 rerun ====="
echo "Array Task:  ${SLURM_ARRAY_TASK_ID}"
echo "Time limit:  11 days"
echo "Start:       $(date)"
Rscript R/run_scenario.R "${SLURM_ARRAY_TASK_ID}"
echo "End: $(date)"
