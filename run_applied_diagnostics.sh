#!/bin/bash
#SBATCH --job-name=applied_nmed
#SBATCH --output=logs/applied_nmed_%j.out
#SBATCH --error=logs/applied_nmed_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

set -euo pipefail
module purge
module load r/4.5.0
cd "$SLURM_SUBMIT_DIR"

export APPLIED_XLSX="/work/users/o/d/odeaug/BART_CMA/completed_liver_pfas_ALL.xlsx"
export APPLIED_NIMP=5
export APPLIED_CONF="primary"

echo "host    : $(hostname)"
echo "started : $(date)"
Rscript applied_nmed_and_support_diagnostics.R
echo "finished: $(date)"
