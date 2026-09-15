#!/bin/bash
#SBATCH --job-name=bart_cma_render
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --time=06:00:00
#SBATCH --output=logs/render_%j.out
#SBATCH --error=logs/render_%j.err
#SBATCH --partition=general

mkdir -p logs

module purge
module load r/4.5.0

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd "${SLURM_SUBMIT_DIR}"

Rscript -e 'rmarkdown::render("BART_CMA_Mediation.Rmd",
  params = list(profile = "production",
                confounder_set = "primary",
                n_cores = 20),
  output_file = "BART_CMA_Mediation_production_primary.html")'
