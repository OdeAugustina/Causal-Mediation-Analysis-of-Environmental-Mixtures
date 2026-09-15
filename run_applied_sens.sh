#!/bin/bash
#SBATCH --job-name=bartapp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=240g
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --output=applied_sens_%j.out
#SBATCH --error=applied_sens_%j.err

cd /work/users/o/d/odeaug/BART_CMA
module load r/4.5.0

Rscript -e 'rmarkdown::render("BART_CMA_Mediation.Rmd",
  params = list(profile = "production", n_cores = 60,
                confounder_set = "sensitivity"),
  output_file = "applied_sensitivity_paired.html")'
