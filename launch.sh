#!/bin/bash
###########################################################
## launch.sh — Submit the full pipeline with dependency chain
##
## Usage:   ./launch.sh
##
## What it does:
##   1. Truth precompute   (array 1-12, CPU, ~1hr)
##   2. BART + BKMR pass   (array 1-12, CPU, ~12hr)   depends on (1)
##   3. DeepMed pass       (array 1-12, GPU, ~4hr)    depends on (2)
##   4. Aggregate metrics  (single job, CPU, ~1hr)    depends on (3)
###########################################################

set -euo pipefail
cd "$(dirname "$0")"

mkdir -p logs

echo "Submitting truth job..."
J1=$(sbatch --parsable slurm/01_truth.sh)
echo "  -> ${J1}"

echo "Submitting BART+BKMR job (depends on ${J1})..."
J2=$(sbatch --parsable --dependency=afterok:${J1} slurm/02_bart_bkmr.sh)
echo "  -> ${J2}"

echo "Submitting DeepMed job (depends on ${J2})..."
J3=$(sbatch --parsable --dependency=afterok:${J2} slurm/03_deepmed.sh)
echo "  -> ${J3}"

echo "Submitting aggregate job (depends on ${J3})..."
J4=$(sbatch --parsable --dependency=afterok:${J3} slurm/04_aggregate.sh)
echo "  -> ${J4}"

echo
echo "Pipeline launched. Monitor with:"
echo "  squeue -u \$USER"
echo "  sacct -j ${J1},${J2},${J3},${J4} --format=JobID,JobName%20,State,Elapsed,MaxRSS"
