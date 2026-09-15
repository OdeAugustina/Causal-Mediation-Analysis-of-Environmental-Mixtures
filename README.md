# Dissertation Simulation — Longleaf Setup

Causal mediation analysis simulation study: BKMR-CMA vs BART-CMA vs DeepMed across 12 scenarios × 100 reps.

This directory is a refactor of `simulation_FULL_PRODUCTION.Rmd` into a SLURM-array-parallel pipeline for Longleaf.

---

## Directory layout

```
sim_proj/
├── R/
│   ├── lib_setup.R       # Config + design grid (sourced by all)
│   ├── lib_dgp.R         # Data generating process functions
│   ├── lib_truth.R       # Monte Carlo truth computation
│   ├── lib_methods.R     # BKMR + BART runners, run_one_rep, compute_metrics
│   ├── lib_deepmed.R     # DeepMed runner (isolated from R-native methods)
│   ├── run_truth.R       # Script: compute truth for ONE scenario
│   ├── run_scenario.R    # Script: BART+BKMR for ONE scenario
│   ├── run_deepmed.R     # Script: DeepMed for ONE scenario
│   └── aggregate.R       # Script: merge + metrics + charts + Excel
├── slurm/
│   ├── 01_truth.sh       # Array 1-12, general partition
│   ├── 02_bart_bkmr.sh   # Array 1-12, general partition, 60 cores
│   ├── 03_deepmed.sh     # Array 1-12, GPU partition
│   └── 04_aggregate.sh   # Single job, general partition
├── logs/                 # SLURM stdout/stderr (created automatically)
├── launch.sh             # One-line full-pipeline submission
└── README.md
```

---

## First-time setup (do once)

### 1. Place the project on `/work`

Your home directory has a small quota (~50 GB) and isn't optimized for high-I/O. Run on `/work`:

```bash
mkdir -p /work/users/o/d/odeaug/sim_proj
cd /work/users/o/d/odeaug/sim_proj
# Upload all files in this directory to here (via OnDemand file browser
# or scp/sftp from your laptop).
```

### 2. Install R packages

Start an **interactive RStudio session via OnDemand** (don't do this on the login node):

- Go to https://ondemand.rc.unc.edu
- Interactive Apps → RStudio Server
- Allocation: 4 CPUs, 8 GB RAM, 2 hours, general partition

In the RStudio console:

```r
install.packages(c("dplyr", "tibble", "MASS", "tidyr", "ggplot2",
                   "dbarts", "bkmr", "reticulate", "keras",
                   "writexl"))

# causalbkmr is GitHub-only:
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("zorabian/causalbkmr")

# DeepMed is also GitHub-only:
remotes::install_github("siqixu/DeepMed")
```

Some packages take a while to compile — be patient.

### 3. Python environment setup (for DeepMed)

DeepMed uses TensorFlow via reticulate. We build it in an isolated venv so it never conflicts with the BART/BKMR R sessions.

On the Longleaf shell (SSH or OnDemand → Clusters → Shell Access):

```bash
module load python/3.11
mkdir -p /work/users/o/d/odeaug/envs
python -m venv /work/users/o/d/odeaug/envs/deepmed_py
source /work/users/o/d/odeaug/envs/deepmed_py/bin/activate
pip install --upgrade pip
pip install "tensorflow==2.15" keras numpy pandas
deactivate
```

This venv path is referenced in `slurm/03_deepmed.sh`. If you put it elsewhere, edit `DEEPMED_VENV` in that file.

---

## Launch the full pipeline

```bash
cd /work/users/o/d/odeaug/sim_proj
chmod +x launch.sh
./launch.sh
```

This submits four jobs with a dependency chain:

```
01_truth (1hr) → 02_bart_bkmr (12hr) → 03_deepmed (4hr) → 04_aggregate (1hr)
```

Each step waits for the previous one to succeed (`--dependency=afterok`).
You'll receive emails when each completes or fails.

---

## Monitoring

```bash
squeue -u $USER                          # what's running / queued
sacct --format=JobID,JobName%20,State,Elapsed,MaxRSS,ExitCode  # job history
seff <jobid>                             # efficiency report for a finished job
tail -f logs/bart_bkmr_<jobid>_<task>.out  # live log
```

If a single scenario fails, you can resubmit just that one:

```bash
sbatch --array=7 slurm/02_bart_bkmr.sh   # rerun scenario 7
```

---

## Outputs

All results go to `results/` (which lives on `/work`):

- `results/truth/truth_<sid>.rds` — ground truth for each scenario
- `results/raw_scenario_<sid>.rds` — BART+BKMR rep results (light, no $data)
- `results/data_for_deepmed/scenario_<sid>.rds` — full data needed for DeepMed
- `results/deepmed/deepmed_<sid>.rds` — DeepMed results per scenario
- `results/errors/scenario_<sid>/rep_<r>_<method>.txt` — per-rep error messages
- `results/all_metrics.rds` and `all_metrics.xlsx` — final aggregated metrics
- `results/figures/*.pdf` — performance plots

---

## Testing before full launch

For a quick smoke test (~10 minutes), edit `R/lib_setup.R` to set `R = 5` in
`SIM_CONFIG`, then submit only scenarios 1 and 11 (the smallest):

```bash
sbatch --array=1,11 slurm/01_truth.sh
# wait for it, then
sbatch --array=1,11 slurm/02_bart_bkmr.sh
```

Once that works, restore `R = 100` and run the full `./launch.sh`.

---

## Why this is faster than the Lambda Rmd

The Rmd ran scenarios **sequentially**: scenario 12 had to wait for scenario 1 to finish.

This refactor runs all 12 scenarios as **independent SLURM array tasks**, so they execute in parallel across multiple cluster nodes. Total wall time becomes "longest single scenario" instead of "sum of all scenarios."

Inside each scenario, the 100 reps still parallelize across cores (now via `SLURM_CPUS_PER_TASK = 60` instead of hardcoded 60), so the per-scenario speed is unchanged. Improvements stack: the overall pipeline finishes roughly **10-12× faster** than serial scenario execution.

The TF/reticulate environment conflict is eliminated structurally: the BART+BKMR jobs never load Python at all, and DeepMed runs in its own isolated venv. There's no shared state for the conflict to occur in.
