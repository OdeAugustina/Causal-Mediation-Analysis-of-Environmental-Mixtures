# Causal Mediation Analysis of Environmental Mixtures

Code accompanying *Estimand Alignment and Implementation Sensitivity in Causal Mediation
Analysis of Environmental Mixtures: Simulation, Plasmode, and Applied Study*.

The repository contains three components:

1. A **simulation study** benchmarking BKMR-CMA, BART-CMA, and DeepMed across twelve
   synthetic scenarios (2 × 2 × 2 core factorial in linearity × exposure–mediator
   interaction × sample size, plus binary-outcome and mixture-size extension arms).
2. A **plasmode arm** (scenarios 13–14) built on the empirical exposure and covariate
   structure of NHANES, preserving the observed correlation and support geometry.
3. An **applied analysis** of a five-component metal/PFAS mixture, the Dietary
   Inflammatory Index as mediator, and six liver outcomes in NHANES.

All estimands are defined on a common componentwise Q25 → Q75 mixture contrast.

---

## Requirements

Developed and run on the UNC Longleaf HPC cluster (SLURM) under **R 4.5.0**.


### R packages

```r
install.packages(c("dplyr", "tibble", "MASS", "tidyr", "ggplot2",
                   "dbarts", "BART", "bkmr", "reticulate",
                   "writexl", "readxl", "rmarkdown"))
```

Both `dbarts` and `BART` are required. The simulation's aligned BART runner uses
`BART::wbart` / `BART::pbart` so that the mediator-draw scheme matches the applied
analysis; `lib_setup.R` controls this through `bart_engine` and `bart_nmed`.

### causalbkmr — install from the correct source, then patch

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("zc2326/causalbkmr")
```

> The maintainer is **zc2326** (Zilan Chai). Katrina Devick authored the original
> method paper. The package has no tagged releases, so record the commit SHA you
> installed if you need an exact pin.

**The installed package must be patched before use.** See
[Patched dependency](#patched-dependency) below — without the patch, the estimated
natural indirect effect is degenerate under the configuration used here.

### DeepMed and its Python environment

```r
remotes::install_github("siqixu/DeepMed")
```

DeepMed calls TensorFlow through `reticulate`. Build it in an isolated virtual
environment so it cannot collide with the R-native BART/BKMR sessions:

```bash
module load python/3.12.4      # 3.11 is NOT available on Longleaf;
                               # available: 3.9.6, 3.12.1, 3.12.2, 3.12.4
mkdir -p /work/users/o/d/odeaug/envs
python -m venv /work/users/o/d/odeaug/envs/deepmed_py
source /work/users/o/d/odeaug/envs/deepmed_py/bin/activate
pip install --upgrade pip
pip install tensorflow keras numpy pandas
deactivate
```

The venv path is referenced in `03_deepmed.sh`; edit `DEEPMED_VENV` there if you
place it elsewhere.

---

## Patched dependency

The analysis depends on a locally patched build of `causalbkmr`. In the version used
here, the internal cross-world prediction evaluates the mediator model at the
comparative exposure level rather than the reference level. Because the algorithm
derives the indirect effect by subtraction, this makes the estimated natural direct
effect equal the total effect and the estimated natural indirect effect degenerate.

- `fix_causalbkmr_astar.R` — applies the correction to the installed package.
- `verify_causalbkmr_fix.R` — runs the same data through the unpatched and patched
  versions and reports both sets of estimates side by side.

Source the fix after loading the package and before any BKMR-CMA fit. To confirm the
behaviour for yourself, run `verify_causalbkmr_fix.R` and compare the two outputs.

> **Scope.** The defect was observed under the configuration used in this study.
> The package's published Quick Start vignette runs a different configuration and
> does not exhibit it. Anyone reproducing this should run
> `verify_causalbkmr_fix.R` under their own settings rather than assuming the
> behaviour is universal.

---

## Repository layout

All files are at the repository root. They group by function as follows.

### Shared libraries

| File | Purpose |
| --- | --- |
| `lib_setup.R` | Configuration, scenario grid, method hyperparameters. Sourced by everything. |
| `lib_dgp.R` | Synthetic data-generating processes (scenarios 1–12). |
| `lib_dgp_plasmode.R` | Plasmode DGP (scenarios 13–14): mediator centering, direct-solve calibration, variance-budget assertion. |
| `lib_truth.R` | Monte Carlo computation of true TE / NDE / NIE. |
| `lib_methods.R` | BKMR and BART runners, `run_one_rep`, `compute_metrics`, per-replicate seeding. |
| `lib_deepmed.R` | DeepMed runner, isolated from the R-native methods. |
| `knots_helper.R` | Knot selection for the Gaussian predictive process approximation in BKMR. |

### Runners

| File | Purpose |
| --- | --- |
| `run_truth.R` | Truth for one scenario. |
| `run_scenario.R` | BART + BKMR for one scenario. |
| `run_method.R` | Method-parameterized driver (`METHOD` environment variable). |
| `run_bart_cma_aligned.R` | BART-CMA with the corrected mediator-draw pairing. |
| `run_plasmode.R` | Plasmode cells 13–14. |
| `run_deepmed.R` | DeepMed, median-split binarization. |
| `run_deepmed_v2.R` | DeepMed, Q25-vs-Q75 binarization. |
| `merge_split.R` | Merges per-replicate output into per-scenario files. |
| `aggregate.R` | Metrics, Excel workbook, performance plots. |

### Corrections and verification

| File | Purpose |
| --- | --- |
| `fix_causalbkmr_astar.R` | Patches the installed `causalbkmr`. |
| `verify_causalbkmr_fix.R` | Before/after comparison for the patch. |
| `patch_pairing.py` | One-time source edit producing the paired mediator-draw scheme. |
| `smoke_paired_bart.R`, `smoke_paired_bart_v2.R` | Confirm interval widths before and after the pairing correction. |
| `mediation_prob_scale.R` | Converts binary-outcome estimands from the probit latent scale to the probability scale. |
| `patch_applied.py` | One-time source edit to the applied analysis. |
| `check_bkmr.R` | BKMR diagnostic checks. |

> The two `patch_*.py` scripts edit source files **in place** and have already been
> applied to the code as committed. They are retained for provenance. Do not re-run
> them against the current files.

### SLURM

| File | Purpose |
| --- | --- |
| `01_truth.sh` | Truth, array 1–12. |
| `02a_bkmr.sh` | BKMR, array 1–12. |
| `02b_bart.sh` | BART, array 1–12. |
| `02c_bkmr_bign.sh` | BKMR rerun for the n = 1000 cells (extended walltime; these took 28–36 h each). |
| `02d_bart_aligned.sh` | BART rerun with the aligned mediator-draw scheme. |
| `02e_bkmr_pkg.sh` | BKMR via the packaged implementation. |
| `03_deepmed.sh` | DeepMed, median split. |
| `03_deepmed_quartile.sh` | DeepMed, quartile contrast. |
| `04_aggregate.sh` | Aggregation. |
| `launch_all.sh` | Submits the pipeline. |

### Applied analysis

| File | Purpose |
| --- | --- |
| `BART_CMA_Mediation.Rmd` | Parameterized applied analysis (`test` and `production` profiles). |
| `render_production.slurm`, `render_sensitivity.slurm` | Render under each confounder set. |
| `run_applied.sh`, `run_applied_sens.sh` | Primary and sensitivity confounder specifications. |
| `run_applied_diagnostics.sh` | Support-geometry and overlap diagnostics. |

The applied analysis reads 20 multiply imputed datasets held as sheets in a single
workbook and applies NHANES PFAS subsample weights. The imputed workbook is not
distributed here; see [Data](#data).

---

## Running the pipeline

Paths in the SLURM scripts assume the project root. Adjust them if you clone
elsewhere. Strip carriage returns after any edit on Windows, or bash and R will
fail with misleading parse errors:

```bash
sed -i 's/\r$//' *.sh *.R
chmod +x *.sh
```

Truth first, then the methods, then aggregation:

```bash
sbatch 01_truth.sh
sbatch 02a_bkmr.sh
sbatch 02b_bart.sh
sbatch 04_aggregate.sh
```

To rerun one scenario:

```bash
sbatch --array=7 02b_bart.sh
```

For a smoke test, set `R = 5` in `SIM_CONFIG` inside `lib_setup.R` and submit only
scenarios 1 and 11, then restore `R = 100`.

---

## Outputs

| Path | Contents |
| --- | --- |
| `results/truth/truth_<sid>.rds` | Ground truth per scenario. |
| `results/raw_scenario_<sid>.rds` | Per-replicate BART and BKMR results. |
| `results/deepmed/deepmed_<sid>.rds` | DeepMed results per scenario. |
| `results/errors/scenario_<sid>/` | Per-replicate error messages. |
| `results/all_metrics.rds` | Aggregated performance metrics. |
| `results/figures/` | Performance plots. |

---

## Reproducing the figures


`clean_results_comments_2_3.R` reads `all_metrics_fullT.csv` and writes the bias,
coverage, and interval-width figures together with a three-sheet workbook separating
the main two-method comparison, the labelled DeepMed exclusion, and the binary-arm
failures. Requires `dplyr`, `ggplot2`, and `writexl`.

---

## Data

The simulation and plasmode arms require no external data; the plasmode DGP draws on
NHANES structure through fitted coefficients supplied in `lib_dgp_plasmode.R`.

The applied analysis uses NHANES public-use files. The imputed analytic workbook is
not redistributed here.

<!-- TODO: state the NHANES cycle(s) used and either commit the download and
     imputation script or point to it. -->

---

## Runtime

BKMR-CMA is substantially more expensive than BART-CMA, by roughly an order of
magnitude per fit, and the n = 1000 cells required extended walltime allocations.
Budget accordingly before launching the full grid.

---

## Known issues

- `run_scenario.R` and `run_method.R` overlap; `run_deepmed.R` and `run_deepmed_v2.R`
  differ only in binarization; `smoke_paired_bart.R` and `smoke_paired_bart_v2.R`
  are successive versions. <!-- TODO: state which of each pair is current, or
  remove the superseded file. -->
- Two different per-replicate seeding formulas appeared across drivers during
  development. <!-- TODO: confirm a single formula is in force. -->
- BKMR componentwise variable selection produces a non-positive-definite kernel under
  the probit model and is disabled in the binary arm.
- One binary-outcome BKMR cell completed 97 of 100 replicates.

---


## Citation


