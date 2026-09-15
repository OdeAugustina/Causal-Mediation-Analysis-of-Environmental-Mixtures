# Running All Three Methods on Longleaf — Step-by-Step Runbook

**Goal:** produce clean per-scenario results for **BKMR-CMA**, **BART-CMA** (run in parallel), and **DeepMed**, across all 12 scenarios × 100 reps, then aggregate. Every stage below has a *gate* — a cheap check that stops you before an expensive job fails.

**The golden rule of this runbook:** never launch a multi-day array job until the exact same code has produced one good `.rds` in a 3-rep smoke test. Smoke tests cost minutes; failed production jobs cost days.

---

## The mental model

| Stage | Script | Where it runs | Output |
|---|---|---|---|
| 0. Files in place | — | login node | `R/`, `SLURM/` populated |
| 1. Verify wiring | (R console) | login node | confidence, no output |
| 2. Interactive node | `srun` | → compute node | a shell you can test in |
| 3. Smoke BKMR + BART | `run_method.R` | compute node | tiny test `.rds` |
| 4. Smoke DeepMed | `run_deepmed.R` | compute node | tiny test `.rds` |
| 5. (Truth precompute) | `00_truth.sh` | array job | `results/truth/` |
| 6. Production launch | `01_bkmr.sh` + `02_bart.sh` + `03_deepmed.sh` | array jobs | `results/{bkmr,bart,deepmed}/` |
| 7. Monitor | `squeue` / logs | login node | — |
| 8. Aggregate | `04_aggregate.sh` | array/single job | final metrics |

BKMR and BART run **in parallel** because they're two independent `sbatch` submissions — SLURM schedules them at the same time on different nodes.

---

## Stage 0 — Put the files in place (login node)

Five files belong in your project:

```
sim_proj/R/run_method.R        # NEW — BKMR + BART driver
sim_proj/R/run_deepmed.R       # NEW — DeepMed driver
sim_proj/SLURM/01_bkmr.sh      # NEW
sim_proj/SLURM/02_bart.sh      # NEW
sim_proj/SLURM/03_deepmed.sh   # REPLACES the old one
```

Upload them via **OnDemand → Files**, then in MobaXterm:

```bash
cd /work/users/o/d/odeaug/sim_proj
```

**⚠️ The single most common Windows→Longleaf error: carriage returns.** Files edited or uploaded from Windows carry `\r` line endings, which make bash scripts fail with `bad interpreter` and make R scripts throw cryptic parse errors. Strip them on every file, every time:

```bash
sed -i 's/\r$//' R/run_method.R R/run_deepmed.R SLURM/01_bkmr.sh SLURM/02_bart.sh SLURM/03_deepmed.sh
chmod +x SLURM/*.sh
```

**Don't overwrite your truth script.** Your original scaffold numbered truth as `01_*`. My method scripts are also `01`/`02`/`03`, so rename truth out of the way to avoid confusion, and retire the old combined scenario script (it's now split into bkmr + bart):

```bash
# only if these exist from the original scaffold:
mv SLURM/01_truth.sh    SLURM/00_truth.sh    2>/dev/null || true
mv SLURM/02_scenario.sh SLURM/_retired_02_scenario.sh 2>/dev/null || true
```

---

## Stage 1 — Verify your project wiring (login node, 2 minutes)

This is the gate that catches *every* name-mismatch error before any compute happens. Start R on the login node (this is read-only inspection, no heavy compute, so it's fine here):

```bash
module load r/4.5.0     # if this errors: `module avail r` and use what's listed
R
```

Then paste, one block, and read the output:

```r
source("R/lib_setup.R")

# (a) scenario grid — must exist and have 12 rows
exists("SIM_GRID"); nrow(SIM_GRID)
names(SIM_GRID)                       # expect columns incl. n, outcome_type

# (b) method hyperparameters — must exist
exists("METHOD_CONFIG"); str(METHOD_CONFIG)

# (c) DGP — must return a list with Y, D, M, Z, X
source("R/lib_dgp.R")
d1 <- generate_data(scenario = SIM_GRID[1, ], n = SIM_GRID$n[1])
str(d1, max.level = 1)

# (d) method wrappers + the mock flag
source("R/lib_methods.R")
exists("run_bkmr_cma"); exists("run_bart_cma"); exists("USE_MOCK"); USE_MOCK
```

**Decision rules:**
- `SIM_GRID` doesn't exist → find the real name (`ls()`), then change `SIM_GRID` in `run_method.R` and `run_deepmed.R` to match.
- The hyperparameters are in `SIM_CONFIG`, not `METHOD_CONFIG` → change `METHOD_CONFIG` → `SIM_CONFIG` in `run_method.R`.
- `generate_data` returns different field names (e.g. `C` instead of `X`, or `A` instead of `D`) → adjust the `dat$...` references in both drivers. **This is the most likely place the drivers need a tweak** — your wrappers evolved through a few naming conventions.
- `SIM_GRID` has no `n` column (n is fixed or stored elsewhere) → change `scen$n` to wherever n lives.
- `USE_MOCK` is `TRUE` → fine, the drivers force it to `FALSE` at runtime. (If `USE_MOCK` doesn't exist at all, delete the `USE_MOCK <- FALSE` line from the drivers.)

Quit R without saving:

```r
q(save = "no")
```

Don't proceed until block (c) prints `Y`, `D`, `M`, `Z`, `X` (or you've mapped the real names into the drivers).

---

## Stage 2 — Grab an interactive compute node

**Never run methods on the login node** — it's shared and jobs there get killed. Get a small interactive node for smoke testing:

```bash
srun --ntasks=1 --cpus-per-task=4 --mem=16g --time=2:00:00 \
     --partition=general --pty bash
```

Your prompt changes to a compute node (e.g. `[odeaug@c0518 ~]$`). Then:

```bash
cd /work/users/o/d/odeaug/sim_proj
module load r/4.5.0
```

---

## Stage 3 — Smoke test BKMR and BART (≈5–15 min)

Run the **real production driver** at tiny scale (3 reps, 3 cores). Same code path as production — that's what makes this a valid test.

**BART first** (it's fast, so failures surface quickly):

```bash
SIM_PROJ=/work/users/o/d/odeaug/sim_proj METHOD=bart DM_REPS=3 N_CORES=3 \
  Rscript R/run_method.R 1
```

**Then BKMR** (slower — give it a few minutes):

```bash
SIM_PROJ=/work/users/o/d/odeaug/sim_proj METHOD=bkmr DM_REPS=3 N_CORES=3 \
  Rscript R/run_method.R 1
```

Each should end with a `Done ... 3/3 reps ok` line. Now **inspect the output** — this is the real gate:

```bash
Rscript -e 'x <- readRDS("results/bart/bart_scenario_01.rds"); print(x); cat("\nNA effects:", sum(is.na(x$effect)), "of", nrow(x), "\n")'
Rscript -e 'x <- readRDS("results/bkmr/bkmr_scenario_01.rds"); print(x); cat("\nNA effects:", sum(is.na(x$effect)), "of", nrow(x), "\n")'
```

**What "good" looks like:** 15 rows each (3 reps × 5 estimands), `effect` values that are finite and in a sane range, `ci_lower < effect < ci_upper`, and the `error` column all `NA`. If every `effect` is `NA`, open the `error` column — it carries the exact message from the wrapper.

> If `effect` values look like random noise around 0.5 with tidy ±0.2 intervals, `USE_MOCK` is still on — the mock estimator is running. Confirm Stage 1 (d) and that the driver's `USE_MOCK <- FALSE` line is present.

---

## Stage 4 — Smoke test DeepMed (the historical blocker)

Two sub-steps, because DeepMed has the most moving parts.

### 4a — Confirm the raw `res$results` shape (the one thing still unconfirmed)

Still on the interactive node:

```bash
Rscript -e '
library(reticulate); use_condaenv("/work/users/o/d/odeaug/envs/deepmed", required = TRUE)
library(tensorflow); library(DeepMed); library(doParallel); registerDoParallel(cores = 2)
set.seed(1); n <- 300
x <- matrix(rnorm(n*2), n, 2); d <- rbinom(n, 1, 0.5)
m <- 0.5*d + x[,1] + rnorm(n); y <- 0.3*m + 0.2*d + x[,2] + rnorm(n)
res <- DeepMed(y, d, m, x, method="DNN",
               hyper_grid = expand.grid(c(0,0.05), 1, 10),
               epochs = 50, batch_size = 50, trim = 0.05)
cat("--- str(res) ---\n"); str(res)
cat("--- res$results ---\n"); print(res$results)
'
```

You want to see `res$results` as a 5-row table with `effect`, `se`, `pval`. **Note the row order** — that's the estimand mapping `aggregate.R` needs. If it's not 5 rows, or columns differ, tell me and I'll adjust `run_deepmed.R` before you scale.

> Expect a wall of CUDA "could not find drivers" / TF-TRT warnings — those are harmless on CPU nodes. Only the `res$results` table matters.

### 4b — Smoke the DeepMed driver itself

```bash
SIM_PROJ=/work/users/o/d/odeaug/sim_proj DEEPMED_ENV=/work/users/o/d/odeaug/envs/deepmed \
  DM_OUTER=2 DM_INNER=2 DM_REPS=2 DM_EPOCHS=50 \
  Rscript R/run_deepmed.R 1
```

Check it:

```bash
Rscript -e 'x <- readRDS("results/deepmed/dm_scenario_01.rds"); print(x)'
```

**If it hangs with no epoch output for several minutes:** that's the TensorFlow-in-a-fork-pool issue I flagged. Don't fight it — drop the inner parallelism by re-running with `DM_INNER=1`, and if you still see hangs at production scale, ping me for the job-array variant of DeepMed (one rep per task, which sidesteps nesting entirely). BKMR and BART have no such risk.

When you're done smoke testing, release the node:

```bash
exit
```

(You're back on the login node.)

---

## Stage 5 — Precompute truth (if not already done)

The method drivers produce *estimates*; bias and coverage need the *true* per-scenario effects. If `results/truth/` isn't already populated from earlier work:

```bash
sbatch SLURM/00_truth.sh
```

This is independent of the method runs and can run concurrently — it doesn't block Stage 6. (If you don't have `00_truth.sh` handy, say so and I'll regenerate it against your `lib_truth.R`.)

---

## Stage 6 — Launch production: all three methods

Now the payoff. BKMR and BART go out together and run in parallel:

```bash
cd /work/users/o/d/odeaug/sim_proj
sbatch SLURM/01_bkmr.sh      # 12 × 50 cores, up to 3 days (MCMC)
sbatch SLURM/02_bart.sh      # 12 × 40 cores, up to 12h (fast)
sbatch SLURM/03_deepmed.sh   # 12 × 80 cores — ONLY after Stage 4 passed
```

Each `sbatch` prints a job ID like `Submitted batch job 12345678`. You now have three independent array jobs (12 tasks each) running concurrently as the scheduler finds room.

> Hold `03_deepmed.sh` until Stage 4 is green. BKMR and BART you can fire immediately once Stage 3 passes.

---

## Stage 7 — Monitor

```bash
squeue -u odeaug                                  # everything queued/running
squeue -u odeaug -t RUNNING | wc -l               # how many tasks actually running
```

Watch a live log (replace with a real scenario number):

```bash
tail -f results/bart_logs/bart_1_*.out
```

Check finished tasks for exit codes and walltime:

```bash
sacct -u odeaug --format=JobID,JobName%14,State,Elapsed,MaxRSS,ExitCode -X
```

**Reading it:** `COMPLETED` is good. `FAILED` or `OUT_OF_MEMORY` → see the table below. `TIMEOUT` → the walltime was too short; bump `--time` and resubmit (Stage 8 resume makes this cheap).

---

## Stage 8 — Partial failures and resume (this is why you won't lose work)

Every rep writes its own part file the moment it finishes. So if a scenario dies at rep 73, the first 72 are safe on disk. **To resume, just resubmit the same script** — the driver skips reps that already have a part file and a scenario that already has its final `.rds`:

```bash
sbatch SLURM/01_bkmr.sh      # picks up exactly where it stopped
```

To see which scenarios actually finished:

```bash
ls -1 results/bkmr/bkmr_scenario_*.rds 2>/dev/null | wc -l    # want 12
ls -1 results/bart/bart_scenario_*.rds 2>/dev/null | wc -l    # want 12
ls -1 results/deepmed/dm_scenario_*.rds 2>/dev/null | wc -l   # want 12
```

To audit failed reps inside a finished scenario:

```bash
Rscript -e 'x <- readRDS("results/bkmr/bkmr_scenario_03.rds"); print(table(is.na(x$error))); print(unique(na.omit(x$error)))'
```

---

## Stage 9 — Aggregate

Once all 36 `.rds` files exist (12 per method):

```bash
sbatch SLURM/04_aggregate.sh
```

**One thing to verify in `aggregate.R` first:** it must read from the three separate directories `results/bkmr/`, `results/bart/`, `results/deepmed/`, and read coverage from `ci_lower`/`ci_upper` for BKMR/BART while deriving DeepMed's CIs as `effect ± 1.96·se`. If your original `aggregate.R` expected a single combined `results/scenario/` folder, it needs that adjustment — tell me and I'll rebuild it to match this layout.

Final sanity check on the aggregated metrics: bias near 0 for well-specified scenarios, coverage near 0.95, RMSE finite, and no method showing all-`NA` for any estimand.

---

## Quick error → fix reference

| Symptom | Cause | Fix |
|---|---|---|
| `bad interpreter: No such file or directory` | Windows `\r` line endings | `sed -i 's/\r$//' <file>` |
| R `unexpected ... in "..."` parse error on source | Same `\r` problem in `.R` files | `sed -i 's/\r$//' R/*.R` |
| `object 'SIM_GRID' not found` | Grid named differently | `ls()` after sourcing `lib_setup.R`; rename in drivers |
| `argument "config" ... ` / wrong field | `METHOD_CONFIG` vs `SIM_CONFIG` | Swap the name in `run_method.R` |
| All `effect` = tidy noise ~0.5 | `USE_MOCK` still TRUE | Confirm driver's `USE_MOCK <- FALSE` line |
| `argument "e.y" is missing` | Old BKMR wrapper | Already fixed — ensure `lib_methods.R` passes `e.y = NULL` |
| BART predictions error / can't predict | Missing `keeptrees` | Ensure `bart(..., keeptrees = TRUE)` in wrapper |
| `incorrect number of dimensions` (DeepMed) | 1-row `hyper_grid` | Use the 8-row grid (already in `run_deepmed.R`) |
| DeepMed hangs, no epochs print | TF inside fork pool | `DM_INNER=1`, or switch to the job-array DeepMed variant |
| `ModuleNotFoundError` / TF not found in job | conda env not reached | Confirm `DEEPMED_ENV` path; env is loaded via `use_condaenv()` inside R, not `conda activate` |
| `OUT_OF_MEMORY` | `--mem` too low | Raise `--mem` in the `.sh`; resubmit (resumes) |
| `TIMEOUT` | `--time` too short | Raise `--time`; resubmit (resumes) |
| Job stuck `PENDING` for hours | 50 cpus/node is a big ask | Lower `--cpus-per-task`; the driver auto-reads the new count |

---

### The three gates that prevent ~90% of failures
1. **Stage 1** — `generate_data()` returns `Y/D/M/Z/X` and `SIM_GRID`/`METHOD_CONFIG` exist by those names.
2. **Stage 3/4** — a 3-rep `.rds` exists with finite effects and an empty `error` column, for each method.
3. **Stage 4a** — `res$results` is a 5-row `effect/se/pval` table.

Clear those three and the production arrays are almost always smooth sailing.
