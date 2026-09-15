#!/usr/bin/env Rscript
###########################################################
## verify_causalbkmr_fix.R
##
## Controlled comparison on ONE replicate of scenario 1:
##   (1) mediation.bkmr with the stock package
##   (2) mediation.bkmr with YaMastar.SamplePred patched (astar fix)
##   (3) truth
##   (4) the archived manual-implementation result, same seed
##
## Only ONE line differs between (1) and (2), so any change in NIE is
## attributable to the astar bug and nothing else.
##
## Usage:  Rscript R/verify_causalbkmr_fix.R
###########################################################

suppressPackageStartupMessages({ library(bkmr); library(causalbkmr) })
source("R/lib_setup.R"); source("R/lib_dgp.R"); source("R/lib_methods.R")

SID <- 1L; REP <- 1L
sc <- as.list(design_grid[design_grid$scenario_id == SID, ])
tr <- readRDS(file.path(TRUTH_DIR, sprintf("truth_%03d.rds", SID)))

set.seed(as.integer(SIM_CONFIG$seed + SID * 1000L + REP))
d <- generate_data(sc)
Xp <- matrix(colMeans(d$C), nrow = 1)

cat("fitting three kmbayes models (iter = 3000)...\n")
fm <- kmbayes(y = d$M, Z = d$Z, X = d$C, iter = 3000, verbose = FALSE, varsel = TRUE)
fy <- kmbayes(y = d$Y, Z = cbind(d$Z, M = d$M), X = d$C, iter = 3000,
              verbose = FALSE, varsel = TRUE)
ft <- kmbayes(y = d$Y, Z = d$Z, X = d$C, iter = 3000, verbose = FALSE, varsel = TRUE)

sel <- seq(METHOD_CONFIG$bkmr_burnin + 1, METHOD_CONFIG$bkmr_iter, by = 5)
K   <- 30L

run <- function(tag) {
  m <- causalbkmr::mediation.bkmr(
    a = tr$a, astar = tr$astar, e.m = NULL, e.y = NULL,
    fit.m = fm, fit.y = fy, fit.y.TE = ft,
    X.predict.M = Xp, X.predict.Y = Xp,
    alpha = 0.05, sel = sel, seed = 123, K = K)
  e <- as.data.frame(m$est)
  g <- function(r) e[match(r, trimws(rownames(e))), "mean"]
  c(TE = g("TE"), NDE = g("NDE"), NIE = g("NIE"))
}

cat("\n--- (1) STOCK PACKAGE ---\n"); before <- run()
print(round(before, 4))

cat("\n--- applying patch ---\n")
source("R/fix_causalbkmr_astar.R")

cat("\n--- (2) PATCHED ---\n"); after <- run()
print(round(after, 4))

cat("\n=================== SUMMARY (scenario 1, rep 1) ===================\n")
cat(sprintf("%-22s %9s %9s %9s\n", "", "TE", "NDE", "NIE"))
cat(sprintf("%-22s %+9.4f %+9.4f %+9.4f\n", "truth", tr$TE, tr$NDE, tr$NIE))
cat(sprintf("%-22s %+9.4f %+9.4f %+9.4f\n", "stock package",
            before["TE"], before["NDE"], before["NIE"]))
cat(sprintf("%-22s %+9.4f %+9.4f %+9.4f\n", "patched (astar fix)",
            after["TE"], after["NDE"], after["NIE"]))

of <- file.path(RESULTS_DIR, sprintf("raw_scenario_%03d_bkmr.rds", SID))
if (file.exists(of)) {
  o <- readRDS(of)[[REP]]$results[["BKMR-CMA"]]
  cat(sprintf("%-22s %+9.4f %+9.4f %+9.4f\n", "manual impl (archived)",
              o$TE_est, o$NDE_est, o$NIE_est))
}
cat("\nDiagnostic: stock NDE - stock TE = %+.4f (0 => bug present)\n")
cat(sprintf("            stock NDE - stock TE = %+.6f\n", before["NDE"] - before["TE"]))
cat(sprintf("            patched NIE vs truth  = %+.4f vs %+.4f\n",
            after["NIE"], tr$NIE))
cat("==================================================================\n")
