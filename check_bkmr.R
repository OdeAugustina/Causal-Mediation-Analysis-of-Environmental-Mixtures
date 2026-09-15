#!/usr/bin/env Rscript
## check_bkmr.R -- did the BKMR run work?
cat(sprintf("\n%-5s %6s %9s %9s %9s %9s %9s %9s  %s\n",
    "cell","n_ok","TE","(true)","NDE","(true)","NIE","(true)","scale"))
cat(strrep("-", 92), "\n")
for (s in sprintf("%03d", 1:12)) {
  f <- sprintf("results/raw_scenario_%s_bkmr_pkg.rds", s)
  if (!file.exists(f)) { cat(sprintf("%-5s  MISSING\n", s)); next }
  x  <- readRDS(f)
  tr <- readRDS(sprintf("results/truth/truth_%s.rds", s))
  g <- function(k) sapply(x, function(r) {
    v <- r$results[["BKMR-CMA"]][[k]]
    if (is.null(v) || length(v) != 1) NA_real_ else as.numeric(v) })
  nie <- g("NIE_est"); ok <- sum(is.finite(nie))
  if (ok == 0) {
    e <- x[[1]]$results[["BKMR-CMA"]]$error
    cat(sprintf("%-5s %6d  ALL FAILED: %s\n", s, 0, substr(e, 1, 55))); next }
  b <- x[[1]]$results[["BKMR-CMA"]]
  cat(sprintf("%-5s %6d %+9.3f %+9.3f %+9.3f %+9.3f %+9.3f %+9.3f  %s\n",
      s, ok, mean(g("TE_est"),na.rm=TRUE), tr$TE,
      mean(g("NDE_est"),na.rm=TRUE), tr$NDE,
      mean(nie,na.rm=TRUE), tr$NIE,
      if (is.null(b$scale)) "?" else b$scale))
}
cat("\n--- the four checks ---\n")
bad <- 0
for (s in sprintf("%03d", 1:12)) {
  f <- sprintf("results/raw_scenario_%s_bkmr_pkg.rds", s)
  if (!file.exists(f)) next
  x <- readRDS(f); tr <- readRDS(sprintf("results/truth/truth_%s.rds", s))
  g <- function(k) mean(sapply(x, function(r) {
    v <- r$results[["BKMR-CMA"]][[k]]
    if (is.null(v) || length(v) != 1) NA_real_ else as.numeric(v) }), na.rm=TRUE)
  nie <- g("NIE_est"); nde <- g("NDE_est"); te <- g("TE_est")
  if (!is.finite(nie)) next
  msg <- character(0)
  if (abs(nie) < 0.01)              msg <- c(msg, "indirect effect is ~zero (old bug)")
  if (abs(nde - te) < 0.01)         msg <- c(msg, "direct effect equals total (old bug)")
  r <- nie / tr$NIE
  if (r > 2.5)                      msg <- c(msg, sprintf("ratio %.1f - wrong units", r))
  if (length(msg)) { cat(sprintf("  s%s: %s\n", s, paste(msg, collapse="; "))); bad <- bad + 1 }
}
if (bad == 0) cat("  all clear\n")
cat("\ns001 indirect effect should be near +1.29 (your hand-written version gave +1.27)\n")
