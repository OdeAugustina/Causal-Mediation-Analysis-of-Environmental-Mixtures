#!/usr/bin/env Rscript
## merge_split.R — recombine method-split raw files into the single
## raw_scenario_%03d.rds aggregate.R expects (both methods per rep).
source("R/lib_setup.R")
sids <- sort(unique(design_grid$scenario_id))
for (sid in sids) {
  fb <- file.path(RESULTS_DIR, sprintf("raw_scenario_%03d_bkmr.rds", sid))
  fa <- file.path(RESULTS_DIR, sprintf("raw_scenario_%03d_bart.rds", sid))
  fo <- file.path(RESULTS_DIR, sprintf("raw_scenario_%03d.rds",      sid))
  miss <- c(bkmr = !file.exists(fb), bart = !file.exists(fa))
  if (any(miss)) {
    cat(sprintf("[%03d] skip - missing: %s\n", sid, paste(names(miss)[miss], collapse=", ")))
    next
  }
  rb <- readRDS(fb); ra <- readRDS(fa)
  if (length(rb) != length(ra)) {
    cat(sprintf("[%03d] skip - length mismatch (bkmr=%d bart=%d)\n", sid, length(rb), length(ra)))
    next
  }
  ra_by  <- setNames(ra, vapply(ra, function(x) x$rep_id, integer(1)))
  merged <- vector("list", length(rb))
  for (i in seq_along(rb)) {
    rid <- rb[[i]]$rep_id
    a   <- ra_by[[as.character(rid)]]
    merged[[i]] <- list(
      rep_id  = rid,
      results = list("BKMR-CMA" = rb[[i]]$results[["BKMR-CMA"]],
                     "BART-CMA" = a$results[["BART-CMA"]]),
      data = NULL)
  }
  saveRDS(merged, fo)
  cat(sprintf("[%03d] merged -> %s\n", sid, basename(fo)))
}
cat("Merge complete.\n")
