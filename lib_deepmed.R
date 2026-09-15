###########################################################
## lib_deepmed.R — DeepMed runner (TensorFlow via reticulate)
## Faithful to Chunks 6 and 11 of original Rmd.
## PATCHED: clears the TF/Keras session + garbage-collects
## after every DeepMed() call, in both the success and error
## paths. This stops the cross-rep memory accumulation that
## caused OUT_OF_MEMORY once ~40-50 reps had run.
###########################################################

## Release TensorFlow/Keras memory between reps.
## Uses Python-side clear_session so it targets whichever
## Keras backend DeepMed actually used (tf-keras here,
## because TF_USE_LEGACY_KERAS=1 is set below).
.clear_tf <- function() {
  try(reticulate::py_run_string(paste(
    "import gc",
    "try:",
    "    import tensorflow as tf",
    "    tf.keras.backend.clear_session()",
    "except Exception:",
    "    pass",
    "gc.collect()",
    sep = "\n"
  )), silent = TRUE)
  invisible(gc(FALSE))
}

run_deepmed <- function(dat, truth) {
  tryCatch({
    Sys.setenv(TF_USE_LEGACY_KERAS = "1")
    suppressPackageStartupMessages({
      library(reticulate); library(keras); library(DeepMed)
    })
    Z <- dat$Z; M <- dat$M; Y <- dat$Y; C <- dat$C
    scenario <- dat$scenario
    mixture_index <- rowSums(scale(Z))
    D    <- as.integer(mixture_index > median(mixture_index))
    X_dm <- cbind(C, Z)
    hyper_grid <- expand.grid(
      METHOD_CONFIG$dm_l1,
      METHOD_CONFIG$dm_layers,
      METHOD_CONFIG$dm_units
    )
    dm <- DeepMed::DeepMed(
      y          = Y,
      d          = D,
      m          = M,
      x          = X_dm,
      method     = METHOD_CONFIG$dm_method,
      hyper_grid = hyper_grid,
      epochs     = METHOD_CONFIG$dm_epochs,
      batch_size = METHOD_CONFIG$dm_batch_size,
      trim       = METHOD_CONFIG$dm_trim
    )
    res     <- dm$results
    nde_est <- res["dir.treat",   "effect"]; nde_se <- res["dir.treat",   "se"]
    nie_est <- res["indir.treat", "effect"]; nie_se <- res["indir.treat", "se"]
    te_est  <- res["total",       "effect"]; te_se  <- res["total",       "se"]

    out <- list(
      method = "DeepMed",
      NIE_est = nie_est, NIE_lo = nie_est - 1.96 * nie_se, NIE_hi = nie_est + 1.96 * nie_se,
      NDE_est = nde_est, NDE_lo = nde_est - 1.96 * nde_se, NDE_hi = nde_est + 1.96 * nde_se,
      TE_est  = te_est,  TE_lo  = te_est  - 1.96 * te_se,  TE_hi  = te_est  + 1.96 * te_se,
      success = TRUE
    )
    .clear_tf()
    out
  }, error = function(e) {
    .clear_tf()
    list(method = "DeepMed", success = FALSE, error = conditionMessage(e))
  })
}
