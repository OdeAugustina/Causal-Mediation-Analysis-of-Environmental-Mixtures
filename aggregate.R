#!/usr/bin/env Rscript
###########################################################
## aggregate.R — Merge all scenarios + DeepMed, compute metrics
## Faithful refactor of Chunks 12, 13, 14, 15.
###########################################################

source("R/lib_setup.R")
source("R/lib_methods.R")  # for compute_metrics

suppressPackageStartupMessages({
  library(ggplot2)
})

## ---- Merge DeepMed results into rep_results for a scenario ----
merge_deepmed <- function(sid, rep_results) {
  dm_file <- file.path(DEEPMED_DIR, sprintf("deepmed_%03d.rds", sid))
  if (!file.exists(dm_file)) {
    cat(sprintf("  [merge] no DeepMed file for scenario %d (continuing without)\n", sid))
    return(rep_results)
  }
  dm_reps <- readRDS(dm_file)
  for (r_str in names(dm_reps)) {
    r <- as.integer(r_str)
    if (r <= length(rep_results)) {
      rep_results[[r]]$results[["DeepMed"]] <- dm_reps[[r_str]]
    }
  }
  rep_results
}

## ---- Loop all scenarios and compute metrics ----
all_metrics <- list(); row_idx <- 0
methods_to_eval <- c("BKMR-CMA", "BART-CMA", "DeepMed")

for (i in seq_len(nrow(design_grid))) {
  sc  <- as.list(design_grid[i, ])
  sid <- sc$scenario_id

  truth_file <- file.path(TRUTH_DIR, sprintf("truth_%03d.rds", sid))
  raw_file   <- file.path(RESULTS_DIR, sprintf("raw_scenario_%03d.rds", sid))
  if (!file.exists(raw_file)) {
    cat(sprintf("Scenario %d: no BART+BKMR results. Skipping.\n", sid))
    next
  }
  if (!file.exists(truth_file)) {
    cat(sprintf("Scenario %d: no truth file. Skipping.\n", sid))
    next
  }
  truth <- readRDS(truth_file)
  reps  <- readRDS(raw_file)
  reps  <- merge_deepmed(sid, reps)

  for (method in methods_to_eval) {
    for (est in c("NIE", "NDE", "TE")) {
      row_idx <- row_idx + 1
      m <- compute_metrics(reps, truth, method, estimand = est)
      m$scenario_id <- sid
      m$linearity    <- sc$linearity
      m$interaction  <- sc$interaction
      m$outcome_type <- sc$outcome_type
      m$n <- sc$n; m$L <- sc$L; m$group <- sc$group
      all_metrics[[row_idx]] <- m
    }
  }
}

metrics_df <- bind_rows(lapply(all_metrics, as.data.frame))
saveRDS(metrics_df, file.path(RESULTS_DIR, "all_metrics.rds"))

cat(sprintf("Total metric rows: %d\n", nrow(metrics_df)))

if (nrow(metrics_df) == 0) {
  cat("No metrics produced; nothing to plot.\n")
  quit(save = "no", status = 0)
}

cat("\n===== NIE Performance =====\n")
nie <- metrics_df %>%
  dplyr::filter(estimand == "NIE") %>%
  dplyr::select(scenario_id, method, n_success, valid, bias, rmse, coverage, ci_width, mean_runtime) %>%
  dplyr::arrange(scenario_id, method)
print(as.data.frame(nie), digits = 3, row.names = FALSE)

## ---- Charts (Chunk 13) ----
plot_df <- metrics_df %>%
  dplyr::filter(valid == TRUE, estimand == "NIE") %>%
  dplyr::mutate(
    scenario_label = paste0(scenario_id, "\n", linearity, "\n", interaction, "\nn=", n),
    method = factor(method, levels = c("BKMR-CMA", "BART-CMA", "DeepMed"))
  )

method_colors <- c("BKMR-CMA" = "#534AB7", "BART-CMA" = "#1D9E75", "DeepMed" = "#D85A30")

mk_bar <- function(df, y, title, sub, ylab) {
  ggplot(df, aes(x = factor(scenario_id), y = .data[[y]], fill = method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = method_colors) +
    labs(title = title, subtitle = sub, x = "Scenario", y = ylab, fill = "Method") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top", axis.text.x = element_text(size = 8))
}

p1 <- mk_bar(plot_df, "bias", "Bias (NIE)", "Closer to 0 is better", "Bias") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")
p2 <- mk_bar(plot_df, "rmse", "RMSE (NIE)", "Lower is better", "RMSE")
p3 <- mk_bar(plot_df, "coverage", "Coverage Probability (95% CI for NIE)",
             "Target = 0.95 (red dashed line)", "Coverage") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_y_continuous(limits = c(0, 1.05))
p4 <- mk_bar(plot_df, "ci_width", "CI Width (NIE)",
             "Narrower is better (given adequate coverage)", "CI Width")
p5 <- mk_bar(plot_df, "mean_runtime", "Runtime per Replication",
             "Seconds (log scale)", "Runtime (seconds)") +
  scale_y_log10()

ggsave(file.path(FIGURES_DIR, "01_bias_nie.pdf"),    p1, width = 10, height = 6)
ggsave(file.path(FIGURES_DIR, "02_rmse_nie.pdf"),    p2, width = 10, height = 6)
ggsave(file.path(FIGURES_DIR, "03_coverage_nie.pdf"),p3, width = 10, height = 6)
ggsave(file.path(FIGURES_DIR, "04_ciwidth_nie.pdf"), p4, width = 10, height = 6)
ggsave(file.path(FIGURES_DIR, "05_runtime.pdf"),     p5, width = 10, height = 6)

cat(sprintf("\nFigures saved to: %s\n", FIGURES_DIR))

## ---- Excel export (Chunk 15) ----
if (!requireNamespace("writexl", quietly = TRUE)) {
  cat("writexl not installed; skipping xlsx export.\n")
} else {
  writexl::write_xlsx(metrics_df, file.path(RESULTS_DIR, "all_metrics.xlsx"))
  cat(sprintf("Excel saved to: %s\n", file.path(RESULTS_DIR, "all_metrics.xlsx")))
}

## ---- Ranking summary (Chunk 14) ----
cat("\n===== Method Ranking (NIE, averaged across scenarios) =====\n")
ranking <- metrics_df %>%
  dplyr::filter(estimand == "NIE", valid == TRUE) %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    scenarios    = dplyr::n(),
    avg_abs_bias = mean(abs(bias)),
    avg_rmse     = mean(rmse),
    avg_coverage = mean(coverage, na.rm = TRUE),
    avg_ci_width = mean(ci_width, na.rm = TRUE),
    avg_runtime  = mean(mean_runtime)
  ) %>%
  dplyr::arrange(avg_rmse)
print(as.data.frame(ranking), digits = 3, row.names = FALSE)

cat("\n[aggregate] complete.\n")
