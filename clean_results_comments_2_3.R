#!/usr/bin/env Rscript
## =====================================================================
## clean_results_comments_2_3.R  (Longleaf-ready)
## Applies committee comments 2 & 3 to the aggregated metrics and
## writes sendable tables + figures. Post-processing only -- it does
## NOT re-run the simulation.
##
##  Comment 2 -- DeepMed is not a valid comparator on the continuous
##    scenarios: it estimates a binary-treatment (median-split) contrast,
##    not the continuous 25th->75th percentile mixture effect. Removed there
##    and retained separately with that reason (N/A-by-design), not deleted.
##
##  Comment 3 -- the binary-outcome scenarios (9, 10) are pulled from the
##    headline comparison (all three methods fail). Retained as a documented-
##    failure table so the negative result is on the record.
##
## Run: interactively in RStudio (source this) or via sbatch. Needs dplyr,
##      ggplot2, writexl (and readxl only if reading an .xlsx input).
## ====================================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(writexl)
})

## --- locate the aggregated metrics (don't assume a fixed name) ---------
find_input <- function() {
  cands <- c("results/all_metrics_fullT.csv", "all_metrics_fullT.csv",
             "results/all_metrics.csv",      "all_metrics.csv",
             "results/all_metrics.rds",      "all_metrics.rds",
             "results/all_metrics.xlsx",     "all_metrics.xlsx")
  hit <- cands[file.exists(cands)][1]
  if (is.na(hit)) stop("Could not find the aggregated metrics. Put all_metrics_fullT.csv ",
                        "(or all_metrics.rds / .xlsx) in the working dir or results/.")
  hit
}
read_metrics <- function(path) {
  cat(sprintf("Reading metrics from: %s\n", path))
  if (grepl("\\.rds$", path)) return(as.data.frame(readRDS(path)))
  if (grepl("\\.xlsx$", path)) { suppressPackageStartupMessages(library(readxl))
    return(as.data.frame(readxl::read_excel(path))) }
  read.csv(path, stringsAsFactors = FALSE)
}

m <- read_metrics(find_input())
m$outcome_type <- tolower(as.character(m$outcome_type))

req <- c("method", "estimand", "scenario_id", "bias", "coverage", "ci_width", "outcome_type")
miss <- setdiff(req, names(m))
if (length(miss)) stop(paste("Input is missing columns: ", paste(miss, collapse = ", "),
  "Expected the aggregate.R metrics schema. Got: ", paste(names(m), collapse = ", ")))

is_bin  <- m$outcome_type == "binary"   # scenarios 9, 10
is_cont <- !is_bin                     # 1-8, 11, 12

## --- Comment 2: DeepMed excluded on continuous, labeled ------------
deepmed_excluded <- m %>%
  filter(method == "DeepMed", is_cont) %>%
  mutate(exclusion_reason =
           "N/A by design: DeepMed estimates a binary-treatment (median-split) contrast, not the continuous 25th-75th percentile mixture effect")

## --- Comment 3: binary arm pulled, kept as failure log ---------------
binary_failures <- m %>%
  filter(is_bin) %>%
  mutate(failure_reason = case_when(
    method == "BART-CMA" ~ "Overconfident: near-zero CI width on the probability scale; coverage 0",
    method == "BKMR-CMA" ~ "Non-convergent probit MCMC (25/100 then 0/100 fits)",
    method == "DeepMed"  ~ "Estimand mismatch (median-split) + binary outcome; TE coverage 0",
    TRUE ~ NA_character_))

## --- Main comparison: continuous scenarios, BKMR vs BART only --------
main_comparison <- m %>%
  filter(is_cont, method %in% c("BKMR-CMA", "BART-CMA")) %>%
  arrange(match(estimand, c("NIE", "NDE", "TE")), scenario_id, method)

## --- output folder --------------------------------------------------
OUTDIR <- "results/comments_2_3"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
p <- function(f) file.path(OUTDIR, f)

rnd <- function(df) {
  keep_int <- c("scenario_id", "n", "L", "n_success")
  for (col in names(df)[sapply(df, is.numeric)])
    if (!col %in% keep_int) df[[col]] <- round(df[[col]], 3)
  df
}

write_xlsx(list(
  Main_Comparison       = rnd(main_comparison),
  DeepMed_Excluded_Cont = rnd(deepmed_excluded),
  Binary_Arm_Failures   = rnd(binary_failures)
), p("results_clean.xlsx"))

write.csv(main_comparison, p("metrics_main_comparison.csv"), row.names = FALSE)
write.csv(binary_failures, p("metrics_binary_failures.csv"), row.names = FALSE)

## --- cleaned figures: BKMR vs BART on continuous scenarios ----------
plot_df <- main_comparison
plot_df$scenario_id <- factor(plot_df$scenario_id, levels = sort(unique(plot_df$scenario_id)))
plot_df$estimand    <- factor(plot_df$estimand, levels = c("NIE", "NDE", "TE"))
pal <- c("BKMR-CMA" = "#5B4FC7", "BART-CMA" = "#1AA179")

gg_bias <- ggplot(plot_df, aes(scenario_id, bias, fill = method)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  facet_wrap(~estimand, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = pal) +
  labs(title = "Bias by scenario (continuous outcomes; DeepMed excluded)",
       subtitle = "Closer to 0 is better", x = "Scenario", y = "Bias", fill = "Method") +
  theme_minimal(base_size = 12)
ggsave(p("fig_bias_clean.pdf"), gg_bias, width = 8, height = 9)

gg_cov <- ggplot(plot_df, aes(scenario_id, coverage, fill = method)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_hline(yintercept = 0.95, linetype = 2, color = "red") +
  facet_wrap(~estimand, ncol = 1) +
  scale_fill_manual(values = pal) + coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Coverage by scenario (continuous outcomes; DeepMed excluded)",
       subtitle = "Target = 0.95 (red line)", x = "Scenario", y = "Coverage", fill = "Method") +
  theme_minimal(base_size = 12)
ggsave(p("fig_coverage_clean.pdf"), gg_cov, width = 8, height = 9)

  geom_col(position = position_dodge(0.8), width = 0.75) +

  facet_wrap(~estimand, ncol = 1, scales = "free_y") +

  scale_fill_manual(values = pal) +

  labs(title = "CI width by scenario (continuous outcomes; DeepMed excluded)",

       subtitle = "Narrower is better given adequate coverage",

       x = "Scenario", y = "CI width", fill = "Method") +

  theme_minimal(base_size = 12)

ggsave(p("fig_ci_width_clean.pdf"), gg_width, width = 8, height = 9)



cat("\nWrote to", OUTDIR, ":\n")

cat("  results_clean.xlsx, metrics_main_comparison.csv, metrics_binary_failures.csv,\n")

cat("  fig_bias_clean.pdf, fig_coverage_clean.pdf, fig_ci_width_clean.pdf\n")

cat("Main comparison scenarios:", paste(sort(unique(main_comparison$scenario_id)), collapse = ", "),

    "| methods: BKMR-CMA, BART-CMA\n")

