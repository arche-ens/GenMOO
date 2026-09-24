#!/usr/bin/env Rscript
# Plot Pareto fronts (SAscore vs My_QED) across multiple series.
#
# To add a new series, just add one entry to `series` below:
#   name = list(dir = "results/<new_dir>", label = "<display name>")

library(ggplot2)

# ---- Configuration (extend here for more series) ----
series <- list(
  test_new  = list(dir = "results/test_new_properties_pareto_results",
                      label = "test_new"),
  test_new1 = list(dir = "results/test_new1_properties_pareto_results",
                      label = "test_new1"),
  test_new2 = list(dir = "results/test_new2_properties_pareto_results",
                   label = "test_new2"),
  test_new3 = list(dir = "results/test_new3_properties_pareto_results",
                   label = "test_new3")
)

# ---- 1. Compare the Pareto optimum (rank 1) across all series ----
read_optimum <- function(s) {
  d <- read.csv(file.path(s$dir, "pareto_front_001.csv"))
  data.frame(series  = s$label,
             SAscore = d$SAscore,
             My_QED  = d$My_QED)
}

optima <- do.call(rbind, lapply(series, read_optimum))

p_compare <- ggplot(optima, aes(x = SAscore, y = My_QED, color = series)) +
  geom_point(size = 2.5, ) +
  labs(title = "Pareto Optimum (Rank 1): SAscore vs My_QED",
       x = "SAscore", y = "My_QED", color = "Series") +
  theme_minimal()

ggsave(file.path("results", "pareto_optima_comparison.png"),
       p_compare, width = 7, height = 5)

# # ---- 2. All Pareto fronts (colored by pareto_rank) for each series ----
# for (s in series) {
#   d <- read.csv(file.path(s$dir, "all_solutions.csv"))
#   d$pareto_rank <- factor(d$pareto_rank)
# 
#   p_all <- ggplot(d, aes(x = SAscore, y = My_QED, color = pareto_rank)) +
#     geom_point(size = 1.5) +
#     labs(title = paste0(s$label, " - All Pareto Fronts"),
#          x = "SAscore", y = "My_QED", color = "Pareto rank") +
#     theme_minimal()
# 
#   ggsave(file.path(s$dir, "pareto_fronts_all.png"),
#          p_all, width = 7, height = 5)
# }

cat("Saved:", file.path("results", "pareto_optima_comparison.png"), "\n")
# for (s in series) {
#   cat("Saved:", file.path(s$dir, "pareto_fronts_all.png"), "\n")
# }
