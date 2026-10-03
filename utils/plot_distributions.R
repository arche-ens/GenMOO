#!/usr/bin/env Rscript
# Plot distributions of `avg`, `SAscore`, `QED_w` across multiple files.
#
# To add a new file, add one entry to `files` below:
#   name = "path/to/properties.csv"

library(ggplot2)

# ---- Configuration (extend here for more files) ----
files <- list(
  Gen_0 = "results/new_molecules_active_properties.csv",
  Gen_1 = "results/new_molecules1_active_properties.csv",
  Gen_2 = "results/new_molecules2_active_properties.csv",
  Gen_3 = "results/new_molecules3_active_properties.csv",
  Gen_4 = "results/new_molecules4_active_properties.csv",
  Gen_5 = "results/new_molecules5_active_properties.csv"
)


factors <- c("avg", "SAscore", "QED_w")

# ---- Read & merge ----
df <- do.call(rbind, lapply(names(files), function(s) {
  d <- read.csv(files[[s]])
  missing <- setdiff(factors, names(d))
  if (length(missing) > 0) {
    stop(sprintf("Missing columns in %s: %s",
                 files[[s]], paste(missing, collapse = ", ")))
  }
  d <- d[, factors, drop = FALSE]
  d$series <- s
  d
}))
df$series <- factor(df$series, levels = names(files))

# ---- Plot per factor ----
series_names <- names(files)

for (f in factors) {
  p_violin <- ggplot(df, aes(x = series, y = .data[[f]], fill = series)) +
    geom_violin() +
    labs(title = paste0(f, " - violin"),
         x = "Series", y = f, fill = "Series") +
    theme_minimal()

  p_box <- ggplot(df, aes(x = series, y = .data[[f]], fill = series)) +
    geom_boxplot() +
    labs(title = paste0(f, " - boxplot"),
         x = "Series", y = f, fill = "Series") +
    theme_minimal()

  # Significance between consecutive series (latter vs former), two-sided Wilcoxon.
  if (length(series_names) >= 2) {
    ymax <- max(df[[f]], na.rm = TRUE)
    ymin <- min(df[[f]], na.rm = TRUE)
    step <- (ymax - ymin) * 0.05
    sig <- data.frame()
    for (i in 2:length(series_names)) {
      former <- df[[f]][df$series == series_names[i - 1]]
      latter <- df[[f]][df$series == series_names[i]]
      p <- wilcox.test(latter, former)$p.value
      stars <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
      sig <- rbind(sig, data.frame(series = series_names[i],
                                   label = stars,
                                   y = ymax + step))
      cat(sprintf("  %s: %s vs %s, p = %.4g (%s)\n",
                  f, series_names[i - 1], series_names[i], p, stars))
    }
    p_box <- p_box +
      geom_text(data = sig,
                aes(x = series, y = y, label = label),
                inherit.aes = FALSE,
                position = position_nudge(x = -0.5))
  }

  ggsave(file.path("plot", paste0("dist_", f, "_violin.png")),
         p_violin, width = 7, height = 5)
  ggsave(file.path("plot", paste0("dist_", f, "_box.png")),
         p_box, width = 7, height = 5)

  cat("Saved:", file.path("plot", paste0("dist_", f, "_violin.png")), "\n")
  cat("Saved:", file.path("plot", paste0("dist_", f, "_box.png")), "\n")
}
