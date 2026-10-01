#!/usr/bin/env Rscript
# Plot distributions of `avg`, `SAscore`, `QED_w` across multiple files.
#
# To add a new file, add one entry to `files` below:
#   name = "path/to/properties.csv"

library(ggplot2)

# ---- Configuration (extend here for more files) ----
files <- list(
  test_new = "results0/test_new_properties.csv",
  test_new1 = "results0/test_new1_properties.csv",
  test_new2 = "results0/test_new2_properties.csv",
  test_new3 = "results0/test_new3_properties.csv"
)
# factors <- c("avg", "SAscore", "QED_w")
factors <- c("SAscore", "My_QED")

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
  p_hist <- ggplot(df, aes(x = .data[[f]], fill = series)) +
    geom_histogram(alpha = 0.4, position = "identity", bins = 50) +
    labs(title = paste0(f, " - distribution"),
         x = f, y = "Count", fill = "Series") +
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

  ggsave(file.path("plot", paste0("dist_", f, "_hist.png")),
         p_hist, width = 7, height = 5)
  ggsave(file.path("plot", paste0("dist_", f, "_box.png")),
         p_box, width = 7, height = 5)

  cat("Saved:", file.path("plot", paste0("dist_", f, "_hist.png")), "\n")
  cat("Saved:", file.path("plot", paste0("dist_", f, "_box.png")), "\n")
}
