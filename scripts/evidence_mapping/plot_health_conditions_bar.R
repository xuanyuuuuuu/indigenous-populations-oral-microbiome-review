#!/usr/bin/env Rscript

# Horizontal bar chart of health conditions represented in the review.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg) > 0L) {
  normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)
} else {
  normalizePath("scripts/evidence_mapping/plot_health_conditions_bar.R", mustWork = FALSE)
}

project_dir <- dirname(dirname(dirname(script_file)))
project_library <- file.path(project_dir, "R_library")
if (dir.exists(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}

required_packages <- c("ggplot2", "dplyr", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install missing packages first: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

health_conditions <- read.csv(
  file.path(project_dir, "data", "evidence_mapping", "health_conditions.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Preserve the supplied order for tied values and display the largest at the top.
health_conditions <- health_conditions |>
  mutate(
    display_order = row_number(),
    Health_Condition = factor(
      Health_Condition,
      levels = rev(Health_Condition)
    )
  )

donut_blue <- "#4F83B6"

bar_plot <- ggplot(
  health_conditions,
  aes(x = n, y = Health_Condition)
) +
  geom_col(
    width = 0.68,
    fill = donut_blue
  ) +
  geom_text(
    aes(label = n),
    hjust = -0.45,
    size = 3.8,
    fontface = "bold",
    color = "#222222"
  ) +
  scale_x_continuous(
    breaks = seq(0, 10, by = 2),
    limits = c(0, 11),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Number of studies",
    y = NULL
  ) +
  theme_minimal(base_family = "sans", base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#D9E1E7", linewidth = 0.45),
    axis.text.x = element_text(size = 9, color = "#333333"),
    axis.text.y = element_text(size = 9.5, color = "#222222"),
    axis.title.x = element_text(size = 10, color = "#222222"),
    axis.ticks = element_blank(),
    plot.margin = margin(8, 18, 8, 8)
  )

output_dir <- file.path(project_dir, "figures", "evidence_mapping")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
base_file <- file.path(output_dir, "health_conditions_bar")

ggsave(
  paste0(base_file, ".png"), bar_plot,
  width = 180, height = 88, units = "mm", dpi = 600, bg = "white"
)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    paste0(base_file, ".svg"), bar_plot,
    width = 180, height = 88, units = "mm",
    device = svglite::svglite, bg = "white"
  )
} else {
  message("Package 'svglite' is unavailable; skipped SVG output.")
}

message("Saved bar chart files to: ", normalizePath(output_dir, winslash = "/"))
