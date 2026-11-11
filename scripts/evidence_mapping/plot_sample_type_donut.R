#!/usr/bin/env Rscript

# Donut chart of sample-type occurrences in the included studies.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg) > 0L) {
  normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)
} else {
  normalizePath("scripts/evidence_mapping/plot_sample_type_donut.R", mustWork = FALSE)
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

sample_types <- read.csv(
  file.path(project_dir, "data", "evidence_mapping", "sample_types.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sample_types <- sample_types |>
  mutate(
    Sample_Type = factor(Sample_Type, levels = Sample_Type),
    proportion = n / sum(n),
    legend_label = paste0(Sample_Type, "  (", n, ")"),
    count_label = if_else(n >= 2, as.character(n), "")
  )

# Original qualitative palette.
sample_colours <- c(
  "#4F83B6",
  "#7FCDBB",
  "#FEC44F",
  "#F03B20",
  "#756BB1",
  "#31A354",
  "#DD1C77",
  "#8C6D31"
)

donut_plot <- ggplot(
  sample_types,
  aes(x = 2, y = n, fill = Sample_Type)
) +
  geom_col(
    width = 0.98,
    color = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = count_label),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 4,
    fontface = "bold",
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = 0.55,
    y = 0,
    label = "Sample-type\noccurrences",
    size = 4.5,
    fontface = "bold",
    lineheight = 0.95,
    color = "#222222"
  ) +
  coord_polar(theta = "y", start = -pi / 2) +
  xlim(0.55, 2.50) +
  scale_fill_manual(
    values = sample_colours,
    breaks = sample_types$Sample_Type,
    labels = sample_types$legend_label,
    name = NULL
  ) +
  guides(
    fill = guide_legend(
      ncol = 1,
      byrow = TRUE,
      override.aes = list(color = NA)
    )
  ) +
  theme_void(base_family = "sans") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9, color = "#222222"),
    legend.key.size = grid::unit(4.5, "mm"),
    legend.spacing.y = grid::unit(1.2, "mm"),
    plot.margin = margin(10, 10, 10, 10)
  )

output_dir <- file.path(project_dir, "figures", "evidence_mapping")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
base_file <- file.path(output_dir, "sample_type_donut")

ggsave(
  paste0(base_file, ".png"), donut_plot,
  width = 175, height = 105, units = "mm", dpi = 600, bg = "white"
)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    paste0(base_file, ".svg"), donut_plot,
    width = 175, height = 105, units = "mm",
    device = svglite::svglite, bg = "white"
  )
} else {
  message("Package 'svglite' is unavailable; skipped SVG output.")
}

message("Saved donut chart files to: ", normalizePath(output_dir, winslash = "/"))
