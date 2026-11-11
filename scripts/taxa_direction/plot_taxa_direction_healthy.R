# Direction synthesis for taxa reported in generally healthy adults and children.
# Taxa are selected separately within each age group.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg) > 0L) {
  normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)
} else {
  normalizePath("scripts/taxa_direction/plot_taxa_direction_healthy.R", mustWork = FALSE)
}

project_dir <- dirname(dirname(dirname(script_file)))
project_library <- file.path(project_dir, "R_library")
if (dir.exists(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}

required_packages <- c("dplyr", "ggplot2", "tidyr", "patchwork", "svglite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install missing packages first: ", paste(missing_packages, collapse = ", "))
}

library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

data_file <- file.path(project_dir, "data", "taxa_direction", "taxa_direction_healthy.csv")
raw <- read.csv(
  data_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
names(raw)[1] <- "Comparison_order"

required_columns <- c(
  "Comparison_order", "Study_ID", "Comparison", "Age_group",
  "Health_status", "Taxon", "Direction", "Plot_label", "Significance"
)
missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

dat <- raw %>%
  filter(
    Health_status == "Generally healthy",
    Age_group %in% c("Adults", "Children")
  )

direction_levels <- c("Increased", "Decreased", "Unchanged", "Not mentioned")
direction_palette <- c(
  "Increased" = "#A2B8DA",
  "Decreased" = "#E9ABB0",
  "Unchanged" = "#8CCECC",
  "Not mentioned" = "#F3DAC9"
)

rank_taxa <- function(age_data) {
  age_data %>%
    group_by(Taxon) %>%
    summarise(
      N_studies = n_distinct(Study_ID),
      N_comparisons = n_distinct(Comparison_order),
      N_significant = sum(Significance == "Significant", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(
      desc(N_studies), desc(N_comparisons), desc(N_significant), Taxon
    )
}

adult_ranking <- rank_taxa(filter(dat, Age_group == "Adults"))
child_ranking <- rank_taxa(filter(dat, Age_group == "Children"))

adult_taxa <- adult_ranking %>% slice_head(n = 10L) %>% pull(Taxon)
# The children subset contains only eight taxa, so retain all of them.
child_taxa <- child_ranking %>% pull(Taxon)

normalise_direction <- function(x) {
  case_when(
    is.na(x) ~ "Not mentioned",
    x == "Higher" ~ "Increased",
    x == "Lower" ~ "Decreased",
    x %in% c("Equal", "Equal (reported)", "Equal (descriptive)", "Similar") ~
      "Unchanged",
    TRUE ~ "Not mentioned"
  )
}

wrap_plot_label <- function(x, width = 26) {
  vapply(
    x,
    function(z) {
      z <- gsub("\\\\n", "\n", z)
      parts <- strsplit(z, "\n", fixed = TRUE)[[1]]
      paste(
        vapply(parts, function(p) paste(strwrap(p, width), collapse = "\n"), ""),
        collapse = "\n"
      )
    },
    character(1)
  )
}

prepare_panel <- function(age, selected_taxa) {
  age_data <- dat %>% filter(Age_group == age)

  comparisons <- age_data %>%
    distinct(
      Comparison_order, Study_ID, Comparison, Plot_label,
      Age_group, Health_status
    ) %>%
    arrange(Comparison_order) %>%
    mutate(
      y = rev(seq_len(n())),
      display_label = wrap_plot_label(Plot_label)
    )

  taxa_table <- tibble(
    Taxon = selected_taxa,
    x = seq_along(selected_taxa)
  )

  matrix_data <- expand_grid(
    Comparison_order = comparisons$Comparison_order,
    Taxon = selected_taxa
  ) %>%
    left_join(
      age_data %>%
        filter(Taxon %in% selected_taxa) %>%
        select(Comparison_order, Taxon, Direction, Significance) %>%
        distinct(Comparison_order, Taxon, .keep_all = TRUE),
      by = c("Comparison_order", "Taxon")
    ) %>%
    left_join(
      comparisons %>% select(Comparison_order, y),
      by = "Comparison_order"
    ) %>%
    left_join(taxa_table, by = "Taxon") %>%
    mutate(
      Direction_plot = factor(
        normalise_direction(Direction),
        levels = direction_levels
      ),
      is_significant = !is.na(Significance) & Significance == "Significant"
    )

  list(
    comparisons = comparisons,
    taxa = taxa_table,
    matrix = matrix_data
  )
}

adult_panel <- prepare_panel("Adults", adult_taxa)
child_panel <- prepare_panel("Children", child_taxa)

make_panel <- function(panel_data, panel_letter, age_label, show_legend = TRUE) {
  comparisons <- panel_data$comparisons
  taxa_table <- panel_data$taxa
  matrix_data <- panel_data$matrix

  taxa_labels <- setNames(
    lapply(taxa_table$Taxon, function(z) bquote(italic(.(z)))),
    taxa_table$x
  )
  y_labels <- setNames(comparisons$display_label, comparisons$y)
  legend_data <- tibble(
    x = rep(1, length(direction_levels)),
    y = rep(1, length(direction_levels)),
    Direction_plot = factor(direction_levels, levels = direction_levels)
  )

  ggplot(matrix_data, aes(x = x, y = y)) +
    geom_point(
      data = legend_data,
      aes(x = x, y = y, fill = Direction_plot),
      inherit.aes = FALSE, shape = 22, size = 7,
      alpha = 0, color = "transparent", show.legend = TRUE
    ) +
    geom_tile(
      aes(fill = Direction_plot),
      width = 0.98, height = 0.98,
      color = "white", linewidth = 0.36,
      show.legend = FALSE
    ) +
    geom_text(
      data = matrix_data %>% filter(is_significant),
      aes(x = x + 0.35, y = y + 0.34, label = "*"),
      hjust = 1, vjust = 1,
      size = 4.0, fontface = "bold", color = "#202830"
    ) +
    scale_fill_manual(
      values = direction_palette,
      limits = direction_levels,
      breaks = direction_levels,
      drop = FALSE,
      name = "Direction"
    ) +
    scale_x_continuous(
      breaks = taxa_table$x,
      labels = taxa_labels,
      limits = c(0.5, nrow(taxa_table) + 0.5),
      position = "top",
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = comparisons$y,
      labels = y_labels[as.character(comparisons$y)],
      limits = c(0.5, nrow(comparisons) + 0.5),
      expand = c(0, 0)
    ) +
    coord_fixed(ratio = 1, clip = "off") +
    labs(
      title = paste0(panel_letter, ". ", age_label, " (generally healthy)"),
      x = NULL,
      y = NULL
    ) +
    guides(
      fill = guide_legend(
        title.position = "top",
        keywidth = grid::unit(7, "mm"),
        keyheight = grid::unit(5.8, "mm"),
        override.aes = list(
          shape = 22, size = 7,
          fill = unname(direction_palette), color = "transparent", alpha = 1
        )
      )
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(
        color = "#C4C9CE", fill = NA, linewidth = 0.55
      ),
      axis.ticks = element_blank(),
      axis.text.x = element_text(
        angle = 32, hjust = 0, vjust = 0.5,
        size = 8.7, color = "#202020", margin = margin(b = 8)
      ),
      axis.text.y = element_text(
        size = 8.3, color = "#202020", lineheight = 0.92,
        margin = margin(r = 7)
      ),
      plot.title = element_text(
        size = 12.5, face = "bold", color = "#202020",
        margin = margin(b = 8)
      ),
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(size = 9.3, face = "bold"),
      legend.text = element_text(size = 8.8),
      plot.margin = margin(8, 8, 8, 8)
    )
}

p_adults <- make_panel(adult_panel, "A", "Adults", show_legend = TRUE)
p_children <- make_panel(child_panel, "B", "Children", show_legend = TRUE)

# The design keeps the same approximate cell width in both panels: the child
# panel uses eight columns and leaves the equivalent of two columns blank.
design <- "
AAAAAAAAAALL
BBBBBBBBSSLL
"

figure <- p_adults + p_children + plot_spacer() + guide_area() +
  plot_layout(
    design = design,
    heights = c(11, 4),
    guides = "collect"
  ) +
  plot_annotation(
    title = "Direction of reported taxa in generally healthy populations",
    caption = "* Statistically significant as reported in the source study",
    theme = theme(
      plot.title = element_text(
        size = 15, face = "bold", color = "#202020", margin = margin(b = 8)
      ),
      plot.caption = element_text(
        size = 8.5, color = "#333333", hjust = 0, margin = margin(t = 8)
      ),
      plot.margin = margin(10, 12, 10, 10)
    )
  )

output_dir <- file.path(project_dir, "figures", "taxa_direction")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
base_file <- file.path(
  output_dir,
  "taxa_direction_healthy"
)

ggsave(
  paste0(base_file, ".svg"), figure,
  width = 13.5, height = 12.4, units = "in",
  device = svglite::svglite, bg = "white"
)
ggsave(
  paste0(base_file, ".png"), figure,
  width = 13.5, height = 12.4, units = "in",
  dpi = 350, bg = "white"
)

message("Adult taxa: ", paste(adult_taxa, collapse = ", "))
message("Child taxa: ", paste(child_taxa, collapse = ", "))
message("Adult comparisons: ", nrow(adult_panel$comparisons))
message("Child comparisons: ", nrow(child_panel$comparisons))
message("Saved healthy-population taxa heatmap to: ", normalizePath(output_dir, winslash = "/"))
