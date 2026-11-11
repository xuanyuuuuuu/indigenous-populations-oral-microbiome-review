# Heatmap of direction and reported statistical significance for the most
# frequently reported taxa across study comparisons.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg) > 0L) {
  normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)
} else {
  normalizePath("scripts/taxa_direction/plot_taxa_direction_top10.R", mustWork = FALSE)
}

project_dir <- dirname(dirname(dirname(script_file)))
project_library <- file.path(project_dir, "R_library")
if (dir.exists(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}

required_packages <- c("ggplot2", "dplyr", "tidyr", "svglite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install missing packages first: ", paste(missing_packages, collapse = ", "))
}

library(ggplot2)
library(dplyr)
library(tidyr)

data_file <- file.path(project_dir, "data", "taxa_direction", "taxa_direction_all.csv")
raw <- read.csv(
  data_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
# The first column is known and may retain a UTF-8 byte-order mark on Windows.
names(raw)[1] <- "Comparison_order"

required_columns <- c(
  "Comparison_order", "Study_ID", "Comparison", "Age_group",
  "Health_status", "Taxon", "Direction", "Significance", "Plot_label"
)
missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

# Number of taxa displayed in the heatmap.
n_top_taxa <- 10L

# Rank taxa by the number of unique comparisons in which they were reported.
# The number of contributing studies is used as the first tie-breaker.
taxa_frequency <- raw %>%
  group_by(Taxon) %>%
  summarise(
    N_comparisons = n_distinct(Comparison_order),
    N_studies = n_distinct(Study_ID),
    .groups = "drop"
  ) %>%
  arrange(desc(N_comparisons), desc(N_studies), Taxon)

top_taxa <- taxa_frequency %>%
  slice_head(n = n_top_taxa) %>%
  pull(Taxon)

# Keep comparisons that contain at least one of the selected taxa, mirroring
# the alpha-direction figure's removal of entirely uninformative columns.
selected_comparisons <- raw %>%
  filter(Taxon %in% top_taxa) %>%
  distinct(Comparison_order)

age_order <- c("Adults", "Children", "Mixed", "Mixed ages")
health_order <- c(
  "Generally healthy", "Dental caries", "Periodontal",
  "Systemic", "Systemic conditions", "Mixed", "Not assessed/reported"
)

comparisons <- raw %>%
  distinct(
    Comparison_order, Study_ID, Comparison, Age_group,
    Health_status, Plot_label
  ) %>%
  semi_join(selected_comparisons, by = "Comparison_order") %>%
  mutate(
    age_rank = match(Age_group, age_order),
    health_rank = match(Health_status, health_order),
    age_rank = if_else(is.na(age_rank), length(age_order) + 1L, age_rank),
    health_rank = if_else(
      is.na(health_rank), length(health_order) + 1L, health_rank
    )
  ) %>%
  arrange(age_rank, health_rank, Comparison_order) %>%
  select(-age_rank, -health_rank) %>%
  mutate(
    display_order = row_number(),
    y = n() - display_order + 1
  )

direction_levels <- c("Increased", "Decreased", "Unchanged", "Not mentioned")
direction_palette <- c(
  "Increased" = "#A2B8DA",
  "Decreased" = "#E9ABB0",
  "Unchanged" = "#8CCECC",
  "Not mentioned" = "#F3DAC9"
)

taxon_x <- setNames(seq_along(top_taxa) + 3, top_taxa)

heatmap_data <- expand_grid(
  Comparison_order = comparisons$Comparison_order,
  Taxon = top_taxa
) %>%
  left_join(
    raw %>%
      filter(Taxon %in% top_taxa) %>%
      select(Comparison_order, Taxon, Direction, Significance),
    by = c("Comparison_order", "Taxon")
  ) %>%
  left_join(
    comparisons %>% select(Comparison_order, y),
    by = "Comparison_order"
  ) %>%
  mutate(
    Direction_plot = case_when(
      is.na(Direction) ~ "Not mentioned",
      Direction == "Higher" ~ "Increased",
      Direction == "Lower" ~ "Decreased",
      Direction %in% c("Equal", "Equal (reported)",
                       "Equal (descriptive)", "Similar") ~ "Unchanged",
      TRUE ~ "Not mentioned"
    ),
    Direction_plot = factor(Direction_plot, levels = direction_levels),
    is_significant = !is.na(Significance) & Significance == "Significant",
    x = unname(taxon_x[Taxon])
  )

make_blocks <- function(values, y_values) {
  run_id <- cumsum(values != dplyr::lag(values, default = values[1]))
  tibble(value = values, y = y_values, run_id = run_id) %>%
    group_by(run_id, value) %>%
    summarise(
      ymin = min(y) - 0.5,
      ymax = max(y) + 0.5,
      y = mean(y),
      .groups = "drop"
    )
}

wrap_label <- function(x, width) {
  vapply(
    x,
    function(z) paste(strwrap(z, width = width), collapse = "\n"),
    character(1)
  )
}

age_light <- c(
  "Adults" = "#E9E8FA",
  "Children" = "#DFF3E5",
  "Mixed" = "#FAEECB",
  "Mixed ages" = "#FAEECB"
)
health_light <- c(
  "Generally healthy" = "#E3EEFA",
  "Dental caries" = "#FBE5D7",
  "Periodontal" = "#F3DDF0",
  "Systemic" = "#DDF3F4",
  "Systemic conditions" = "#DDF3F4",
  "Mixed" = "#F9EDC8",
  "Not assessed/reported" = "#E7EAEE"
)

age_blocks <- make_blocks(comparisons$Age_group, comparisons$y) %>%
  mutate(
    fill = unname(age_light[value]),
    fill = if_else(is.na(fill), "#EEEEEE", fill),
    label = if_else(value %in% c("Mixed", "Mixed ages"), "Mixed\nages", value)
  )

health_blocks <- make_blocks(comparisons$Health_status, comparisons$y) %>%
  mutate(
    fill = unname(health_light[value]),
    fill = if_else(is.na(fill), "#EEEEEE", fill),
    label = wrap_label(value, 13)
  )

age_boundaries <- head(age_blocks$ymin, -1)
health_boundaries <- head(health_blocks$ymin, -1)

comparison_labels <- setNames(
  wrap_label(comparisons$Plot_label, 25),
  comparisons$y
)

x_breaks <- c(0.85, 2.35, unname(taxon_x))
x_labels <- parse(text = c(
  "'Age group'",
  "'Health status'",
  paste0("italic('", top_taxa, "')")
))

heatmap_plot <- ggplot() +
  geom_rect(
    data = age_blocks,
    aes(xmin = 0.18, xmax = 1.52, ymin = ymin, ymax = ymax, fill = I(fill)),
    color = "white", linewidth = 0.6
  ) +
  geom_text(
    data = age_blocks,
    aes(x = 0.85, y = y, label = label),
    size = 2.45, lineheight = 0.88, fontface = "bold", color = "#222222"
  ) +
  geom_rect(
    data = health_blocks,
    aes(xmin = 1.58, xmax = 3.12, ymin = ymin, ymax = ymax, fill = I(fill)),
    color = "white", linewidth = 0.6
  ) +
  geom_text(
    data = health_blocks,
    aes(x = 2.35, y = y, label = label),
    size = 2.2, lineheight = 0.86, fontface = "bold", color = "#222222"
  ) +
  geom_tile(
    data = heatmap_data,
    aes(x = x, y = y, fill = Direction_plot),
    width = 0.98, height = 0.98, color = "white", linewidth = 0.32
  ) +
  geom_text(
    data = heatmap_data %>% filter(is_significant),
    aes(x = x + 0.38, y = y + 0.36, label = "*"),
    hjust = 1, vjust = 1, size = 4.3, fontface = "bold", color = "#1F2933"
  ) +
  geom_hline(
    yintercept = health_boundaries,
    color = "#FFFFFF", linewidth = 1.0
  ) +
  geom_hline(
    yintercept = age_boundaries,
    color = "#FFFFFF", linewidth = 2.0
  ) +
  scale_fill_manual(
    values = direction_palette,
    breaks = direction_levels,
    drop = FALSE,
    name = "Direction"
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    limits = c(0.15, max(unname(taxon_x)) + 0.52),
    position = "top",
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = comparisons$y,
    labels = comparison_labels[as.character(comparisons$y)],
    limits = c(0.48, nrow(comparisons) + 0.52),
    expand = c(0, 0)
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      keywidth = grid::unit(7, "mm"),
      keyheight = grid::unit(5.5, "mm")
    )
  ) +
  labs(
    title = "Direction of reported taxa across comparisons",
    x = NULL,
    y = NULL,
    caption = "* Statistically significant as reported in the source study"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      size = 8.6, face = "bold", color = "#202020",
      margin = margin(b = 8)
    ),
    axis.text.y = element_text(
      size = 7.8, color = "#202020", lineheight = 0.92,
      margin = margin(r = 7)
    ),
    plot.title = element_text(
      size = 13, face = "bold", color = "#202020", margin = margin(b = 8)
    ),
    plot.caption = element_text(
      size = 8, color = "#333333", hjust = 0, margin = margin(t = 8)
    ),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8.5),
    legend.spacing.y = grid::unit(2, "mm"),
    plot.margin = margin(12, 14, 10, 10)
  )

output_dir <- file.path(project_dir, "figures", "taxa_direction")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
base_file <- file.path(output_dir, "taxa_direction_top10")
figure_height <- max(8.2, 2.6 + 0.36 * nrow(comparisons))
figure_width <- max(10.8, 5.8 + 0.92 * length(top_taxa))

ggsave(
  paste0(base_file, ".svg"), heatmap_plot,
  width = figure_width, height = figure_height, units = "in",
  device = svglite::svglite, bg = "white"
)
ggsave(
  paste0(base_file, ".png"), heatmap_plot,
  width = figure_width, height = figure_height, units = "in",
  dpi = 350, bg = "white"
)

message("Top taxa: ", paste(top_taxa, collapse = ", "))
message("Comparisons shown: ", nrow(comparisons))
message("Saved taxa-direction heatmap to: ", normalizePath(output_dir, winslash = "/"))
