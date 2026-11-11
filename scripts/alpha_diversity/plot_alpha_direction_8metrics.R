# Direction matrix for reported alpha-diversity indices.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg) > 0L) {
  normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)
} else {
  normalizePath("scripts/alpha_diversity/plot_alpha_direction_8metrics.R", mustWork = FALSE)
}

project_dir <- dirname(dirname(dirname(script_file)))
project_library <- file.path(project_dir, "R_library")
if (dir.exists(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}

required_packages <- c("dplyr", "ggplot2", "patchwork", "svglite")
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

library(dplyr)
library(ggplot2)
library(patchwork)
library(svglite)

data_file <- file.path(
  project_dir,
  "data", "alpha_diversity",
  "alpha_direction_8metrics.csv"
)
if (!file.exists(data_file)) {
  stop("Input file does not exist: ", data_file, call. = FALSE)
}

raw <- read.csv(
  data_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
# Remove a UTF-8 byte-order mark retained by base R on some Windows locales.
names(raw)[1] <- "Comparison_order"

# Display all eight indices in the established figure order.
metrics <- c(
  "Shannon", "Simpson", "Chao1", "Observed richness",
  "Phylogenetic diversity (PD)", "ACE", "BWPD (0.25)",
  "Pielou evenness"
)

dat <- raw %>%
  filter(Metric %in% metrics) %>%
  group_by(Comparison_order) %>%
  filter(any(Direction != "Not reported")) %>%
  ungroup() %>%
  mutate(Metric = factor(Metric, levels = rev(metrics)))

cols <- dat %>%
  distinct(Comparison_order, Comparison, Plot_label, Age_group, Health_status) %>%
  mutate(
    age_rank = match(Age_group, c("Adults", "Children", "Mixed", "Mixed ages")),
    health_rank = match(
      Health_status,
      c(
        "Generally healthy", "Dental caries", "Periodontal conditions",
        "Systemic conditions", "Mixed", "Not assessed/reported"
      )
    )
  ) %>%
  arrange(age_rank, health_rank, Comparison_order) %>%
  select(-age_rank, -health_rank) %>%
  mutate(x = row_number())

dat <- dat %>%
  left_join(cols %>% select(Comparison_order, x), by = "Comparison_order") %>%
  mutate(
    Direction = factor(
      Direction,
      levels = c("Higher", "Lower", "Equal", "Not reported")
    )
  )

# Numeric x positions prevent repeated population labels from merging.
age_pal <- c(
  "Adults" = "#8871B5",
  "Children" = "#50A45C",
  "Mixed" = "#E0B83F",
  "Mixed ages" = "#E0B83F"
)
age_light <- c(
  "Adults" = "#E9E8FA",
  "Children" = "#DFF3E5",
  "Mixed" = "#FAEECB",
  "Mixed ages" = "#FAEECB"
)
health_pal <- c(
  "Generally healthy" = "#4F83B6",
  "Dental caries" = "#F07858",
  "Periodontal conditions" = "#D76BB1",
  "Systemic conditions" = "#37A9B2",
  "Mixed" = "#E7B840",
  "Not assessed/reported" = "#9EA7B2"
)
health_light <- c(
  "Generally healthy" = "#E3EEFA",
  "Dental caries" = "#FBE5D7",
  "Periodontal conditions" = "#F3DDF0",
  "Systemic conditions" = "#DDF3F4",
  "Mixed" = "#F9EDC8",
  "Not assessed/reported" = "#E7EAEE"
)

make_blocks <- function(values) {
  block_id <- cumsum(values != dplyr::lag(values, default = values[1]))
  data.frame(x = cols$x, value = values, block_id = block_id) %>%
    group_by(block_id, value) %>%
    summarise(
      xmin = min(x) - 0.5,
      xmax = max(x) + 0.5,
      xmid = mean(range(x)),
      n = n(),
      .groups = "drop"
    )
}

age_blocks <- make_blocks(cols$Age_group) %>%
  mutate(
    label = if_else(value %in% c("Mixed", "Mixed ages"), "Mixed\nages", value),
    text_size = ifelse(n == 1, 3.1, 3.6)
  )

health_blocks <- make_blocks(cols$Health_status) %>%
  mutate(
    short_name = recode(
      value,
      "Periodontal conditions" = "Periodontal",
      "Systemic conditions" = "Systemic",
      "Not assessed/reported" = "Not\nassessed"
    ),
    label = short_name,
    text_size = case_when(n == 1 ~ 2.45, n == 2 ~ 2.8, TRUE ~ 3.25)
  )

x_limits <- c(0.5, nrow(cols) + 0.5)
shared_margin <- margin(0, 18, 0, 8)

make_header <- function(blocks, palette, header_label) {
  ggplot(blocks) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = value),
      color = "white",
      linewidth = 0.7
    ) +
    geom_text(
      aes(x = xmid, y = 0.5, label = label, size = text_size),
      lineheight = 0.92,
      fontface = "bold"
    ) +
    scale_fill_manual(values = palette, guide = "none") +
    scale_size_identity() +
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    scale_y_continuous(
      limits = c(0, 1), breaks = 0.5, labels = header_label,
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    theme_void(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 9.5, face = "bold", color = "#222222"),
      axis.text.y.left = element_text(margin = margin(r = 8)),
      plot.margin = shared_margin
    )
}

age_header <- make_header(age_blocks, age_light, "Age group")
health_header <- make_header(health_blocks, health_light, "Health status")

age_boundaries <- head(age_blocks$xmax, -1)
health_boundaries <- head(health_blocks$xmax, -1)
comparison_labels <- vapply(
  cols$Plot_label,
  function(z) {
    z <- gsub("\\\\n", "\n", z)
    paste(strwrap(z, width = 22), collapse = "\n")
  },
  character(1)
)

p <- ggplot(dat, aes(x = x, y = Metric)) +
  geom_vline(
    xintercept = seq(1.5, nrow(cols) - 0.5, 1),
    color = "#E0E4E8",
    linewidth = 0.35
  ) +
  geom_hline(
    yintercept = seq(1.5, length(metrics) - 0.5, 1),
    color = "#E0E4E8",
    linewidth = 0.35
  ) +
  geom_vline(
    xintercept = health_boundaries,
    color = "#BCC4CB",
    linewidth = 0.7
  ) +
  geom_vline(
    xintercept = age_boundaries,
    color = "#A6AFB8",
    linewidth = 1.05
  ) +
  geom_point(
    aes(fill = Direction, size = Direction),
    shape = 21,
    stroke = 0.25,
    color = "white"
  ) +
  scale_fill_manual(
    values = c(
      "Higher" = "#4C82C3",
      "Lower" = "#F1714C",
      "Equal" = "#C4CBD1",
      "Not reported" = "#B8BEC3"
    ),
    breaks = c("Higher", "Lower", "Equal", "Not reported"),
    drop = FALSE,
    name = "Direction (relative\nto comparator)"
  ) +
  scale_size_manual(
    values = c(
      "Higher" = 12.8,
      "Lower" = 12.8,
      "Equal" = 9.0,
      "Not reported" = 2.8
    ),
    breaks = c("Higher", "Lower", "Equal", "Not reported"),
    drop = FALSE,
    name = "Direction (relative\nto comparator)"
  ) +
  scale_x_continuous(
    breaks = cols$x,
    labels = comparison_labels,
    limits = x_limits,
    expand = c(0, 0)
  ) +
  coord_fixed(ratio = 1, clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "#AEB6BD", fill = NA, linewidth = 0.65),
    axis.text.x = element_text(
      angle = 65,
      hjust = 1,
      vjust = 1,
      size = 7.0,
      color = "#202020",
      lineheight = 0.93
    ),
    axis.text.y = element_text(size = 10.7, color = "#202020"),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 9.5, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.height = grid::unit(6, "mm"),
    plot.margin = margin(0, 18, 12, 8)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        size = c(8.5, 8.5, 6.5, 2.5),
        color = "white"
      )
    ),
    size = "none"
  )

fig <- age_header / health_header / p +
  plot_layout(heights = c(0.72, 0.88, 4.8))

output_dir <- file.path(project_dir, "figures", "alpha_diversity")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
base_file <- file.path(output_dir, "alpha_direction_8metrics")

ggsave(
  paste0(base_file, ".png"), fig,
  width=20.5, height=9.4, units="in", dpi=400, bg="white"
)
ggsave(
  paste0(base_file, ".svg"), fig,
  width=20.5, height=9.4, units="in",
  device=svglite::svglite, bg="white"
)

message("Saved alpha-direction figure to: ", normalizePath(output_dir, winslash = "/"))
