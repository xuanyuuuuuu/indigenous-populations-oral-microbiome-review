#!/usr/bin/env Rscript

# Draw the geographic coverage of studies by continent.
# Hue identifies the continent; dark/light shading identifies review inclusion.
# The original count-gradient map remains in plot_review_country_distribution.R.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg) > 0L) {
  normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)
} else {
  normalizePath(
    "scripts/evidence_mapping/plot_review_country_distribution.R",
    mustWork = FALSE
  )
}

project_dir <- dirname(dirname(dirname(script_file)))
project_library <- file.path(project_dir, "R_library")
if (dir.exists(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}

required_packages <- c(
  "ggplot2", "sf", "rnaturalearth", "dplyr", "units"
)
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
  library(sf)
  library(rnaturalearth)
  library(dplyr)
  library(units)
})

review_data <- read.csv(
  file.path(project_dir, "data", "evidence_mapping", "country_studies.csv"),
  stringsAsFactors = FALSE
)

# The saturated versions of these hues mark countries represented in the review.
continent_colours <- c(
  "Africa" = "#D89B00",
  "Asia" = "#2166AC",
  "Europe" = "#7B3294",
  "North America" = "#1B7837",
  "South America" = "#D6604D",
  "Oceania" = "#C51B7D",
  "Antarctica" = "#737373"
)

sf::sf_use_s2(FALSE)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
  st_make_valid() |>
  mutate(
    map_iso3 = if_else(iso_a3 == "-99", adm0_a3, iso_a3),
    Continent = case_when(
      continent == "Africa" ~ "Africa",
      continent == "Asia" ~ "Asia",
      continent == "Europe" ~ "Europe",
      continent == "North America" ~ "North America",
      continent == "South America" ~ "South America",
      continent == "Oceania" ~ "Oceania",
      continent == "Antarctica" ~ "Antarctica",
      TRUE ~ "Oceania"
    ),
    Continent = factor(Continent, levels = names(continent_colours))
  ) |>
  left_join(review_data, by = c("map_iso3" = "ISO3")) |>
  mutate(
    Review_status = if_else(
      !is.na(N_studies) & N_studies > 0,
      "Included in review",
      "Not represented"
    ),
    Review_status = factor(
      Review_status,
      levels = c("Included in review", "Not represented")
    )
  ) |>
  st_transform("+proj=eqearth")

unmatched_iso3 <- setdiff(review_data$ISO3, unique(world$map_iso3))
if (length(unmatched_iso3) > 0L) {
  warning(
    "ISO3 code(s) not found in the map data: ",
    paste(unmatched_iso3, collapse = ", "),
    call. = FALSE
  )
}

graticules <- st_graticule(
  lon = seq(-180, 180, by = 20),
  lat = seq(-80, 80, by = 20)
) |>
  st_transform("+proj=eqearth")

map_plot <- ggplot() +
  geom_sf(data = graticules, color = "grey88", linewidth = 0.2) +
  geom_sf(
    data = world,
    aes(fill = Continent, alpha = Review_status),
    color = "white",
    linewidth = 0.22
  ) +
  scale_fill_manual(
    values = continent_colours,
    breaks = setdiff(names(continent_colours), "Antarctica"),
    drop = FALSE,
    name = "Continent"
  ) +
  scale_alpha_manual(
    values = c("Included in review" = 1, "Not represented" = 0.20),
    drop = FALSE,
    name = "Review coverage"
  ) +
  coord_sf(datum = NA) +
  labs(title = "Geographic distribution of included studies") +
  theme_void(base_family = "sans") +
  guides(
    fill = guide_legend(
      order = 1,
      nrow = 2,
      title.position = "top",
      override.aes = list(alpha = 1, color = NA)
    ),
    alpha = guide_legend(
      order = 2,
      title.position = "top",
      override.aes = list(fill = "#444444", color = NA)
    )
  ) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.width = unit(8, "mm"),
    legend.key.height = unit(4, "mm"),
    legend.spacing.y = unit(2, "mm"),
    plot.margin = margin(8, 8, 5, 8)
  )

output_dir <- file.path(project_dir, "figures", "evidence_mapping")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
base_file <- file.path(output_dir, "country_distribution")

ggsave(
  paste0(base_file, ".png"), map_plot,
  width = 180, height = 125, units = "mm", dpi = 600, bg = "white"
)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    paste0(base_file, ".svg"), map_plot,
    width = 180, height = 125, units = "mm",
    device = svglite::svglite, bg = "white"
  )
} else {
  message("Package 'svglite' is unavailable; skipped SVG output.")
}

message("Saved continent-colour map files to: ", normalizePath(output_dir, winslash = "/"))
