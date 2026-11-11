# CARE-related practices across studies
# Input: data/care_principles/care_indicators.csv
# Outputs: publication-ready SVG and 600 dpi PNG

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)
} else {
  normalizePath("scripts/care_principles/plot_care_reporting_heatmap.R", mustWork = FALSE)
}

project_dir <- dirname(dirname(dirname(script_file)))
project_library <- file.path(project_dir, "R_library")
if (dir.exists(project_library)) .libPaths(c(project_library, .libPaths()))

required <- c("dplyr", "tidyr", "ggplot2", "svglite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

library(dplyr)
library(tidyr)
library(ggplot2)

if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}

input_file <- file.path(project_dir, "data", "care_principles", "care_indicators.csv")
figure_dir <- file.path(project_dir, "figures", "care_principles")
result_dir <- file.path(project_dir, "results", "care_principles")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

codes <- c("C1", "C2", "C3", "A1", "A2", "A3", "R1", "R2", "R3", "E1", "E2", "E3")

# This mapping follows the order and wording supplied with the current CARE dataset.
definitions <- data.frame(
  code = codes,
  group = c(rep("Collective Benefit (C)", 3), rep("Authority to Control (A)", 3),
            rep("Responsibility (R)", 3), rep("Ethics (E)", 3)),
  description = c(
    "For inclusive development and innovation",
    "For improved governance and citizen engagement",
    "For equitable outcomes",
    "Recognizing rights and interests",
    "Data for governance",
    "Governance of data",
    "For positive relationships",
    "For expanding capability and capacity",
    "For Indigenous languages and worldviews",
    "For minimizing harm and maximizing benefit",
    "For justice",
    "For future use"
  ),
  x = seq_along(codes),
  stringsAsFactors = FALSE
)

group_bands <- data.frame(
  group = c("Collective Benefit (C)", "Authority to Control (A)", "Responsibility (R)", "Ethics (E)"),
  xmin = c(0.52, 3.52, 6.52, 9.52),
  xmax = c(3.48, 6.48, 9.48, 12.48),
  fill = c("#CFE3F4", "#D8EEE1", "#F3D7D0", "#F8E7B7"),
  stringsAsFactors = FALSE
)

wrap_text <- function(x, width) paste(strwrap(x, width = width), collapse = "\n")

clean_author <- function(x) {
  match <- regexpr("^.*?et\\s*\\.?\\s*al\\.?", x, ignore.case = TRUE, perl = TRUE)
  author <- if (match[1] > 0) regmatches(x, match) else trimws(strsplit(x, ";", fixed = TRUE)[[1]][1])
  author <- gsub("et\\s*\\.\\s*al", "et al", author, ignore.case = TRUE)
  author <- gsub("et\\s+al", "et al", author, ignore.case = TRUE)
  author <- sub("[.;]+$", "", author)
  if (grepl("et al$", author, ignore.case = TRUE)) author <- paste0(author, ".")
  author
}

extract_year <- function(x) {
  hits <- regmatches(x, gregexpr("(19|20)[0-9]{2}", x, perl = TRUE))[[1]]
  if (!length(hits)) return(NA_integer_)
  as.integer(tail(hits, 1))
}

dat <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
needed <- c("Study(Author; Country; Year)", codes)
if (!all(needed %in% names(dat))) {
  stop("care_indicators.csv is missing required columns: ", paste(setdiff(needed, names(dat)), collapse = ", "))
}

dat <- dat %>%
  mutate(
    author = vapply(`Study(Author; Country; Year)`, clean_author, character(1)),
    year = vapply(`Study(Author; Country; Year)`, extract_year, integer(1)),
    base_label = sprintf("%s (%d)", author, year)
  )

# Add a/b only when author-year labels are duplicated.
dup_n <- ave(seq_len(nrow(dat)), dat$base_label, FUN = length)
dup_i <- ave(seq_len(nrow(dat)), dat$base_label, FUN = seq_along)
dat$study_label <- ifelse(
  dup_n > 1,
  paste0(sub("\\)$", "", dat$base_label), letters[dup_i], ")"),
  dat$base_label
)

invalid <- setdiff(unique(unlist(dat[codes], use.names = FALSE)), c("Y", "NR"))
invalid <- invalid[!is.na(invalid)]
if (length(invalid)) stop("Unexpected CARE coding: ", paste(invalid, collapse = ", "))

n_studies <- nrow(dat)
dat$y <- rev(seq_len(n_studies))

heat <- dat %>%
  select(study_label, y, all_of(codes)) %>%
  pivot_longer(all_of(codes), names_to = "code", values_to = "status") %>%
  mutate(
    x = match(code, codes),
    status = factor(status, levels = c("Y", "NR"), labels = c("Reported", "Not reported"))
  )

summary_dat <- heat %>%
  group_by(code, x) %>%
  summarise(reported = sum(status == "Reported"), .groups = "drop") %>%
  mutate(
    percent = 100 * reported / n_studies,
    label = sprintf("%d\n(%.1f)", reported, percent)
  )

study_summary <- heat %>%
  group_by(study_label, y) %>%
  summarise(indicators_reported = sum(status == "Reported"), .groups = "drop") %>%
  mutate(label = sprintf("%d/12", indicators_reported))

reported_blue <- "#4F83B6"
not_reported_grey <- "#E3E8EC"
summary_fill <- "#EDF2F5"
header_ymin <- n_studies + 0.75
header_ymax <- n_studies + 2.55
group_ymin <- n_studies + 2.75
group_ymax <- n_studies + 4.15
summary_y <- -1.15

p <- ggplot() +
  # Study-by-subprinciple matrix
  geom_tile(
    data = heat,
    aes(x = x, y = y, fill = status),
    width = 0.98, height = 0.92, color = "white", linewidth = 0.35
  ) +
  # Per-study total across the 12 CARE indicators
  geom_tile(
    data = study_summary, aes(x = 13, y = y),
    width = 0.98, height = 0.92, fill = "#F1F4F6", color = "white", linewidth = 0.35
  ) +
  geom_text(
    data = study_summary, aes(x = 13, y = y, label = label),
    family = "Arial", fontface = "bold", size = 3.25, color = "#1A1A1A"
  ) +
  # Study labels
  geom_text(
    data = dat, aes(x = 0.38, y = y, label = study_label),
    hjust = 1, family = "Arial", size = 3.25, color = "#1A1A1A"
  ) +
  # Subprinciple header cells
  geom_rect(
    data = definitions,
    aes(xmin = x - 0.48, xmax = x + 0.48, ymin = header_ymin, ymax = header_ymax),
    fill = "#F4F7F8", color = "white", linewidth = 0.6
  ) +
  geom_text(
    data = definitions, aes(x = x, y = header_ymax - 0.72, label = code),
    family = "Arial", fontface = "bold", size = 4.6
  ) +
  # CARE group bands
  geom_rect(
    data = group_bands,
    aes(xmin = xmin, xmax = xmax, ymin = group_ymin, ymax = group_ymax),
    fill = group_bands$fill, color = "white", linewidth = 0.7
  ) +
  geom_text(
    data = group_bands,
    aes(x = (xmin + xmax) / 2, y = (group_ymin + group_ymax) / 2, label = group),
    family = "Arial", fontface = "bold", size = 4.7
  ) +
  annotate("rect", xmin = 12.52, xmax = 13.48, ymin = header_ymin, ymax = group_ymax,
           fill = "#E5EAED", color = "white", linewidth = 0.7) +
  annotate("text", x = 13, y = (header_ymin + group_ymax) / 2,
           label = "Indicators\nReported\n(n/12)", family = "Arial",
           fontface = "bold", size = 3.45, lineheight = 0.95) +
  # Left header block
  annotate("rect", xmin = -3.9, xmax = 0.45, ymin = header_ymin, ymax = group_ymax,
           fill = "#EDF2F5", color = "white", linewidth = 0.7) +
  annotate("text", x = -3.62, y = (header_ymin + group_ymax) / 2,
           label = sprintf("Study (n = %d)", n_studies), hjust = 0,
           family = "Arial", fontface = "bold", size = 4.6) +
  # Summary row
  annotate("rect", xmin = -3.9, xmax = 0.45, ymin = summary_y - 0.70, ymax = summary_y + 0.70,
           fill = summary_fill, color = "white", linewidth = 0.6) +
  annotate("text", x = -3.62, y = summary_y, label = "Reported, n (%)",
           hjust = 0, family = "Arial", fontface = "bold", size = 4.0) +
  geom_tile(
    data = summary_dat, aes(x = x, y = summary_y),
    width = 0.98, height = 1.4, fill = summary_fill, color = "white", linewidth = 0.6
  ) +
  geom_text(
    data = summary_dat, aes(x = x, y = summary_y, label = label),
    family = "Arial", fontface = "bold", size = 3.65, lineheight = 1.05
  ) +
  annotate("rect", xmin = 12.52, xmax = 13.48, ymin = summary_y - 0.70, ymax = summary_y + 0.70,
           fill = summary_fill, color = "white", linewidth = 0.6) +
  scale_fill_manual(
    values = c("Reported" = reported_blue, "Not reported" = not_reported_grey),
    breaks = c("Reported", "Not reported"),
    labels = c("Reported (Y)", "Not reported (NR)"),
    name = NULL
  ) +
  scale_x_continuous(limits = c(-3.9, 13.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(summary_y - 0.8, group_ymax + 0.25), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(title = sprintf("CARE-related practices across studies (n = %d)", n_studies)) +
  guides(fill = guide_legend(override.aes = list(color = NA), nrow = 1)) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(family = "Arial", face = "bold", size = 21, hjust = 0,
                              margin = margin(b = 14)),
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(family = "Arial", size = 10),
    legend.key.width = grid::unit(0.55, "in"),
    plot.margin = margin(12, 18, 12, 12)
  )

svg_file <- file.path(figure_dir, "care_reporting_heatmap.svg")
png_file <- file.path(figure_dir, "care_reporting_heatmap.png")

ggsave(svg_file, p, width = 16.8, height = 13.2, units = "in",
       device = svglite::svglite, bg = "white")
ggsave(png_file, p, width = 16.8, height = 13.2, units = "in",
       dpi = 600, bg = "white", limitsize = FALSE)

write.csv(summary_dat[, c("code", "reported", "percent")],
          file.path(result_dir, "care_reporting_summary.csv"), row.names = FALSE)

message("Saved CARE heatmap to: ", normalizePath(figure_dir, winslash = "/"))
