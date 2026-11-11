# Generally healthy adult alpha-diversity meta-analyses.
# Primary model: random effects with REML heterogeneity estimation and the
# standard Hartung-Knapp adjustment. Also performs leave-one-study-out analyses.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[1]), mustWork = FALSE)
} else {
  normalizePath("scripts/meta_analysis/run_adult_healthy_alpha_meta_analysis.R", mustWork = FALSE)
}

project_dir <- dirname(dirname(dirname(script_file)))
project_library <- file.path(project_dir, "R_library")
if (dir.exists(project_library)) .libPaths(c(project_library, .libPaths()))

required <- c("metafor", "ggplot2", "patchwork", "svglite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

library(metafor)
library(ggplot2)
library(patchwork)

analyses <- data.frame(
  id = c("faith_pd", "shannon"),
  file = c("adult_healthy_faith_pd.csv", "adult_healthy_shannon.csv"),
  title = c(
    "Faith's phylogenetic diversity: generally healthy adults",
    "Shannon diversity: generally healthy adults"
  ),
  stringsAsFactors = FALSE
)

figure_dir <- file.path(project_dir, "figures", "meta_analysis")
result_dir <- file.path(project_dir, "results", "meta_analysis")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

fmt_num <- function(x) formatC(x, format = "fg", digits = 4, flag = "#")
fmt_p <- function(x) {
  if (x < 0.001) "< 0.001" else paste0("= ", formatC(x, format = "f", digits = 3))
}
fmt_effect <- function(est, lb, ub) sprintf("%.2f [%.2f, %.2f]", est, lb, ub)

fit_hk <- function(es) {
  rma(yi, vi, data = es, method = "REML", test = "knha")
}

extract_fit <- function(fit, analysis_id, omitted_study) {
  data.frame(
    analysis = analysis_id,
    omitted_study = omitted_study,
    k = fit$k,
    estimate = as.numeric(fit$b),
    se = fit$se,
    tval = fit$zval,
    df = fit$ddf,
    pval = fit$pval,
    ci.lb = fit$ci.lb,
    ci.ub = fit$ci.ub,
    tau2 = fit$tau2,
    Q = fit$QE,
    Qp = fit$QEp,
    I2 = fit$I2,
    H2 = fit$H2,
    stringsAsFactors = FALSE
  )
}

read_effect_data <- function(input_file) {
  dat <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
  needed <- c("study", "comparison", "n1", "m1", "sd1", "n0", "m0", "sd0")
  if (!all(needed %in% names(dat))) {
    stop("Missing columns in ", input_file)
  }
  if (anyDuplicated(dat$study)) {
    stop("Expected one aggregated row per study in ", input_file)
  }
  numeric_columns <- c("n1", "m1", "sd1", "n0", "m0", "sd0")
  dat[numeric_columns] <- lapply(dat[numeric_columns], as.numeric)
  if (any(!complete.cases(dat[c("study", "comparison", numeric_columns)]))) {
    stop("Missing study labels or meta-analysis inputs in ", input_file)
  }
  if (any(dat$n1 < 2 | dat$n0 < 2 | dat$sd1 <= 0 | dat$sd0 <= 0)) {
    stop("Each group must have n >= 2 and SD > 0 in ", input_file)
  }

  es <- escalc(
    measure = "SMD",
    m1i = m1, sd1i = sd1, n1i = n1,
    m2i = m0, sd2i = sd0, n2i = n0,
    data = dat
  )
  dat$yi <- es$yi
  dat$vi <- es$vi
  dat$sei <- sqrt(dat$vi)
  dat$ci_lb <- dat$yi - qnorm(0.975) * dat$sei
  dat$ci_ub <- dat$yi + qnorm(0.975) * dat$sei
  list(dat = dat, es = es)
}

make_forest_plot <- function(dat, fit, title, output_stem) {
  dat$weight_raw <- 1 / (dat$vi + fit$tau2)
  dat$weight <- 100 * dat$weight_raw / sum(dat$weight_raw)

  k <- nrow(dat)
  dat$y <- rev(seq_len(k)) + 2
  total_y <- 1.55
  header_y <- k + 3.35
  subheader_y <- k + 2.55

  table_theme <- theme_void(base_size = 6.2) +
    theme(plot.margin = margin(2, 0, 2, 3))

  table_plot <- ggplot() +
    annotate(
      "text", x = 1.5, y = header_y + 0.55,
      label = title, hjust = 0, fontface = "bold", size = 2.35
    ) +
    annotate(
      "text", x = 1.5, y = subheader_y,
      label = "Study", hjust = 0, fontface = "bold", size = 2.0
    ) +
    annotate(
      "text", x = 53, y = header_y,
      label = "First-listed group", fontface = "bold", size = 2.0
    ) +
    annotate(
      "text", x = 79, y = header_y,
      label = "Comparator", fontface = "bold", size = 2.0
    ) +
    annotate(
      "text", x = c(48, 56, 64, 74, 82, 90), y = subheader_y,
      label = rep(c("Mean", "SD", "Total"), 2),
      fontface = "bold", size = 1.85
    ) +
    annotate(
      "text", x = 98, y = subheader_y,
      label = "Weight", fontface = "bold", size = 1.85
    ) +
    annotate(
      "text", x = 114, y = header_y,
      label = "Hedges' g", fontface = "bold", size = 2.0
    ) +
    annotate(
      "text", x = 114, y = subheader_y,
      label = "REML + HK, 95% CI", fontface = "bold", size = 1.85
    ) +
    geom_segment(
      aes(x = 0, xend = 125, y = subheader_y - 0.38, yend = subheader_y - 0.38),
      linewidth = 0.45, color = "#333333"
    ) +
    geom_text(
      data = dat, aes(x = 1.5, y = y, label = study),
      hjust = 0, fontface = "bold", size = 2.1
    ) +
    geom_text(data = dat, aes(x = 48, y = y, label = fmt_num(m1)), size = 1.45) +
    geom_text(data = dat, aes(x = 56, y = y, label = fmt_num(sd1)), size = 1.45) +
    geom_text(data = dat, aes(x = 64, y = y, label = n1), size = 1.45) +
    geom_text(data = dat, aes(x = 74, y = y, label = fmt_num(m0)), size = 1.45) +
    geom_text(data = dat, aes(x = 82, y = y, label = fmt_num(sd0)), size = 1.45) +
    geom_text(data = dat, aes(x = 90, y = y, label = n0), size = 1.45) +
    geom_text(
      data = dat, aes(x = 98, y = y, label = sprintf("%.1f%%", weight)),
      size = 1.45
    ) +
    geom_text(
      data = dat, aes(x = 114, y = y, label = fmt_effect(yi, ci_lb, ci_ub)),
      size = 1.45
    ) +
    annotate(
      "text", x = 0, y = total_y, label = "Total (95% CI)",
      hjust = 0, fontface = "bold", size = 1.7
    ) +
    annotate("text", x = 64, y = total_y, label = sum(dat$n1),
             fontface = "bold", size = 1.6) +
    annotate("text", x = 90, y = total_y, label = sum(dat$n0),
             fontface = "bold", size = 1.6) +
    annotate("text", x = 98, y = total_y, label = "100.0%",
             fontface = "bold", size = 1.6) +
    annotate(
      "text", x = 114, y = total_y,
      label = fmt_effect(as.numeric(fit$b), fit$ci.lb, fit$ci.ub),
      fontface = "bold", size = 1.6
    ) +
    annotate(
      "text", x = 0, y = 0.78, hjust = 0, size = 1.45,
      label = sprintf(
        "Heterogeneity: tau^2 = %.2f; Q = %.2f, df = %d (P %s); I^2 = %.0f%%",
        fit$tau2, fit$QE, fit$k - 1, fmt_p(fit$QEp), fit$I2
      )
    ) +
    annotate(
      "text", x = 0, y = 0.25, hjust = 0, size = 1.45,
      label = sprintf(
        "Hartung-Knapp test: t(%d) = %.2f (P %s)",
        fit$ddf, fit$zval, fmt_p(fit$pval)
      )
    ) +
    scale_x_continuous(limits = c(0, 125), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, header_y + 1.0), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    table_theme

  x_min <- min(dat$ci_lb, fit$ci.lb, 0)
  x_max <- max(dat$ci_ub, fit$ci.ub, 0)
  x_pad <- max(0.15, 0.12 * (x_max - x_min))
  x_limits <- c(x_min - x_pad, x_max + x_pad)

  forest_plot <- ggplot(dat, aes(x = yi, y = y)) +
    annotate(
      "segment", x = 0, xend = 0, y = 0, yend = subheader_y - 0.38,
      color = "#555555", linewidth = 0.55
    ) +
    geom_errorbar(
      aes(xmin = ci_lb, xmax = ci_ub), orientation = "y", width = 0,
      linewidth = 0.65, color = "#343434"
    ) +
    geom_point(aes(size = weight), shape = 15, color = "#15965B") +
    annotate(
      "polygon",
      x = c(fit$ci.lb, as.numeric(fit$b), fit$ci.ub, as.numeric(fit$b)),
      y = c(total_y, total_y + 0.28, total_y, total_y - 0.28),
      fill = "black"
    ) +
    annotate(
      "text", x = mean(x_limits), y = header_y,
      label = "Hedges' g", fontface = "bold", size = 2.0
    ) +
    annotate(
      "text", x = mean(x_limits), y = subheader_y,
      label = "REML + HK, 95% CI", fontface = "bold", size = 1.85
    ) +
    annotate(
      "segment", x = x_limits[1], xend = x_limits[2],
      y = subheader_y - 0.38, yend = subheader_y - 0.38,
      linewidth = 0.45, color = "#333333"
    ) +
    scale_size_area(max_size = 3.6, guide = "none") +
    scale_x_continuous(
      limits = x_limits, breaks = scales::breaks_pretty(n = 5), expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, header_y + 1.0), breaks = NULL, expand = c(0, 0)
    ) +
    labs(x = "Comparator higher        First-listed group higher", y = NULL) +
    theme_classic(base_size = 6.2) +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 5.8),
      axis.text.x = element_text(size = 5.5),
      plot.margin = margin(2, 3, 2, 0)
    )

  fig <- table_plot + forest_plot + plot_layout(widths = c(2.9, 1.0))
  figure_height <- max(2.4, 1.15 + 0.25 * nrow(dat))
  ggsave(
    file.path(figure_dir, paste0(output_stem, "_forest.svg")), fig,
    width = 6.75, height = figure_height, units = "in",
    device = svglite::svglite, bg = "white"
  )
  ggsave(
    file.path(figure_dir, paste0(output_stem, "_forest.png")), fig,
    width = 6.75, height = figure_height, units = "in",
    dpi = 350, bg = "white"
  )
}

make_loo_outputs <- function(es, fit, dat, id, title) {
  rows <- list(extract_fit(fit, id, "None (full model)"))
  for (i in seq_len(nrow(dat))) {
    reduced_fit <- fit_hk(es[-i, , drop = FALSE])
    rows[[i + 1L]] <- extract_fit(reduced_fit, id, dat$study[i])
  }
  result <- do.call(rbind, rows)

  write.csv(
    result,
    file.path(result_dir, paste0(id, "_leave_one_out.csv")),
    row.names = FALSE
  )
  report <- capture.output({
    cat(title, "\n")
    cat("Random-effects model: REML + standard Hartung-Knapp\n\n")
    cat("FULL MODEL\n")
    print(fit)
    cat("\nLEAVE-ONE-STUDY-OUT MODELS\n")
    print(result, row.names = FALSE, digits = 5)
  })
  writeLines(report, file.path(result_dir, paste0(id, "_leave_one_out.txt")))

  plot_dat <- result
  plot_dat$label <- ifelse(
    plot_dat$omitted_study == "None (full model)",
    "Full model",
    paste0("Omit ", plot_dat$omitted_study)
  )
  plot_dat$label <- factor(plot_dat$label, levels = rev(plot_dat$label))

  p <- ggplot(plot_dat, aes(x = estimate, y = label)) +
    geom_vline(xintercept = 0, color = "#777777", linewidth = 0.5) +
    geom_errorbar(
      aes(xmin = ci.lb, xmax = ci.ub), orientation = "y",
      width = 0, linewidth = 0.75, color = "#333333"
    ) +
    geom_point(
      aes(fill = omitted_study == "None (full model)"),
      shape = 21, size = 4.2, color = "white", stroke = 0.35
    ) +
    scale_fill_manual(
      values = c("TRUE" = "#111111", "FALSE" = "#15965B"),
      guide = "none"
    ) +
    labs(
      title = title,
      subtitle = "Leave-one-study-out sensitivity analysis (REML + Hartung-Knapp)",
      x = "Hedges' g (Hartung-Knapp 95% CI)", y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(10, 14, 10, 10)
    )

  height <- max(4.2, 1.0 + 0.55 * nrow(plot_dat))
  ggsave(
    file.path(figure_dir, paste0(id, "_leave_one_out.svg")), p,
    width = 8.5, height = height, units = "in",
    device = svglite::svglite, bg = "white"
  )
  ggsave(
    file.path(figure_dir, paste0(id, "_leave_one_out.png")), p,
    width = 8.5, height = height, units = "in", dpi = 350, bg = "white"
  )

  # Compact appendix panel: results table and forest plot share the same rows.
  panel_dat <- result
  panel_dat$row_label <- ifelse(
    panel_dat$omitted_study == "None (full model)",
    "Full model", panel_dat$omitted_study
  )
  panel_dat$y <- rev(seq_len(nrow(panel_dat))) + 1
  header_y <- nrow(panel_dat) + 2.45
  subheader_y <- nrow(panel_dat) + 1.82
  rule_y <- nrow(panel_dat) + 1.48

  table_part <- ggplot() +
    annotate("text", x = 0, y = header_y + 0.48,
             label = paste0(ifelse(id == "faith_pd", "A. ", "B. "), title),
             hjust = 0, fontface = "bold", size = 3.0) +
    annotate("text", x = 0, y = subheader_y, label = "Study omitted",
             hjust = 0, fontface = "bold", size = 2.35) +
    annotate("text", x = 48, y = header_y, label = "Studies",
             fontface = "bold", size = 2.25) +
    annotate("text", x = 48, y = subheader_y, label = "retained",
             fontface = "bold", size = 2.25) +
    annotate("text", x = 72, y = header_y, label = "Hedges' g",
             fontface = "bold", size = 2.25) +
    annotate("text", x = 72, y = subheader_y,
             label = "REML + HK, 95% CI", fontface = "bold", size = 2.15) +
    annotate("text", x = 94, y = subheader_y, label = "P value",
             fontface = "bold", size = 2.25) +
    annotate("text", x = 108, y = subheader_y, label = "I-squared",
             fontface = "bold", size = 2.25) +
    geom_segment(aes(x = 0, xend = 114, y = rule_y, yend = rule_y),
                 linewidth = 0.45, color = "#222222") +
    geom_text(data = panel_dat, aes(x = 2, y = y, label = row_label),
              hjust = 0, size = 2.35,
              fontface = ifelse(panel_dat$row_label == "Full model", "bold", "plain")) +
    geom_text(data = panel_dat, aes(x = 48, y = y, label = k), size = 2.15) +
    geom_text(data = panel_dat,
              aes(x = 72, y = y, label = fmt_effect(estimate, ci.lb, ci.ub)),
              size = 2.15) +
    geom_text(data = panel_dat,
              aes(x = 94, y = y,
                  label = ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))),
              size = 2.15) +
    geom_text(data = panel_dat, aes(x = 108, y = y, label = sprintf("%.0f%%", I2)),
              size = 2.15) +
    scale_x_continuous(limits = c(0, 114), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.7, header_y + 0.8), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void(base_size = 8) +
    theme(plot.margin = margin(2, 0, 2, 3))

  x_min <- min(panel_dat$ci.lb, 0)
  x_max <- max(panel_dat$ci.ub, 0)
  x_pad <- max(0.12, 0.08 * (x_max - x_min))
  x_limits <- c(x_min - x_pad, x_max + x_pad)

  forest_part <- ggplot(panel_dat, aes(x = estimate, y = y)) +
    annotate("segment", x = 0, xend = 0, y = 0.7, yend = rule_y,
             linewidth = 0.45, color = "#555555") +
    geom_errorbar(aes(xmin = ci.lb, xmax = ci.ub), orientation = "y",
                  width = 0, linewidth = 0.55, color = "#333333") +
    geom_point(aes(fill = row_label == "Full model"), shape = 22,
               size = 3.0, stroke = 0.35, color = "#333333") +
    scale_fill_manual(values = c("TRUE" = "#111111", "FALSE" = "#15965B"),
                      guide = "none") +
    annotate("text", x = mean(x_limits), y = header_y, label = "Hedges' g",
             fontface = "bold", size = 2.25) +
    annotate("text", x = mean(x_limits), y = subheader_y,
             label = "REML + HK, 95% CI", fontface = "bold", size = 2.15) +
    annotate("segment", x = x_limits[1], xend = x_limits[2],
             y = rule_y, yend = rule_y, linewidth = 0.45, color = "#222222") +
    scale_x_continuous(limits = x_limits, breaks = scales::breaks_pretty(n = 5),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.7, header_y + 0.8), breaks = NULL,
                       expand = c(0, 0)) +
    labs(x = "Non-Indigenous higher        Indigenous higher", y = NULL) +
    theme_classic(base_size = 8) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 7), axis.title.x = element_text(size = 7.2),
          plot.margin = margin(2, 3, 2, 0))

  compact_panel <- table_part + forest_part + plot_layout(widths = c(3.05, 1.0))
  list(result = result, panel = compact_panel)
}

all_results <- vector("list", nrow(analyses))
loo_panels <- vector("list", nrow(analyses))
for (i in seq_len(nrow(analyses))) {
  loaded <- read_effect_data(file.path(project_dir, "data", "meta_analysis", analyses$file[i]))
  fit <- fit_hk(loaded$es)
  make_forest_plot(loaded$dat, fit, analyses$title[i], analyses$id[i])
  loo_output <- make_loo_outputs(
    loaded$es, fit, loaded$dat, analyses$id[i], analyses$title[i]
  )
  all_results[[i]] <- loo_output$result
  loo_panels[[i]] <- loo_output$panel
}

combined <- do.call(rbind, all_results)
write.csv(
  combined,
  file.path(result_dir, "adult_healthy_alpha_hk_results.csv"),
  row.names = FALSE
)

appendix_panel <- wrap_plots(loo_panels, ncol = 1)
ggsave(
  file.path(figure_dir, "leave_one_out_panel.svg"),
  appendix_panel, width = 8.25, height = 7.1, units = "in",
  device = svglite::svglite, bg = "white"
)
ggsave(
  file.path(figure_dir, "leave_one_out_panel.png"),
  appendix_panel, width = 8.25, height = 7.1, units = "in",
  dpi = 600, bg = "white"
)

message(
  "Saved two REML + standard Hartung-Knapp analyses to: ",
  normalizePath(figure_dir, winslash = "/")
)
