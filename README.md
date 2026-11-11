# Indigenous populations and the oral microbiome review

This repository contains cleaned study-level data, R scripts, figures, and selected analysis results supporting a systematic review of oral microbiome research involving Indigenous populations.

The repository is intentionally limited to the files needed to reproduce the main visualisations and the final meta-analyses. Internal extraction notes, article files, manuscript drafts, and superseded analyses are not included.

## Repository structure

```text
data/
  evidence_mapping/   Country, sample type, method, and health-condition summaries
  alpha_diversity/    Direction-of-effect synthesis data
  taxa_direction/     Taxonomic direction and significance data
  meta_analysis/      Generally healthy adult Shannon and Faith's PD data
  care_principles/    CARE C1-E3 reporting matrix

scripts/
  evidence_mapping/   Map, donut, and bar-chart scripts
  alpha_diversity/    Four- and eight-metric alpha-diversity direction figures
  taxa_direction/     Overall and healthy-population taxa figures
  meta_analysis/      REML + Hartung-Knapp meta-analysis and sensitivity analysis
  care_principles/    CARE reporting heatmap

figures/              Selected SVG outputs
results/              Meta-analysis and CARE summary tables
```

## Main analyses

### Evidence mapping

- Geographic distribution of included studies
- Sample-type occurrences
- Microbiome-method occurrences
- Health-condition distribution

### Alpha-diversity direction synthesis

The four-metric figure includes Shannon diversity, Chao1 richness, observed richness, and Faith's phylogenetic diversity. The eight-metric version retains the additional measures available in the extracted dataset.

### Taxa direction synthesis

The direction plots distinguish increased, decreased, unchanged, and not mentioned taxa. An asterisk denotes a statistically significant reported difference. Separate figures are provided for the overall synthesis and generally healthy adult/child comparisons.

### Meta-analysis

The final meta-analyses are restricted to generally healthy adults and include:

- Shannon diversity
- Faith's phylogenetic diversity
- random-effects models with REML heterogeneity estimation
- standard Hartung-Knapp adjustment
- leave-one-study-out sensitivity analyses

### CARE principles

The CARE heatmap reports whether each study described practices corresponding to 12 subprinciples:

- C1: for inclusive development and innovation
- C2: for improved governance and citizen engagement
- C3: for equitable outcomes
- A1: recognizing rights and interests
- A2: data for governance
- A3: governance of data
- R1: for positive relationships
- R2: for expanding capability and capacity
- R3: for Indigenous languages and worldviews
- E1: for minimizing harm and maximizing benefit
- E2: for justice
- E3: for future use

`Y` means the indicator was reported and `NR` means it was not reported.

## Running the scripts

Run scripts from the repository root. Each script derives the repository location from its own file path and creates its output directory when needed.

Examples:

```bash
Rscript scripts/alpha_diversity/plot_alpha_direction_4metrics.R
Rscript scripts/meta_analysis/run_adult_healthy_alpha_meta_analysis.R
Rscript scripts/care_principles/plot_care_reporting_heatmap.R
```

### Input-to-output map

| Analysis | Main input | Script | Main output |
|---|---|---|---|
| Geographic evidence map | `data/evidence_mapping/country_studies.csv` | `scripts/evidence_mapping/plot_review_country_distribution.R` | `figures/evidence_mapping/country_distribution.svg` |
| Descriptive summaries | `data/evidence_mapping/*.csv` | `scripts/evidence_mapping/plot_*` | `figures/evidence_mapping/` |
| Alpha-diversity direction | `data/alpha_diversity/alpha_direction_8metrics.csv` | `scripts/alpha_diversity/plot_alpha_direction_*.R` | `figures/alpha_diversity/` |
| Taxa direction | `data/taxa_direction/*.csv` | `scripts/taxa_direction/plot_taxa_direction_*.R` | `figures/taxa_direction/` |
| Adult healthy meta-analysis | `data/meta_analysis/*.csv` | `scripts/meta_analysis/run_adult_healthy_alpha_meta_analysis.R` | `figures/meta_analysis/` and `results/meta_analysis/` |
| CARE reporting | `data/care_principles/care_indicators.csv` | `scripts/care_principles/plot_care_reporting_heatmap.R` | `figures/care_principles/care_reporting_heatmap.svg` |

## Software

Analyses and visualisations were developed in R 4.5.1. Key packages include `ggplot2`, `dplyr`, `tidyr`, `patchwork`, `svglite`, `sf`, `rnaturalearth`, and `metafor`.

Package requirements are checked at the beginning of each script. Local package libraries and machine-specific paths are not included in the repository.

## Data scope

The deposited files contain aggregated, study-level information extracted from publicly available reports. They do not contain participant-level data. CARE internal notes and other working annotations have been excluded from the public dataset.
