# =============================================================================
# Complexity of Infection (COI) Analysis: Fws and COI Comparison
# =============================================================================
#
# Description:
#   Compares genetic complexity metrics (Fws and COI) between two temporal
#   populations of P. falciparum from Gambia (1966-71 vs 2015).
#   Creates publication-quality box plots with statistical comparisons.
#
# Input:
#   - source_data/fws_coi_comparison.rds
#   Expected columns: Fws, COI (Complexity of Infection)
#
# Output:
#   - Combined panel plot (Panel A: Fws, Panel B: COI)
#   - Statistical test results
#
# Author: Mouhamadou Fadel DIOP
# Date: [Date]
# =============================================================================

# Load required libraries
library(tidyverse)    # Data manipulation and visualization
library(readxl)       # Excel file reading
library(ggpubr)       # Publication-ready plots and statistics
library(patchwork)    # Combining plots

# -----------------------------------------------------------------------------
# Configuration and Parameters
# -----------------------------------------------------------------------------

# File paths
data_dir <- "source_data"

# Output
output_file <- "results/figures/fws_coi_comparison.pdf"

# Plot parameters
panel_labels <- c("A", "B")
y_label_fws <- "Fws"
y_label_coi <- "Complexity of Infection (COI)"
x_label <- "Year of Sampling"

# Color scheme
color_1966_71 <- "#E15759"  # Red
color_2015 <- "#4E79A7"     # Blue

# Point parameters
point_size <- 2
point_alpha <- 0.7
jitter_width <- 0.2

# Box plot parameters
box_width <- 0.5
box_alpha <- 0.3

# Statistical test parameters
stat_test_method <- "wilcox.test"  # Wilcoxon rank-sum test (Mann-Whitney U)
stat_label_y_fws <- 1.05
stat_label_y_coi <- 3.5

# Plot dimensions
plot_width_mm <- 180
plot_height_mm <- 80
plot_dpi <- 600

# -----------------------------------------------------------------------------
# Create Output Directory
# -----------------------------------------------------------------------------

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Data Import and Preparation
# -----------------------------------------------------------------------------

cat("Loading data...\n")

# Load data
data = readRDS("source_data/fws_coi_comparison.rds")

data_1966_71 <- data$data_1966_71
data_2015 <- data$data_2015

# Rename column and add temporal group identifier
data_1966_71 <- data_1966_71 %>%
   rename("sample_id" = `WSI sample ID`) %>%
   mutate(period = "pop1966-71")

data_2015 <- data_2015 %>%
   rename("sample_id" = Sample_ID) %>%
   mutate(period = "pop2015")

# Combine datasets
combined_data <- bind_rows(data_1966_71, data_2015) %>%
  mutate(
    period = factor(period, levels = c("pop1966-71", "pop2015"))
  )

# Display sample sizes
cat("\nSample sizes:\n")
combined_data %>%
  group_by(period) %>%
  summarise(
    n_samples = n(),
    .groups = "drop"
  ) %>%
  print()

# Display data summary
cat("\nData summary:\n")
combined_data %>%
  group_by(period) %>%
  summarise(
    across(
      c(Fws, COI),
      list(
        mean = ~mean(., na.rm = TRUE),
        median = ~median(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  print()

# -----------------------------------------------------------------------------
# Statistical Tests
# -----------------------------------------------------------------------------

cat("\nPerforming statistical tests...\n")

# Fws comparison
fws_test <- wilcox.test(
  Fws ~ period,
  data = combined_data,
  exact = FALSE
)

cat("\nFws comparison (Wilcoxon rank-sum test):\n")
print(fws_test)

# COI comparison
coi_test <- wilcox.test(
  COI ~ period,
  data = combined_data,
  exact = FALSE
)

cat("\nCOI comparison (Wilcoxon rank-sum test):\n")
print(coi_test)

# -----------------------------------------------------------------------------
# Panel A: Fws Plot
# -----------------------------------------------------------------------------

cat("\nGenerating Fws plot...\n")

plot_fws <- ggplot(combined_data, aes(x = period, y = Fws, fill = period)) +

  # Box plot
  geom_boxplot(
    width = box_width,
    alpha = box_alpha,
    outlier.shape = NA  # Don't show outliers (points will show all data)
  ) +

  # Individual data points with jitter
  geom_jitter(
    aes(color = period),
    width = jitter_width,
    size = point_size,
    alpha = point_alpha
  ) +

  # Statistical comparison
  stat_compare_means(
    method = stat_test_method,
    label = "p.signif",  # Show significance stars
    label.y = stat_label_y_fws,
    size = 6
  ) +

  # Color scheme
  scale_fill_manual(values = c(color_1966_71, color_2015)) +
  scale_color_manual(values = c(color_1966_71, color_2015)) +

  # Axis labels
  labs(
    x = x_label,
    y = y_label_fws
  ) +

  # Y-axis limits
  coord_cartesian(ylim = c(0.4, 1.1)) +

  # Theme
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10, color = "#000000"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.line = element_line(color = "#000000"),
    panel.border = element_rect(color = "#000000", fill = NA, linewidth = 1)
  )

# -----------------------------------------------------------------------------
# Panel B: COI Plot
# -----------------------------------------------------------------------------

cat("Generating COI plot...\n")

plot_coi <- ggplot(combined_data, aes(x = period, y = COI, fill = period)) +

  # Box plot
  geom_boxplot(
    width = box_width,
    alpha = box_alpha,
    outlier.shape = NA
  ) +

  # Individual data points with jitter
  geom_jitter(
    aes(color = period),
    width = jitter_width,
    size = point_size,
    alpha = point_alpha
  ) +

  # Statistical comparison
  stat_compare_means(
    method = stat_test_method,
    label = "p.signif",
    label.y = stat_label_y_coi,
    size = 6
  ) +

  # Color scheme
  scale_fill_manual(values = c(color_1966_71, color_2015)) +
  scale_color_manual(values = c(color_1966_71, color_2015)) +

  # Axis labels
  labs(
    x = x_label,
    y = y_label_coi
  ) +

  # Y-axis limits
  coord_cartesian(ylim = c(0, 3.8)) +

  # Theme
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10, color = "#000000"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.line = element_line(color = "#000000"),
    panel.border = element_rect(color = "#000000", fill = NA, linewidth = 1)
  )

# -----------------------------------------------------------------------------
# Combine Plots
# -----------------------------------------------------------------------------

cat("Combining plots...\n")

# Combine panels side by side
combined_plot <- plot_fws + plot_coi +
  plot_annotation(
    tag_levels = list(panel_labels)
  ) &
  theme(
    plot.tag = element_text(size = 14, face = "bold")
  )

# Display combined plot
print(combined_plot)

# -----------------------------------------------------------------------------
# Save Plot
# -----------------------------------------------------------------------------

cat("Saving plot...\n")

ggsave(
  filename = output_file,
  plot = combined_plot,
  width = plot_width_mm,
  height = plot_height_mm,
  units = "mm",
  dpi = plot_dpi
)

cat("Plot saved to:", output_file, "\n")

# -----------------------------------------------------------------------------
# Save Statistical Results
# -----------------------------------------------------------------------------

# Create results summary
results_summary <- data.frame(
  Comparison = c("Fws", "COI"),
  Test = c("Wilcoxon rank-sum", "Wilcoxon rank-sum"),
  P_value = c(fws_test$p.value, coi_test$p.value),
  Statistic = c(fws_test$statistic, coi_test$statistic)
)

# Save to CSV
results_file <- "results/statistical_tests_fws_coi.csv"
write_csv(results_summary, results_file)

cat("\nStatistical results saved to:", results_file, "\n")

# Display results
cat("\nStatistical Test Results:\n")
print(results_summary)

cat("\nAnalysis complete!\n")

# -----------------------------------------------------------------------------
# Session Info (for reproducibility)
# -----------------------------------------------------------------------------

# Uncomment to print session information
# sessionInfo()

# =============================================================================
# Notes:
# =============================================================================
#
# Fws (Within-sample F-statistic):
#   - Measure of inbreeding/clonality within individual infections
#   - Range: 0 to 1
#   - Higher values indicate lower diversity (more clonal)
#   - Lower values indicate higher diversity (polyclonal infections)
#
# COI (Complexity of Infection):
#   - Number of distinct parasite clones within a single infection
#   - Integer values (typically 1-5+)
#   - Higher values indicate multiple clone infections (polyclonal)
#
# Statistical Test:
#   - Wilcoxon rank-sum test (Mann-Whitney U test)
#   - Non-parametric test for comparing two independent groups
#   - Appropriate when data may not be normally distributed
#
# Significance levels:
#   ns: p > 0.05
#   *: p <= 0.05
#   **: p <= 0.01
#   ***: p <= 0.001
#   ****: p <= 0.0001
#
# =============================================================================
