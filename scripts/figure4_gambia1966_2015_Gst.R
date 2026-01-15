# =============================================================================
# Manhattan Plot: Genetic Differentiation (Gst) Across Chromosomes
# =============================================================================
#
# Description:
#   Creates a Manhattan plot visualizing Gst values across chromosomes for
#   Gambia samples (1966-2015). This plot shows genetic differentiation
#   patterns genome-wide.
#
# Input:
#   - File: results/tables/gambia1966_2015_mmod.txt
#   - Expected columns: chr, pos, Gst
#
# Output:
#   - PDF figure: results/figures/figure4_gambia1966_2015_Gst.pdf
#
# Author: Mouhamadou Fadel DIOP
# Date: [Date]
# =============================================================================

# Load required libraries
library(tidyverse)

# -----------------------------------------------------------------------------
# Configuration and Parameters
# -----------------------------------------------------------------------------

# Define color scheme for alternating chromosomes
chromosome_colors <- c("#000000", "#FFCC00")  # Black and gold

# Input/output file paths
input_file <- "results/tables/gambia1966_2015_mmod.txt"
output_file <- "results/figures/figure4_gambia1966_2015_Gst.pdf"

# Plot parameters
y_axis_limits <- c(0, 1)
point_size <- 2
plot_width_mm <- 190
plot_height_mm <- 150
plot_dpi <- 600

# -----------------------------------------------------------------------------
# Data Import
# -----------------------------------------------------------------------------

# Read Fst/Gst data
# Expected columns: chr (chromosome), pos (position), Gst (genetic differentiation)
fst_data <- read_delim(
  input_file,
  show_col_types = FALSE
)

# -----------------------------------------------------------------------------
# Data Transformation: Calculate Cumulative Positions
# -----------------------------------------------------------------------------

# Calculate cumulative positions for Manhattan plot
# This transforms chromosome-specific positions into genome-wide positions
fst_transformed <- fst_data %>%
  # Step 1: Compute the length of each chromosome
  group_by(chr) %>%
  summarise(chr_len = max(pos), .groups = "drop") %>%

  # Step 2: Calculate cumulative position offset for each chromosome
  # (sum of all previous chromosome lengths)
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%

  # Step 3: Join cumulative offsets back to original data
  left_join(fst_data, ., by = "chr") %>%

  # Step 4: Calculate genome-wide cumulative position for each SNP
  # and convert chromosome to numeric factor for coloring
  arrange(chr, pos) %>%
  mutate(
    BPcum = pos + tot,                    # Cumulative base pair position
    chr = as.integer(factor(chr))         # Numeric chromosome ID
  ) %>%
  select(-tot)

# -----------------------------------------------------------------------------
# Calculate Chromosome Centers for X-axis Labels
# -----------------------------------------------------------------------------

# Calculate the center position of each chromosome for axis labeling
chromosome_axis <- fst_transformed %>%
  group_by(chr) %>%
  summarize(center = mean(BPcum), .groups = "drop")

# -----------------------------------------------------------------------------
# Create Manhattan Plot
# -----------------------------------------------------------------------------

manhattan_plot <- fst_transformed %>%
  ggplot(aes(x = BPcum, y = Gst, color = as_factor(chr))) +

  # Plot all SNP points
  geom_point(size = point_size) +

  # X-axis: chromosome numbers at chromosome centers
  scale_x_continuous(
    label = chromosome_axis$chr,
    breaks = chromosome_axis$center) +

  # Y-axis: Gst values from 0 to 1
  scale_y_continuous(
    expand = c(0, 0),
    limits = y_axis_limits) +

  # Alternate chromosome colors
  scale_color_manual(
    values = rep(chromosome_colors, length.out = nrow(chromosome_axis))) +

  # Axis labels
  labs(
    x = "Chromosome",
    y = "Gst") +

  # Theme customization
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 9, color = "#000000"),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(color = "#000000"))

# -----------------------------------------------------------------------------
# Save Plot
# -----------------------------------------------------------------------------

# Create output directory if it doesn't exist
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# Save the plot
ggsave(
  filename = output_file,
  plot = manhattan_plot,
  units = "mm",
  width = plot_width_mm,
  height = plot_height_mm,
  dpi = plot_dpi
)

# Print confirmation message
cat("Manhattan plot saved to:", output_file, "\n")

# -----------------------------------------------------------------------------
# Session Info (for reproducibility)
# -----------------------------------------------------------------------------

# Uncomment to print session info
# sessionInfo()
