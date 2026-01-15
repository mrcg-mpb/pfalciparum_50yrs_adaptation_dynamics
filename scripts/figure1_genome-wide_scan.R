# =============================================================================
# Genome-wide Selection Scan: Beta Score and iHS Analysis
# =============================================================================
#
# Description:
#   Creates a two-panel Manhattan plot showing signatures of natural selection
#   in P. falciparum from The Gambia (1966-1971). Panel A displays beta scores
#   (balancing selection) and Panel B shows iHS scores (directional selection).
#
# Methods:
#   - Beta scores: Identify balancing selection via correlated SNP allele frequencies
#   - iHS (integrated Haplotype Score): Detects recent directional selection
#     based on extended haplotype homozygosity
#
# Input:
#   - Beta score data: data frame with columns Chr, Start, betaScore
#   - iHS data: data frame with columns Chr, Start, iHS (or similar p-value metric)
#   - Gene annotations: data frame with Chr, Start, End, Gene for labeling
#
# Output:
#   - Two-panel Manhattan plot (PDF/PNG)
#   - Panel A: Beta scores across genome
#   - Panel B: -log10(p) for iHS scores
#
# Author: [Your Name]
# Date: [Date]
# =============================================================================

# Load required libraries
library(tidyverse)
library(ggrepel)      # For gene label positioning
library(patchwork)    # For combining plots

# -----------------------------------------------------------------------------
# Configuration and Parameters
# -----------------------------------------------------------------------------

# File paths (adjust according to your data structure)
betascore_file <- "source_data/gambia1966_71_betaScores.xlsx"
ihs_file <- "source_data/gampf1966_71_iHS.xlsx"
gene_annotations_file <- "source_data/metadata/PF_GeneName.tsv"

# Output
output_file <- "results/figures/figure1_selection_scan.pdf"
output_format <- "pdf"  # Options: "pdf", "png"

# Plot parameters
chromosome_colors <- c("gray20", "orange")  # Alternating chromosome colors
point_size_beta <- 1.2
point_size_ihs <- 1.0

# Significance thresholds
ihs_significance_threshold <- -log10(0.05)  # Adjust based on your threshold
beta_score_highlight_threshold <- 10  # Adjust based on your data

# Gene labeling parameters
n_top_genes_beta <- 26      # Number of genes to label in beta score panel
n_top_genes_ihs <- 20       # Number of genes to label in iHS panel
label_text_size <- 2.5
label_box_padding <- 0.25

# Genes with significant iHS scores (bold labels in boxes)
# Update this list based on Table 1 from your manuscript
genes_with_significant_ihs <- c(
  "Clag3.1", "Clag3.2", "Surf 1.2", "Surf 4.2",
  "eba-165", "Msp1", "AMA-1", "TRAP"
  # Add other genes from your analysis
)

# Plot dimensions
plot_width_mm <- 180
plot_height_mm <- 140
plot_dpi <- 600

# -----------------------------------------------------------------------------
# Helper Function: Manhattan Plot for Beta Scores
# -----------------------------------------------------------------------------

plot_betaScore <- function(data,
                          highlight_genes = NULL,
                          threshold = NULL,
                          point_size = 1.2) {

  # Calculate cumulative positions
  data_transformed <- data %>%
    group_by(Chr) %>%
    summarise(chr_len = max(Start), .groups = "drop") %>%
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    left_join(data, ., by = "Chr") %>%
    arrange(Chr, Start) %>%
    mutate(
      BPcum = Start + tot,
      Chr = as.integer(factor(Chr, levels = unique(Chr)))
    )

  # Calculate chromosome centers for x-axis
  axis_set <- data_transformed %>%
    group_by(Chr) %>%
    summarize(center = mean(BPcum), .groups = "drop")

  # Set y-axis limits
  ylimits <- c(
    floor(min(data_transformed$betaScore, na.rm = TRUE)),
    ceiling(max(data_transformed$betaScore, na.rm = TRUE)) + 6
  )

  # Create base plot
  p <- data_transformed %>%
    ggplot(aes(x = BPcum, y = betaScore, color = as_factor(Chr))) +

    # All points
    geom_point(size = point_size) +

    # Baseline at y = 0
    geom_hline(yintercept = 0, color = "grey40", linetype = "dashed") +

    # Highlight threshold if provided
    {if (!is.null(threshold))
      geom_hline(yintercept = threshold, color = "red", linetype = "dashed", alpha = 0.5)
    } +

    # Color scheme
    scale_color_manual(values = rep(chromosome_colors, length.out = nrow(axis_set))) +

    # X-axis: chromosomes
    scale_x_continuous(
      label = axis_set$Chr,
      breaks = axis_set$center,
      expand = c(0.01, 0.01)
    ) +

    # Y-axis
    scale_y_continuous(expand = c(0, 0), limits = ylimits) +

    # Labels
    labs(
      x = "Chromosomes",
      y = expression("Beta Scores (" * beta * "1)")) +

    # Theme
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 8, color = "#000000"),
      axis.title = element_text(size = 9, face = "bold"),
      axis.line = element_line(color = "#000000", linewidth = 1),
      axis.ticks = element_line(color = 'black', linewidth = .7),
      axis.ticks.length = unit(.20, "cm")
    )

  # Add gene labels if provided
  if (!is.null(highlight_genes)) {

    # Merge gene annotations with data
    labeled_data <- highlight_genes %>% 
      group_by(Chr) %>%
      summarise(chr_len = max(Start), .groups = "drop") %>%
      mutate(tot = cumsum(chr_len) - chr_len) %>%
      select(-chr_len) %>%
      left_join(highlight_genes, ., by = "Chr") %>%
      select(-c(NewGeneName, GeneName)) %>% 
      arrange(Chr, Start) %>%
      mutate(
        BPcum = Start + tot,
        Chr = as.integer(factor(Chr, levels = unique(Chr)))
      ) %>%
      arrange(betaScore) %>% 
      rename(Gene = `Gene Symbol or description`) %>% 
      mutate(Gene = tolower(gsub("_", " ", Gene))) %>% 
      distinct(Gene, .keep_all = TRUE)

    # Add gene labels with text repel
    p <- p +
      geom_text_repel(
        data = labeled_data,
        aes(label = Gene),
        size = label_text_size,
        box.padding = label_box_padding,
        max.overlaps = Inf,
        segment.size = 0.2,
        segment.color = "grey50",
        color = "black",
        fontface = "bold"
        # fontface = ifelse(
        #   labeled_data$Gene %in% genes_with_significant_ihs,
        #   "bold",
        #   "plain"
        # )
      )

    # Add boxes around genes with significant iHS
    genes_with_boxes <- labeled_data %>%
      filter(Gene %in% genes_with_significant_ihs)

    if (nrow(genes_with_boxes) > 0) {
      p <- p +
        geom_point(
          data = genes_with_boxes,
          aes(x = BPcum, y = betaScore),
          shape = 0,  # Square outline
          size = 3,
          color = "black",
          stroke = 0.5
        )
    }
  }

  return(p)
}

# -----------------------------------------------------------------------------
# Helper Function: Manhattan Plot for iHS Scores
# -----------------------------------------------------------------------------

plot_iHS <- function(data,
                    highlight_genes = NULL,
                    threshold = NULL,
                    point_size = 1.0) {

  # Transform p-values to -log10 scale if needed
  # Adjust column name based on your data
  if (!"neg_log10_p" %in% names(data)) {
    data <- data %>%
      mutate(neg_log10_p = -log10(p_value))  # Adjust 'p_value' to your column name
  }

  # Calculate cumulative positions
  data_transformed <- data %>%
    group_by(Chr) %>%
    summarise(chr_len = max(Start), .groups = "drop") %>%
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    left_join(data, ., by = "Chr") %>%
    arrange(Chr, Start) %>%
    mutate(
      BPcum = Start + tot,
      Chr = as.integer(factor(Chr, levels = unique(Chr)))
    )

  # Calculate chromosome centers
  axis_set <- data_transformed %>%
    group_by(Chr) %>%
    summarize(center = mean(BPcum), .groups = "drop")

  # Set y-axis limits
  ylimits <- c(
    0,
    ceiling(max(data_transformed$neg_log10_p, na.rm = TRUE)) + 1
  )

  # Create base plot
  p <- data_transformed %>%
    ggplot(aes(x = BPcum, y = neg_log10_p, color = as_factor(Chr))) +

    # All points
    geom_point(size = point_size, alpha = 0.7) +

    # Significance threshold
    {if (!is.null(threshold))
      geom_hline(yintercept = threshold, color = "red", linetype = "dashed", alpha = 0.5)
    } +

    # Color scheme
    scale_color_manual(values = rep(chromosome_colors, length.out = nrow(axis_set))) +

    # X-axis
    scale_x_continuous(
      label = axis_set$Chr,
      breaks = axis_set$center,
      expand = c(0.01, 0.01)
    ) +

    # Y-axis
    scale_y_continuous(expand = c(0, 0), limits = ylimits) +

    # Labels
    labs(
      x = "Chromosomes",
      y = expression("-log"[10] * "(p) for iHS")
    ) +

    # Theme
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 8, color = "#000000"),
      axis.title = element_text(size = 9, face = "bold"),
      axis.line = element_line(color = "#000000"),
      panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.5)
    )

  # Add gene labels if provided
  if (!is.null(highlight_genes)) {

    labeled_data <- data_transformed %>%
      inner_join(highlight_genes, by = c("Chr", "Start")) %>%
      arrange(desc(neg_log10_p))

    p <- p +
      geom_text_repel(
        data = labeled_data,
        aes(label = Gene),
        size = label_text_size,
        box.padding = label_box_padding,
        max.overlaps = Inf,
        segment.size = 0.2,
        segment.color = "grey50",
        color = "black"
      )
  }

  return(p)
}

# -----------------------------------------------------------------------------
# Data Import
# -----------------------------------------------------------------------------

cat("Loading data...\n")

# Load beta score data
# Expected columns: Chr, Start, betaScore
beta_data <- read_xlsx(betascore_file, skip = 1)

# Load iHS data
# Expected columns: Chr, Start, p_value (or iHS score)
ihs_data <- read_delim(
  ihs_file,
  show_col_types = FALSE
)

# Load gene annotations
# Expected columns: Chr, Start, End, Gene
gene_annotations <- read_delim(
  gene_annotations_file,
  show_col_types = FALSE
)

# -----------------------------------------------------------------------------
# Identify Top Genes for Labeling
# -----------------------------------------------------------------------------

cat("Identifying genes for labeling...\n")

# Top genes by beta score
top_beta_genes <- beta_data %>%
  rename(Gene_ID = `Gene ID`, Start = `Start (Bp)`, betaScore = BetaScore) %>%  # Chr = Chromosome, 
  arrange(desc(betaScore)) %>%
  head(n_top_genes_beta) %>%
  left_join(gene_annotations, by = c("Gene_ID", "Chromosome", "Start"))

# Top genes by iHS
top_ihs_genes <- ihs_data %>%
  mutate(neg_log10_p = -log10(p_value)) %>%  # Adjust column name
  arrange(desc(neg_log10_p)) %>%
  head(n_top_genes_ihs) %>%
  left_join(gene_annotations, by = c("Chr", "Start"))

# -----------------------------------------------------------------------------
# Create Plots
# -----------------------------------------------------------------------------

cat("Generating Panel A (Beta Scores)...\n")
beta_data <- beta_data %>% 
  mutate(Chr = as.numeric(factor(Chromosome))) %>% 
  rename(Start = `Start (Bp)`, betaScore = BetaScore)

panel_a <- plot_betaScore(
  data = beta_data,
  highlight_genes = top_beta_genes,
  threshold = beta_score_highlight_threshold,
  point_size = point_size_beta
)

cat("Generating Panel B (iHS Scores)...\n")

panel_b <- plot_iHS(
  data = ihs_data,
  highlight_genes = top_ihs_genes,
  threshold = ihs_significance_threshold,
  point_size = point_size_ihs
)

# -----------------------------------------------------------------------------
# Combine Panels
# -----------------------------------------------------------------------------

cat("Combining panels...\n")

combined_plot <- panel_a / panel_b +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(size = 12, face = "bold")
  )

# Display plot
print(combined_plot)

# -----------------------------------------------------------------------------
# Save Plot
# -----------------------------------------------------------------------------

cat("Saving plot...\n")

# Create output directory
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# Save based on format
if (output_format == "pdf") {
  ggsave(
    filename = output_file,
    plot = combined_plot,
    width = plot_width_mm,
    height = plot_height_mm,
    units = "mm",
    dpi = plot_dpi
  )
} else if (output_format == "png") {
  output_file <- sub("\\.pdf$", ".png", output_file)
  ggsave(
    filename = output_file,
    plot = combined_plot,
    width = plot_width_mm,
    height = plot_height_mm,
    units = "mm",
    dpi = plot_dpi
  )
}

cat("Plot saved to:", output_file, "\n")

# -----------------------------------------------------------------------------
# Summary Statistics
# -----------------------------------------------------------------------------

cat("\nSummary Statistics:\n")

cat("\nBeta Scores:\n")
summary(beta_data$betaScore)

cat("\niHS p-values (top 10):\n")
ihs_data %>%
  arrange(p_value) %>%
  head(10) %>%
  print()

cat("\nGenes with significant iHS (labeled in bold):\n")
cat(paste(genes_with_significant_ihs, collapse = ", "), "\n")

cat("\nAnalysis complete!\n")

# -----------------------------------------------------------------------------
# Session Info (for reproducibility)
# -----------------------------------------------------------------------------

# Uncomment to print session information
# sessionInfo()

# =============================================================================
# Data Format Requirements
# =============================================================================
#
# Beta Score Data (betascore_file):
#   - Chr: Chromosome number (1-14)
#   - Start: SNP position (base pairs)
#   - betaScore: Beta score value
#
# iHS Data (ihs_file):
#   - Chr: Chromosome number (1-14)
#   - Start: SNP position (base pairs)
#   - p_value: p-value for iHS test (or iHS score that will be converted)
#
# Gene Annotations (gene_annotations_file):
#   - Chr: Chromosome number
#   - Start: Gene start position
#   - End: Gene end position
#   - Gene: Gene name/symbol
#
# =============================================================================
