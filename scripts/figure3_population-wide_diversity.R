# =============================================================================
# Population Genomic Diversity Analysis: MDS and IBS Heatmap
# =============================================================================
# 
# Description:
#   Analyzes population-wide genomic diversity of P. falciparum infections
#   from The Gambia across two time periods (1966-71 vs 2015) using:
#   - Multidimensional Scaling (MDS) of SNP profiles
#   - Identity-by-State (IBS) heatmap with hierarchical clustering
#
# Methods:
#   - MDS based on genetic distances from SNP profiles (predominant SNP per locus)
#   - SNPs filtered for linkage disequilibrium (LD-pruned)
#   - Neighbor-Joining clustering of IBS matrix
#
# Input:
#   - VCF file: Filtered SNP data (LD-pruned, MAF filtered)
#   - Sample metadata: Sample IDs, Year/Period, Population info
#
# Output:
#   - Panel A: MDS scatter plot (Axis 1 vs Axis 2)
#   - Panel B: IBS heatmap with dendrogram clustering
#   - Combined two-panel figure
#
# Author: [Your Name]
# Date: [Date]
# =============================================================================

# Load required libraries
library(tidyverse)    # Data manipulation and ggplot2
library(vcfR)         # VCF file handling
library(adegenet)     # Population genetics analysis
library(ape)          # Phylogenetics (Neighbor-Joining)
library(SNPRelate)    # SNP-based analyses and IBS
library(pheatmap)     # Heatmap visualization
library(ggpubr)       # Publication-ready plots
library(patchwork)    # Combine plots

# -----------------------------------------------------------------------------
# Configuration and Parameters
# -----------------------------------------------------------------------------

# File paths
vcf_input <- "data/gampf1966_72vs2015_ldpruned_maf03.recode.vcf"
metadata_file <- "data/sample_metadata.txt"
gds_output <- "data/gampf1966_72vs2015_ldpruned.gds"

# Output files
output_figure <- "results/figures/figure_population_diversity.pdf"
output_mds_data <- "results/mds_coordinates.csv"
output_ibs_matrix <- "results/ibs_matrix.csv"

# Sample information
n_samples_1966_71 <- 54
n_samples_2015 <- 89

# Color scheme
color_1966_71 <- "#E15759"  # Red
color_2015 <- "#4E79A7"     # Blue
period_colors <- c("1966-71" = color_1966_71, "2015" = color_2015)

# MDS parameters
n_mds_axes <- 3  # Number of MDS axes to calculate
mds_point_size <- 3
mds_point_alpha <- 0.7

# IBS heatmap parameters
heatmap_colors <- colorRampPalette(c("yellow", "orange", "red"))(256)
dendrogram_height <- 20
clustering_method <- "average"  # For Neighbor-Joining-like clustering

# Analysis parameters
n_threads <- 2

# Plot dimensions
plot_width_mm <- 180
plot_height_mm <- 90
plot_dpi <- 600

# -----------------------------------------------------------------------------
# Create Output Directories
# -----------------------------------------------------------------------------

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Data Import
# -----------------------------------------------------------------------------

cat("Loading data...\n")

# Read VCF file
cat("  Reading VCF file...\n")
vcf_data <- read.vcfR(vcf_input)

# Load sample metadata
# Expected columns: SampleID, Year, Period (or similar)
cat("  Reading sample metadata...\n")
sample_metadata <- read_delim(metadata_file, show_col_types = FALSE)

# Ensure Period is a factor with correct order
sample_metadata <- sample_metadata %>%
   mutate(Period = factor(Period, levels = c("1966-71", "2015")))

# Display sample summary
cat("\nSample Summary:\n")
sample_metadata %>%
   group_by(Period) %>%
   summarise(n = n(), .groups = "drop") %>%
   print()

# -----------------------------------------------------------------------------
# Convert VCF to Genind Object
# -----------------------------------------------------------------------------

cat("\nConverting VCF to genind object...\n")

# Convert to genind for adegenet analysis
genind_data <- vcfR2genind(vcf_data)

# Assign population information
pop(genind_data) <- sample_metadata$Period

# Display genind summary
cat("\nGenind Object Summary:\n")
print(genind_data)

cat(sprintf("  Total SNPs analyzed: %d\n", nLoc(genind_data)))
cat(sprintf("  Total samples: %d\n", nInd(genind_data)))

# -----------------------------------------------------------------------------
# Multidimensional Scaling (MDS) Analysis
# -----------------------------------------------------------------------------

cat("\nPerforming MDS analysis...\n")

# Scale genind data (center but don't scale to preserve genetic distance)
genind_scaled <- scaleGen(
   x = genind_data,
   center = TRUE,
   scale = FALSE,
   truenames = TRUE
)

# Calculate genetic distances (pairwise distances)
cat("  Calculating genetic distances...\n")
genetic_distances <- dist.gene(
   x = genind_scaled,
   method = "pairwise"
)

# Perform MDS (Principal Coordinate Analysis)
cat("  Performing PCoA/MDS...\n")
mds_result <- cmdscale(
   d = genetic_distances,
   k = n_mds_axes,
   eig = TRUE
)

# Extract MDS coordinates
mds_coords <- as.data.frame(mds_result$points)
colnames(mds_coords) <- paste0("MDS", 1:n_mds_axes)

# Add sample information
mds_coords$SampleID <- rownames(mds_coords)
mds_coords <- mds_coords %>%
   left_join(sample_metadata, by = "SampleID")

# Calculate variance explained by each axis
variance_explained <- mds_result$eig / sum(mds_result$eig) * 100

cat(sprintf("\nVariance explained by MDS axes:\n"))
cat(sprintf("  MDS1: %.2f%%\n", variance_explained[1]))
cat(sprintf("  MDS2: %.2f%%\n", variance_explained[2]))
if (n_mds_axes >= 3) {
   cat(sprintf("  MDS3: %.2f%%\n", variance_explained[3]))
}

# Save MDS coordinates
write_csv(mds_coords, output_mds_data)
cat(sprintf("\nMDS coordinates saved to: %s\n", output_mds_data))

# -----------------------------------------------------------------------------
# Panel A: MDS Scatter Plot
# -----------------------------------------------------------------------------

cat("\nGenerating Panel A (MDS plot)...\n")

# Create axis labels with variance explained
x_label <- sprintf("MDS Axis 1 (%.1f%%)", variance_explained[1])
y_label <- sprintf("MDS Axis 2 (%.1f%%)", variance_explained[2])

panel_a <- ggplot(mds_coords, aes(x = MDS1, y = MDS2, color = Period)) +
   
   # Data points
   geom_point(
      size = mds_point_size,
      alpha = mds_point_alpha
   ) +
   
   # Color scheme
   scale_color_manual(
      values = period_colors,
      name = "Sampling Period",
      labels = c(
         sprintf("1966-71 (N=%d)", n_samples_1966_71),
         sprintf("2015 (N=%d)", n_samples_2015)
      )
   ) +
   
   # Labels
   labs(
      x = x_label,
      y = y_label
   ) +
   
   # Theme
   theme_classic() +
   theme(
      legend.position = "top",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      axis.text = element_text(size = 8, color = "#000000"),
      axis.title = element_text(size = 9, face = "bold"),
      axis.line = element_line(color = "#000000"),
      panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.5)
   )

# Display plot
print(panel_a)

# -----------------------------------------------------------------------------
# Convert VCF to GDS and Calculate IBS Matrix
# -----------------------------------------------------------------------------

cat("\nCalculating Identity-by-State (IBS) matrix...\n")

# Close any existing GDS connection
if (exists("gds_connection")) {
   tryCatch(
      snpgdsClose(gds_connection),
      error = function(e) NULL
   )
}

# Convert VCF to GDS format
cat("  Converting VCF to GDS format...\n")
snpgdsVCF2GDS(
   vcf.fn = vcf_input,
   out.fn = gds_output,
   method = "biallelic.only"
)

# Open GDS file
gds_connection <- snpgdsOpen(gds_output)

# Calculate IBS matrix
cat("  Computing IBS matrix...\n")
ibs_result <- snpgdsIBS(
   gdsobj = gds_connection,
   num.thread = n_threads,
   autosome.only = FALSE
)

# Extract IBS matrix
ibs_matrix <- ibs_result$ibs

# Add sample names to matrix
rownames(ibs_matrix) <- sample_metadata$SampleID
colnames(ibs_matrix) <- sample_metadata$SampleID

# Save IBS matrix
write.csv(ibs_matrix, output_ibs_matrix)
cat(sprintf("IBS matrix saved to: %s\n", output_ibs_matrix))

# Close GDS connection
snpgdsClose(gds_connection)

# -----------------------------------------------------------------------------
# Hierarchical Clustering Using Neighbor-Joining Approach
# -----------------------------------------------------------------------------

cat("\nPerforming hierarchical clustering...\n")

# Convert IBS to distance (distance = 1 - similarity)
ibs_distance <- as.dist(1 - ibs_matrix)

# Perform hierarchical clustering (using average linkage for NJ-like results)
# Note: For true Neighbor-Joining, use ape::nj()
# Here we use hierarchical clustering for pheatmap compatibility

# Option 1: Hierarchical clustering (compatible with pheatmap)
hclust_result <- hclust(ibs_distance, method = clustering_method)

# Option 2: Neighbor-Joining tree (convert to hclust for pheatmap)
# Uncomment to use true NJ clustering
# nj_tree <- nj(ibs_distance)
# hclust_result <- as.hclust(nj_tree)

# -----------------------------------------------------------------------------
# Prepare Annotations for Heatmap
# -----------------------------------------------------------------------------

cat("Preparing heatmap annotations...\n")

# Create annotation dataframe
annotation_df <- sample_metadata %>%
   select(SampleID, Period) %>%
   column_to_rownames("SampleID")

# Define annotation colors
annotation_colors <- list(
   Period = period_colors
)

# -----------------------------------------------------------------------------
# Panel B: IBS Heatmap with Clustering
# -----------------------------------------------------------------------------

cat("\nGenerating Panel B (IBS heatmap)...\n")

# Create heatmap
panel_b <- pheatmap(
   mat = ibs_matrix,
   
   # Clustering
   cluster_rows = hclust_result,
   cluster_cols = hclust_result,
   clustering_distance_rows = ibs_distance,
   clustering_distance_cols = ibs_distance,
   
   # Visual styling
   color = heatmap_colors,
   show_rownames = FALSE,
   show_colnames = FALSE,
   
   # Dendrogram
   treeheight_row = dendrogram_height,
   treeheight_col = dendrogram_height,
   
   # Annotations
   annotation_row = annotation_df,
   annotation_col = annotation_df,
   annotation_colors = annotation_colors,
   annotation_names_row = FALSE,
   annotation_names_col = FALSE,
   annotation_legend = TRUE,
   
   # Legend
   legend = TRUE,
   legend_labels = c("Low", "High"),
   
   # Sizing
   fontsize = 8,
   
   # Don't display immediately (we'll save directly)
   silent = TRUE
)

# -----------------------------------------------------------------------------
# Combine Panels and Save
# -----------------------------------------------------------------------------

cat("\nCombining panels and saving figure...\n")

# Save the combined figure
pdf(output_figure, width = plot_width_mm / 25.4, height = plot_height_mm / 25.4)

# Create layout: Panel A on left, Panel B on right
layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1.2))

# Panel A
par(mar = c(4, 4, 2, 1))
print(panel_a)

# Panel B
par(mar = c(4, 1, 2, 4))
print(panel_b[[4]])  # Extract the plot from pheatmap object

dev.off()

cat(sprintf("Combined figure saved to: %s\n", output_figure))

# Alternative: Save panels separately
cat("\nSaving individual panels...\n")

# Save Panel A
ggsave(
   filename = "results/figures/panel_a_mds.pdf",
   plot = panel_a,
   width = 90,
   height = 90,
   units = "mm",
   dpi = plot_dpi
)

# Save Panel B
pdf("results/figures/panel_b_ibs_heatmap.pdf", 
    width = 90 / 25.4, 
    height = 90 / 25.4)
print(panel_b)
dev.off()

# -----------------------------------------------------------------------------
# Summary Statistics
# -----------------------------------------------------------------------------

cat("\n=== Analysis Summary ===\n")

# MDS summary
cat("\nMDS Analysis:\n")
cat(sprintf("  Total SNPs: %d\n", nLoc(genind_data)))
cat(sprintf("  Total samples: %d\n", nInd(genind_data)))
cat(sprintf("  Cumulative variance (MDS1-2): %.2f%%\n", 
            sum(variance_explained[1:2])))

# IBS summary
cat("\nIBS Analysis:\n")
cat(sprintf("  Mean IBS: %.4f\n", mean(ibs_matrix[lower.tri(ibs_matrix)])))
cat(sprintf("  Median IBS: %.4f\n", median(ibs_matrix[lower.tri(ibs_matrix)])))
cat(sprintf("  Range: %.4f - %.4f\n", 
            min(ibs_matrix[lower.tri(ibs_matrix)]),
            max(ibs_matrix[lower.tri(ibs_matrix)])))

# Identify highly related pairs (IBS > 0.95)
high_relatedness_threshold <- 0.95
highly_related <- which(
   ibs_matrix > high_relatedness_threshold & 
      ibs_matrix < 1,  # Exclude self-comparisons
   arr.ind = TRUE
)

if (nrow(highly_related) > 0) {
   cat(sprintf("\nHighly related pairs (IBS > %.2f): %d pairs\n", 
               high_relatedness_threshold, 
               nrow(highly_related) / 2))  # Divide by 2 as matrix is symmetric
}

# Identify 2015 cluster with high MDS2 values
mds2_cluster <- mds_coords %>%
   filter(Period == "2015") %>%
   arrange(desc(MDS2)) %>%
   head(11)

cat("\nTop 11 samples from 2015 with high MDS2 values:\n")
print(mds2_cluster %>% select(SampleID, MDS1, MDS2))

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
# VCF File (vcf_input):
#   - LD-pruned SNPs (no significant linkage disequilibrium)
#   - MAF filtered (e.g., MAF > 0.03)
#   - Biallelic SNPs only
#   - Missing data handled appropriately
#
# Sample Metadata (metadata_file):
#   Columns:
#   - SampleID: Unique sample identifier (must match VCF)
#   - Year: Sampling year (numeric or character)
#   - Period: Sampling period ("1966-71" or "2015")
#   - Optional: Additional metadata (location, etc.)
#
# Expected format (tab or comma delimited):
#   SampleID  Year  Period
#   Sample1   1966  1966-71
#   Sample2   1967  1966-71
#   Sample3   2015  2015
#   ...
#
# =============================================================================