# =============================================================================
# Population Genetic Structure Analysis: PCoA and IBS Heatmap
# =============================================================================
#
# Description:
#   Performs Principal Coordinate Analysis (PCoA) on genetic data from
#   Gambia P. falciparum samples (1966-72 vs 2015) and creates an Identity
#   By State (IBS) heatmap with hierarchical clustering.
#
# Workflow:
#   1. Load VCF file and convert to genind object
#   2. Perform PCoA on genetic distances
#   3. Test different clustering methods
#   4. Calculate IBS matrix using SNPRelate
#   5. Generate clustered heatmap
#
# Input:
#   - VCF file: gampf1966_72vs2015BK_maf03nomiss.recode.vcf
#   - Population metadata: gampf1966_71pops (must be pre-loaded)
#   - Sample year annotation: gamPf6671_2015_SampleYear (must be pre-loaded)
#   - Annotation colors: annotation_colors (must be pre-loaded)
#
# Output:
#   - GDS file: gampf1966_72vs2015BK_maf03nomiss.recode.gds
#   - PCoA plots (interactive)
#   - IBS heatmap (interactive)
#
# Author: Mouhamadou Fadel DIOP
# Date: [Date]
# =============================================================================

# Load required libraries
library(adegenet)   # Genetic data analysis
library(ape)        # Phylogenetics and evolution
library(vcfR)       # VCF file handling
library(SNPRelate)  # SNP-based relationship analysis
library(pheatmap)   # Heatmap visualization

# -----------------------------------------------------------------------------
# Configuration and Parameters
# -----------------------------------------------------------------------------

# File paths
vcf_input <- "gampf1966_72vs2015BK_maf03nomiss.recode.vcf"
gds_output <- "gampf1966_72vs2015BK_maf03nomiss.recode.gds"

# Analysis parameters
n_pcoa_axes <- 3              # Number of PCoA axes to retain
n_threads <- 2                # Number of threads for parallel processing
heatmap_dendrogram_height <- 15  # Height of dendrograms in heatmap

# Heatmap color scheme
heatmap_colors <- colorRampPalette(c("yellow", "orange", "red"))(256)

# -----------------------------------------------------------------------------
# Data Import and Conversion
# -----------------------------------------------------------------------------

# Read VCF file
cat("Reading VCF file...\n")
vcf_data <- read.vcfR(vcf_input)

# Convert VCF to genind object for population genetics analysis
cat("Converting to genind object...\n")
genind_data <- vcfR2genind(vcf_data)

# Assign population information
# NOTE: 'gampf1966_71pops' must be loaded before running this script
# Expected structure: data frame with 'pop' column matching sample order
if (!exists("gampf1966_71pops")) {
  stop("Error: 'gampf1966_71pops' object not found. Please load population metadata first.")
}
pop(genind_data) <- gampf1966_71pops$pop

# Display genind object summary
print(genind_data)

# -----------------------------------------------------------------------------
# Principal Coordinate Analysis (PCoA)
# -----------------------------------------------------------------------------

cat("Performing PCoA...\n")

# Calculate genetic distances and perform PCoA
# Steps:
#   1. Scale genetic data (center but don't scale)
#   2. Calculate pairwise genetic distances
#   3. Perform PCoA on distance matrix
pcoa_result <- dudi.pco(
  d = dist.gene(
    x = scaleGen(
      x = genind_data,
      center = TRUE,
      scale = FALSE,
      truenames = TRUE
    ),
    method = "pairwise"
  ),
  scannf = FALSE,  # Don't display scree plot
  nf = n_pcoa_axes  # Number of axes to retain
)

# Visualize PCoA results
cat("Generating PCoA plots...\n")

# Plot first two principal coordinates
plot(
  pcoa_result$li[, 1:2],
  main = "PCoA: PC1 vs PC2",
  xlab = "PC1",
  ylab = "PC2"
)

# Plot first three principal coordinates (if available)
if (n_pcoa_axes >= 3) {
  plot(
    pcoa_result$li[, 1:3],
    main = "PCoA: First 3 Principal Coordinates"
  )
}

# -----------------------------------------------------------------------------
# Hierarchical Clustering Exploration
# -----------------------------------------------------------------------------

cat("Testing hierarchical clustering methods...\n")

# Calculate distance matrix from PCoA coordinates
pcoa_distance <- dist(pcoa_result$li[, 1:n_pcoa_axes])

# Test different linkage methods
# Note: The last method assigned will be used for the heatmap
clustering_methods <- c("ward.D2", "complete", "average", "single")

for (method in clustering_methods) {
  cat("  Testing method:", method, "\n")
  pcoa_clustering <- hclust(pcoa_distance, method = method)
  # Note: Uncomment to visualize each clustering
  # plot(pcoa_clustering, main = paste("Hierarchical Clustering -", method))
}

# Use the final clustering method for subsequent analysis
# (Currently set to 'single' based on original code)
pcoa_clustering_final <- hclust(pcoa_distance, method = "single")

# -----------------------------------------------------------------------------
# Identity By State (IBS) Analysis using SNPRelate
# -----------------------------------------------------------------------------

cat("Converting VCF to GDS format...\n")

# Close any existing GDS file connection (if open)
if (exists("gds_connection")) {
  tryCatch(
    snpgdsClose(gds_connection),
    error = function(e) NULL
  )
}

# Convert VCF to GDS format for efficient SNP analysis
snpgdsVCF2GDS(
  vcf.fn = vcf_input,
  out.fn = gds_output,
  method = "biallelic.only"
)

# Open GDS file
cat("Opening GDS file and calculating IBS...\n")
gds_connection <- snpgdsOpen(gds_output)

# Calculate Identity By State (IBS) matrix
# IBS measures genetic similarity between samples
ibs_result <- snpgdsIBS(
  gdsobj = gds_connection,
  num.thread = n_threads,
  autosome.only = FALSE  # Include all chromosomes
)

# Convert IBS similarity to distance matrix
# Distance = 1 - similarity
ibs_distance <- dist(1 - ibs_result$ibs)

# -----------------------------------------------------------------------------
# Generate IBS Heatmap with Clustering
# -----------------------------------------------------------------------------

cat("Generating IBS heatmap...\n")

# Check required annotation objects
if (!exists("gamPf6671_2015_SampleYear")) {
  warning("'gamPf6671_2015_SampleYear' not found. Heatmap will be generated without sample year annotations.")
  gamPf6671_2015_SampleYear <- NULL
}

if (!exists("annotation_colors")) {
  warning("'annotation_colors' not found. Using default colors.")
  annotation_colors <- NULL
}

# Create heatmap
pheatmap(
  mat = ibs_result$ibs,

  # Clustering configuration
  cluster_rows = pcoa_clustering_final,
  cluster_cols = pcoa_clustering_final,
  clustering_distance_rows = ibs_distance,
  clustering_distance_cols = ibs_distance,

  # Visual styling
  color = heatmap_colors,
  show_rownames = FALSE,  # Hide row names to reduce clutter
  show_colnames = FALSE,  # Hide column names to reduce clutter
  main = "Identity By State (IBS) Heatmap",

  # Dendrogram height
  treeheight_row = heatmap_dendrogram_height,
  treeheight_col = heatmap_dendrogram_height,

  # Sample annotations
  annotation_row = gamPf6671_2015_SampleYear,
  annotation_col = gamPf6671_2015_SampleYear,
  annotation_colors = annotation_colors,
  annotation_names_row = FALSE,
  annotation_names_col = FALSE
)

# -----------------------------------------------------------------------------
# Cleanup
# -----------------------------------------------------------------------------

# Close GDS file connection
snpgdsClose(gds_connection)

cat("Analysis complete!\n")

# -----------------------------------------------------------------------------
# Session Info (for reproducibility)
# -----------------------------------------------------------------------------

# Uncomment to print session information
# sessionInfo()

# =============================================================================
# Expected Data Objects (must be loaded before running this script)
# =============================================================================
#
# gampf1966_71pops:
#   - Data frame with population assignments
#   - Must have 'pop' column
#   - Order must match samples in VCF file
#
# gamPf6671_2015_SampleYear:
#   - Data frame with sample year annotations
#   - Used for heatmap row/column annotations
#
# annotation_colors:
#   - Named list of colors for annotations
#   - Format: list(ColumnName = c(level1 = "color1", level2 = "color2", ...))
#
# =============================================================================
