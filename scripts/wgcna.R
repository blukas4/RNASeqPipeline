library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads() # Skip this line if using R Studio

# Get arguments
args <- commandArgs(trailingOnly = TRUE)

counts_file <- args[1]
output_file <- args[2]

# Read in normalized count data
dat_expr <- t(read.csv(counts_file, row.names = 1))

# ====================
# Choosing the soft-thresholding power: analysis of network topology
# ====================
if (FALSE) {
  # Traditional approach of choosing the soft-thresholding power
  # Based on human judgement of visual inspection of these plots
  powers <- c(1:10, seq(from = 12, to = 20, by = 2))
  sft <- pickSoftThreshold(dat_expr, powerVector = powers, verbose = 5)

  sizeGrWindow(9, 5)
  par(mfrow = c(1, 2))
  cex1 <- 0.9
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit,signed R^2",
    type = "n", main = "Scale independence"
  )
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers, cex = cex1, col = "red"
  )
  abline(h = 0.90, col = "red")
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
    type = "n", main = "Mean connectivity"
  )
  text(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    labels = powers,
    cex = cex1,
    col = "red"
  )
}

# Automated function to select soft thresholding power
select_soft_power <- function(data, fit_threshold = 0.8) {
  # Calculate the correlation matrix
  powers <- c(1:10, seq(from = 12, to = 20, by = 2))
  sft <- pickSoftThreshold(data, powerVector = powers, verbose = 5)

  # Find the optimal power
  optimal_idx <- which(
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2] > fit_threshold
  )[1]

  # Check if any power meets the threshold
  if (is.na(optimal_idx)) {
    warning(
      "No power reached the desired scale-free fit threshold.",
      " Consider choosing the power manually based on the fit index plot."
    )
    return(NULL)
  } else {
    return(powers[optimal_idx])
  }
}

# ====================
# Co-expression similarity and adjacency
# ====================
soft_power <- select_soft_power(dat_expr)
print(paste("Selected soft-thresholding power:", soft_power))
adjacency <- adjacency(dat_expr, power = soft_power)

# ====================
# Topological Overlap Matrix (TOM)
# ====================
tom <- TOMsimilarity(adjacency)
diss_tom <- 1 - tom

# ====================
# Clustering using TOM
# ====================
gene_tree <- hclust(as.dist(diss_tom), method = "average")
if (FALSE) {
  # Visualize gene tree
  sizeGrWindow(12, 9)
  plot(gene_tree,
    xlab = "", sub = "",
    main = "Gene clustering on TOM-based dissimilarity",
    labels = FALSE, hang = 0.04
  )
}

min_module_size <- 30
dynamic_mods <- cutreeDynamic(
  dendro = gene_tree, distM = diss_tom,
  deepSplit = 2, pamRespectsDendro = FALSE,
  minClusterSize = min_module_size
)

table(dynamic_mods)
dynamic_colors <- labels2colors(dynamic_mods)
table(dynamic_colors)

if (FALSE) {
  # Visualize gene tree with colors
  sizeGrWindow(8, 6)
  plotDendroAndColors(gene_tree, dynamic_colors, "Dynamic Tree Cut",
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05,
    main = "Gene dendrogram and module colors"
  )
}

# ====================
# Merging of modules whose expression profiles are very similar
# ====================
me_list <- moduleEigengenes(dat_expr, colors = dynamic_colors)
mes <- me_list$eigengenes
me_diss <- 1 - cor(mes)
me_tree <- hclust(as.dist(me_diss), method = "average")
me_diss_thres <- 0.25

if (FALSE) {
  # Visualize clustering of module eigengenes
  sizeGrWindow(7, 6)
  plot(me_tree, main = "Clustering of module eigengenes", xlab = "", sub = "")

  abline(h = me_diss_thres, col = "red")
}

merge_data <- mergeCloseModules(dat_expr, dynamic_colors,
  cutHeight = me_diss_thres, verbose = 3
)
merged_colors <- merge_data$colors
merged_mes <- merge_data$newMEs

if (FALSE) {
  # Visualize gene tree with original and merged colors
  sizeGrWindow(12, 9)
  plotDendroAndColors(gene_tree, cbind(dynamic_colors, merged_colors),
    c("Dynamic Tree Cut", "Merged dynamic"),
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05
  )
}

module_colors <- merged_colors
color_order <- c("grey", standardColors(50))
module_labels <- match(module_colors, color_order) - 1
mes <- merged_mes
table(merged_colors)

# ====================
# Create a dataframe mapping genes to their merged_colors
# ====================
# Set Gene names as row names and create the dataframe
gene_to_color_df <- data.frame(
  Dynamic_Color = dynamic_colors,
  Merged_Color = merged_colors,
  row.names = colnames(dat_expr)
)

# Write dataframe to CSV
write.csv(gene_to_color_df, file = output_file)
