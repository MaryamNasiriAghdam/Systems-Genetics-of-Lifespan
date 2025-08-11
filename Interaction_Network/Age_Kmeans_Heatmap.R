
# Clear workspace and set library paths
rm(list = ls())
.libPaths(c("/data/mackanholt_lab/maryamn/r_packages", .libPaths()))

# Load required libraries
library(ComplexHeatmap)
library(circlize)

# Set working directory for input files
setwd("/data/mackanholt_lab/maryamn/Lifespan/SAS_April/Anova_Omni_GLM/")

# Read and inspect the gene expression data
data_matrix <- read.csv("Lsmeans_Age_Sig.csv", row.names = 1)
head(data_matrix)

# Scale the data by rows (genes)
scaled_matrix <- t(apply(data_matrix, 1, scale))
colnames(scaled_matrix) <- colnames(data_matrix)
head(scaled_matrix)


# Determine the optimal number of clusters using the WSS (Elbow Method)
max_clusters <- 10  # Set the maximum number of clusters to test
wss <- numeric(max_clusters)

# Calculate total within-cluster sum of squares for each k
for (k in 1:max_clusters) {
  set.seed(123)  # Ensure reproducibility
  km_result <- kmeans(scaled_matrix, centers = k, nstart = 25)
  wss[k] <- km_result$tot.withinss
}

# Plot the WSS values to visualize the "elbow"
plot(1:max_clusters, wss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters",
     ylab = "Total Within-Cluster Sum of Squares",
     main = "Elbow Method for Determining Optimal Clusters")


# Perform K-means clustering on the scaled data
num_clusters <- 6  # Adjust the number of clusters as needed
kmeans_result <- kmeans(scaled_matrix, centers = num_clusters)
cluster_assignments <- kmeans_result$cluster

# Set seed for reproducibility
set.seed(1)


# Define a custom color palette
my_palette <- colorRamp2(c(-2, 0, 2), c("#522D80", "#FDE8DB", "#F66733"))


# Create a heatmap with k-means clustering and annotations
heatmap_object <- Heatmap(
  scaled_matrix,
  name = "Expression",
  col = my_palette,
  row_names_gp = gpar(fontsize = 8),
  column_title = "OMNIBUS_SIG_AGE",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Z-scaled Normalized Expression Level",
    at = seq(-2, 2, length.out = 5),
    labels = c("Low", "", "", "", "High")
  ),
  show_row_names = FALSE,
  row_km = num_clusters,
  use_raster = FALSE,
  clustering_distance_rows = "euclidean",
  cluster_rows = TRUE,
  row_dend_side = "left",
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 8),
  row_gap = unit(1, "mm"),  # Adds separation between clusters
  border = TRUE,  # Adds black border around clusters
  column_order = order(as.numeric(gsub("column", "", colnames(data_matrix))))
)

# Display the heatmap
print(heatmap_object)

# Change working directory for output files
setwd("/data/mackanholt_lab/maryamn/Lifespan/Kmeans_Paper_Feb2025/")

# extract clustering details
set.seed(1)
drawn_heatmap <- draw(heatmap_object)
row_dendrogram <- row_dend(drawn_heatmap)  # Extract row dendrogram
cluster_list <- row_order(drawn_heatmap)    # Extract ordered clusters (as a list)
lapply(cluster_list, length)  # Optionally print the size of each cluster

# Create a data frame with gene cluster assignments
gene_clusters <- NULL
for (i in seq_along(cluster_list)) {
  genes <- rownames(scaled_matrix)[cluster_list[[i]]]
  cluster_data <- cbind(genes, paste("cluster", i, sep = ""))
  gene_clusters <- if (is.null(gene_clusters)) {
    cluster_data
  } else {
    rbind(gene_clusters, cluster_data)
  }
}

gene_clusters <- as.data.frame(gene_clusters, stringsAsFactors = FALSE)
colnames(gene_clusters) <- c("GeneID", "Cluster")

# Write the gene cluster assignments to a file
write.table(gene_clusters, file = "gene_clusters.txt", sep = "\t", quote = FALSE, row.names = FALSE)
