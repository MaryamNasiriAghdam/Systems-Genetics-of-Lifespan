
# Load required libraries
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(tidyr)

# Read input datasets

FB_Lsmeans <- read.csv("FB_Lsmeans.csv", header = TRUE)
FH_Lsmeans <- read.csv("FH_Lsmeans.csv", header = TRUE)
MB_Lsmeans <- read.csv("MB_Lsmeans.csv", header = TRUE)
MH_Lsmeans <- read.csv("MH_Lsmeans.csv", header = TRUE)

# Organize data into a list
group_data <- list(FB = FB_Lsmeans, FH = FH_Lsmeans, MB = MB_Lsmeans, MH = MH_Lsmeans)

gene_symbols <- FB_Lsmeans$Symbol

# Define comparisons
comparisons <- list(
  c("FB", "FH"),
  c("FH", "MH"),
  c("MH", "MB"),
  c("FB", "MB"),
  c("FB", "MH"),
  c("FH", "MB")
)

# Initialize results dataframe
cor_results <- data.frame(Symbol = gene_symbols)

# Compute Spearman correlation for each gene across comparisons
for (comp in comparisons) {
  group1 <- comp[1]
  group2 <- comp[2]
  cor_values <- sapply(gene_symbols, function(gene) {
    vec1 <- as.numeric(group_data[[group1]][group_data[[group1]]$Symbol == gene, -1])
    vec2 <- as.numeric(group_data[[group2]][group_data[[group2]]$Symbol == gene, -1])
    if (length(vec1) > 0 && length(vec2) > 0) {
      cor(vec1, vec2, method = "spearman", use = "complete.obs")
    } else {
      NA
    }
  })
  cor_results[[paste(group1, group2, sep = "_vs_")]] <- cor_values
}

# Save results
write.csv(cor_results, "Spearman_Correlation.csv", row.names = FALSE)

# Display first few rows
head(cor_results)

# Prepare heatmap
cor_matrix <- cor_results[,-1]
row.names(cor_matrix) <- cor_results$Symbol

gene_matrix <- as.matrix(cor_matrix)


# Function to calculate WSS (Within-Cluster Sum of Squares)
wss <- function(k) {
  kmeans(gene_matrix, centers = k, nstart = 25)$tot.withinss
}

# Set range of k values to test
k_values <- 1:10  # Adjust if needed

# Compute WSS for different k values
wss_values <- sapply(k_values, wss)

# Create the Elbow Plot
elbow_plot <- data.frame(K = k_values, WSS = wss_values)


library(ggplot2)
ggplot(elbow_plot, aes(x = K, y = WSS)) +
  geom_line() +
  geom_point(size = 2) +
  ggtitle("Elbow Method for Optimal Clusters") +
  xlab("Number of Clusters (K)") +
  ylab("Total Within-Cluster Sum of Squares") +
  theme_minimal()


# Perform K-means clustering
set.seed(123)
kmeans_result <- kmeans(gene_matrix, centers = 6)

# Define color function
col_fun <- colorRamp2(c(min(gene_matrix, na.rm = TRUE), 0, max(gene_matrix, na.rm = TRUE)),
                      c("#F5BD1E", "white", "#5F9C3F"))


Heatmap(gene_matrix, 
        name = "Spearman Correlation", 
        col = col_fun, 
        show_row_names = FALSE,  # Hide gene names
        show_column_names = TRUE, 
        column_names_rot = 45,  # Rotate column names by 45 degrees
        column_names_gp = gpar(fontsize = 14, fontfamily = "Arial"),  # Set column font to Arial 14
        row_names_gp = gpar(fontsize = 14, fontfamily = "Arial"),  # Ensure row fonts match (even if hidden)
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 14, fontfamily = "Arial"),  # Legend title font
          labels_gp = gpar(fontsize = 14, fontfamily = "Arial")  # Legend labels font
        ),
        row_split = kmeans_result$cluster,  # Split rows by clusters
        cluster_rows = FALSE,  # Use k-means clusters instead of hierarchical clustering
        cluster_columns = FALSE,  # Remove column clustering to avoid the dashed line
        show_column_dend = FALSE,  # Remove the column dendrogram numbers
        row_gap = unit(1, "mm"),  # Adds separation between clusters
        column_gap = unit(1, "mm"),  # Adds space between columns
        border = TRUE,  # Adds black border around clusters
        column_split = factor(1:ncol(gene_matrix)))  # Ensure column separation with outlines

## Save the genes in each cluster in separate files

for (i in 1:6) {
  cluster_genes <- gene_symbols[kmeans_result$cluster == i]
  write.csv(cluster_genes, paste0("CorrelationCluster_", i, ".csv"), row.names = FALSE)
}


##### Heatmap for each cluster (Repeated for each of 1-6 correlation clusters)

# Read data
FB_Lsmeans <- read.csv("FB_Lsmeans.csv", header = TRUE)
FH_Lsmeans <- read.csv("FH_Lsmeans.csv", header = TRUE)
MB_Lsmeans <- read.csv("MB_Lsmeans.csv", header = TRUE)
MH_Lsmeans <- read.csv("MH_Lsmeans.csv", header = TRUE)
Cluster1 <- read.csv("CorrelationCluster_1.csv", header = TRUE)

# Filter data for Cluster1 genes
FB_cluster <- FB_Lsmeans %>% filter(Symbol %in% Cluster1$x)
FH_cluster <- FH_Lsmeans %>% filter(Symbol %in% Cluster1$x)
MB_cluster <- MB_Lsmeans %>% filter(Symbol %in% Cluster1$x)
MH_cluster <- MH_Lsmeans %>% filter(Symbol %in% Cluster1$x)

# Define custom colors using -2 and 2 as extremes
custom_colors <- colorRamp2(c(-2, 0, 2), c("#522D80", "white", "#F56600"))

# Function to create a heatmap
create_heatmap <- function(data, title) {
  rownames(data) <- data$Symbol  # Assign row names before removing the Symbol column
  data <- data[,-1]  # Remove the Symbol column
  data <- as.matrix(data)  # Ensure it's a matrix
  data_scaled <- t(apply(data, 1, scale))  # Scale data
  
  Heatmap(
    data_scaled,
    name = title,
    col = custom_colors,  # Use pre-defined custom colors
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 4),  # Reduce font size for row names
    column_title = title,
    row_title = NULL,
    border = TRUE  # Add a black outline
  )
}

# Generate heatmaps for FB, FH, MB, MH
FB_heatmap <- create_heatmap(FB_cluster, "FB")
FH_heatmap <- create_heatmap(FH_cluster, "FH")
MB_heatmap <- create_heatmap(MB_cluster, "MB")
MH_heatmap <- create_heatmap(MH_cluster, "MH")

# Combine the heatmaps
combined_heatmap <- FB_heatmap + FH_heatmap + MB_heatmap + MH_heatmap

# Draw the combined heatmap
draw(combined_heatmap)



####### Draw line plots for representative genes from each cluster ########

########### For Cluster 1: Select genes with positive values 0.8 >= 3 
# Load required data
Spearman <- read.csv("Spearman_Correlation.csv", header = TRUE)
Cluster1 <- read.csv("CorrelationCluster_1.csv", header = TRUE)

# Filter Spearman table to keep only genes present in Cluster
Cluster1_genes <- Cluster1$x
Filtered_Spearman <- Spearman %>% filter(Symbol %in% Cluster1_genes)

# Filter genes where at least three correlation values are greater than 0.8
filtered_genes <- Filtered_Spearman[apply(Filtered_Spearman[, -1], 1, function(row) sum(row > 0.8) >= 3), ]

# Load Lsmeans data
Lsmeans <- read.csv("Lsmeans_WOR.csv", header = TRUE)

# Filter Lsmeans for overlapping genes
Lsmeans_overlapping_genes <- Lsmeans %>%
  filter(Symbol %in% filtered_genes$Symbol)

# Reshape data to long format for easier plotting
Lsmeans_overlapping_long <- Lsmeans_overlapping_genes %>%
  select(FlyBase.ID, FB_Day5:MH_Day70) %>%
  pivot_longer(cols = -FlyBase.ID, 
               names_to = "Time_Point", 
               values_to = "Lsmeans") %>%
  mutate(
    Group = sub("_.*", "", Time_Point),
    Time_Point = sub(".*_", "", Time_Point),
    Time_Point = factor(Time_Point, levels = c("Day5", "Day7", "Day14", "Day21", "Day28", 
                                               "Day35", "Day42", "Day49", "Day56", 
                                               "Day63", "Day70"))
  )

# Define visualization properties
group_colors <- c("FB" = scales::alpha("#F56600", 0.5),
                  "FH" = scales::alpha("#522D80", 0.5),
                  "MB" = scales::alpha("#546223", 0.5),
                  "MH" = scales::alpha("#005EB8", 0.5))

group_shapes <- c("FB" = 16, "FH" = 15, "MB" = 17, "MH" = 18)

# Loop through all unique genes to create and display plots
for (gene in unique(Lsmeans_overlapping_long$FlyBase.ID)) {
  processed_data <- Lsmeans_overlapping_long %>%
    filter(FlyBase.ID == gene)
  
  gene_plot <- ggplot(processed_data, aes(x = Time_Point, y = Lsmeans, group = Group, color = Group, shape = Group)) +
    geom_line(linewidth = 1.2, alpha = 0.9) + 
    geom_point(size = 4, alpha = 0.7) +
    scale_color_manual(values = group_colors) +
    scale_shape_manual(values = group_shapes) +
    scale_x_discrete(labels = function(x) gsub("Day", "", x)) +
    labs(title = paste("Expression Profile -", gene),
         x = "Time (Days)", 
         y = "Normalized Expression Level", 
         color = "Groups", 
         shape = "Groups") +
    theme_minimal(base_size = 14) +
    theme(
      text = element_text(family = "Arial", color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 14),
      axis.text.y = element_text(margin = margin(r = 10), color = "black", size = 14),
      axis.title.x = element_text(face = "bold", color = "black", size = 14, family = "Arial"),
      axis.title.y = element_text(face = "bold", color = "black", size = 14, family = "Arial"),
      panel.background = element_rect(fill = "white", color = "gray", linewidth = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "black", family = "Arial")
    )  
  
  
  print(gene_plot)
  
  readline(prompt = "Press [Enter] to continue to the next plot or manually stop execution if done.")
}


########### For Cluster 2: Select genes with the highest negative correlation in any comparison
# Load required data
Spearman <- read.csv("Spearman_Correlation_Feb15.csv", header = TRUE)
Cluster2 <- read.csv("CorrelationCluster_2.csv", header = TRUE)


# Filter Spearman table to keep only genes present in Cluster
Cluster2_genes <- Cluster2$x
Filtered_Spearman <- Spearman %>% filter(Symbol %in% Cluster2_genes)

# Find genes with the highest negative correlation in any comparison
high_neg_corr_genes <- Filtered_Spearman %>%
  rowwise() %>%
  mutate(min_correlation = min(c_across(-Symbol))) %>%
  ungroup() %>%
  arrange(min_correlation) %>%
  slice_head(n = 50)


# Load Lsmeans data
Lsmeans <- read.csv("Lsmeans_WOR.csv", header = TRUE)

# Filter Lsmeans for overlapping genes
Lsmeans_overlapping_genes <- Lsmeans %>%
  filter(Symbol %in% high_neg_corr_genes$Symbol)


# Reshape data to long format for easier plotting
Lsmeans_overlapping_long <- Lsmeans_overlapping_genes %>%
  select(FlyBase.ID, FB_Day5:MH_Day70) %>%
  pivot_longer(cols = -FlyBase.ID, 
               names_to = "Time_Point", 
               values_to = "Lsmeans") %>%
  mutate(
    Group = sub("_.*", "", Time_Point),
    Time_Point = sub(".*_", "", Time_Point),
    Time_Point = factor(Time_Point, levels = c("Day5", "Day7", "Day14", "Day21", "Day28", 
                                               "Day35", "Day42", "Day49", "Day56", 
                                               "Day63", "Day70"))
  )

# Define visualization properties
group_colors <- c("FB" = scales::alpha("#F56600", 0.5),
                  "FH" = scales::alpha("#522D80", 0.5),
                  "MB" = scales::alpha("#546223", 0.5),
                  "MH" = scales::alpha("#005EB8", 0.5))

group_shapes <- c("FB" = 16, "FH" = 15, "MB" = 17, "MH" = 18)

# Loop through all unique genes to create and display plots
for (gene in unique(Lsmeans_overlapping_long$FlyBase.ID)) {
  processed_data <- Lsmeans_overlapping_long %>%
    filter(FlyBase.ID == gene)
  
  gene_plot <- ggplot(processed_data, aes(x = Time_Point, y = Lsmeans, group = Group, color = Group, shape = Group)) +
    geom_line(linewidth = 1.2, alpha = 0.9) + 
    geom_point(size = 4, alpha = 0.7) + 
    scale_color_manual(values = group_colors) + 
    scale_shape_manual(values = group_shapes) +
    scale_x_discrete(labels = function(x) gsub("Day", "", x)) + 
    labs(title = paste("Expression Profile -", gene),
         x = "Time (Days)", 
         y = "Normalized Expression Level", 
         color = "Groups", 
         shape = "Groups") + 
    theme_minimal(base_size = 14) + 
    theme(
      text = element_text(family = "Arial", color = "black"), 
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 14), 
      axis.text.y = element_text(margin = margin(r = 10), color = "black", size = 14), 
      axis.title.x = element_text(face = "bold", color = "black", size = 14, family = "Arial"),  
      axis.title.y = element_text(face = "bold", color = "black", size = 14, family = "Arial"), 
      panel.background = element_rect(fill = "white", color = "gray", linewidth = 1), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "black", family = "Arial")
    )
  
  print(gene_plot)
  
  readline(prompt = "Press [Enter] to continue to the next plot or manually stop execution if done.")
}


########### For Cluster 1: Select genes with positive values 0.8 >= 6 
# Load required data
Spearman <- read.csv("Spearman_Correlation.csv", header = TRUE)
Cluster3 <- read.csv("CorrelationCluster_3.csv", header = TRUE)


# Filter Spearman table to keep only genes present in Cluster
Cluster3_genes <- Cluster3$x
Filtered_Spearman <- Spearman %>% filter(Symbol %in% Cluster3_genes)


# Filter genes where at least six correlation values are greater than 0.8
filtered_genes <- Filtered_Spearman[apply(Filtered_Spearman[, -1], 1, function(row) sum(row > 0.8) >= 6), ]

# Load Lsmeans data
Lsmeans <- read.csv("Lsmeans_WOR.csv", header = TRUE)

# Filter Lsmeans for overlapping genes
Lsmeans_overlapping_genes <- Lsmeans %>%
  filter(Symbol %in% filtered_genes$Symbol)

# Reshape data to long format for easier plotting
Lsmeans_overlapping_long <- Lsmeans_overlapping_genes %>%
  select(FlyBase.ID, FB_Day5:MH_Day70) %>%
  pivot_longer(cols = -FlyBase.ID, 
               names_to = "Time_Point", 
               values_to = "Lsmeans") %>%
  mutate(
    Group = sub("_.*", "", Time_Point),
    Time_Point = sub(".*_", "", Time_Point),
    Time_Point = factor(Time_Point, levels = c("Day5", "Day7", "Day14", "Day21", "Day28", 
                                               "Day35", "Day42", "Day49", "Day56", 
                                               "Day63", "Day70"))
  )

# Define visualization properties
group_colors <- c("FB" = scales::alpha("#F56600", 0.5),
                  "FH" = scales::alpha("#522D80", 0.5),
                  "MB" = scales::alpha("#546223", 0.5),
                  "MH" = scales::alpha("#005EB8", 0.5))

group_shapes <- c("FB" = 16, "FH" = 15, "MB" = 17, "MH" = 18)

# Loop through all unique genes to create and display plots
for (gene in unique(Lsmeans_overlapping_long$FlyBase.ID)) {
  processed_data <- Lsmeans_overlapping_long %>%
    filter(FlyBase.ID == gene)
  
  gene_plot <- ggplot(processed_data, aes(x = Time_Point, y = Lsmeans, group = Group, color = Group, shape = Group)) +
    geom_line(linewidth = 1.2, alpha = 0.9) + 
    geom_point(size = 4, alpha = 0.7) +
    scale_color_manual(values = group_colors) +
    scale_shape_manual(values = group_shapes) +
    scale_x_discrete(labels = function(x) gsub("Day", "", x)) +
    labs(title = paste("Expression Profile -", gene),
         x = "Time (Days)", 
         y = "Normalized Expression Level", 
         color = "Groups", 
         shape = "Groups") +
    theme_minimal(base_size = 14) +
    theme(
      text = element_text(family = "Arial", color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 14),
      axis.text.y = element_text(margin = margin(r = 10), color = "black", size = 14),
      axis.title.x = element_text(face = "bold", color = "black", size = 14, family = "Arial"),
      axis.title.y = element_text(face = "bold", color = "black", size = 14, family = "Arial"),
      panel.background = element_rect(fill = "white", color = "gray", linewidth = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "black", family = "Arial")
    )  
  
  
  print(gene_plot)
  
  readline(prompt = "Press [Enter] to continue to the next plot or manually stop execution if done.")
}


########### For Cluster 4: Select genes with the highest negative correlation in any comparison

# Load required data
Spearman <- read.csv("Spearman_Correlation_Feb15.csv", header = TRUE)
Cluster4 <- read.csv("CorrelationCluster_4.csv", header = TRUE)


# Filter Spearman table to keep only genes present in Cluster
Cluster4_genes <- Cluster4$x
Filtered_Spearman <- Spearman %>% filter(Symbol %in% Cluster4_genes)


# Find genes with the highest negative correlation in any comparison
high_neg_corr_genes <- Filtered_Spearman %>%
  rowwise() %>%
  mutate(min_correlation = min(c_across(-Symbol))) %>%
  ungroup() %>%
  arrange(min_correlation) %>%
  slice_head(n = 50)


# Load Lsmeans data
Lsmeans <- read.csv("Lsmeans_WOR.csv", header = TRUE)

# Filter Lsmeans for overlapping genes
Lsmeans_overlapping_genes <- Lsmeans %>%
  filter(Symbol %in% high_neg_corr_genes$Symbol)

# Reshape data to long format for easier plotting
Lsmeans_overlapping_long <- Lsmeans_overlapping_genes %>%
  select(FlyBase.ID, FB_Day5:MH_Day70) %>%
  pivot_longer(cols = -FlyBase.ID, 
               names_to = "Time_Point", 
               values_to = "Lsmeans") %>%
  mutate(
    Group = sub("_.*", "", Time_Point),
    Time_Point = sub(".*_", "", Time_Point),
    Time_Point = factor(Time_Point, levels = c("Day5", "Day7", "Day14", "Day21", "Day28", 
                                               "Day35", "Day42", "Day49", "Day56", 
                                               "Day63", "Day70"))
  )

# Define visualization properties
group_colors <- c("FB" = scales::alpha("#F56600", 0.5),
                  "FH" = scales::alpha("#522D80", 0.5),
                  "MB" = scales::alpha("#546223", 0.5),
                  "MH" = scales::alpha("#005EB8", 0.5))

group_shapes <- c("FB" = 16, "FH" = 15, "MB" = 17,"MH" = 18)

# Loop through all unique genes to create and display plots
for (gene in unique(Lsmeans_overlapping_long$FlyBase.ID)) {
  processed_data <- Lsmeans_overlapping_long %>%
    filter(FlyBase.ID == gene)
  
  gene_plot <- ggplot(processed_data, aes(x = Time_Point, y = Lsmeans, group = Group, color = Group, shape = Group)) +
    geom_line(linewidth = 1.2, alpha = 0.9) +
    geom_point(size = 4, alpha = 0.7) +
    scale_color_manual(values = group_colors) +
    scale_shape_manual(values = group_shapes) +
    scale_x_discrete(labels = function(x) gsub("Day", "", x)) +
    labs(title = paste("Expression Profile -", gene),
         x = "Time (Days)", 
         y = "Normalized Expression Level", 
         color = "Groups", 
         shape = "Groups") +
    theme_minimal(base_size = 14) +
    theme(
      text = element_text(family = "Arial", color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 14),
      axis.text.y = element_text(margin = margin(r = 10), color = "black", size = 14),
      axis.title.x = element_text(face = "bold", color = "black", size = 14, family = "Arial"),
      axis.title.y = element_text(face = "bold", color = "black", size = 14, family = "Arial"),
      panel.background = element_rect(fill = "white", color = "gray", linewidth = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "black", family = "Arial")
    )
  
  print(gene_plot)
  
  readline(prompt = "Press [Enter] to continue to the next plot or manually stop execution if done.")
}



########### For Cluster 5: Select genes with the highest negative correlation in any comparison
# Load required data
Spearman <- read.csv("Spearman_Correlation.csv", header = TRUE)
Cluster5 <- read.csv("CorrelationCluster_5.csv", header = TRUE)

# Filter Spearman table to keep only genes present in Cluster
Cluster5_genes <- Cluster5$x
Filtered_Spearman <- Spearman %>% filter(Symbol %in% Cluster5_genes)


# Find genes with the highest negative correlation in any comparison
high_neg_corr_genes <- Filtered_Spearman %>%
  rowwise() %>%
  mutate(min_correlation = min(c_across(-Symbol))) %>%
  ungroup() %>%
  arrange(min_correlation) %>%
  slice_head(n = 50)

# Load Lsmeans data
Lsmeans <- read.csv("Lsmeans_WOR.csv", header = TRUE)

# Filter Lsmeans for overlapping genes
Lsmeans_overlapping_genes <- Lsmeans %>%
  filter(Symbol %in% high_neg_corr_genes$Symbol)


# Reshape data to long format for easier plotting
Lsmeans_overlapping_long <- Lsmeans_overlapping_genes %>%
  select(FlyBase.ID, FB_Day5:MH_Day70) %>%
  pivot_longer(cols = -FlyBase.ID, 
               names_to = "Time_Point", 
               values_to = "Lsmeans") %>%
  mutate(
    Group = sub("_.*", "", Time_Point),
    Time_Point = sub(".*_", "", Time_Point),
    Time_Point = factor(Time_Point, levels = c("Day5", "Day7", "Day14", "Day21", "Day28", 
                                               "Day35", "Day42", "Day49", "Day56", 
                                               "Day63", "Day70"))
  )

# Define visualization properties
group_colors <- c("FB" = scales::alpha("#F56600", 0.5),
                  "FH" = scales::alpha("#522D80", 0.5),
                  "MB" = scales::alpha("#546223", 0.5),
                  "MH" = scales::alpha("#005EB8", 0.5))

group_shapes <- c("FB" = 16, "FH" = 15, "MB" = 17, "MH" = 18)

# Loop through all unique genes to create and display plots
for (gene in unique(Lsmeans_overlapping_long$FlyBase.ID)) {
  processed_data <- Lsmeans_overlapping_long %>%
    filter(FlyBase.ID == gene)
  
  gene_plot <- ggplot(processed_data, aes(x = Time_Point, y = Lsmeans, group = Group, color = Group, shape = Group)) +
    geom_line(linewidth = 1.2, alpha = 0.9) +
    geom_point(size = 4, alpha = 0.7) +
    scale_color_manual(values = group_colors) +
    scale_shape_manual(values = group_shapes) +
    scale_x_discrete(labels = function(x) gsub("Day", "", x)) +
    labs(title = paste("Expression Profile -", gene),
         x = "Time (Days)", 
         y = "Normalized Expression Level", 
         color = "Groups", 
         shape = "Groups") +
    theme_minimal(base_size = 14) +
    theme(
      text = element_text(family = "Arial", color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 14),
      axis.text.y = element_text(margin = margin(r = 10), color = "black", size = 14),
      axis.title.x = element_text(face = "bold", color = "black", size = 14, family = "Arial"),
      axis.title.y = element_text(face = "bold", color = "black", size = 14, family = "Arial"),
      panel.background = element_rect(fill = "white", color = "gray", linewidth = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "black", family = "Arial")
    )
  
  print(gene_plot)
  
  readline(prompt = "Press [Enter] to continue to the next plot or manually stop execution if done.")
}


########### For Cluster 6: Select genes with the highest negative correlation in any comparison
# Load required data
Spearman <- read.csv("Spearman_Correlation.csv", header = TRUE)
Cluster6 <- read.csv("CorrelationCluster_6.csv", header = TRUE)


# Filter Spearman table to keep only genes present in Cluster
Cluster6_genes <- Cluster6$x
Filtered_Spearman <- Spearman %>% filter(Symbol %in% Cluster6_genes)


# Find genes with the highest negative correlation in any comparison
high_neg_corr_genes <- Filtered_Spearman %>%
  rowwise() %>%
  mutate(min_correlation = min(c_across(-Symbol))) %>%
  ungroup() %>%
  arrange(min_correlation) %>%
  slice_head(n = 50)

# Load Lsmeans data
Lsmeans <- read.csv("Lsmeans_WOR.csv", header = TRUE)

# Filter Lsmeans for overlapping genes
Lsmeans_overlapping_genes <- Lsmeans %>%
  filter(Symbol %in% high_neg_corr_genes$Symbol)

# Reshape data to long format for easier plotting
Lsmeans_overlapping_long <- Lsmeans_overlapping_genes %>%
  select(FlyBase.ID, FB_Day5:MH_Day70) %>%
  pivot_longer(cols = -FlyBase.ID, 
               names_to = "Time_Point", 
               values_to = "Lsmeans") %>%
  mutate(
    Group = sub("_.*", "", Time_Point),
    Time_Point = sub(".*_", "", Time_Point),
    Time_Point = factor(Time_Point, levels = c("Day5", "Day7", "Day14", "Day21", "Day28", 
                                               "Day35", "Day42", "Day49", "Day56", 
                                               "Day63", "Day70"))
  )

# Define visualization properties
group_colors <- c("FB" = scales::alpha("#F56600", 0.5),
                  "FH" = scales::alpha("#522D80", 0.5),
                  "MB" = scales::alpha("#546223", 0.5),
                  "MH" = scales::alpha("#005EB8", 0.5))

group_shapes <- c("FB" = 16, "FH" = 15, "MB" = 17, "MH" = 18)

# Loop through all unique genes to create and display plots
for (gene in unique(Lsmeans_overlapping_long$FlyBase.ID)) {
  processed_data <- Lsmeans_overlapping_long %>%
    filter(FlyBase.ID == gene)
  
  gene_plot <- ggplot(processed_data, aes(x = Time_Point, y = Lsmeans, group = Group, color = Group, shape = Group)) +
    geom_line(linewidth = 1.2, alpha = 0.9) +
    geom_point(size = 4, alpha = 0.7) +
    scale_color_manual(values = group_colors) +
    scale_shape_manual(values = group_shapes) +
    scale_x_discrete(labels = function(x) gsub("Day", "", x)) +
    labs(title = paste("Expression Profile -", gene),
         x = "Time (Days)", 
         y = "Normalized Expression Level", 
         color = "Groups", 
         shape = "Groups") +
    theme_minimal(base_size = 14) +
    theme(
      text = element_text(family = "Arial", color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 14),
      axis.text.y = element_text(margin = margin(r = 10), color = "black", size = 14),
      axis.title.x = element_text(face = "bold", color = "black", size = 14, family = "Arial"),
      axis.title.y = element_text(face = "bold", color = "black", size = 14, family = "Arial"),
      panel.background = element_rect(fill = "white", color = "gray", linewidth = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "black", family = "Arial")
    )
  
  print(gene_plot)
  
  readline(prompt = "Press [Enter] to continue to the next plot or manually stop execution if done.")
}
