# ============================================================
# Title: Figure 1. Genetic variation for lifespan in the AIP founder lines.
# Manuscript: Systems Genetics of Lifespan and Senescence in Drosophila melanogaster (Nasiri Aghdam et al., 2025)
# Author: Maria E. Adonay (@amalgamaria)
# Date: 2025-08-20
# Description: 
#   This script generates:
#     (A) Frequency distribution of lifespan by sex
#     (B) Mean lifespan per DGRP line (with SE), ordered by males
#   Outputs a combined figure in TIFF format
# R Version: 4.1.2
# ============================================================

# -----------------------------
# Prepare environment
# -----------------------------
rm(list = ls())
set.seed(123)

library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(patchwork)

setwd("/path/to/input")
infile <- "filename.csv"
outfile <- "lifespan_frequency-dist_line-means_12w5h_300dpi.tif"

# -----------------------------
# Read, clean data
# -----------------------------
data <- read.csv(infile, na.strings = ".", header = TRUE)

data_clean <- data %>%
  filter(!is.na(lifespan)) %>%
  mutate(sex = str_to_title(sex))

# -----------------------------
# Define common theme
# -----------------------------
common_theme <- theme_minimal(base_family = "Arial") +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "none",
    
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.3),

    axis.text = element_text(size = 12, family = "Arial"),
    axis.title = element_text(size = 14, family = "Arial"),
    
    plot.tag = element_text(size = 16, face = "bold", family = "Arial"))

# -----------------------------
# A: Frequency histogram
# -----------------------------
freq <- ggplot(data_clean, aes(x = lifespan, fill = sex)) +
  geom_histogram(bins = 10, alpha = 0.7, position = "dodge") +  
  scale_fill_manual(values = c("Female" = "#522D80", "Male" = "#F56600")) +
  labs(
    x = "\nLifespan (Days)",
    y = "Frequency\n") +
  scale_y_continuous(
    breaks = seq(0, 600, 200),
    minor_breaks = seq(0, 600, 100),
    expand = expansion(mult = c(0, 0.05))) +
  common_theme +
  theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, family = "Arial"))

# -----------------------------
# Compute line-level stats
# -----------------------------
line_stats <- data_clean %>%
  mutate(line = paste0("DGRP_", line)) %>%
  group_by(line, sex) %>%
  summarise(
    mean_lifespan = mean(lifespan),
    n = n(),
    SE = sqrt(var(lifespan) / n),
    .groups = "drop")

# -----------------------------
# Order lines by male mean lifespan
# -----------------------------
line_order_m <- line_stats %>%
  filter(sex == "Male") %>%
  arrange(mean_lifespan) %>%
  pull(line)

line_stats <- line_stats %>%
  mutate(line = factor(line, levels = line_order_m))

# -----------------------------
# B: Line means ± SE
# -----------------------------
means <- ggplot(line_stats, aes(x = line, y = mean_lifespan, color = sex)) +
  geom_point(size = 1.5, alpha = 0.7) +  
  geom_errorbar(aes(ymin = mean_lifespan - SE, ymax = mean_lifespan + SE),
                width = 0.2, alpha = 0.7, linewidth = 0.75) +
  scale_color_manual(values = c("Female" = "#522D80", "Male" = "#F56600")) +
  labs(
    x = "\nLine (Ordered by Increasing Male Mean Lifespan)",
    y = "Mean Lifespan (Days)\n") +
  scale_y_continuous(
    limits = c(0, 80),
    breaks = seq(0, 80, 20),
    minor_breaks = seq(0, 80, 10)) +
  common_theme +
  theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1, family = "Arial"))

# -----------------------------
# Combine plots
# -----------------------------
freq_means <- freq + means + 
  plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = 'A')

# -----------------------------
# Save
# -----------------------------
ggsave(outfile, 
       plot = freq_means, 
       width = 12, height = 5, dpi = 300, units = "in", 
       device = "tiff", bg = "white")
