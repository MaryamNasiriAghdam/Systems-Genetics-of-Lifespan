# ============================================================
# Title: Figure 2. Figure 2. xQTL mapping results. 
# Manuscript: Systems Genetics of Lifespan and Senescence in Drosophila melanogaster (Nasiri Aghdam et al., 2025)
# Author: Maria E. Adonay (@amalgamaria)
# Date: 2025-08-20
# Description: 
#   This script generates:
#     (A) Female xQTL Manhattan plot
#     (B) Male xQTL Manhattan plot
#   Outputs a combined figure in TIFF format
# R Version: 4.1.2
# ============================================================

# -----------------------------
# Prepare environment
# -----------------------------
rm(list = ls())
set.seed(123)

library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)

setwd("/path/to/input")
infile_f <- "filename_f.txt"
infile_m <- "filename_m.txt"
outfile <- "lifespan_manhattan_7w5h_300dpi.tif"

# -----------------------------
# Function: Read, preprocess data
# -----------------------------
process_data <- function(file, sex_label) {
  df <- fread(file, header = TRUE) %>% 
    setnames(c("chr", "pos", "ref", "major", "minor",
               "o_1", "o_2", "o_3", "o_ave",
               "r_1", "r_2", "r_3", "r_ave",
               "diff", "p")) %>% 
    na.omit() %>% 
    mutate(
      sex = sex_label,
      p_adj = p.adjust(p, method = "bonferroni"),
      sig = if_else(p.adjust(p, method = "bonferroni") <= 0.05, 1, 0))
  return(df)
}

# -----------------------------
# Read, preprocess female + male data separately
# -----------------------------
df_f <- process_data(infile_f, "Female")
df_m <- process_data(infile_m, "Male")

# -----------------------------
# Calculate thresholds per sex
# -----------------------------
calc_threshold <- function(df_sex) {
  df_sex <- arrange(df_sex, p_adj)
  
  p_below <- max(df_sex$p[df_sex$p_adj <= 0.05], na.rm = TRUE)
  p_above <- min(df_sex$p[df_sex$p_adj > 0.05], na.rm = TRUE)
  
  raw_p_thresh <- mean(c(p_below, p_above))
  log10_raw_p_thresh <- -log10(raw_p_thresh)
  
  data.frame(
    sex = unique(df_sex$sex),
    raw_p_thresh = raw_p_thresh,
    log10_raw_p_thresh = log10_raw_p_thresh)
}

thresholds <- bind_rows(calc_threshold(df_f), calc_threshold(df_m))

# -----------------------------
# Combine data
# -----------------------------
df <- bind_rows(df_f, df_m)

# -----------------------------
# Define chromosome levels, lengths
# -----------------------------
chr_levels <- c("2L", "2R", "3L", "3R", "4", "X", "Y")
chr_lengths <- data.frame(
  chr = factor(chr_levels, levels = chr_levels),
  chr_len = c(23513712, 25286936, 28110227, 32079331, 1348131, 23542271, 3667352))

# -----------------------------
# Add Female Y placeholder
# -----------------------------
fake_Y_female <- tibble(
  chr = factor("Y", levels = chr_levels),
  pos = 1,
  sig = 0,
  sex = "Female",
  p = 1,
  p_adj = 1,
  fill_color = "#FFFFFF",
  plot_alpha = 0)         

df <- bind_rows(df, fake_Y_female)

# -----------------------------
# Assign colors, alpha for plotting
# -----------------------------
df <- df %>%
  mutate(
    fill_color = case_when(
      sex == "Female" & chr == "Y" ~ "#FFFFFF",
      sex == "Female" & sig == 1 & chr != "Y" ~ "#522D80",
      sex == "Female" & sig == 0 & chr != "Y" ~ "#A186C2",
      sex == "Male" & sig == 1 ~ "#F56600",
      sex == "Male" & sig == 0 ~ "#FFB380"
    ),
    plot_alpha = if_else(sex == "Female" & chr == "Y", 0, 0.7))

# -----------------------------
# Define common theme
# -----------------------------
common_theme <- theme_minimal(base_family = "Arial") +
  theme(
    text = element_text(family = "Arial"),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(family = "Arial", size = 8),
    axis.title = element_text(family = "Arial", size = 10),
    strip.text = element_text(family = "Arial", size = 8, face = "italic"),
    plot.tag = element_text(family = "Arial", size = 12, face = "bold"),
    plot.margin = margin(10, 25, 10, 25))

# -----------------------------
# Function: Manhattan plot
# -----------------------------
plot_manhattan <- function(df, sex_label, thresholds_df, y_max = NULL) {
  df_sex <- filter(df, sex == sex_label)
  thresh <- filter(thresholds_df, sex == sex_label) %>% pull(log10_raw_p_thresh)
  
  if (is.null(y_max)) {y_max <- max(-log10(df_sex$p), na.rm = TRUE)}
  
  p <- ggplot(df_sex, aes(x = pos, y = -log10(p))) +
    geom_point(aes(color = fill_color, alpha = plot_alpha), size = 1.2) +
    scale_color_identity() +
    scale_alpha_identity() +
    geom_hline(yintercept = thresh, linetype = "dashed",
               color = "black", linewidth = 0.7, alpha = 0.75) +
    facet_grid(cols = vars(chr), scales = "free_x", space = "free_x", switch = "x") +
    geom_blank(data = chr_lengths, aes(x = chr_len), inherit.aes = FALSE) +
    scale_y_continuous(
      breaks = seq(0, ceiling(y_max/10)*10, by = 10),
      minor_breaks = seq(0, ceiling(y_max/5)*5, by = 5),
      expand = c(0,0)
    ) +
    coord_cartesian(ylim = c(0, y_max), clip = "off") +
    labs(x = ifelse(sex_label == "Male", "Chromosome Arm", ""), 
         y = expression(-log[10](italic("p")))) +
    common_theme
  
  if (sex_label == "Female") {p <- p + theme(strip.text.x = element_blank())}
  
  return(p)
}

# -----------------------------
# Generate, combine female + male plots
# -----------------------------
manhat_f <- plot_manhattan(df, "Female", thresholds, y_max = 35)
manhat_m <- plot_manhattan(df, "Male", thresholds, y_max = 20)

manhat_combined <- manhat_f / manhat_m +
  plot_layout(heights = c(3, 2)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag.position = c(-0.02, 1))

# -----------------------------
# Save
# -----------------------------
ggsave(outfile, 
       plot = manhat_combined, 
       width = 7, height = 5, dpi = 300, units = "in", 
       device = "tiff", bg = "white")
