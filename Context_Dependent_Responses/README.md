# Context-Dependent Responses in Gene Expression

This repository contains R scripts and datasets for analyzing **context-dependent gene expression patterns**.  
The analysis performs:
- **Spearman correlation** of gene expression profiles.
- Clustering of genes based on correlation patterns.
- Visualization of expression profiles over time.

---

## R Script
- **`Context_Dependent_Responses.R`** — Main analysis script for:
  - Reading and processing correlation data.
  - Filtering genes.
  - Performing clustering.
  - Generating visualizations.

---

## Input Files
- **`Lsmeans_WOR.csv`** — Lsmeans values for all groups, excluding **Reproductive** tissue.
- **`FB_Lsmeans.csv`** — Lsmeans values for **Female Body** samples.
- **`FH_Lsmeans.csv`** — Lsmeans values for **Female Head** samples.
- **`MB_Lsmeans.csv`** — Lsmeans values for **Male Body** samples.
- **`MH_Lsmeans.csv`** — Lsmeans values for **Male Head** samples.

---

## Output Files
All results are generated in the `Output_files/` directory:
- **`Spearman_Correlation.csv`** — Spearman correlation coefficients for six tissue/condition comparisons.
- **`WSS_Plot.svg`** — Within Sum of Squares (WSS) plot used to determine the optimal number of clusters.
- **`Correlation_Heatmap.svg`** — Heatmap of Spearman correlation coefficients.
- **`CorrelationCluster_1.csv`** — Genes in Cluster 1.
- **`CorrelationCluster_2.csv`** — Genes in Cluster 2.
- **`CorrelationCluster_3.csv`** — Genes in Cluster 3.
- **`CorrelationCluster_4.csv`** — Genes in Cluster 4.
- **`CorrelationCluster_5.csv`** — Genes in Cluster 5.
- **`CorrelationCluster_6.csv`** — Genes in Cluster 6.

---

## Requirements
- **R** (version ≥ 4.0 recommended)
- R packages:
  - `tidyverse`
  - `ggplot2`
  - `dplyr`
  - `tidyr`

**Install packages in R:**
```r
install.packages(c("tidyverse", "ggplot2", "dplyr", "tidyr"))
