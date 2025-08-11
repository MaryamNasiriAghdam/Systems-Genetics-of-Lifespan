# Variation Partition Analysis

This repository contains an R script to visualize the **percentage contribution** of different sources of variation using a **bar plot**.  
The visualization is based on **sum of squares (SS)** calculations performed externally (e.g., in Excel).

---

## R Script
- **`Variation_Partition.R`** — Generates an **inverted bar plot** showing the percentage contribution of each source of variation.

---

## Input Files
- **`Variation_Partition.xlsx`** — Excel file containing sum of squares (SS) values for each source of variation.

---

## Output Files
- **`Variation_Partition.svg`** — Inverted bar plot illustrating the percentage contribution of each variation source.

---

## Requirements
- **R** (version ≥ 4.0 recommended)
- R package:
  - `ggplot2`

**Install in R:**
```r
install.packages("ggplot2")
library(ggplot2)

