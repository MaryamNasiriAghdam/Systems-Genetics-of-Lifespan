# Pre-processing Count Table and RNA-seq Normalization

This repository contains R scripts and input files for the following workflow:  
1. **Pre-processing Count Table**  
2. **RNA-seq Normalization**  

---

## Part 1 — Pre-processing Count Table

**Script:** `RNA_Seq_Normalization.R`  

**Input Files**  
- `combined_counts.txt` — Raw combined count table from RNA-seq analysis.

**Output Files**  
- `New_Combined_Cleaned_Raw_Counts.txt` — Cleaned and formatted count table suitable for downstream normalization.

**Requirements**  
- **R**  
- R packages:  
  - `dplyr`  

---

## Part 2 — RNA-seq Normalization

**Script:** `RNA_Seq_Normalization.R`  

**Input Files**  
- `New_Combined_Cleaned_Raw_Counts.txt` — Pre-processed count table.

**Output Files**  
- `New_TMM_Normalized.csv` — TMM-normalized expression matrix.  

**Requirements**  
- **R**  
- R packages:  
  - `dplyr`  
  - `matrixStats`  
  - `data.table`  
  - `tidyr`  
  - `edgeR`  
  - `limma`  

---

**Acknowledgments**  
Part of this workflow was developed in collaboration with [Colleague Name](https://github.com/username).  

