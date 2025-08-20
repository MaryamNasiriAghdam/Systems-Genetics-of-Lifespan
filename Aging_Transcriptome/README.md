# Pre-SAS Table Preparation, ANOVA in SAS, and FDR Correction (Benjamini–Hochberg)

This repository contains R scripts and input files for the following workflow:
1. **Pre-SAS Table Preparation**
2. **ANOVA in SAS**
3. **FDR Correction (Benjamini–Hochberg)**

---

## Part 1 — Pre-SAS Table Preparation

**Script:** `ANOVA.R`  

**Input Files**  
- `TMM_Transposed.csv` — Rows: 132 samples (11 time points); columns: expression values for each gene.

**Output Files**  
- `TMM_PRE_SAS.csv` — CSV file with columns for Sample, Age, Replicate, Tissue, Sex, FBgn, and Expression.

---

## Part 2 — ANOVA in SAS

**Script:** `ANOVA.R`  

**Input Files**  
- `TMM_PRE_SAS.csv`

**Output Files**  
- `CLASSLEVELS_OMNI.csv`  
- `FS_OMNI.csv`  
- `LSMEANS_OMNI.csv`  
- `MODELANOVA_OMNI.csv`  
- `NOBS_OMNI.csv`  
- `OVERALLANOVA_OMNI.csv`  
- `OMNI_Type3.csv`  

**Requirements**  
- **SAS** (`PROC GLM`)

---

## Part 3 — FDR Correction (Benjamini–Hochberg)

**Script:** `ANOVA.R`  

**Input Files**  
- `OMNI_Type3.csv`

**Output Files**  
- `GLM_Type3_Padj.txt` — Table of p-values adjusted using the Benjamini–Hochberg method.

**Requirements**  
- **R**
- R packages:  
  - `rstatix`  
  - `data.table`  
  - `tidyr`  
  - `dplyr`  
  - `bit64`
