# Gene Expression Clustering, Heatmap Visualization and Interactive Network 

This repository contains R scripts and input files for:
1. **K-means clustering** of age-significant genes, followed by heatmap visualization using the [`ComplexHeatmap`](https://jokergoo.github.io/ComplexHeatmap-reference/book/) package.
2. **Interactive network analysis** and **gene set enrichment** of clusters using Cytoscape.

The goal is to identify gene expression patterns across time points and explore network-level relationships.

---

## Part 1 — K-means Clustering & Heatmap

**Script:** `Age_Kmeans_Heatmap.R`  
Performs **K-means clustering** on gene expression data and generates a publication-quality heatmap.

### Input Files
- **`Lsmeans_Age_Sig.csv`** — Rows: genes; Columns: expression values for 11 time points (Day 5 to Day 70).

### Output Files
- **`gene_clusters.txt`** — Tab-delimited file mapping each gene to its assigned cluster.
- **`NumberofClusters_WSS.svg`** — Elbow plot used to determine the optimal number of clusters.
- **`Age_Heatmap.svg`** — Heatmap of age-significant genes in 6 K-means clusters.

### Requirements
- **R**
- R packages:
  - `ComplexHeatmap`
  - `circlize`

---

## Part 2 — Interactive Network Analysis

Generates an **interactive network** for each of the K-means clusters from Part 1.

### Input Files
- **`OMNI_Kcluster1.csv`**  
- **`OMNI_Kcluster2.csv`**  
- **`OMNI_Kcluster3.csv`**  
- **`OMNI_Kcluster4.csv`**  
- **`OMNI_Kcluster5.csv`**  
- **`OMNI_Kcluster6.csv`**  
- **`All_xQTL.csv`**  
- **`Fly_orthologs.csv`** — From [DIOPT](https://www.flyrnai.org/diopt).  
- **`flybase_interaction_database_2024_02.txt`** — Compiled from:  
  - [Gene genetic interactions (FB2022_03)](https://ftp.flybase.net/releases/FB2022_03/precomputed_files/genes/gene_genetic_interactions_fb_2024_02.tsv.gz)  
  - [Physical interactions (FB2022_03)](https://ftp.flybase.net/releases/FB2022_03/precomputed_files/genes/physical_interactions_mitab_fb_2022_03.tsv.gz)

### Cytoscape Parameters

**Preprocessing:**
- Imported `flybase_interaction_database_2024_02.txt` into Cytoscape.
- Removed duplicate edges and self-loops.

**Topology Filtering:**
- Selected nodes with at least 4 neighbors within a distance of 1.
- Match **any (OR)** condition enabled.

**MCODE Settings:**
- Haircut mode enabled
- Node score cutoff: `0.4`

### Output Files
- **`Network_Cluster1.csv`**
- **`Network_Cluster2.csv`**
- **`Network_Cluster3.csv`**
- **`Network_Cluster4.csv`**
- **`Network_Cluster5.csv`**
- **`Network_Cluster6.csv`**

---

## Workflow Overview
```mermaid
graph TD
  A[Lsmeans_Age_Sig.csv] --> B[K-means clustering]
  B --> C[Heatmap visualization]
  C --> D[Cluster-specific CSV files]
  D --> E[Network construction in Cytoscape]
  E --> F[Topology filtering & MCODE]
  F --> G[Network_Cluster*.csv outputs]
