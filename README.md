# Systems Genetics of Lifespan and Senescence in *Drosophila melanogaster*

**Authors:**  
Maryam Nasiri Aghdam<sup>1</sup>, Desireé Unselt<sup>2,3</sup>, Maria Adonay<sup>1</sup>, Tatania V. Morozova<sup>2,4</sup>, Gunjan H. Arya<sup>2,5</sup>, Lavanya Turlapati<sup>2</sup>, Vijay Shankar<sup>1</sup>, Robert R. H. Anholt<sup>1,2,6</sup>, Trudy F. C. Mackay<sup>1,2,6</sup>  

**Affiliations:**  
<sup>1</sup> Center for Human Genetics, Clemson University, Greenwood, SC, USA  
<sup>2</sup> Program in Genetics, W. M. Keck Center for Behavioral Biology and Department of Biological Sciences, North Carolina State University, Raleigh, NC 27695, USA  
<sup>3</sup> Labcorp, 3601 Davis Drive, Morrisville, NC 27560, USA  
<sup>4</sup> BioSkryb Genomics, 2810 Meridian Parkway, Suite 110, Durham, NC 27713, USA  
<sup>5</sup> Life Edit Therapeutics, 300 Morrise Street, Suite 300, Durham, NC 27701, USA  
<sup>6</sup> Department of Genetics and Biochemistry, Clemson University, Clemson, SC, USA  



---

## Study Overview
We investigated the genetic architecture of lifespan and senescence in *Drosophila melanogaster* using an outbred advanced intercross population derived from 37 inbred, fully sequenced lines of the *Drosophila* Genetic Reference Panel (DGRP). Whole-genome sequencing was performed on young and old flies, separately for males and females, to identify extreme quantitative trait loci (xQTLs) and genes associated with increased lifespan. We also performed bulk RNA sequencing of male and female heads, bodies, and reproductive tissues weekly to 10 weeks of age.

---

## **Analysis Workflow and Code Resources**

This repository and the linked resources contain all scripts and workflows used in the project. The components are organized as follows:

- **xQTL Mapping**  
  Workflow for mapping quantitative trait loci in the Advanced Intercross Population (AIP), adapted from  
  [xQTL mapping in the AIP – Wen Huang](https://github.com/qgg-lab/dgrp-lifespan).

- **AIP Parental Population**
  Scripts for (i) combined frequency distribution and line means plot and (ii) full and reduced ANOVAs are available in the `Parental_Population` directory by Maria E. Adonay.
  
- **RNA-seq Quality Control, Alignment, and Read Counting**  
  RNA-seq preprocessing pipeline adapted from  
  [RNA-seq pipeline – Vijay Shankar](https://github.com/vshanka23/snakemake_rnaseq).

- **RNA-seq Normalization**  
  Scripts for normalization and transformation of RNA-seq data are available in the `RNA_seq_Normalization` directory.

- **ANOVA for Transcriptome Analysis**  
  Transcriptome-wide ANOVA scripts are located in the `Aging_Transriptome` directory.

- **Variance Partitioning**  
  Code for partitioning variance components is provided in the `Variation_Partition` directory.

- **Context-Dependent Responses**  
  Analysis scripts for context-dependent effects are available in the `Context_Dependent_Responses` directory.

- **Interaction Network Construction and Analysis**  
  Scripts for building and analyzing interaction networks are contained in the `Interaction_Network` directory.



 