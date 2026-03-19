# Longitudinal single-cell RNA-seq analysis of immune responses

This repository contains code for longitudinal single-cell transcriptomic analysis to investigate dynamic immune responses over time.

---

## Overview

The analysis focuses on identifying temporal changes in gene expression and uncovering regulatory programs underlying immune heterogeneity.

All analyses are implemented in **R**, primarily using **Seurat** and downstream statistical and functional analysis frameworks.

---

## Analysis Structure

The scripts are organized according to the main analysis steps:

- **01_longitudinal_DE_screening**  
  Paired differential expression analysis across time points using **MAST**, performed across multiple cell types to identify the most responsive populations.

- **02_monocyte_module_analysis**  
  Downstream analysis focused on monocytes, including:
  - Construction of effect-size matrices across time points  
  - Gene module identification using hierarchical clustering  
  - Heatmap visualization of temporal expression patterns  

- **03_monocyte_module_enrichment**  
  Functional interpretation of gene modules using:
  - Gene Ontology (GO) enrichment  
  - Reactome pathway analysis  

---

## Methods and Tools

- **Single-cell analysis**: Seurat  
- **Differential expression**: MAST (paired design)  
- **Clustering**: hierarchical clustering (Ward.D2)  
- **Visualization**: ComplexHeatmap  
- **Functional analysis**: clusterProfiler, ReactomePA  

---

## Key Features

- Longitudinal analysis of transcriptomic changes  
- Integration of multiple time-point comparisons  
- Identification of gene modules capturing dynamic regulatory programs  
- Functional interpretation of immune-related pathways  

---

## Notes

- Data are not included due to privacy and data protection restrictions  
- Scripts are provided as representative examples of the computational workflow  
