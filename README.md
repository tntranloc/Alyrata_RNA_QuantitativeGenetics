# PhD Project Scripts Repository  
**Quantitative Genetics in Plant Transcriptomics & Fitness**

**Author:** Nhu L. T. Tran  
**Email:** ntran5@uni-koeln.de  
**Last updated:** 06 April 2026  

---

## Overview
This repository contains a collection of scripts developed during my PhD project focusing on **quantitative genetics**, **plant transcriptomics**, and **fitness-related analyses**.  

The workflows cover a wide range of analyses, including:
- Animal models for genetic variance estimation  
- RNA-seq processing and variant calling  
- Population genetics statistics (pN/pS, Ka/Ks)  
- gene expression-QTL (eQTL) and methylation-QTL (meQTL)
- Transcriptome-wide association studies (TWAS), traditionally as well as using Machine learning approaches
- Machine learning approaches for associating genetic architecture with modes of inheritance  
- Functional enrichment and clustering  
- High-throughput phenotyping via computer vision  

---

## Main Workflows

### 🧬 Animal Models & Genetic Variance
- **`AnimalModel_collected_scripts.Rmd`**  
  Overview of animal model implementations over time, from `pedigreemm` to modern `brms`.

- **`complete_AnimalModel_BRMS_scripts`**  
  Current, fully developed pipeline for running large-scale animal models on HPC clusters.

- **`genetic_variance_analysis_brms.R`**  
  Script for large-scale genetic variance analysis (e.g., gene expression matrices) using `brms`.

---

### RNA-seq & Variant Calling
- **`complete_RNAseq_workflow_from_fastq_to_counttable.sh`**  
  End-to-end RNA-seq pipeline: FASTQ → processed reads → count table.

- **`complete_VariantCalling_workflow_for_RNAseq_gatk.sh`**  
  Variant calling workflow using RNA-seq data with GATK best practices.

- **`RNA_integrity_check.sh`**  
  Script for assessing RNA integrity post-mapping.

---

### Population Genetics Analyses
- **`complete_pNpS_ratios_for_gene_groups_workflow.txt`**  
  Workflow to compute pN/pS ratios using bootstrap resampling in parental populations.

- **`KaKs/`**  
  Scripts and instructions for calculating Ka/Ks (aka dN/dS) ratios in *Arabidopsis lyrata*.

---

### Functional Genomics & Association Studies
- **`eQTL_analysis_with_qtltools.sh`**  
  Pipeline for expression QTL analysis using QTLtools.

- **`transcriptome_wide_association_study_withGEMMA.sh`**  
  TWAS workflow using GEMMA.

- **`twas-like-ML/`**  
  Python-based TWAS-like pipeline including:
  - PCA-based regression (PCA-OLS with back-transformed gene effects)  
  - Gradient Boosting (regression & classification)  
  - Hypergeometric overlap enrichment  
  - YAML-configured reproducible runs  
  - Outputs ranked gene lists and performance metrics  

---

### Machine Learning & Statistical Methods
- **`randomForest_for_geneticArchitecture.R`**  
  Random forest analysis for exploring genetic architecture.

- **`enrichment_and_clustering.R`**  
  Functional enrichment and hierarchical clustering workflows.

---

### Phenotyping & Image Analysis
- **`plantCV_computervision_...md`**  
  Python-based workflow for high-throughput plant image analysis.

---

### Variant Annotation
- **`snpeff_annotation_and_how_to_build_own_refgenome_database.sh`**  
  Workflow for:
  - Annotating variants (synonymous vs non-synonymous) using SnpEff  
  - Building a custom reference genome database  

---

## 🧩 Supporting Scripts

- **`plot_quantGene_concept_model.md`**  
  Visualization of genetic concepts:
  - Additive effects  
  - Dominance  
  - Overdominance / underdominance  
  - Genetic variance components  

---

## Requirements
Most workflows rely on a combination of:
- **R** (e.g., `brms`, `pedigreemm`, ML packages)  
- **Python** (for ML and image analysis pipelines)  
- **Bash** (HPC workflows)  

External tools are mentioned in the scripts. 

---

## Notes
- Many scripts are designed for **high-performance computing (HPC)** environments.  
- Workflows are modular and can be adapted to different datasets or species.  
- Configuration files (e.g., YAML) are used where possible for reproducibility.  
