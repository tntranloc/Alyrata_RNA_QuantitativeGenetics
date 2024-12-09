#### Authors: Nhu L.T.Tran
#### Email: ntran5@uni-koeln.de
#### Last updated: 30.08.2024

#### Collected scripts for my PhD project in Quantitative Genetics in plant transcriptome and fitness. 

Main scripts
- AnimalModel_collected_scripts.Rmd shows methods to run Animal Model over the past decade, starting with pedigreemm and ending with current brms package
- complete_RNAseq_workflow_from_fastq_to_counttable.sh shows how to perform RNA-seq analysis from scratch (fastq) to final count table ready for any downstream analysis
- genetic_variance_analysis_brms.R is an Rscript to run Animal Model on large scale data such as gene count matrix on cluster
- eQTL_analysis_with_qtltools.sh shows how to perform expression QTL analysis with qtltools
- plantCV_computervision_...md is a Python workflow to analyse plant images in a higher throughput way
- calculate pips_ratio.sh: calculate nucleotide diversity in a population using variant calling input
- snpeff.sh: using snpeff to extract synonymous and non.synonymous variants using variant calling input
- enrichment_and_clustering.R: functional enrichment analysis and hierachial clustering
- randomForest_for_geneticArchitecture.R: random forest analysis
- transcriptome_wide_association_study_withGEMMA.sh: running TWAS using simple workflow with GEMMA


Supporting scripts:
- plot_quantGene_concept_model.md shows how to visualise allele substitution effect and genetic variance components in different scenarios of dominance, additive, underdominance, overdominance, etc
- RNA_integrity_check.sh shows how to check RNA integrity after mapping
