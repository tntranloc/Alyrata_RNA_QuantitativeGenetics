# ------------------------------------------------------------
# Script Generator for Parallelized BRMS Gene Expression Models
# Nhu L.T. Tran
# March 2025
#
# This script generates multiple R and SLURM job files for submitting
# brms-based animal models in parallel on a computing cluster.
#
# Each job processes a batch of 200 genes and skips genes already processed.
# The main brms model is wrapped in: animal_model_RNAseq_brms_March2025_wrapped_function.r
# ------------------------------------------------------------

# Number of genes to process per job script
genes_per_script = 200
# Total number of genes to model
total_genes = 6270
# Number of scripts/jobs needed
num_scripts = ceiling(total_genes / genes_per_script)

# Set working directory where scripts will be written to
setwd("/projects/ag-demeaux/ltntran/scripts/AnimalModel_F2/scripts_6k")

# Generate R and SLURM script for each chunk of genes
for (i in 1:num_scripts) {
  st = (i - 1) * genes_per_script + 1
  en = min(i * genes_per_script, total_genes)
  
# ---------- Create the R script to be executed by SLURM ----------
  script_file = paste0("chunk_", st, "_", en, ".R")
  writeLines(c(
  'Sys.setenv(TMPDIR = "/scratch/ntran5/tmp/")',
  "library(brms)",
  "library(posterior)",
  "library(bayestestR)",
  "library(Matrix)",
  "library(dplyr)",
  # Define the gene indices to process
  paste0("st = ", st),
  paste0("end = ", en),
  '',
  '# Load input files',
  'load("/projects/ag-demeaux/ltntran/scripts/AnimalModel_F2/input_6k/counts.RData")',
  'load("/projects/ag-demeaux/ltntran/scripts/AnimalModel_F2/input_6k/dam_info.RData")',
  'load("/projects/ag-demeaux/ltntran/scripts/AnimalModel_F2/input_6k/Amat.RData")',
  'load("/projects/ag-demeaux/ltntran/scripts/AnimalModel_F2/input_6k/Dmat.RData")',
  '',
  '# Set output directory',
  'outdir = "/projects/ag-demeaux/ltntran/scripts/AnimalModel_F2/output_6k/"',
  'dir.create(outdir, showWarnings = FALSE)',
  '',
  '# Run gene processing loop',
  'individual_files = TRUE',
   # Run the actual model function
  'source("/projects/ag-demeaux/ltntran/scripts/AnimalModel_F2/animal_model_RNAseq_brms_March2025_wrapped_function_annotated.r")'
), con = script_file)

# ---------- SLURM script generation follows here ----------
  job_file = paste0("job_", st, "_", en, ".sh")
  writeLines(c(
    "#!/bin/bash",
    "#SBATCH --nodes=1",
    "#SBATCH --time=24:00:00",
    "#SBATCH --mem=14G",
    "#SBATCH --cpus-per-task=16",
    "#SBATCH --account=ag-demeaux",
    "#SBATCH --error=/projects/ag-demeaux/ltntran/scripts/AnimalModel_F2/errors/bayes_batch.%J.err",
    "#SBATCH --out=/projects/ag-demeaux/ltntran/scripts/AnimalModel_F2/errors/bayes_batch.%J.out",
    "module load lang/R/4.4.1-gfbf-2023b",
    "export TMPDIR=/scratch/ntran5/tmp/",
    paste0("Rscript ", script_file)
  ), con = job_file)
}

