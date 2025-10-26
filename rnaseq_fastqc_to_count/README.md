# RNA-seq ETL: FastQC → Trimmomatic → STAR (+ GeneCounts)

**Author:** Nhu L.T. Tran · **Email:** ntran5@uni-koeln.de  
SLURM-ready pipeline for RNA-seq preprocessing: raw QC → trimming → (optional) post-trim QC → alignment with STAR and per-gene counts.

## Requirements
- SLURM scheduler
- STAR ≥ 2.7, FastQC ≥ 0.11.9, Trimmomatic 0.39, samtools ≥ 1.13
- Java (for Trimmomatic)
- Optional: `conda env create -f envs/rnaseq.yaml`

## Inputs
- `config/samples.tsv`: tab-separated: `SAMPLE_ID  R1.fastq.gz  R2.fastq.gz`
- `config/config.yaml`: paths & parameters (adapters, index, GTF/GFF, threads)

## Quick start
```bash
# 0) (once) prepare reference (STAR index, GFF->GTF if needed)
sbatch scripts/00_prep_reference.sh

# 1) Raw FastQC
sbatch --array=1-$(($(wc -l < config/samples.tsv))) scripts/10_fastqc_raw.sh

# 2) Trimming
sbatch --array=1-$(($(wc -l < config/samples.tsv))) scripts/20_trimming.sh

# 3) (optional) FastQC on trimmed reads
sbatch --array=1-$(($(wc -l < config/samples.tsv))) scripts/30_fastqc_trimmed.sh

# 4) STAR mapping + GeneCounts
sbatch --array=1-$(($(wc -l < config/samples.tsv))) scripts/40_star_map_counts.sh
