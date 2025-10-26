# GATK RNA-seq Variant Calling ETL Pipeline

**Last updated:** 12 Nov 2024  

This pipeline performs variant calling from RNA-seq BAMs using GATK 4.6.1.0, 
following an ETL (Extract–Transform–Load) pattern:

1. **Extract:** per-sample preprocessing and GVCF generation  
2. **Transform:** combine all GVCFs and genotype across samples  
3. **Load:** filter, subset, and generate QC statistics  

---

## Requirements

| Tool | Version (tested) |
|------|------------------|
| GATK | 4.6.1.0 |
| Java | 17.0.6 |
| Picard | ≥ 2.27 |
| bcftools | ≥ 1.17 |
| samtools | ≥ 1.17 |
| SLURM | (for HPC array jobs) |

Optional: `conda env create -f envs/gatk.yaml`

## Basic usage

```bash
# 1️⃣ submit per-sample extraction (array)
sbatch --array=1-$(wc -l < config/samples.tsv) scripts/gatk_rnaseq_etl.sh

# 2️⃣ after all samples complete, run cohort merge/filter
sbatch scripts/gatk_rnaseq_cohort.sh

