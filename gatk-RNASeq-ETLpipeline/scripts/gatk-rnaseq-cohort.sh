#!/usr/bin/env bash
#SBATCH -J rnaseq_gatk_cohort
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 24:00:00
#SBATCH -o logs/%x.%j.log
#SBATCH -e logs/%x.%j.err

set -euo pipefail
REF=$(grep '^ref:' config/config.yaml | awk '{print $2}')
OUTDIR=$(grep '^outdir:' config/config.yaml | awk '{print $2}')

mkdir -p "$OUTDIR/cohort" logs
cd "$OUTDIR/cohort"

# === combine and genotype ===
VCFS=$(find ../per_sample -name "*.g.vcf.gz" | sort)
gatk GenomicsDBImport \
  --genomicsdb-workspace-path genomicsdb \
  $(for v in $VCFS; do echo -V "$v"; done) \
  --reader-threads 4

gatk GenotypeGVCFs \
  -R "$REF" \
  -V gendb://genomicsdb \
  -O genotyped.vcf.gz \
  --max-alternate-alleles 2

# === SNPs only ===
gatk SelectVariants \
  -R "$REF" -V genotyped.vcf.gz \
  --select-type-to-include SNP \
  -O genotyped.snps.vcf.gz

# === bcftools filters ===
bcftools view -i 'F_MISSING < 0.1 && FORMAT/DP > 10 && FORMAT/DP < 100 && MAF >= 0.05' \
  genotyped.snps.vcf.gz -Oz -o genotyped.snps.filtered.vcf.gz
tabix -p vcf genotyped.snps.filtered.vcf.gz

# === stats ===
bcftools stats genotyped.snps.filtered.vcf.gz > stats.txt

echo "Final filtered VCF: $OUTDIR/cohort/genotyped.snps.filtered.vcf.gz"
