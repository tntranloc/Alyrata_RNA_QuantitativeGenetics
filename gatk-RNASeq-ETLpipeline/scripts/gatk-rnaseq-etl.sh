#!/usr/bin/env bash
#SBATCH -J rnaseq_gatk
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 24:00:00
#SBATCH -o logs/%x.%A_%a.log
#SBATCH -e logs/%x.%A_%a.err

set -euo pipefail
IFS=$'\n\t'

# === load config ===
REF=$(grep '^ref:' config/config.yaml | awk '{print $2}')
OUTDIR=$(grep '^outdir:' config/config.yaml | awk '{print $2}')
PICARD=$(grep '^picard:' config/config.yaml | awk '{print $2}')
SAMPLES=config/samples.tsv

mkdir -p "$OUTDIR/per_sample" logs

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES")
SAMPLE=$(awk '{print $1}' <<< "$LINE")
BAM=$(awk '{print $2}' <<< "$LINE")
SAMP_OUT="$OUTDIR/per_sample/$SAMPLE"
mkdir -p "$SAMP_OUT"

echo "Processing $SAMPLE"

# 0. Add read groups
java -Xmx4g -jar "$PICARD" AddOrReplaceReadGroups \
  I="$BAM" \
  O="$SAMP_OUT/${SAMPLE}.rg.bam" \
  RGID="$SAMPLE" RGLB=lib RGPL=ILLUMINA RGPU=unit RGSM="$SAMPLE" \
  CREATE_INDEX=true

# 1. Mark duplicates
gatk MarkDuplicatesSpark \
  -I "$SAMP_OUT/${SAMPLE}.rg.bam" \
  -O "$SAMP_OUT/${SAMPLE}.dedup.bam" \
  -M "$SAMP_OUT/${SAMPLE}.dup_metrics.txt" \
  --create-output-bam-index true

# 2. SplitNCigarReads
gatk SplitNCigarReads \
  -R "$REF" \
  -I "$SAMP_OUT/${SAMPLE}.dedup.bam" \
  -O "$SAMP_OUT/${SAMPLE}.split.bam"

# 3. HaplotypeCaller
gatk --java-options "-XX:+UseSerialGC" HaplotypeCaller \
  -R "$REF" \
  -I "$SAMP_OUT/${SAMPLE}.split.bam" \
  -O "$SAMP_OUT/${SAMPLE}.g.vcf.gz" \
  -ERC GVCF \
  --dont-use-soft-clipped-bases \
  --disable-read-filter MappingQualityAvailableReadFilter \
  --min-pruning 1

echo "Done $SAMPLE → $SAMP_OUT/${SAMPLE}.g.vcf.gz"
