#!/usr/bin/env bash
#SBATCH --job-name=fastqc_trim
#SBATCH --partition=smp
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err

module load fastqc/0.11.9 || true
source scripts/_lib.sh

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" config/samples.tsv)
SID=$(awk '{print $1}' <<< "$LINE")

for fq in "${TRIM}/${SID}.R1P.fastq.gz" "${TRIM}/${SID}.R2P.fastq.gz"; do
  echo "[fastqc_trim] ${SID}: $(basename "$fq")"
  fastqc -t "${THREADS}" -o "${TRIM_QC}" "$fq"
done
