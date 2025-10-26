#!/usr/bin/env bash
#SBATCH --job-name=fastqc_raw
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
R1=$(awk '{print $2}' <<< "$LINE")
R2=$(awk '{print $3}' <<< "$LINE")

echo "[fastqc_raw] $SID"
fastqc -t "${THREADS}" -o "${RAW_QC}" "${R1}"
fastqc -t "${THREADS}" -o "${RAW_QC}" "${R2}"
