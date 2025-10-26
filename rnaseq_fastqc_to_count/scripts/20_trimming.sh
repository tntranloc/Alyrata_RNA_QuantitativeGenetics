#!/usr/bin/env bash
#SBATCH --job-name=trimmo
#SBATCH --partition=smp
#SBATCH --cpus-per-task=4
#SBATCH --mem=42gb
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err

module load trimmomatic/0.39 || true
module load openjdk/17 || module load openjdk/1.8.0_60 || true

source scripts/_lib.sh

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" config/samples.tsv)
SID=$(awk '{print $1}' <<< "$LINE")
R1=$(awk '{print $2}' <<< "$LINE")
R2=$(awk '{print $3}' <<< "$LINE")

base="${TRIM}/${SID}"
echo "[trim] ${SID}"

# Output names (paired/unpaired)
R1P="${base}.R1P.fastq.gz"; R1U="${base}.R1U.fastq.gz"
R2P="${base}.R2P.fastq.gz"; R2U="${base}.R2U.fastq.gz"

trimmomatic PE -threads "${THREADS}" -phred33 \
  "${R1}" "${R2}" \
  "${R1P}" "${R1U}" "${R2P}" "${R2U}" \
  ILLUMINACLIP:${ADAPT}:2:30:10 \
  LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50

echo "[trim] Done ${SID}"
