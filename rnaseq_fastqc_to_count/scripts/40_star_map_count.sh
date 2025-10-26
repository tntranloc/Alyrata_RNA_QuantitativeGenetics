#!/usr/bin/env bash
#SBATCH --job-name=star_map
#SBATCH --partition=smp
#SBATCH --cpus-per-task=4
#SBATCH --mem=42gb
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err

module load star/2.7.8a || true
module load samtools/1.13 || true

source scripts/_lib.sh

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" config/samples.tsv)
SID=$(awk '{print $1}' <<< "$LINE")

R1P="${TRIM}/${SID}.R1P.fastq.gz"
R2P="${TRIM}/${SID}.R2P.fastq.gz"

echo "[star] ${SID}"
STAR --runThreadN "${THREADS}" \
     --genomeDir "${STAR_INDEX}" \
     --sjdbGTFfile "${GTF}" \
     --readFilesIn "${R1P}" "${R2P}" \
     --readFilesCommand zcat \
     --outFileNamePrefix "${MAP}/${SID}." \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts

# index BAM
samtools index -@ 2 "${MAP}/${SID}.Aligned.sortedByCoord.out.bam"
echo "[star] Done ${SID}"
