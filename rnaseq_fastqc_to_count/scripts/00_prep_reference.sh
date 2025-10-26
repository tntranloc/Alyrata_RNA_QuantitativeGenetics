#!/usr/bin/env bash
#SBATCH --job-name=prep_ref
#SBATCH --partition=smp
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err

module load star/2.7.8a || true
module load gffread || true

source scripts/_lib.sh

# 1) Convert GFF -> GTF if GFF provided
if [[ -n "${GFF:-}" && "${GFF}" != "" ]]; then
  echo "[prep] Converting GFF to GTF..."
  gffread "${GFF}" -T -o "${GTF}"
fi

# 2) STAR genome index
echo "[prep] Building STAR index at ${STAR_INDEX} ..."
mkdir -p "${STAR_INDEX}"
STAR --runThreadN "${THREADS}" --runMode genomeGenerate \
     --genomeDir "${STAR_INDEX}" \
     --genomeFastaFiles "${REF}" \
     --sjdbGTFfile "${GTF}" \
     --genomeSAindexNbases $(cfg_get star_genomeSAindexNbases)

echo "[prep] Done."
