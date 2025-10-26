#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

cfg_get() { grep -E "^$1:" config/config.yaml | awk '{print $2}'; }

REF=$(cfg_get ref_fasta)
GTF=$(cfg_get gtf || true)
GFF=$(cfg_get gff || true)
STAR_INDEX=$(cfg_get star_index_dir)

RES=$(cfg_get results_dir)
RAW_QC=$(cfg_get fastqc_raw_dir)
TRIM=$(cfg_get trim_dir)
TRIM_QC=$(cfg_get fastqc_trim_dir)
MAP=$(cfg_get map_dir)

ADAPT=$(cfg_get adapters_fa)
THREADS=$(cfg_get threads)

mkdir -p "$RES" "$RAW_QC" "$TRIM" "$TRIM_QC" "$MAP" logs
