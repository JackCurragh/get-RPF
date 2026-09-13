#!/usr/bin/env bash
# Fetch the leading reads of each real-panel run from ENA
# (validation/panel_v1/truth.yaml). Only the first N reads are transferred:
# the gzip stream is closed once N records have been read.
#
# Usage: scripts/fetch_panel_subsets.sh [N_READS] [OUT_DIR]
# Leading reads come from the first flow-cell tiles; the library structure is
# the same, but quality at the tile edges can be slightly lower.
set -euo pipefail

N_READS=${1:-1000000}
OUT_DIR=${2:-validation/panel_v1/data}
mkdir -p "$OUT_DIR"

while read -r run path; do
  [ -z "$run" ] && continue
  out="$OUT_DIR/$run.fastq.gz"
  if [ -s "$out" ]; then
    echo "$run: already present"
    continue
  fi
  # head exits after N records; curl and gzip then stop on SIGPIPE, which is
  # expected here, so their exit codes are ignored and the result is checked.
  { curl -sS --fail "https://$path" || true; } \
    | { gzip -dc 2>/dev/null || true; } \
    | head -n $((N_READS * 4)) \
    | gzip -c > "$out.part"
  lines=$(gzip -dc "$out.part" | wc -l | tr -d ' ')
  if [ "$lines" -eq 0 ] || [ $((lines % 4)) -ne 0 ]; then
    echo "$run: incomplete download ($lines lines)" >&2
    exit 1
  fi
  mv "$out.part" "$out"
  echo "$run: $((lines / 4)) reads -> $out"
done <<'EOF'
SRR3945920 ftp.sra.ebi.ac.uk/vol1/fastq/SRR394/000/SRR3945920/SRR3945920.fastq.gz
SRR3945930 ftp.sra.ebi.ac.uk/vol1/fastq/SRR394/000/SRR3945930/SRR3945930.fastq.gz
SRR12693498 ftp.sra.ebi.ac.uk/vol1/fastq/SRR126/098/SRR12693498/SRR12693498.fastq.gz
SRR23242345 ftp.sra.ebi.ac.uk/vol1/fastq/SRR232/045/SRR23242345/SRR23242345.fastq.gz
EOF
