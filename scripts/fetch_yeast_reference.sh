#!/usr/bin/env bash
# Fetch the S. cerevisiae R64-1-1 genome (Ensembl release 108, matching
# bench/sacCer3's GTF) and build a STAR index for the P7 alignment check.
#
# Usage: scripts/fetch_yeast_reference.sh [OUT_DIR]
# Then:  GETRPF_YEAST_STAR_INDEX=OUT_DIR/star pytest tests/test_structure/test_real_alignment.py
set -euo pipefail

OUT_DIR=${1:-validation/panel_v1/data/reference}
URL=https://ftp.ensembl.org/pub/release-108/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
FASTA="$OUT_DIR/$(basename "$URL")"

mkdir -p "$OUT_DIR"
[ -s "$FASTA" ] || curl -sS --fail -o "$FASTA" "$URL"
python -c "
from pathlib import Path
from getRPF.core.structure.align import build_star_index
build_star_index(Path('$FASTA'), Path('$OUT_DIR/star'))
print('STAR index: $OUT_DIR/star')
"
