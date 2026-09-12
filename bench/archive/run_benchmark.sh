#!/usr/bin/env bash
# run_benchmark.sh: Standardized benchmark for multi-accession path walking.
set -e
cd "$(dirname "$0")/.."
ROOT=$(pwd)

PYTHON="${PYTHON:-python3}"
WALKER="$ROOT/python/logan_path_walker.py"

echo "=== Logan Multi-Accession Path Walker Benchmark ==="
echo "Host: $(hostname) | CPUs: $(nproc) | Memory: $(free -h | awk '/Mem:/ {print $2}')"
echo "Target: HIV-1 pol (50 accessions)"

OUT_DIR="$ROOT/bench/output_hiv1"
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

START=$(date +%s.%N)
$PYTHON "$WALKER" \
  --accession-list "$ROOT/bench/hiv1_50.txt" \
  --query "$ROOT/queries/query_hiv1_pol.fasta" \
  --out-dir "$OUT_DIR/work" \
  --out-fasta "$OUT_DIR/reconstructed.fasta" \
  --workers 16 \
  --threads 2

END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
N_RECS=$(grep -c '^>' "$OUT_DIR/reconstructed.fasta" 2>/dev/null || echo 0)

echo ""
echo "=== Benchmark Results ==="
echo "Total Wall Time:  ${DIFF}s"
echo "Reconstructed:    ${N_RECS} / 50 accessions"
THROUGHPUT=$(echo "scale=2; 50 / $DIFF" | bc)
echo "Throughput:       ${THROUGHPUT} accessions/second"
echo "Output Directory: $OUT_DIR"
