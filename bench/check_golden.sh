#!/usr/bin/env bash
# Golden concordance gate for logan-walker.
#
# Runs the walker under old-behaviour flags (--prune-frac 0 --no-prefer-majority --no-emit-alts)
# and checks the primary haplotype against bench/golden/*.primary.fasta.
#
# Requires GraphAligner on PATH and network access to the Logan S3 bucket, so this is NOT part of
# the offline every-push CI. Run it locally, or via the manual "golden" CI job.
#
#   bench/check_golden.sh          # exit 0 = concordant, 1 = regression
#
# ERR2028820 (PfCRT) is the byte-exact gate: a clonal isolate, a small subgraph, and a
# deterministic GraphAligner run, so its primary is reproducible byte for byte.
#
# SRR8209082 (HIV-1 pol) is reported but NOT gated. Its 26k-node subgraph is ~a third of the
# whole graph; the Rust k-mer anchoring builds a different subgraph than the Python BLAST run
# the fixture was captured with, and GraphAligner is not deterministic on a graph that size when
# unpruned. It is therefore not byte-reproducible across implementations. See
# docs/rust_walker_as_built.md section 7 and docs/issue_intra_accession_variance_lost.md
# section 15. Its primary is stabilised by the default pruning, which this old-behaviour gate
# turns off on purpose.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.."; pwd)"
WALKER="${LOGAN_WALKER_BIN:-$ROOT/walker/target/release/logan-walker}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

if [ ! -x "$WALKER" ]; then
  echo "building walker..."; ( cd "$ROOT/walker" && cargo build --release ) || exit 1
fi
command -v GraphAligner >/dev/null || { echo "FAIL: GraphAligner not on PATH"; exit 1; }

run_primary() {  # acc query -> writes $WORK/$acc.got (bare sequence)
  local acc="$1" query="$2"
  "$WALKER" -a "$acc" -q "$ROOT/$query" -o "$WORK/$acc" -f "$WORK/$acc/recon.fasta" \
    --prune-frac 0.0 --no-prefer-majority --no-emit-alts >/dev/null 2>&1
  grep -v '^>' "$WORK/$acc/recon.fasta" 2>/dev/null | tr -d '\n ' > "$WORK/$acc.got"
}

fails=0

# --- byte-exact gate: ERR2028820 / PfCRT ---
run_primary ERR2028820 queries/pfcrt_only.fasta
grep -v '^>' "$ROOT/bench/golden/ERR2028820.pfcrt.primary.fasta" | tr -d '\n ' > "$WORK/exp"
if [ -s "$WORK/ERR2028820.got" ] && cmp -s "$WORK/ERR2028820.got" "$WORK/exp"; then
  echo "PASS ERR2028820  ($(wc -c < "$WORK/exp") nt, byte-identical)"
else
  echo "FAIL ERR2028820  (got $(wc -c < "$WORK/ERR2028820.got" 2>/dev/null || echo 0) nt vs golden $(wc -c < "$WORK/exp") nt)"
  fails=$((fails+1))
fi

# --- informational: SRR8209082 / HIV-1 pol (not gated) ---
run_primary SRR8209082 queries/query_hiv1_pol.fasta
got_len=$(wc -c < "$WORK/SRR8209082.got" 2>/dev/null || echo 0)
exp_len=$(grep -v '^>' "$ROOT/bench/golden/SRR8209082.hiv1_pol.primary.fasta" | tr -d '\n ' | wc -c)
echo "INFO SRR8209082  (walker $got_len nt, golden $exp_len nt; not byte-gated, see header)"

if [ "$fails" -eq 0 ]; then echo "golden: concordant"; exit 0; else echo "golden: $fails regression(s)"; exit 1; fi
