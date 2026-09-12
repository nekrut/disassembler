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
# ERR2028820 (PfCRT) and SRR12363134 (C. auris ERG11) are the byte-exact gates: clonal
# isolates with small subgraphs and a deterministic GraphAligner run, so their primaries are
# reproducible byte for byte. SRR12363134's fixture was cross-validated: the Rust walker and
# the (now removed) Python walker produced the identical 1,574 nt primary.
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

# --- byte-exact gates: clonal, deterministic accessions ---
#     accession   query                                golden primary
GATES=(
  "ERR2028820  queries/pfcrt_only.fasta               bench/golden/ERR2028820.pfcrt.primary.fasta"
  "SRR12363134 queries/markers/candida_auris_erg11.fasta bench/golden/SRR12363134.candida_auris_erg11.primary.fasta"
)
for g in "${GATES[@]}"; do
  set -- $g; acc="$1"; query="$2"; golden="$3"
  run_primary "$acc" "$query"
  grep -v '^>' "$ROOT/$golden" | tr -d '\n ' > "$WORK/$acc.exp"
  if [ -s "$WORK/$acc.got" ] && cmp -s "$WORK/$acc.got" "$WORK/$acc.exp"; then
    echo "PASS $acc  ($(wc -c < "$WORK/$acc.exp") nt, byte-identical)"
  else
    echo "FAIL $acc  (got $(wc -c < "$WORK/$acc.got" 2>/dev/null || echo 0) nt vs golden $(wc -c < "$WORK/$acc.exp") nt)"
    fails=$((fails+1))
  fi
done

# --- informational: SRR8209082 / HIV-1 pol (not gated) ---
run_primary SRR8209082 queries/query_hiv1_pol.fasta
got_len=$(wc -c < "$WORK/SRR8209082.got" 2>/dev/null || echo 0)
exp_len=$(grep -v '^>' "$ROOT/bench/golden/SRR8209082.hiv1_pol.primary.fasta" | tr -d '\n ' | wc -c)
echo "INFO SRR8209082  (walker $got_len nt, golden $exp_len nt; not byte-gated, see header)"

if [ "$fails" -eq 0 ]; then echo "golden: concordant"; exit 0; else echo "golden: $fails regression(s)"; exit 1; fi
