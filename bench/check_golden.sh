#!/usr/bin/env bash
# Golden concordance gate for logan-walker. Two sets, both run here:
#
#   1. old-behaviour: the walker under --prune-frac 0 --no-prefer-majority --no-emit-alts, with
#      the primary haplotype checked against bench/golden/*.primary.fasta. These fixtures
#      predate 0.1.0's defaults, which is why the flags turn those defaults off.
#   2. default-flag: the walk with no flags at all, with every output it writes checked against
#      bench/golden/default/. Added 2026-11-13; see the second block below for what it covers
#      and, importantly, what it does not.
#
# Requires GraphAligner on PATH and network access to the Logan S3 bucket, so this is NOT part of
# the offline every-push CI (.github/workflows/ci.yml). It runs in its own workflow,
# .github/workflows/golden.yml, on pushes to main, on v* and logan-walker-v* tags, on pull
# requests that touch walker/, weekly, and on manual dispatch.
#
#   bench/check_golden.sh          # exit 0 = concordant, 1 = regression
#
# Set LOGAN_UNITIGS_CACHE to a directory to keep the downloaded unitig objects between runs.
#
# ERR2028820 (PfCRT) and SRR12363134 (C. auris ERG11) are the byte-exact gates: clonal
# isolates with small subgraphs and a deterministic GraphAligner run, so their primaries are
# reproducible byte for byte. SRR12363134's fixture was cross-validated: the Rust walker and
# the (now removed) Python walker produced the identical 1,574 nt primary under GraphAligner
# 1.0.20. All three fixtures were recaptured from the Rust walker under 1.0.17b on 2026-10-18,
# which gives the two gated ones their final base (bench/golden/README.md).
#
# SRR8209082 (HIV-1 pol) is reported but NOT gated. Its 26k-node subgraph is ~a third of the
# whole graph, so its unpruned route depends on the anchoring implementation and on the
# GraphAligner version (1.0.17b and 1.0.20 take different routes), and GraphAligner is not
# deterministic on a graph that size with more than one thread. It is therefore not
# byte-reproducible across implementations. See the Reproducibility section of docs/walker.md
# and docs/issue_intra_accession_variance_lost.md. Its primary is stabilised by the default
# pruning, which this old-behaviour gate turns off on purpose.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.."; pwd)"
WALKER="${LOGAN_WALKER_BIN:-$ROOT/walker/target/release/logan-walker}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

if [ ! -x "$WALKER" ]; then
  echo "building walker..."; ( cd "$ROOT/walker" && cargo build --release ) || exit 1
fi
command -v GraphAligner >/dev/null || { echo "FAIL: GraphAligner not on PATH"; exit 1; }

# Optional shared store for the downloaded Logan objects: a re-run then costs no download, and
# a transient S3 failure does not turn the gate red. Logan's `u/` objects are immutable, so a
# cached object is never stale. The CI workflow points this at a cached directory; left unset,
# every run downloads afresh, which is the old behaviour.
CACHE_ARG=()
if [ -n "${LOGAN_UNITIGS_CACHE:-}" ]; then
  mkdir -p "$LOGAN_UNITIGS_CACHE" || { echo "FAIL: cannot create $LOGAN_UNITIGS_CACHE"; exit 1; }
  CACHE_ARG=(--unitigs-cache "$LOGAN_UNITIGS_CACHE")
  echo "unitig cache: $LOGAN_UNITIGS_CACHE"
fi

# Record the aligner version with the result. The fixtures are version-dependent: 1.0.19+
# clips the query's final base, so these goldens (captured under 1.0.17b, which aligns end to
# end) FAIL on 1.0.19 and 1.0.20 by exactly that base. Without this line a green gate does not
# say which behaviour it green-lit. See docs/issue_graphaligner_end_clipping.md.
running="$(GraphAligner --version 2>&1 | head -n 1)"
echo "using $running at $(command -v GraphAligner)"
captured="$(cat "$ROOT/bench/golden/GRAPHALIGNER_VERSION" 2>/dev/null)"
if [ -n "$captured" ] && [ "$running" != "$captured" ]; then
  echo "WARN fixtures were captured under '$captured'; a FAIL below may be the aligner, not the walker"
fi

# Writes $WORK/$acc.got (bare sequence) and $WORK/$acc.log, and RETURNS the walker's exit
# status. That status used to be thrown away along with stderr, so an S3 timeout produced an
# absent recon.fasta, an empty .got, and the message "FAIL <acc> (got 0 nt vs golden 1476 nt)":
# an infrastructure outage reported as a walk regression, with every diagnostic discarded.
run_primary() {  # acc query -> rc
  local acc="$1" query="$2"
  # Separate statement on purpose: every word of a simple command is expanded before the
  # command runs, so `local acc="$1" log="$WORK/$acc.log"` would expand $acc to whatever the
  # CALLER last left in it, not to "$1".
  local log="$WORK/$acc.log"
  "$WALKER" -a "$acc" -q "$ROOT/$query" -o "$WORK/$acc" -f "$WORK/$acc/recon.fasta" \
    ${CACHE_ARG[@]+"${CACHE_ARG[@]}"} \
    --prune-frac 0.0 --no-prefer-majority --no-emit-alts >"$log" 2>&1
  local rc=$?
  grep -v '^>' "$WORK/$acc/recon.fasta" 2>/dev/null | tr -d '\n ' > "$WORK/$acc.got"
  return $rc
}

# Print the tail of a failed run's log, indented, so CI output carries the cause.
show_log() {
  local acc="$1"
  echo "    --- last 5 lines of $acc.log ---"
  tail -n 5 "$WORK/$acc.log" 2>/dev/null | sed 's/^/    | /'
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
  if ! run_primary "$acc" "$query"; then
    # Since logan-walker exits non-zero when nothing reconstructed, this catches the fetch
    # and GraphAligner failures that used to masquerade as a golden mismatch.
    echo "ERROR $acc  (walker exited non-zero: a run failure, NOT a walk regression)"
    show_log "$acc"
    fails=$((fails+1))
    continue
  fi
  grep -v '^>' "$ROOT/$golden" | tr -d '\n ' > "$WORK/$acc.exp"
  if [ -s "$WORK/$acc.got" ] && cmp -s "$WORK/$acc.got" "$WORK/$acc.exp"; then
    echo "PASS $acc  ($(wc -c < "$WORK/$acc.exp") nt, byte-identical)"
  else
    echo "FAIL $acc  (got $(wc -c < "$WORK/$acc.got" 2>/dev/null || echo 0) nt vs golden $(wc -c < "$WORK/$acc.exp") nt)"
    fails=$((fails+1))
  fi
done

# --- default-flag gates: the walk every user actually runs ---
#
# The gate above gives up three of the walk's four shaping steps on purpose, to reproduce the
# pre-0.1.0 behaviour its fixtures were captured under: it passes --prune-frac 0
# --no-prefer-majority --no-emit-alts and compares only the bare primary sequence. So the
# prune/realign acceptance, variants::prefer_majority, and the alt, VCF and SNV-table writers
# had no end-to-end check at all, in any test (review TC-02). These fixtures close that: they
# are the default walk's own outputs, captured 2026-11-13 with --emit-vcf --emit-snv-table
# under GraphAligner 1.0.17b at graphaligner_threads=1 (the default).
#
# SRR8209082 IS gated here, unlike in the old-behaviour set above. That is not an inconsistency:
# the header explains that its route is stabilised BY the default pruning, which the gate above
# turns off, and two consecutive default runs reproduced byte for byte across all five outputs
# on 2026-11-13. Note the limit of that evidence: one machine, one thread. GraphAligner is not
# deterministic on a graph this size with more than one thread. If this flaps in CI, move it to
# the informational section below rather than deleting the gate.
#
# No normalisation is needed. The VCF carries no ##fileDate and no command line, and the FASTA
# headers hold only walk statistics, so every file here is byte-comparable as written.

# Compare one produced file against its fixture. A missing fixture means the file must NOT
# exist: ERR2028820 and SRR12363134 recover no SNPs, so an appearing VCF is a real change.
cmp_fixture() {  # got fixture label
  local got="$1" fixture="$2" label="$3"
  if [ ! -e "$fixture" ]; then
    if [ -s "$got" ]; then
      echo "  FAIL $label  (no fixture, but the walker wrote $(wc -c < "$got") B)"
      fails=$((fails+1))
    else
      echo "  PASS $label  (absent, as expected)"
    fi
    return
  fi
  if [ ! -e "$got" ]; then
    echo "  FAIL $label  (fixture exists, but the walker wrote nothing)"
    fails=$((fails+1))
  elif cmp -s "$got" "$fixture"; then
    echo "  PASS $label  ($(wc -c < "$fixture") B, byte-identical)"
  else
    echo "  FAIL $label  (got $(wc -c < "$got") B vs fixture $(wc -c < "$fixture") B)"
    fails=$((fails+1))
  fi
}

# bubbles.json is gated by content where it is small and by sha256 where it is not:
# SRR8209082 has 3,632 bubbles and a 1 MB file, which does not belong in the repository and
# could not be read as a diff anyway. It is regenerable, so a checksum is the right gate.
cmp_bubbles() {  # got tag
  local got="$1" tag="$2"
  local body="$ROOT/bench/golden/default/$tag.bubbles.json"
  local sum="$body.sha256"
  if [ -e "$body" ]; then
    cmp_fixture "$got" "$body" "$tag bubbles.json"
  elif [ -e "$sum" ]; then
    local want got_sum
    want="$(cat "$sum")"
    got_sum="$(sha256sum "$got" | awk '{print $1}')"
    if [ "$want" = "$got_sum" ]; then
      echo "  PASS $tag bubbles.json  (sha256 $(echo "$want" | cut -c1-12)..., $(wc -l < "$got") lines)"
    else
      echo "  FAIL $tag bubbles.json  (sha256 ${got_sum:0:12}... != fixture ${want:0:12}...)"
      fails=$((fails+1))
    fi
  else
    echo "  WARN $tag bubbles.json  (no fixture and no checksum; not gated)"
  fi
}

# summary.json is gated with its volatile fields masked: the clock, the paths, and the
# versions under test (bench/mask_summary.py says exactly which, and why). Everything that
# describes the WALK is gated, including every count in `discards` and `cohort`, which is the
# point: those numbers used to exist only on stderr, so nothing could regress them visibly.
cmp_summary() {  # got fixture label
  local got="$1" fixture="$2" label="$3"
  if [ ! -e "$fixture" ]; then
    echo "  WARN $label  (no fixture; not gated)"
    return
  fi
  if [ ! -e "$got" ]; then
    echo "  FAIL $label  (fixture exists, but the walker wrote no summary)"
    fails=$((fails+1))
    return
  fi
  local a="$WORK/.mask.got.$$" b="$WORK/.mask.want.$$"
  if ! python3 "$ROOT/bench/mask_summary.py" "$got" > "$a" 2>"$a.err"; then
    echo "  FAIL $label  (could not mask the produced summary: $(head -n 1 "$a.err"))"
    fails=$((fails+1)); rm -f "$a" "$a.err"; return
  fi
  if ! python3 "$ROOT/bench/mask_summary.py" "$fixture" > "$b" 2>"$b.err"; then
    echo "  FAIL $label  (could not mask the fixture: $(head -n 1 "$b.err"))"
    fails=$((fails+1)); rm -f "$a" "$b" "$a.err" "$b.err"; return
  fi
  if cmp -s "$a" "$b"; then
    echo "  PASS $label  (masked summary identical)"
  else
    echo "  FAIL $label  (masked summary differs)"
    diff -u "$b" "$a" | head -n 20 | sed 's/^/    | /'
    fails=$((fails+1))
  fi
  rm -f "$a" "$b" "$a.err" "$b.err"
}

DEFAULT_GATES=(
  "ERR2028820  queries/pfcrt_only.fasta                   ERR2028820.pfcrt"
  "SRR12363134 queries/markers/candida_auris_erg11.fasta  SRR12363134.candida_auris_erg11"
  "SRR8209082  queries/query_hiv1_pol.fasta               SRR8209082.hiv1_pol"
)
echo "--- default-flag gates ---"
for g in "${DEFAULT_GATES[@]}"; do
  set -- $g; acc="$1"; query="$2"; tag="$3"
  d="$WORK/def-$tag"
  mkdir -p "$d"
  # Default flags throughout: the only additions are the two writers that are off by default,
  # because an output that is never produced cannot be gated.
  "$WALKER" -a "$acc" -q "$ROOT/$query" -o "$d" -f "$d/primary.fasta" \
    ${CACHE_ARG[@]+"${CACHE_ARG[@]}"} \
    --emit-vcf --emit-snv-table >"$d/run.log" 2>&1
  if [ $? -ne 0 ]; then
    echo "  ERROR $tag  (walker exited non-zero: a run failure, NOT a walk regression)"
    echo "    --- last 5 lines of $tag run.log ---"
    tail -n 5 "$d/run.log" 2>/dev/null | sed 's/^/    | /'
    fails=$((fails+1))
    continue
  fi
  F="$ROOT/bench/golden/default/$tag"
  cmp_fixture "$d/primary.fasta"        "$F.primary.fasta" "$tag primary"
  cmp_fixture "$d/alt_haplotypes.fasta" "$F.alts.fasta"    "$tag alts"
  cmp_fixture "$d/snvs.tsv"             "$F.snvs.tsv"      "$tag snvs"
  cmp_fixture "$d/$acc.vcf"             "$F.vcf"           "$tag vcf"
  cmp_bubbles "$d/bubbles.json"         "$tag"
  # The summary names the walk's own query stem, which is the query file's stem.
  qstem="$(basename "$query")"; qstem="${qstem%.*}"
  cmp_summary "$d/$qstem.summary.json"  "$F.summary.json"  "$tag summary"
done

# --- informational: SRR8209082 / HIV-1 pol (not gated) ---
exp_len=$(grep -v '^>' "$ROOT/bench/golden/SRR8209082.hiv1_pol.primary.fasta" | tr -d '\n ' | wc -c)
if run_primary SRR8209082 queries/query_hiv1_pol.fasta; then
  got_len=$(wc -c < "$WORK/SRR8209082.got" 2>/dev/null || echo 0)
  echo "INFO SRR8209082  (walker $got_len nt, golden $exp_len nt; not byte-gated, see header)"
else
  # Not gated, so this does not fail the run, but it must not read as a 0 nt walk either.
  echo "INFO SRR8209082  (walker exited non-zero; golden $exp_len nt; not byte-gated)"
  show_log SRR8209082
fi

if [ "$fails" -eq 0 ]; then echo "golden: concordant"; exit 0; else echo "golden: $fails regression(s)"; exit 1; fi
