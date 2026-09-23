#!/usr/bin/env bash
#
# test_sieve.sh -- concordance test for the k-mer diversity sieve (l1_pipeline_serial).
#
# Runs the sieve directly on a fixture of 33 REAL reconstructed HIV-pol haplotypes
# (~2737 nt each) and asserts:
#   (a) it runs (exit 0),
#   (b) the promoted count is >= 1 and <= the input count (33),
#   (c) it is deterministic: two runs on the same input promote the same count,
#   (d) the promoted count is SMALL (<= 3).
#
# WHY (c)+(d) look the way they do -- read before "fixing" the bound:
#   The sieve is a HAPLOTYPE/LINEAGE-level diversity reducer built on MinHash/LSH.
#   The fixture haplotypes are INTRA-HOST alternative alleles of one another: they
#   differ by a single SNV on a ~2737 nt backbone. For a 2737 nt sequence differing
#   by one base the Jaccard similarity of the k-mer sets is ~= 1.0 -- far above ANY
#   usable LSH promotion threshold -- so the sieve treats them as the SAME lineage
#   and collapses them to ONE representative. That single-SNV signal is BELOW the
#   MinHash resolution by design. This is exactly why the sieve is meant for LARGE
#   cross-accession cohorts of near-duplicate LINEAGES, and NOT for preserving
#   intra-host SNVs (disassembler's --emit-vcf / --emit-snv-table do that instead).
#   So we assert "small" (<= 3), not a brittle exact number, while still requiring
#   >= 1 (it must promote a representative) and determinism (same count both runs).
#
# Dependency-free: needs only a C compiler (cc) + bash. No conda, no python.

set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
SIEVE_DIR="$REPO/extern/kmer-diversity-sieve"
SIEVE_BIN="$SIEVE_DIR/l1_pipeline_serial"
FIXTURE="$HERE/fixture_reconstructed.fasta"

# sieve parameters (per task spec)
MIN_LEN=2000
MAX_N_FRAC=0.05
THRESHOLD=0.15
WARMUP=10

EXPECTED_INPUT_N=33     # MAX_PROMOTED is calibrated against a fixture of this size
MAX_PROMOTED=3          # intra-host single-SNV alts collapse hard (see header)

fail() { echo "FAIL: $*" >&2; exit 1; }

skip() { echo "SKIP: $*"; exit 0; }

# --- 0. the sieve source lives in a submodule; skip cleanly if it isn't -----
# checked out (e.g. CI without submodules: true). A missing submodule is not a
# test failure -- it just means there's nothing to test here.
[[ -f "$SIEVE_DIR/Makefile" ]] || skip "submodule $SIEVE_DIR not checked out (run: git submodule update --init)"
command -v cc >/dev/null 2>&1 || skip "no C compiler (cc) available"

# --- 1. build the sieve if the binary is missing ---------------------------
if [[ ! -x "$SIEVE_BIN" ]]; then
    echo "[test_sieve] binary missing, building: (cd $SIEVE_DIR && make serial)"
    ( cd "$SIEVE_DIR" && make serial ) || fail "sieve build (make serial) failed"
fi
[[ -x "$SIEVE_BIN" ]] || fail "sieve binary not present/executable after build: $SIEVE_BIN"
[[ -f "$FIXTURE"   ]] || fail "fixture not found: $FIXTURE"

INPUT_N="$(grep -c '^>' "$FIXTURE")"
echo "[test_sieve] fixture: $FIXTURE ($INPUT_N haplotypes)"
# Hard failure, not a note. MAX_PROMOTED below is a property of THIS fixture, so a fixture
# that grew or shrank silently changes what the bound means while the bound stays enforced.
(( INPUT_N == EXPECTED_INPUT_N )) || fail \
    "fixture has $INPUT_N haplotypes, expected $EXPECTED_INPUT_N. MAX_PROMOTED=$MAX_PROMOTED is \
calibrated against that fixture: if the fixture changed deliberately, re-derive the bound and \
update EXPECTED_INPUT_N and MAX_PROMOTED together."

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# --- helper: run the sieve once, echo the promoted count -------------------
# Counts headers in the emitted PREFIX.rank*.fasta (concat = promoted subset),
# which is the authoritative promoted set the pipeline actually consumes.
run_sieve() {
    local tag="$1" prefix="$WORK/$1"
    local log="$WORK/$1.log"
    "$SIEVE_BIN" --input "$FIXTURE" \
        --min-len "$MIN_LEN" --max-n-frac "$MAX_N_FRAC" \
        --threshold "$THRESHOLD" --warmup "$WARMUP" \
        --emit-promoted "$prefix" >"$log" 2>&1
    local rc=$?
    [[ $rc -eq 0 ]] || { echo "[test_sieve] run '$tag' exit=$rc"; cat "$log" >&2; return 1; }
    # Count promoted sequences from the emitted shards.
    #
    # `grep -c` exits 1 when the count is 0, and it used to be the last command here, so a
    # sieve that ran fine and promoted nothing was reported as "did not exit 0" and the
    # dedicated `N1 >= 1` assertion below could never run. An unexpanded glob (if
    # --emit-promoted ever changes its shard naming) hits the same path. Keep the count and
    # the exit status separate.
    local n
    n="$(cat "$prefix".rank*.fasta 2>/dev/null | grep -c '^>')" || true
    printf '%s\n' "${n:-0}"
    return 0
}

# --- 2/3. run twice; capture counts and exit status ------------------------
N1="$(run_sieve run1)" || fail "sieve run 1 did not exit 0"
N2="$(run_sieve run2)" || fail "sieve run 2 did not exit 0"

echo "[test_sieve] promoted: run1=$N1  run2=$N2  (of $INPUT_N input)"

# --- assertions ------------------------------------------------------------
# (a) ran: covered above (run_sieve fails on non-zero exit).
# (b) 1 <= promoted <= input
[[ "$N1" =~ ^[0-9]+$ ]] || fail "run1 promoted count not numeric: '$N1'"
[[ "$N2" =~ ^[0-9]+$ ]] || fail "run2 promoted count not numeric: '$N2'"
(( N1 >= 1 ))          || fail "run1 promoted $N1 < 1 (nothing promoted)"
(( N1 <= INPUT_N ))    || fail "run1 promoted $N1 > input $INPUT_N (impossible)"

# (c) determinism: same count across the two runs
(( N1 == N2 ))         || fail "non-deterministic: run1=$N1 != run2=$N2"

# (d) intra-host single-SNV alts must collapse to a SMALL set (see header)
(( N1 <= MAX_PROMOTED )) \
    || fail "promoted $N1 > $MAX_PROMOTED: intra-host single-SNV alts should collapse to ~1 representative"

echo "PASS: sieve ran, promoted $N1/$INPUT_N (1<=n<=$INPUT_N), deterministic, collapsed (<=$MAX_PROMOTED)"
exit 0
