"""sieve: thin near-duplicate haplotypes with the k-mer diversity sieve.

Wraps `extern/kmer-diversity-sieve/l1_pipeline_serial`. The sieve is a MinHash/LSH
lineage-level reducer: on defaults tuned to gene scale (tau=0.005, sketch=256) it
collapses a cohort of ~identical lineages while preserving topology; the sieve
binary's own defaults (tau=0.15, sketch=64) target archive-scale collapse and reduce
a gene-scale cohort to ~1 representative. See docs/issue_sieve_calibration.md.

Lifted from the sieve block in the old `python/disassemble.py`; behaviour is
unchanged. Both the CLI and the composite pipelines call `sieve_run(...)`.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from disassembler.resolve import resolve_sieve
from disassembler.validators import fraction, positive_int


def first_record_len(path: Path):
    """Residues in the FIRST record of a FASTA, plus the number of records.

    The walker aligns only record 1, so sizing --min-len against every record
    would turn a multi-record query into a rejection floor above every reconstruction.
    """
    n = 0
    records = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                records += 1
                continue
            if records == 1:
                n += len(line.strip())
    return n, records


def sieve_run(
    input_fasta: Path,
    output_fasta: Path,
    *,
    ref: Path | None = None,
    min_len: int | None = None,
    threshold: float = 0.005,
    sketch: int = 256,
    sieve_bin: str | None = None,
) -> int:
    """Run the sieve on `input_fasta`, concatenate rank shards to `output_fasta`.

    If `min_len` is None: derived as int(0.95 * first_record_len(ref)). One of
    `min_len` or `ref` must be given.

    Returns process exit code (0 on success).
    """
    resolved = resolve_sieve(sieve_bin)
    if not resolved:
        print(
            "sieve: kmer-diversity-sieve not found. Build it "
            "(`cd extern/kmer-diversity-sieve && make serial`), put it on PATH, "
            "or pass --sieve-bin.",
            file=sys.stderr,
        )
        return 127

    if min_len is None:
        if ref is None:
            print("sieve: pass --min-len, or --ref to derive it from the query.",
                  file=sys.stderr)
            return 2
        cds_len, n_records = first_record_len(ref)
        if n_records > 1:
            print(f"[sieve] warn: {ref} has {n_records} records; sizing --min-len on "
                  f"the first one ({cds_len} nt), which is the record the walker "
                  f"aligns. One gene per query.", file=sys.stderr)
        if cds_len == 0:
            print(f"sieve: {ref}: first record has no sequence; cannot derive --min-len.",
                  file=sys.stderr)
            return 2
        # Reconstructions are typically a few bases shorter than the reference (the walk
        # may not reach the exact termini). Defaulting to the exact query length would
        # then drop every haplotype, so the auto-default takes a 5% haircut.
        min_len = int(cds_len * 0.95)

    input_fasta = Path(input_fasta)
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    emit_prefix = output_fasta.with_suffix("")
    # Clear any stale shards from a previous run so the concat glob can't pick them up.
    for stale in output_fasta.parent.glob(f"{emit_prefix.name}.rank*.fasta"):
        stale.unlink()

    cmd = [
        resolved,
        "--input", str(input_fasta),
        "--min-len", str(min_len),
        "--max-n-frac", "0.05",
        "--threshold", str(threshold),
        "--warmup", "2000",
        "--sketch", str(sketch),
        "--emit-promoted", str(emit_prefix),
    ]
    print(f"[sieve] {' '.join(cmd)}", flush=True)
    rc = subprocess.run(cmd).returncode
    if rc != 0:
        print(f"sieve: kmer-diversity-sieve exited {rc}", file=sys.stderr)
        return rc

    rank_files = sorted(output_fasta.parent.glob(f"{emit_prefix.name}.rank*.fasta"))
    with output_fasta.open("w") as fh:
        for rf in rank_files:
            fh.write(rf.read_text())

    n_in = sum(1 for line in input_fasta.open() if line.startswith(">"))
    n_out = sum(1 for line in output_fasta.open() if line.startswith(">"))
    print(f"[sieve] {n_in} -> {n_out} haplotypes "
          f"(tau={threshold}, sketch={sketch}) -> {output_fasta}", flush=True)
    if n_out == 0:
        print("sieve: promoted 0 haplotypes; try a lower --threshold or a lower "
              "--min-len (reconstructions may be shorter than the reference).",
              file=sys.stderr)
        return 1
    return 0


def add_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "sieve",
        help="Thin near-duplicate haplotypes with the k-mer diversity sieve.",
        description=(
            "MinHash/LSH lineage-level reducer. Use only for LARGE cross-accession "
            "cohorts of near-duplicate lineages, NOT to preserve intra-host SNVs "
            "(the walker's --emit-vcf / --emit-snv-table do that)."
        ),
    )
    p.add_argument("-i", "--input", required=True, help="Input FASTA")
    p.add_argument("-o", "--output", required=True, help="Output (promoted) FASTA")
    p.add_argument("--min-len", type=positive_int, default=None,
                   help="Minimum sequence length gate. Default: auto = 95%% of the "
                        "first record in --ref.")
    p.add_argument("--ref", default=None,
                   help="Reference CDS FASTA, used only to auto-derive --min-len "
                        "when --min-len is not given.")
    p.add_argument("--threshold", type=fraction, default=0.005,
                   help="LSH promotion threshold tau (default 0.005, gene-scale). "
                        "The sieve binary's own default of 0.15 targets archive-scale.")
    p.add_argument("--sketch", type=positive_int, default=256,
                   help="MinHash sketch size m (default 256, gene-scale). The "
                        "estimator cannot resolve novelty finer than 1/m.")
    p.add_argument("--sieve-bin", default=None,
                   help="Path to l1_pipeline_serial (auto-discovered otherwise).")
    p.set_defaults(func=run)


def run(args) -> int:
    return sieve_run(
        Path(args.input),
        Path(args.output),
        ref=Path(args.ref) if args.ref else None,
        min_len=args.min_len,
        threshold=args.threshold,
        sketch=args.sketch,
        sieve_bin=args.sieve_bin,
    )
