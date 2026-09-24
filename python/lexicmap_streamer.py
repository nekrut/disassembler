#!/usr/bin/env python3
"""lexicmap_streamer.py: Streaming multi-HSP coordinate tiling, codon QC, and haplotype collapsing.

Streams LexicMap tabular search results directly from usegalaxy.org (over HTTP API),
a local TSV file, or standard input (stdin).

Core capabilities:
  1. Multi-HSP coordinate tiling:
     Solves the Logan contig-shattering problem by grouping HSPs by SRA accession
     (`sgenome`), projecting alignment blocks onto the reference CDS coordinate system,
     and resolving overlapping positions by highest percent identity. LexicMap is indexed
     on Logan's *contigs* (`s3://logan-pub/c/`), not its unitig graph. On a diverse sample
     the assembler shatters the gene at every bubble it cannot cleanly pop -- on Zika NS5
     the gene tiles across ~50 short contigs, longest under 100 nt -- so LexicMap emits
     one HSP per surviving contig. See docs/contigs_vs_unitigs.md.
  2. Reading-frame codon QC & Stop-codon screening:
     Translates in-frame CDS (NCBI Code 11 or Standard), detects premature internal stop
     codons (e.g. Sanger am3 amber mutants, sequencing errors, or out-of-frame projections),
     and audits exact codon coordinates.
  3. Four-tier cohort partitioning:
     - Tier 1: Clean Full-Length (coverage >= min_coverage, 0 premature stops)
     - Tier 2: Clean Partial (min_coverage_partial <= coverage < min_coverage, 0 stops)
     - Tier 3: Flagged Premature Stops (stops > max_stops)
     - Tier 4: Low-Coverage Fragments (coverage < min_coverage_partial)
  4. Online haplotype collapser:
     Collapses clean full-length sequences without holding the input, calls nucleotide and
     amino acid mutations, classifies mutation impact (WT, WT_partial, Synonymous, Missense,
     Nonsense), and outputs clean MSAs, mutation frequency tables, and complete cohort ledgers.

     Memory is O(distinct clean full-length haplotypes), NOT O(1): `clean_full_haps`,
     `clean_full_pident_sum` and `clean_full_samples` are keyed by the full merged sequence,
     so on a divergent target where few accessions share a haplotype this approaches one
     entry per accession, each holding a copy of the sequence. It is the ROW count that
     streams, which is the thing that is unbounded; the haplotype count is not.
  5. Regime diagnostics:
     Calculates coverage and identity quantiles (p10, p25, median, p75, p90), detects
     divergent targets (median pid < 85% or cov < 0.60), and alerts when de novo
     graph path-walking with `logan-walker` is required.

Outputs generated:
  - <out-prefix>.clean.msa.fasta             (collapsed clean haplotypes with rank and count)
  - <out-prefix>.clean.haplotypes.tsv        (detailed mutation audit, frequencies, sample accessions)
  - <out-prefix>.clean.accessions.fasta      (uncollapsed full-length clean accessions)
  - <out-prefix>.clean_expanded.accessions.fasta (clean cohort including >=50% partials)
  - <out-prefix>.flagged.accessions.fasta    (flagged sequences: stops, frameshifts, or low cov)
  - <out-prefix>.flagged.tsv                 (audit table of all flagged accessions and reasons)
  - <out-prefix>.cohort_ledger.tsv           (complete audit ledger of all evaluated accessions)
  - <out-prefix>.summary.json                (quantiles, cohort breakdown, multi-HSP rescue counts)
  - <out-prefix>.clean.a2m                   (A2M formatted MSA for protein language models / EVcouplings)

Examples:
  # Stream directly from a Galaxy dataset ID over HTTP:
  python/lexicmap_streamer.py \\
      -q queries/markers/zika_ns5.fasta \\
      -d f9cad7b01a47213500416a19986dca5c \\
      -o out/zika_lexicmap

  # Stream from a local LexicMap tabular file:
  python/lexicmap_streamer.py \\
      -q queries/query_hiv1_pol.fasta \\
      -f /data/hiv1_pol_lexicmap.tsv \\
      -o out/hiv1_pol

  # Pipe directly from a compressed file or stream via STDIN:
  zcat /data/lexicmap_output.tsv.gz | python/lexicmap_streamer.py \\
      -q queries/master/pfcrt.fasta \\
      -f - \\
      -o out/pfcrt
"""

import os
import sys
import gzip
import math
import json
import time
import argparse
import urllib.request
from pathlib import Path
from collections import Counter, defaultdict

# Translation table: NCBI standard / bacterial code 11
CODON_TABLE_11 = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def get_galaxy_api_key(api_key: str = None, key_file: str = None) -> str:
    """Resolves Galaxy API key from CLI, environment, or default credentials file."""
    if api_key:
        return api_key.strip()
    if os.environ.get("GALAXY_API_KEY"):
        return os.environ["GALAXY_API_KEY"].strip()
    if key_file and Path(key_file).exists():
        return Path(key_file).read_text().strip()
    default_tmp = Path("/tmp/gxy.txt")
    if default_tmp.exists():
        return default_tmp.read_text().strip()
    raise RuntimeError(
        "Galaxy API key not found. Provide --api-key, set GALAXY_API_KEY, or write /tmp/gxy.txt."
    )

def translate_nt(seq: str, codon_table: dict = CODON_TABLE_11) -> str:
    """Translates nucleotide sequence to amino acids."""
    seq = seq.upper()
    aa_list = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if '-' in codon or len(codon) < 3:
            aa_list.append('-')
        else:
            aa_list.append(codon_table.get(codon, 'X'))
    return "".join(aa_list)

def find_internal_stops(seq: str, ref_len: int, codon_table: dict = CODON_TABLE_11) -> list:
    """
    Finds 1-based codon positions of premature stop codons in coding sequence body.
    Excludes the terminal natural stop codon (last 3 nt) of the CDS.
    """
    stop_positions = []
    for i in range(0, ref_len - 3, 3):
        codon = seq[i:i + 3]
        if '-' not in codon and len(codon) == 3:
            if codon_table.get(codon, 'X') == '*':
                stop_positions.append(i // 3 + 1)
    return stop_positions

def compute_quantiles(values: list) -> dict:
    """Calculates summary percentiles for a metric."""
    if not values:
        return {"p10": 0.0, "p25": 0.0, "p50_median": 0.0, "p75": 0.0, "p90": 0.0, "mean": 0.0}
    s = sorted(values)
    n = len(s)
    def q(p):
        idx = int(round(p * (n - 1)))
        return round(s[max(0, min(n - 1, idx))], 4)
    return {
        "p10": q(0.10),
        "p25": q(0.25),
        "p50_median": q(0.50),
        "p75": q(0.75),
        "p90": q(0.90),
        "mean": round(sum(s) / n, 4)
    }

# How many malformed rows to name individually before falling back to a count.
MALFORMED_REPORT_LIMIT = 5


def load_reference(fasta_path: str):
    """Loads the FIRST record's sequence, length, and header from a FASTA file.

    Records after the first are ignored, with a warning. Concatenating them instead built a
    pseudo-gene under record 1's name: the 160-record master panel came back as
    ref_len=313311 called 'h5n1_pb2_reference', which sent every accession to tier 4 at a
    coverage around 0.009 and discarded the whole cohort at exit 0.
    """
    path = Path(fasta_path)
    # is_file, not exists: a directory passed `exists()` and then raised IsADirectoryError
    # from open() instead of the intended FileNotFoundError.
    if not path.is_file():
        raise FileNotFoundError(f"Reference FASTA file not found (or not a file): {path}")

    header = ""
    seq_parts = []
    n_records = 0
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                n_records += 1
                if n_records == 1:
                    header = line[1:]
                continue
            if n_records == 1:
                seq_parts.append(line)

    if n_records > 1:
        print(f"[warn] {path}: {n_records} records present, using only the first "
              f"({header.split()[0] if header else path.stem!r}). One gene per reference.",
              file=sys.stderr)

    ref_seq = "".join(seq_parts).upper()
    # Every coverage figure divides by len(ref_seq). An empty or header-only FASTA sailed
    # through to `cov = covered_bases / ref_len` and died there with a bare
    # ZeroDivisionError, naming neither the file nor the reason.
    if not ref_seq:
        raise ValueError(
            f"Reference FASTA has no sequence: {path} "
            f"({n_records} record(s) read, first header {header[:60]!r})"
        )
    ref_name = header.split()[0] if header else path.stem
    return ref_name, header, ref_seq

def get_galaxy_stream(dataset_id: str, base_url: str = "https://usegalaxy.org", api_key: str = None, key_file: str = None):
    """Opens a buffered HTTP streaming connection to Galaxy dataset display endpoint."""
    key = get_galaxy_api_key(api_key, key_file)
    url = f"{base_url.rstrip('/')}/api/datasets/{dataset_id}/display"
    req = urllib.request.Request(url, headers={"x-api-key": key})
    response = urllib.request.urlopen(req, timeout=120)
    return response

def resolve_dataset_state(dataset_id: str, base_url: str = "https://usegalaxy.org", api_key: str = None, key_file: str = None):
    """Checks the live state of a Galaxy dataset."""
    try:
        key = get_galaxy_api_key(api_key, key_file)
        url = f"{base_url.rstrip('/')}/api/datasets/{dataset_id}"
        req = urllib.request.Request(url, headers={"x-api-key": key})
        with urllib.request.urlopen(req, timeout=30) as res:
            data = json.loads(res.read())
            return data.get("state")
    except Exception:
        return None

def merge_accession_hsps(hsps: list, ref_len: int, min_pident: float, min_coverage: float):
    """
    Merges all qualifying HSPs of an accession onto the reference coordinate system.
    Overlapping positions between HSPs are resolved by selecting the base from
    the HSP with higher pident.
    """
    grid = [None] * ref_len
    single_ge_min_cov = False
    single_full = False

    for pident, qcov, qstart, qend, qseq, sseq in hsps:
        if pident < min_pident:
            continue
        if qcov >= (min_coverage * 100.0):
            single_ge_min_cov = True
        if qstart == 1 and qend == ref_len and pident >= min_pident:
            single_full = True

        q_pos = qstart - 1
        for qc, sc in zip(qseq, sseq):
            if qc != '-':
                if 0 <= q_pos < ref_len:
                    ex = grid[q_pos]
                    if ex is None or pident > ex[1]:
                        grid[q_pos] = (sc.upper(), pident)
                q_pos += 1

    # Mean identity over the positions that actually carry a subject base. The numerator must
    # use the same non-gap positions as `covered_bases`; summing pident over gap positions too
    # (subject deletions, where x[0] == '-') inflates the mean and can push it above 100%.
    covered = [x for x in grid if x is not None and x[0] != '-']
    covered_bases = len(covered)
    mean_pid = (sum(x[1] for x in covered) / covered_bases) if covered_bases > 0 else 0.0
    merged_seq = "".join(grid[i][0] if grid[i] is not None else '-' for i in range(ref_len))

    return merged_seq, covered_bases, mean_pid, single_ge_min_cov, single_full

def stream_lexicmap_to_msa(
    ref_path: str,
    dataset_id: str = None,
    file_path: str = None,
    base_url: str = "https://usegalaxy.org",
    api_key: str = None,
    api_key_file: str = "/tmp/gxy.txt",
    min_coverage: float = 0.80,
    min_coverage_partial: float = 0.50,
    min_pident: float = 60.0,
    max_internal_stops: int = 0,
    allow_frameshifts: bool = False,
    full_length_only: bool = False,
    max_records: int = None,
    sample_cap: int = 10,
    out_prefix: str = None,
    write_a2m: bool = True
) -> dict:
    """
    Main streaming pipeline:
    - Streams rows, groups by accession, merges multi-HSP alignments
    - Partitions cohort into:
        1. Clean Full-Length (cov >= min_coverage, 0 internal stops)
        2. Clean Partial (min_coverage_partial <= cov < min_coverage, 0 stops)
        3. Flagged (internal stops, frameshifts, or low coverage)
    - Writes clean MSAs, flagged audit FASTA/TSV, and complete cohort ledger.
    """
    ref_name, ref_header, ref_seq = load_reference(ref_path)
    ref_len = len(ref_seq)
    ref_is_cds = (ref_len % 3 == 0)
    ref_aa = translate_nt(ref_seq) if ref_is_cds else ""

    if dataset_id:
        ds_state = resolve_dataset_state(dataset_id, base_url, api_key, api_key_file)
        if ds_state and ds_state != "ok":
            print(f"[!] Warning: Galaxy dataset {dataset_id} is in state '{ds_state}' (not 'ok').")
            if ds_state == "queued":
                print("    Job is still queued on Galaxy. Skipping for now.")
                return {"status": "queued", "dataset_id": dataset_id}

    print("=" * 80)
    print(" LexicMapStreamer: Multi-HSP Coordinate Tiler & Cohort Partitioner")
    print("=" * 80)
    print(f" Reference Sequence : {ref_name} ({ref_len} bp, {len(ref_aa)} codons)")
    print(f" Reference Header   : {ref_header}")
    print(f" Target QC Settings : full_length_cov={min_coverage*100:.1f}%, partial_cov={min_coverage_partial*100:.1f}%, "
          f"min_pident={min_pident:.1f}%, max_stops={max_internal_stops if ref_is_cds else 'N/A'}")

    if dataset_id:
        print(f" Input Source       : Streaming Galaxy dataset: {dataset_id} ({base_url})")
        stream = get_galaxy_stream(dataset_id, base_url, api_key, api_key_file)
        is_http = True
    elif file_path == "-":
        print(" Input Source       : Reading from STDIN")
        stream = sys.stdin
        is_http = False
    elif file_path:
        print(f" Input Source       : Local file: {file_path}")
        if file_path.endswith(".gz"):
            stream = gzip.open(file_path, "rt", encoding="utf-8")
        else:
            stream = open(file_path, "r", encoding="utf-8")
        is_http = False
    else:
        raise ValueError("Must specify either --dataset or --file.")

    if not out_prefix:
        out_dir = Path("msa")
        out_dir.mkdir(parents=True, exist_ok=True)
        out_prefix = f"msa/{ref_name}"
    else:
        out_path = Path(out_prefix)
        if out_path.is_dir() or str(out_prefix).endswith("/"):
            out_path.mkdir(parents=True, exist_ok=True)
            out_prefix = str(out_path / ref_name)
        else:
            out_path.parent.mkdir(parents=True, exist_ok=True)

    # Read header line (skipping comment lines if any)
    while True:
        header_line = stream.readline()
        if not header_line:
            raise ValueError("Input stream was empty or did not contain a header line.")
        if is_http and isinstance(header_line, bytes):
            header_line = header_line.decode("utf-8")
        header_line = header_line.strip()
        if header_line and not (header_line.startswith("#") and "\t" not in header_line):
            break

    # Strip leading '#' if present on first column
    if header_line.startswith("#"):
        header_line = header_line[1:].strip()

    header_cols = header_line.split("\t")
    col_idx = {c.strip(): i for i, c in enumerate(header_cols)}

    # Schema resolution
    def find_col(candidates):
        for c in candidates:
            if c in col_idx:
                return col_idx[c]
        # Case-insensitive fallback
        lower_map = {k.lower(): v for k, v in col_idx.items()}
        for c in candidates:
            if c.lower() in lower_map:
                return lower_map[c.lower()]
        return None

    idx_sgenome = find_col(["sgenome", "subject", "sseqid", "acc", "accession"])
    idx_pident = find_col(["pident", "pid", "perc_identity", "identity"])
    idx_qcov = find_col(["qcovHSP", "qcov", "qcovhsp", "q_cov"])
    idx_qstart = find_col(["qstart", "q_start", "query_start"])
    idx_qend = find_col(["qend", "q_end", "query_end"])
    idx_qseq = find_col(["qseq", "query_seq", "q_seq"])
    idx_sseq = find_col(["sseq", "subject_seq", "s_seq"])
    # Neither of these was resolved at all, so a multi-target search was silently merged onto
    # one reference and a mismatched reference was never noticed. Both are OPTIONAL: older or
    # hand-made TSVs may not carry them, and requiring them would reject working inputs.
    idx_query = find_col(["query", "qseqid", "qaccver", "qname", "query_id"])
    idx_qlen = find_col(["qlen", "query_length", "qseqlen", "q_len"])

    missing = []
    if idx_sgenome is None: missing.append("sgenome/subject")
    if idx_pident is None: missing.append("pident")
    if idx_qcov is None: missing.append("qcovHSP")
    if idx_qstart is None: missing.append("qstart")
    if idx_qend is None: missing.append("qend")
    if idx_qseq is None: missing.append("qseq")
    if idx_sseq is None: missing.append("sseq")

    if missing:
        raise KeyError(f"Missing required columns {missing} in LexicMap header: {header_cols}")

    if idx_query is None:
        print(f"[warn] no query-name column in the header ({header_cols}); rows cannot be "
              f"filtered to '{ref_name}'. If this file holds hits for more than one target, "
              f"their HSPs will be merged onto this one reference.", file=sys.stderr)
    if idx_qlen is None:
        print(f"[*] no qlen column; skipping the reference length agreement check.")

    # Cohort file handles
    f_clean_fa = open(f"{out_prefix}.clean.accessions.fasta", "w", encoding="utf-8")
    f_clean_exp_fa = open(f"{out_prefix}.clean_expanded.accessions.fasta", "w", encoding="utf-8")
    f_flagged_fa = open(f"{out_prefix}.flagged.accessions.fasta", "w", encoding="utf-8")
    f_flagged_tsv = open(f"{out_prefix}.flagged.tsv", "w", encoding="utf-8")
    f_cohort_tsv = open(f"{out_prefix}.cohort_ledger.tsv", "w", encoding="utf-8")

    f_flagged_tsv.write("accession\tcoverage\tmean_pident\tflag\tn_stops\tstop_codons\thsps_merged\n")
    f_cohort_tsv.write("accession\tcoverage\tmean_pident\tstatus\ttier\tn_stops\tstop_codons\thsps_merged\n")

    # Tracking & State
    clean_full_haps = Counter()
    clean_full_pident_sum = defaultdict(float)
    clean_full_samples = defaultdict(list)

    all_covs = []
    all_pids = []
    cov_ge_50_count = 0
    cov_ge_50_stops = 0

    total_rows = 0
    n_malformed = 0
    n_other_query = 0
    n_qlen_mismatch = 0
    total_accessions = 0
    count_clean_full = 0
    count_clean_partial = 0
    count_flagged_stops = 0
    count_flagged_low_cov = 0
    total_rescued_multihsp = 0
    total_single_hsp_full = 0

    start_time = time.time()
    last_report_time = start_time

    print("\n---> Streaming rows, grouping by accession, and partitioning into Clean & Flagged cohorts...")

    curr_acc = None
    curr_hsps = []
    # The stream is grouped by accession one block at a time, so all of an accession's HSPs must
    # be contiguous. If they are not, an accession is split into several blocks: it is counted
    # more than once and its multi-HSP coverage is undercounted, silently. Rather than corrupt
    # the cohort, detect a reappearance and stop with a clear instruction to sort the input.
    seen_accessions = set()

    def process_accession_block(acc, hsps):
        nonlocal total_accessions, count_clean_full, count_clean_partial
        nonlocal count_flagged_stops, count_flagged_low_cov
        nonlocal total_rescued_multihsp, total_single_hsp_full
        nonlocal cov_ge_50_count, cov_ge_50_stops

        if not acc or not hsps:
            return

        if acc in seen_accessions:
            raise ValueError(
                f"Accession '{acc}' reappears after its block was already processed, so the "
                f"input is not grouped by accession. All HSPs of an accession must be "
                f"contiguous. Sort the LexicMap output by the accession column first, e.g. "
                f"`sort -k<sgenome_col> -s input.tsv`, then re-run."
            )
        seen_accessions.add(acc)

        total_accessions += 1

        merged_seq, covered_bases, mean_pid, single_ge_cov, single_full = merge_accession_hsps(
            hsps, ref_len, min_pident, min_coverage
        )

        cov = covered_bases / ref_len
        all_covs.append(cov)
        if cov > 0:
            all_pids.append(mean_pid)

        # Check premature stop codons in coding body
        stop_codons = find_internal_stops(merged_seq, ref_len) if ref_is_cds else []
        n_stops = len(stop_codons)
        stop_str = ",".join(str(p) for p in stop_codons) if stop_codons else "None"

        if cov >= 0.50:
            cov_ge_50_count += 1
            if n_stops > 0:
                cov_ge_50_stops += 1

        # Stratify accession. With --full-length-only, the clean full-length tier requires
        # complete coverage of the reference (every position tiled), not merely min_coverage.
        full_cov = 1.0 if full_length_only else min_coverage
        if cov >= full_cov and (n_stops <= max_internal_stops or allow_frameshifts):
            status = "CLEAN_FULL_LENGTH"
            tier = 1
            count_clean_full += 1
            if single_full:
                total_single_hsp_full += 1
            if not single_ge_cov:
                total_rescued_multihsp += 1

            clean_full_haps[merged_seq] += 1
            clean_full_pident_sum[merged_seq] += mean_pid
            samples = clean_full_samples[merged_seq]
            if len(samples) < sample_cap and acc not in samples:
                samples.append(acc)

            f_clean_fa.write(f">{acc} cov={cov:.4f} pid={mean_pid:.2f} status=CLEAN_FULL_LENGTH\n{merged_seq}\n")
            f_clean_exp_fa.write(f">{acc} cov={cov:.4f} pid={mean_pid:.2f} status=CLEAN_FULL_LENGTH\n{merged_seq}\n")

        elif cov >= min_coverage_partial and (n_stops <= max_internal_stops or allow_frameshifts):
            status = "CLEAN_PARTIAL"
            tier = 2
            count_clean_partial += 1
            f_clean_exp_fa.write(f">{acc} cov={cov:.4f} pid={mean_pid:.2f} status=CLEAN_PARTIAL\n{merged_seq}\n")

        elif n_stops > max_internal_stops and not allow_frameshifts:
            status = "FLAGGED_PREMATURE_STOPS"
            tier = 3
            count_flagged_stops += 1
            f_flagged_fa.write(f">{acc} cov={cov:.4f} pid={mean_pid:.2f} flag=PREMATURE_STOPS stops={n_stops} codons={stop_str}\n{merged_seq}\n")
            f_flagged_tsv.write(f"{acc}\t{cov:.4f}\t{mean_pid:.2f}\tPREMATURE_STOPS\t{n_stops}\t{stop_str}\t{len(hsps)}\n")

        else:  # cov < min_coverage_partial and no stops
            status = "LOW_COVERAGE_FRAGMENT"
            tier = 4
            count_flagged_low_cov += 1
            f_flagged_fa.write(f">{acc} cov={cov:.4f} pid={mean_pid:.2f} flag=LOW_COVERAGE\n{merged_seq}\n")
            f_flagged_tsv.write(f"{acc}\t{cov:.4f}\t{mean_pid:.2f}\tLOW_COVERAGE\t0\tNone\t{len(hsps)}\n")

        f_cohort_tsv.write(f"{acc}\t{cov:.4f}\t{mean_pid:.2f}\t{status}\t{tier}\t{n_stops}\t{stop_str}\t{len(hsps)}\n")

    try:
        while True:
            line = stream.readline()
            if not line:
                break
            if is_http and isinstance(line, bytes):
                line = line.decode("utf-8")
            line = line.strip()
            if not line:
                continue

            total_rows += 1
            parts = line.split("\t")
            needed = [idx_sgenome, idx_pident, idx_qcov, idx_qstart, idx_qend, idx_qseq, idx_sseq]
            needed += [i for i in (idx_query, idx_qlen) if i is not None]
            if len(parts) <= max(needed):
                continue

            # A LexicMap search against a multi-target panel emits rows for every target in
            # one file, distinguished ONLY by this column. Grouping is by sgenome alone, so
            # without this an accession hit by several genes had HSPs from DIFFERENT genes
            # tiled onto a single grid, and a file sorted by query rather than accession
            # tripped the "reappears after its block was processed" ValueError instead.
            if idx_query is not None and parts[idx_query] != ref_name:
                n_other_query += 1
                continue

            # Every HSP is projected onto a grid of len(ref_seq). If the search used a
            # different version of the gene, every projection is wrong and nothing said so.
            if idx_qlen is not None:
                try:
                    row_qlen = int(parts[idx_qlen])
                except ValueError:
                    row_qlen = None
                if row_qlen is not None and row_qlen != ref_len:
                    n_qlen_mismatch += 1
                    if n_qlen_mismatch == 1:
                        raise ValueError(
                            f"reference length disagreement: the search reports qlen="
                            f"{row_qlen} for query '{parts[idx_query] if idx_query is not None else ref_name}', "
                            f"but {ref_path} is {ref_len} nt. Every HSP is projected onto a "
                            f"grid of the reference length, so this would mis-place all of "
                            f"them. Use the same reference FASTA the search used."
                        )

            acc = parts[idx_sgenome]
            # Bare float()/int() here aborted the whole stream on one bad field, AFTER the
            # cohort files had been partly written and BEFORE .summary.json was written, so
            # the output directory looked half-plausible. Realistic triggers: concatenated
            # LexicMap outputs with an embedded second header row, `NA` or `-` in a numeric
            # column, a qseq containing a tab. Skip the row, count it, report at the end.
            try:
                pident = float(parts[idx_pident])
                qcov = float(parts[idx_qcov])
                qstart = int(parts[idx_qstart])
                qend = int(parts[idx_qend])
            except ValueError:
                n_malformed += 1
                if n_malformed <= MALFORMED_REPORT_LIMIT:
                    print(f"[warn] row {total_rows}: unparseable numeric field, skipping: "
                          f"{line[:120]}", file=sys.stderr)
                continue
            qseq = parts[idx_qseq]
            sseq = parts[idx_sseq]

            hsp = (pident, qcov, qstart, qend, qseq, sseq)

            if acc != curr_acc:
                if curr_acc is not None:
                    process_accession_block(curr_acc, curr_hsps)
                curr_acc = acc
                curr_hsps = [hsp]
            else:
                curr_hsps.append(hsp)

            now = time.time()
            if now - last_report_time >= 5.0 or (max_records and total_rows % 25000 == 0):
                elapsed = now - start_time
                rate = total_rows / elapsed if elapsed > 0 else 0
                print(f"     [Stream] Rows: {total_rows:,} | Accs: {total_accessions:,} | "
                      f"Clean Full: {count_clean_full:,} | Rescued: {total_rescued_multihsp:,} | "
                      f"Flagged Stops: {count_flagged_stops:,} | Speed: {rate:,.0f} rows/s")
                last_report_time = now

            if max_records and total_rows >= max_records:
                print(f"\n[*] Reached row limit (--max-records {max_records:,}).")
                break

        if curr_acc is not None:
            process_accession_block(curr_acc, curr_hsps)

        if n_other_query:
            print(f"[*] {n_other_query:,} row(s) belonged to another query and were skipped; "
                  f"only '{ref_name}' was tiled.")
            # Everything filtered out means the reference name does not match the search's
            # query name, not that there were no hits. An empty cohort would look plausible.
            if n_other_query == total_rows:
                raise ValueError(
                    f"no row matched query '{ref_name}': every one of {total_rows:,} rows "
                    f"names a different query. The reference FASTA's record name must match "
                    f"the query column of the search."
                )

        if n_malformed:
            print(f"[!] {n_malformed:,} of {total_rows:,} rows had an unparseable numeric "
                  f"field and were skipped"
                  f"{' (only the first %d were named)' % MALFORMED_REPORT_LIMIT if n_malformed > MALFORMED_REPORT_LIMIT else ''}.",
                  file=sys.stderr)
            # Every row unusable is a schema mismatch, not dirty data, and the cohort below
            # would be empty rather than incomplete. Fail rather than write a plausible-looking
            # empty output directory.
            if n_malformed == total_rows:
                raise ValueError(
                    f"all {total_rows:,} data rows had unparseable numeric fields; "
                    f"check the column layout of the input against the resolved schema above"
                )

    finally:
        if hasattr(stream, "close"):
            stream.close()
        f_clean_fa.close()
        f_clean_exp_fa.close()
        f_flagged_fa.close()
        f_flagged_tsv.close()
        f_cohort_tsv.close()

    total_time = time.time() - start_time

    # Distributions & Diagnostics
    cov_stats = compute_quantiles(all_covs)
    pid_stats = compute_quantiles(all_pids)
    med_cov = cov_stats["p50_median"]
    med_pid = pid_stats["p50_median"]
    stop_rate_ge_50 = (cov_ge_50_stops / cov_ge_50_count * 100.0) if cov_ge_50_count > 0 else 0.0

    is_divergent = (med_pid < 85.0 or med_cov < 0.60)
    regime = "divergent" if is_divergent else "conserved"

    # A valid stream can carry a header and no data rows; guard every by-total percentage so an
    # empty cohort reports zeros rather than raising ZeroDivisionError.
    def pct_of_total(n):
        return (n / total_accessions * 100.0) if total_accessions > 0 else 0.0

    print(f"\n[✓] Stream processing complete in {total_time:.2f} seconds.")
    print(f"    Total rows streamed           : {total_rows:,}")
    print(f"    Total unique SRA accessions   : {total_accessions:,}")
    print(f"    1. Clean Full-Length Cohort   : {count_clean_full:,} accessions ({pct_of_total(count_clean_full):.2f}% of total)")
    print(f"       - Multi-HSP rescued        : {total_rescued_multihsp:,} accessions")
    print(f"       - Single-HSP full-length   : {total_single_hsp_full:,} accessions")
    print(f"    2. Clean Partials (50%-80% cov): {count_clean_partial:,} accessions")
    print(f"       --> Combined Clean Cohort  : {count_clean_full + count_clean_partial:,} accessions ({pct_of_total(count_clean_full + count_clean_partial):.2f}%)")
    print(f"    3. Flagged for Premature Stops: {count_flagged_stops:,} accessions")
    print(f"    4. Low-Coverage Fragments (<50%): {count_flagged_low_cov:,} accessions")
    print(f"    Coverage Distribution (All)   : p10={cov_stats['p10']:.2f}, median={med_cov:.2f}, p90={cov_stats['p90']:.2f}")
    print(f"    Identity Distribution (All)   : p10={pid_stats['p10']:.1f}%, median={med_pid:.1f}%, p90={pid_stats['p90']:.1f}%")

    warning_msg = None
    if is_divergent:
        warning_msg = (
            f"DIVERGENT REGIME DETECTED: Median identity is {med_pid:.1f}% (< 85%) and/or "
            f"median coverage is {med_cov:.2f} (< 0.60) with {stop_rate_ge_50:.1f}% premature stops in >=50% covered accessions. "
            f"Reference projection is fragmentary. Clean full-length sequences ({count_clean_full:,}) and clean partials "
            f"({count_clean_partial:,}) have been safely partitioned from flagged stop-codon projections ({count_flagged_stops:,}). "
            f"For unbiased whole-cohort haplotype assembly on this target, use disassembler (logan-walker)."
        )
        print("\n" + "!" * 80)
        print(" [!] DIVERGENT TARGET REGIME DETECTED")
        print("!" * 80)
        print(warning_msg)
        print("!" * 80 + "\n")

    # Sort clean full-length haplotypes
    sorted_haplotypes = clean_full_haps.most_common()

    # Compute Shannon entropy
    shannon_entropy = 0.0
    for seq, count in sorted_haplotypes:
        p = count / count_clean_full if count_clean_full > 0 else 0
        if p > 0:
            shannon_entropy -= p * math.log2(p)

    msa_clean_fasta = Path(f"{out_prefix}.clean.msa.fasta")
    msa_compat_fasta = Path(f"{out_prefix}.msa.fasta")
    counts_tsv_path = Path(f"{out_prefix}.clean.haplotypes.tsv")
    counts_compat_tsv = Path(f"{out_prefix}.haplotypes.tsv")
    summary_json_path = Path(f"{out_prefix}.summary.json")

    print(f"\n---> Emitting clean alignments and cohort tables to '{out_prefix}'...")

    # Write clean MSA FASTA
    # One pass, two files. The compat alias is byte-identical, so it is written from the
    # same formatting rather than regenerated.
    with open(msa_clean_fasta, "w", encoding="utf-8") as f_fa, \
         open(msa_compat_fasta, "w", encoding="utf-8") as f_fa_compat:
        def _w(s):
            f_fa.write(s)
            f_fa_compat.write(s)

        _w(f">REFERENCE_{ref_name} len={ref_len} count={count_clean_full} type=reference\n{ref_seq}\n")
        for rank, (seq, count) in enumerate(sorted_haplotypes, 1):
            freq = count / count_clean_full if count_clean_full > 0 else 0
            is_wt = (seq == ref_seq)
            tag = "WT" if is_wt else f"variant_rank_{rank}"
            _w(f">haplotype_{rank:05d} count={count} freq={freq:.5f} tag={tag}\n{seq}\n")

    # Write A2M
    if write_a2m:
        with open(f"{out_prefix}.clean.a2m", "w", encoding="utf-8") as f_a2m, open(f"{out_prefix}.a2m", "w", encoding="utf-8") as f_a2m_compat:
            for f in [f_a2m, f_a2m_compat]:
                f.write(f">{ref_name}/1-{ref_len}\n{ref_seq}\n")
                for rank, (seq, count) in enumerate(sorted_haplotypes, 1):
                    f.write(f">hap_{rank}_count_{count}/1-{ref_len}\n{seq}\n")

    # Write Haplotype TSV
    top_preview = []
    # One computation pass, two files. The loop body translates the haplotype and diffs it
    # against the reference per position; running it again for the byte-identical compat
    # alias was measured at 30% of this stage's runtime on a 2,000-haplotype cohort.
    with open(counts_tsv_path, "w", encoding="utf-8") as f_tsv, \
         open(counts_compat_tsv, "w", encoding="utf-8") as f_tsv_compat:
        def _wt(s):
            f_tsv.write(s)
            f_tsv_compat.write(s)

        headers = [
            "haplotype_id", "rank", "count", "frequency", "mean_pident",
            "coverage_pct", "num_nt_muts", "nt_mutations",
            "num_aa_muts", "aa_mutations", "mutation_class",
            "representative_samples", "sequence"
        ]
        _wt("\t".join(headers) + "\n")

        for rank, (seq, count) in enumerate(sorted_haplotypes, 1):
            freq = count / count_clean_full if count_clean_full > 0 else 0
            mean_pid = clean_full_pident_sum[seq] / count if count > 0 else 0
            cov_pct = sum(1 for c in seq if c != '-') / ref_len * 100.0

            nt_muts = []
            for i, (r_char, h_char) in enumerate(zip(ref_seq, seq)):
                if h_char != '-' and h_char != r_char:
                    nt_muts.append(f"{i+1}{r_char}>{h_char}")

            aa_muts = []
            if seq == ref_seq:
                mut_class = "WT"
            elif ref_is_cds:
                h_aa = translate_nt(seq)
                for r_idx, (r_a, h_a) in enumerate(zip(ref_aa, h_aa)):
                    if h_a != '-' and h_a != r_a:
                        aa_muts.append(f"{r_a}{r_idx+1}{h_a}")

                if not aa_muts and not nt_muts:
                    # NOT "WT": `seq == ref_seq` was ruled out above, so this sequence
                    # matches the reference everywhere it is COVERED and differs only
                    # where it is not. Labelling it WT made this table disagree with
                    # summary.json's wildtype_count, which counts exact matches only.
                    mut_class = "WT_partial"
                elif not aa_muts:
                    mut_class = "Synonymous"
                elif any(a.endswith("*") for a in aa_muts):
                    mut_class = "Nonsense"
                else:
                    # Everything left has amino-acid changes and none is a stop. The old
                    # `elif aa_muts` / `else: "Complex"` pair could not reach Complex:
                    # the two arms above consume every `not aa_muts` case, so `aa_muts`
                    # is always truthy here.
                    mut_class = "Missense"
            elif not nt_muts:
                mut_class = "NonCoding_WT_partial"
            else:
                mut_class = "NonCoding_Variant"

            nt_mut_str = ",".join(nt_muts) if nt_muts else "WT"
            aa_mut_str = ",".join(aa_muts) if aa_muts else (
                "WT" if mut_class in ("WT", "WT_partial", "Synonymous") else "NA")
            samples_str = ",".join(clean_full_samples[seq])

            row_data = [
                f"haplotype_{rank:05d}",
                str(rank),
                str(count),
                f"{freq:.6f}",
                f"{mean_pid:.2f}",
                f"{cov_pct:.2f}",
                str(len(nt_muts)),
                nt_mut_str,
                str(len(aa_muts)),
                aa_mut_str,
                mut_class,
                samples_str,
                seq
            ]
            _wt("\t".join(row_data) + "\n")

            if rank <= 10:
                top_preview.append({
                    "rank": rank,
                    "count": count,
                    "freq": freq,
                    "mut_class": mut_class,
                    "nt_muts": nt_mut_str,
                    "aa_muts": aa_mut_str,
                    "samples": samples_str
                })

    # Summary JSON
    summary_data = {
        "status": "ok",
        "reference_name": ref_name,
        "reference_length_bp": ref_len,
        "is_cds": ref_is_cds,
        "regime": regime,
        "regime_warning": warning_msg,
        "total_streamed_rows": total_rows,
        "total_accessions_evaluated": total_accessions,
        "cohort_breakdown": {
            "clean_full_length": count_clean_full,
            "clean_partial": count_clean_partial,
            "combined_clean_cohort": count_clean_full + count_clean_partial,
            "flagged_premature_stops": count_flagged_stops,
            "flagged_low_coverage": count_flagged_low_cov
        },
        "recovery_metrics": {
            "multi_hsp_rescued_accessions": total_rescued_multihsp,
            "single_hsp_full_length_accessions": total_single_hsp_full,
            "premature_stop_rate_pct_at_50cov": round(stop_rate_ge_50, 2)
        },
        "coverage_quantiles": cov_stats,
        "pident_quantiles": pid_stats,
        "clean_unique_haplotypes": len(sorted_haplotypes),
        "shannon_entropy_bits": round(shannon_entropy, 4),
        "wildtype_count": clean_full_haps.get(ref_seq, 0),
        "wildtype_fraction": round(clean_full_haps.get(ref_seq, 0) / count_clean_full, 6) if count_clean_full > 0 else 0,
        "top_10_haplotypes": top_preview,
        "runtime_seconds": round(total_time, 2)
    }
    summary_json_path.write_text(json.dumps(summary_data, indent=2))

    print(f"     [✓] Clean MSA FASTA       : {msa_clean_fasta} ({len(sorted_haplotypes) + 1} sequences)")
    print(f"     [✓] Clean Full Accessions : {out_prefix}.clean.accessions.fasta ({count_clean_full:,} accessions)")
    print(f"     [✓] Clean Expanded (>=50%): {out_prefix}.clean_expanded.accessions.fasta ({count_clean_full + count_clean_partial:,} accessions)")
    print(f"     [✓] Flagged Accessions    : {out_prefix}.flagged.accessions.fasta ({count_flagged_stops + count_flagged_low_cov:,} accessions)")
    print(f"     [✓] Flagged Audit Table   : {out_prefix}.flagged.tsv")
    print(f"     [✓] Complete Cohort Ledger: {out_prefix}.cohort_ledger.tsv ({total_accessions:,} rows)")
    print(f"     [✓] Summary JSON          : {summary_json_path}")

    # Display Top Table
    print("\n" + "=" * 80)
    print(f" Top Clean Haplotypes ({ref_name})")
    print("=" * 80)
    print(f"{'Rank':<5} {'Count':<9} {'Freq':<8} {'Class':<12} {'AA Mutations':<18} {'NT Mutations':<25}")
    print("-" * 80)
    for row in top_preview:
        aa_disp = row['aa_muts'][:16] + ".." if len(row['aa_muts']) > 16 else row['aa_muts']
        nt_disp = row['nt_muts'][:23] + ".." if len(row['nt_muts']) > 23 else row['nt_muts']
        print(f"{row['rank']:<5} {row['count']:<9} {row['freq']:<8.4f} {row['mut_class']:<12} {aa_disp:<18} {nt_disp:<25}")
    print("=" * 80)

    return summary_data

def parse_args():
    parser = argparse.ArgumentParser(
        description="LexicMapStreamer: Stream LexicMap search results, tile multi-HSP alignments, stratify cohorts, and collapse haplotypes."
    )
    parser.add_argument("-q", "-r", "--query", "--ref", dest="ref", type=str, required=True,
                        help="Path to reference FASTA file (e.g. queries/query_hiv1_pol.fasta).")
    parser.add_argument("-d", "--dataset", type=str, default=None,
                        help="Galaxy dataset ID to stream from over HTTP.")
    parser.add_argument("-f", "-i", "--file", "--input", dest="file", type=str, default=None,
                        help="Path to local LexicMap tabular file (or '-' for stdin).")
    parser.add_argument("--base-url", type=str, default="https://usegalaxy.org",
                        help="Galaxy server base URL (default: https://usegalaxy.org).")
    parser.add_argument("--api-key", type=str, default=None,
                        help="Galaxy API key string.")
    parser.add_argument("--api-key-file", type=str, default="/tmp/gxy.txt",
                        help="Path to file containing Galaxy API key (default: /tmp/gxy.txt).")
    parser.add_argument("--min-coverage", type=float, default=0.80,
                        help="Minimum merged query coverage fraction for full-length clean cohort (default: 0.80 = 80%%).")
    parser.add_argument("--min-coverage-partial", type=float, default=0.50,
                        help="Minimum coverage fraction for clean partial expanded cohort (default: 0.50 = 50%%).")
    parser.add_argument("--min-pident", type=float, default=60.0,
                        help="Minimum percent identity for candidate HSPs (default: 60.0%%).")
    parser.add_argument("--max-internal-stops", type=int, default=0,
                        help="Maximum premature stop codons allowed in CDS body (default: 0).")
    parser.add_argument("--allow-frameshifts", action="store_true",
                        help="Allow indels/frameshifts with premature stop codons into clean cohort.")
    parser.add_argument("--full-length-only", action="store_true",
                        help="Only include accessions that cover 100%% of query length.")
    parser.add_argument("--max-records", type=int, default=None,
                        help="Stop after processing N lines (useful for testing/previews).")
    parser.add_argument("--sample-cap", type=int, default=10,
                        help="Maximum sample accessions to record per haplotype in TSV (default: 10).")
    parser.add_argument("-o", "--out-prefix", "--out-dir", dest="out_prefix", type=str, default=None,
                        help="Output prefix or directory for generated files (default: msa/<ref_name>).")
    parser.add_argument("--no-a2m", action="store_true",
                        help="Disable generating .a2m alignment files.")
    return parser.parse_args()

def main():
    args = parse_args()
    if not args.dataset and not args.file:
        sys.exit("Error: Must specify either -d/--dataset <GALAXY_ID> or -f/-i/--file <PATH> (or '-' for stdin).")

    stream_lexicmap_to_msa(
        ref_path=args.ref,
        dataset_id=args.dataset,
        file_path=args.file,
        base_url=args.base_url,
        api_key=args.api_key,
        api_key_file=args.api_key_file,
        min_coverage=args.min_coverage,
        min_coverage_partial=args.min_coverage_partial,
        min_pident=args.min_pident,
        max_internal_stops=args.max_internal_stops,
        allow_frameshifts=args.allow_frameshifts,
        full_length_only=args.full_length_only,
        max_records=args.max_records,
        sample_cap=args.sample_cap,
        out_prefix=args.out_prefix,
        write_a2m=not args.no_a2m
    )

if __name__ == "__main__":
    main()
