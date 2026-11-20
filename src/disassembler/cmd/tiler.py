#!/usr/bin/env python3
"""tiler (`disassembler tile`): streaming multi-HSP coordinate tiling, codon QC, and haplotype collapsing.

Reads a LexicMap tabular search from a Galaxy dataset, a local TSV or stdin, and writes a
codon-checked cohort. The user documentation is docs/tiler.md; this is the contract.

What it does, in order:

  1. **Tiles multi-HSP accessions.** LexicMap is indexed on Logan's *contigs*, not its unitig
     graph, and on a diverse sample the assembler shatters the gene at every bubble it cannot
     pop, so one accession arrives as many HSPs. Rows are grouped by `sgenome`, projected onto
     the reference CDS coordinates, and each overlapping position is resolved on its own by
     how deep it lies inside each HSP, recording where two HSPs disagree. See
     docs/notes/contigs_vs_unitigs.md.
  2. **Checks the reading frame**, translating in-frame CDS (NCBI code 11 or standard) and
     auditing premature internal stops at exact codon coordinates.
  3. **Sorts accessions into six tiers**, written to the ledger as `status` and `tier`:

       1 CLEAN_FULL_LENGTH       coverage >= min_coverage, no premature stop, no frameshift
       2 CLEAN_PARTIAL           coverage >= min_coverage_partial, same stop/frame rules
       3 FLAGGED_PREMATURE_STOPS covered and in frame, stops > max_stops
       4 LOW_COVERAGE_FRAGMENT   coverage < min_coverage_partial, whatever its frame says
       5 FLAGGED_FRAMESHIFT      covered, frame shifted for >= FRAMESHIFT_MIN_SPAN bases
       6 FLAGGED_MIXED_HSPS      covered, in frame, no stops, but its own HSPs disagree

     The rules apply in the order coverage, frameshift, stops, conflict, which is not tier
     order: too little of the gene on the grid to judge is a fragment whatever else is true.
     Indels come from the aligned qseq/sseq; an in-frame deletion stays in the row as `-` and
     counts as covered, while an insertion has no column and is recorded separately. An indel
     within INDEL_EDGE_WINDOW bases of an HSP end is reported but not trusted.
  4. **Collapses tier 1 into haplotypes** as they stream, calling nucleotide and amino-acid
     mutations and classifying each (WT, WT_partial, Synonymous, Missense, Nonsense,
     InFrameIndel, Frameshift). A haplotype is its row plus its trusted indels.

     Memory is O(distinct clean full-length haplotypes), **not** O(1): `clean_full_haps`,
     `clean_full_pident_sum` and `clean_full_samples` are keyed by the merged sequence, so on
     a divergent target this approaches one entry per accession. The ROW count is what
     streams, and rows are what is unbounded.
  5. **Reports the regime**: coverage and identity quantiles, and a warning when the target is
     divergent enough (median pident < 85% or cov < 0.60) that `logan-walker` is the better
     route.

Outputs, all at `<out-prefix>`:

  .clean.msa.fasta              collapsed clean haplotypes, with rank and count
  .clean.haplotypes.tsv         mutation audit, frequencies, sample accessions
  .clean.accessions.fasta       uncollapsed full-length clean accessions
  .clean_expanded.accessions.fasta   clean cohort including partials
  .flagged.accessions.fasta     stops, frameshifts, low coverage, mixed HSPs
  .cohort_ledger.tsv            --cohort-ledger only: the summary's records[] as a TSV
  .summary.json                 quantiles, cohort breakdown, multi-HSP rescue counts
  .clean.a2m                    NUCLEOTIDE A2M of the reference plus tier 1; translate it
                                before use with protein MSA models

Every row is projected onto the reference grid, so the FASTA and MSA outputs are already
column-aligned at reference length and need no MAFFT. A position no HSP covers, or that the
sample deletes, is `-`, so at the default --min-cov 0.80 a clean row can be a fifth `-` plus
its deletions; under --require-full-length the only `-` are in-frame deletions.
"""

import os
import sys
import math
import json

from disassembler import __version__
from disassembler import galaxy
from disassembler import summary as summary_io
import time
import argparse
import array
import bisect
import http.client
from Bio.Data import CodonTable
import urllib.request
from pathlib import Path
from collections import Counter, defaultdict
from typing import NamedTuple

from disassembler.textio import open_text, text_stream
from disassembler.validators import (
    file_safe_name,
    fraction,
    non_negative_int,
    percentage,
    positive_int,
    validate_api_key_file,
)

# Translation table: NCBI standard / bacterial code 11
def codon_table_for(genetic_code: int = 11) -> dict:
    """`{codon: residue}` with `*` for stops, for an NCBI genetic code id.

    The table used to be a single hand-written code 11 literal with no way to choose another.
    TGA is Trp in code 4 and in the mitochondrial codes, and CTG is Ser in code 12, which the
    shipped `queries/markers/candida_auris_erg11.fasta` uses: it has CTG at codon 74, read as
    Leu, so that residue was mislabelled and TTG to CTG changes came out `Synonymous`. A
    reference using another code also has its own sense codons counted as premature stops, so
    every accession, perfect matches included, lands in FLAGGED_PREMATURE_STOPS and the run
    exits 1 saying the clean cohort is empty, which names the wrong cause.

    Biopython's tables rather than more literals: it is already a dependency, and it is the
    reference the rest of the field checks against. The code 11 table it returns was compared
    against the literal this replaced, entry for entry, and is identical.
    """
    try:
        tbl = CodonTable.unambiguous_dna_by_id[genetic_code]
    except KeyError:
        known = ", ".join(str(i) for i in sorted(CodonTable.unambiguous_dna_by_id))
        raise ValueError(f"unknown genetic code {genetic_code}; NCBI ids are {known}") from None
    out = dict(tbl.forward_table)
    out.update({c: '*' for c in tbl.stop_codons})
    return out


CODON_TABLE_11 = codon_table_for(11)

def get_galaxy_api_key(api_key: str = None, key_file: str = None) -> str:
    """The Galaxy API key. One implementation, in `disassembler.galaxy`; this name
    stays because the tests and `harvest` both reach for it."""
    return galaxy.api_key(api_key, key_file)


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

def find_internal_stops(seq: str, ref_len: int, codon_table: dict = CODON_TABLE_11,
                        skip_last: bool = True) -> list:
    """
    Finds 1-based codon positions of premature stop codons in the coding body.

    `skip_last` leaves the final codon out, which is right only when the reference ends in its
    natural stop. It cannot be unconditional: six shipped references are polyprotein segments
    that end in a sense codon (`dengue_1_ns5` to `dengue_4_ns5`, `westnile_ns5`, `zika_ns5`,
    which ends CTG). Skipping the last codon of those lets a nonsense change there pass the
    tier rule while the haplotype table, which translates the whole row, calls the same
    haplotype `Nonsense`. Dengue NS5 ends TGG, one G>A from
    TAG or TGA.
    """
    stop_positions = []
    for i in range(0, (ref_len - 3) if skip_last else (ref_len - 2), 3):
        codon = seq[i:i + 3]
        if '-' not in codon and len(codon) == 3:
            if codon_table.get(codon, 'X') == '*':
                stop_positions.append(i // 3 + 1)
    return stop_positions

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
    with open_text(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                n_records += 1
                if n_records == 1:
                    header = line[1:]
                continue
            if n_records == 0:
                # `msa` rejects text before the first header and this reader skipped
                # it, so one file had two lengths in one `from-lexicmap` run.
                raise ValueError(f"{path}: text before the first '>' header ({line[:40]!r})")
            if n_records == 1:
                # A space or tab inside a sequence line must not become a reference
                # POSITION. `validators.fasta_record_lengths` and Biopython both drop it, so
                # a 12 nt CDS written `ATGGCT GAATAA` would load as 13 nt here and 12
                # everywhere else: `is_cds` goes false, the stop and frame rules switch off,
                # and an identical sample scores 0.9231 coverage, all at exit 0.
                seq_parts.append("".join(line.split()))

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
    try:
        return urllib.request.urlopen(req, timeout=120)
    except OSError as e:
        raise RuntimeError(f"could not download Galaxy dataset {dataset_id} from "
                           f"{base_url}: {e}") from e

# Galaxy dataset states that mean "not finished yet", as opposed to "failed". These are
# Galaxy's own `Dataset.non_ready_states`.
GALAXY_PENDING_STATES = {"new", "upload", "queued", "running", "setting_metadata"}

# States that will never become `ok` on their own. `deferred` and `paused` used to be
# treated as pending, so the run exited EX_TEMPFAIL saying "rerun once it is 'ok'" and a
# script retrying on 75, which the docs invite, looped until its own limit. A deferred
# dataset is terminal and needs a materialize request; a paused job needs a person to resume
# it. Both are exit 1 with the action to take.
GALAXY_NEEDS_ACTION = {
    "deferred": ("it is a deferred dataset, which is a terminal state: materialize it in "
                 "Galaxy and pass the new dataset's id"),
    "paused": "its job is paused: resume it in Galaxy and rerun once the dataset is 'ok'",
}


class DatasetNotReady(RuntimeError):
    """The Galaxy dataset is still being produced. `run` exits EX_TEMPFAIL (75) for it, so a
    script can tell "try again later" from a failure."""


def resolve_dataset_state(dataset_id: str, base_url: str = "https://usegalaxy.org", api_key: str = None, key_file: str = None):
    """The live state of a Galaxy dataset. Raises if it cannot be read: proceeding without it
    used to stream whatever a running or failed job had left."""
    key = get_galaxy_api_key(api_key, key_file)
    url = f"{base_url.rstrip('/')}/api/datasets/{dataset_id}"
    req = urllib.request.Request(url, headers={"x-api-key": key})
    try:
        with urllib.request.urlopen(req, timeout=30) as res:
            return json.loads(res.read()).get("state")
    except (OSError, ValueError) as e:
        raise RuntimeError(f"could not read the state of Galaxy dataset {dataset_id} at "
                           f"{base_url}: {e}") from e


def check_dataset_ready(dataset_id: str, state) -> None:
    """Raise unless the dataset is `ok`. Only `queued` was caught before, and it returned
    exit 0 having written nothing, so `from-lexicmap` went on to a stale cohort; any other
    state was streamed as if it were finished."""
    if state == "ok":
        return
    if state in GALAXY_PENDING_STATES:
        raise DatasetNotReady(f"Galaxy dataset {dataset_id} is '{state}', not finished yet. "
                              f"Rerun once it is 'ok'.")
    if state in GALAXY_NEEDS_ACTION:
        raise ValueError(f"Galaxy dataset {dataset_id} is '{state}': "
                         f"{GALAXY_NEEDS_ACTION[state]}.")
    raise ValueError(f"Galaxy dataset {dataset_id} is in state '{state}'; only an 'ok' "
                     f"dataset can be read.")


# An indel this close to either end of its HSP is not trusted. LexicMap can close an alignment
# with a forced gap a few bases from where the contig or the extension stops: on the C. auris
# ERG11 panel, three runs from three studies share a 1-nt insertion 3 aligned bases from the start
# of their HSP. Such an indel is still reported (marked "(edge)"), but it never sends an accession
# to the frameshift tier, never counts as covered, and never splits a haplotype.
INDEL_EDGE_WINDOW = 10
# How many conflicting positions a row names before it is elided. The count is always exact;
# this bounds the ledger cell and the summary record (docs/summaries.md, 1 KB per record).
CONFLICT_REPORT_LIMIT = 10

# A disagreement counts only if **at least one** of the two HSPs spans this many reference bases.
# A conserved gene's termini match short fragments all over a Logan assembly, and two such
# fragments disagreeing with each other say nothing about the allele the sample carries. The depth
# rule does not stop them, because 10 bases of depth fit inside a 64-base scrap.
#
# At least one, not both, and the difference matters. An unpopped Logan SNP bubble is a short
# branch by nature (the regression fixture's is 61 bases, over LexicMap's -l 50), so requiring
# both to be long would silence the case this tier exists for. What it must silence is a position
# contested *only* by fragments, with no HSP over it long enough to settle the question.
#
# Measured on C. auris EF1A, where all 6,502 accessions were flagged: the conflicts sat on four
# reference positions out of 1,377, two near each end, 99.9% lay outside the accession's own
# longest HSP, and the contesting pair was 64 and 84 bases. At 100 this leaves **6 of 6,502**
# flagged, while H. pylori babA keeps **all 1,576** of its accessions and 100% of its 15,098
# conflicting positions, where the disagreement is real babA/babB paralog mixing inside a
# 3,873-base alignment. Requiring both to be long instead costs babA 126 accessions and 9% of its
# positions and buys 4 accessions on EF1A, which is the wrong trade. 50 is too low either way: the
# EF1A scraps clear it.
CONFLICT_MIN_HSP_BP = 100
# ...but never more than a tenth of the reference, so a short reference is not silenced outright.
# The floor is CONFLICT_MIN_HSP_BP for any reference of 1 kb or more, which is every coding gene
# in the shipped panels.
CONFLICT_MIN_HSP_REF_FRACTION = 10

# A shift of the reading frame counts as a frameshift only if it lasts this many reference bases,
# or never ends. Indels that cancel each other sooner (net length a multiple of 3) are usually how
# an aligner writes a divergent stretch: on H. pylori babA, sabA and vacA most such clusters close
# within 5-20 bases, while real frameshifts (long-read homopolymer deletions in C. auris, the oipA
# phase-variation repeat) stay shifted for hundreds of bases or to the end of the gene.
FRAMESHIFT_MIN_SPAN = 30


class Indel(NamedTuple):
    """One insertion or deletion in the sample, relative to the reference."""
    kind: str    # "ins" or "del"
    pos: int     # 0-based: a deletion starts here; an insertion sits between pos-1 and pos
    length: int
    bases: str   # inserted bases, or the reference bases deleted
    edge: bool   # within INDEL_EDGE_WINDOW aligned bases of its HSP's end

    def notation(self) -> str:
        """1-based, HGVS-like: `685_686insT`, `2316_2321del`, `709del`."""
        if self.kind == "ins":
            s = f"{self.pos}_{self.pos + 1}ins{self.bases}"
        elif self.length == 1:
            s = f"{self.pos + 1}del"
        else:
            s = f"{self.pos + 1}_{self.pos + self.length}del"
        return s + "(edge)" if self.edge else s


class Tiling(NamedTuple):
    seq: str             # reference-length row; '-' where no base was placed
    covered: int         # positions carrying a base, plus positions of trusted deletions
    mean_pident: float   # over the positions carrying a base
    single_ge_min_cov: bool
    single_full: bool
    indels: list         # Indel, in reference order
    conflicts: tuple = ()  # 0-based reference positions where two HSPs disagreed


def _covered_runs(seq: str):
    """Half-open runs of positions carrying a base. `-` is an untiled position or a trusted
    deletion; which it is does not matter here, because folding compares only keys that carry
    the same trusted indels, so a deletion is `-` on both sides of every comparison."""
    runs, start = [], None
    for i, c in enumerate(seq):
        if c != '-':
            if start is None:
                start = i
        elif start is not None:
            runs.append((start, i))
            start = None
    if start is not None:
        runs.append((start, len(seq)))
    return runs


def coverage_independent_haplotypes(haps, hap_indels, ref_seq: str, ref_len: int) -> dict:
    """Haplotype counts that do not move when an accession's contigs stop short.

    Tier 1 rows are collapsed on the whole reference-length row, which carries `-` at every
    untiled position, and with the default `--min-cov 0.80` a tier 1 row may be a fifth `-`.
    Two accessions with identical bases wherever they are covered but different contig ends
    therefore become different haplotypes, so `clean_unique_haplotypes`,
    `shannon_entropy_bits` and `wildtype_fraction` measure assembly fragmentation as well as
    sequence diversity: reported diversity tracks sequencing depth, and numbers from two
    cohorts or two targets cannot be compared.

    The original fields are left exactly as they were, for comparability with runs already
    done. These are computed beside them, two ways:

    * **full coverage** -- the subset of keys that cover the whole reference (bases placed,
      plus the reference bases a trusted deletion removes). Nothing is inferred; the cost is
      that a cohort of short contigs can leave it empty.
    * **folded** -- each partial key is folded into the *one* full-coverage key that carries
      the same trusted indels and agrees at every position the partial covers. A partial that
      matches none, or more than one, is counted in `partial_rows_unresolved` and never
      guessed, because a partial whose variant falls in its own uncovered region is genuinely
      ambiguous. Unresolved rows are left out of the folded entropy and fraction rather than
      being counted as haplotypes of their own, which is the artefact being removed.

    Folding on variant signatures alone would be wrong for the same reason: it would call a
    partial row wild type whenever the variant it is missing lies outside its coverage.
    """
    ref_key = (ref_seq, ())
    fulls, partials = [], []
    for key, count in haps.items():
        seq = key[0]
        dels = sum(e.length for e in hap_indels.get(key, ()) if e.kind == "del")
        based = sum(1 for c in seq if c != '-')
        (fulls if based + dels >= ref_len else partials).append((key, count))

    full_rows = sum(c for _k, c in fulls)
    out = {
        "full_coverage_rows": full_rows,
        "clean_unique_haplotypes_full_coverage": len(fulls),
        "shannon_entropy_bits_full_coverage": round(shannon_bits([c for _k, c in fulls]), 4),
    }

    # Index the full-coverage keys by trusted-indel signature, and within that by each
    # position where they differ from the reference. A near-reference cohort has few such
    # positions, so a partial's candidate set is found without scanning every full key.
    buckets = defaultdict(list)                    # indel notation tuple -> [full index]
    variants = []                                  # full index -> {pos: base}
    by_variant = defaultdict(set)                  # (indel sig, pos, base) -> {full index}
    for i, (key, _count) in enumerate(fulls):
        seq, sig = key
        v = {j: c for j, c in enumerate(seq) if c != '-' and c != ref_seq[j]}
        variants.append(v)
        buckets[sig].append(i)
        for pos, base in v.items():
            by_variant[(sig, pos, base)].add(i)

    folded = Counter({key: count for key, count in fulls})
    n_folded = n_unresolved = 0
    for key, count in partials:
        seq, sig = key
        bucket = buckets.get(sig, [])
        if not bucket:
            n_unresolved += count
            continue
        pv = {j: c for j, c in enumerate(seq) if c != '-' and c != ref_seq[j]}
        if pv:
            cand = set(bucket)
            for pos, base in pv.items():
                cand &= by_variant.get((sig, pos, base), set())
                if not cand:
                    break
        else:
            cand = set(bucket)
        # A candidate must also carry no variant of its own inside what the partial covers:
        # agreeing on the partial's variants is not enough if the full key differs somewhere
        # the partial calls reference.
        runs = _covered_runs(seq)
        cand = {
            i for i in cand
            if not any(pos not in pv and any(a <= pos < b for a, b in runs)
                       for pos in variants[i])
        }
        if len(cand) == 1:
            folded[fulls[next(iter(cand))][0]] += count
            n_folded += count
        else:
            n_unresolved += count

    folded_rows = sum(folded.values())
    wt = folded.get(ref_key, 0)
    out.update({
        "partial_rows_folded": n_folded,
        "partial_rows_unresolved": n_unresolved,
        "clean_unique_haplotypes_folded": len(folded),
        "shannon_entropy_bits_folded": round(shannon_bits(list(folded.values())), 4),
        "wildtype_compatible_count": wt,
        "wildtype_fraction_folded": round(wt / folded_rows, 6) if folded_rows else 0,
    })
    return out


def shannon_bits(counts) -> float:
    """Shannon entropy in bits over a list of counts."""
    total = sum(counts)
    if total <= 0:
        return 0.0
    h = 0.0
    for c in counts:
        if c > 0:
            pr = c / total
            h -= pr * math.log2(pr)
    return h


class StagedOutputs:
    """Writes every bulk output to `<name>.part` and renames them in only on success.

    The cohort files stream from the first row while the MSA, A2M, haplotype tables and
    summary are written at the end, so opening them under their final names leaves a failed
    rerun's cohort files beside the previous run's MSA, with nothing on disk saying which
    belong together.

    `summary.json` is deliberately NOT staged: a failed run still writes one, saying it
    failed, which is how a reader is meant to tell (docs/summaries.md rule 4). Anything else
    left over from an earlier run is removed on commit, so a success is always a complete,
    self-consistent set.
    """

    # Every bulk output this stage can write, including the compatibility aliases.
    SUFFIXES = (
        ".clean.accessions.fasta", ".clean_expanded.accessions.fasta",
        ".flagged.accessions.fasta", ".cohort_ledger.tsv",
        ".clean.msa.fasta", ".msa.fasta",
        ".clean.haplotypes.tsv", ".haplotypes.tsv",
        ".clean.a2m", ".a2m",
    )

    def __init__(self, prefix: str):
        self.prefix = str(prefix)
        self.staged = []

    def path(self, suffix: str):
        assert suffix in self.SUFFIXES, suffix
        if suffix not in self.staged:
            self.staged.append(suffix)
        return Path(f"{self.prefix}{suffix}.part")

    def open(self, suffix: str):
        return open(self.path(suffix), "w", encoding="utf-8")

    def abort(self):
        self.clear(self.prefix)

    @classmethod
    def clear(cls, prefix):
        """Remove every `.part` under `prefix`, whoever wrote it.

        Called from `_write_failed_summary` as well as from the stream's own handler, so a
        failure anywhere after the stream, in the MSA or haplotype writers, leaves nothing
        half-written either.
        """
        for suffix in cls.SUFFIXES:
            Path(f"{prefix}{suffix}.part").unlink(missing_ok=True)

    def commit(self):
        for suffix in self.SUFFIXES:
            final = Path(f"{self.prefix}{suffix}")
            part = Path(f"{final}.part")
            if part.exists():
                os.replace(part, final)
            else:
                # An output this run did not produce, such as the A2M under --no-a2m or the
                # ledger without --cohort-ledger. A copy from an earlier run would otherwise
                # be read as part of this one.
                final.unlink(missing_ok=True)


def output_stem(ref_name: str, ref_path: str) -> str:
    """The file stem for a run's outputs: `ref_name` made safe as one path component.

    `ref_name` itself is unchanged and keeps driving the query-name filter, the FASTA headers
    and `summary.json`, because that is the name the search used and it has to keep matching.
    """
    return (file_safe_name(ref_name) or file_safe_name(Path(ref_path).stem) or "tiler")


def resolve_out_prefix(out_prefix, ref_name: str) -> str:
    """`<dir>/<ref_name>` for a directory, the path itself for a stem, `msa/<ref_name>` for
    nothing, creating the directory in each case.

    Shared with `_write_failed_summary` deliberately. The two used to resolve the prefix
    differently, so a run that failed wrote its summary to `<dir>/.summary.json`: a dotfile,
    under a name no wrapper looking for `<prefix>.summary.json` would find, and not the name
    the same run would have used had it succeeded. Rule 4 of docs/summaries.md says a summary
    is always written, which is only useful if it is written where it is looked for.
    """
    if not out_prefix:
        Path("msa").mkdir(parents=True, exist_ok=True)
        return f"msa/{ref_name}"
    out_path = Path(out_prefix)
    if out_path.is_dir() or str(out_prefix).endswith("/"):
        out_path.mkdir(parents=True, exist_ok=True)
        return str(out_path / ref_name)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    return str(out_prefix)


def tile_accession_hsps(hsps: list, ref_len: int, min_pident: float, min_coverage: float,
                        min_conflict_hsp_bp: int = CONFLICT_MIN_HSP_BP) -> Tiling:
    """
    Tiles all qualifying HSPs of an accession onto the reference coordinate system.

    Overlapping positions are resolved **per position**, by depth: the base comes from the HSP
    the position lies furthest inside, measured to the nearer of that HSP's two ends, then by
    `pident`, then by input order. Where two HSPs place different symbols at a position that
    is at least INDEL_EDGE_WINDOW inside both, and at least one of which spans
    `min_conflict_hsp_bp` reference bases, the position is recorded in `conflicts`.

    Per position, and not per HSP. Ranking whole HSPs by identity to the reference makes the
    reference branch win every position it touches: Logan leaves an unpopped SNP bubble as two
    contigs, both of which pass LexicMap's defaults, and only the reference branch is 100%
    identical. That overrides the real variant whatever the sample carries, and a deletion on
    the losing HSP goes with it, so a mixed sample comes out CLEAN_FULL_LENGTH wild type.
    Measured at 37-46 wrong bases per accession on emulated HSPs over real Logan contigs.

    Depth is the right tiebreak because a short bubble branch is short: every one of its
    positions is near one of its ends, while the same positions sit deep inside the contig
    that spans the gene. Length and bitscore do not work: branches tie on length, and bitscore
    favours the reference.

    Conflicts are not resolved silently either way, because the tiler cannot know which
    branch the sample carries: they are reported, and `--allow-mixed` is how a caller says to
    tile them anyway. The INDEL_EDGE_WINDOW requirement on *both* sides keeps an ordinary
    ragged overlap at an HSP end, where one aligner simply ran out of sequence, from being
    read as two contigs disagreeing.

    `min_conflict_hsp_bp` is the second half of that guard and the reason the length of an HSP
    matters here and nowhere else in the tiler. Depth scales with the HSP, so 10 bases of it fit
    inside a 64-base fragment, and a conserved gene's termini match such fragments throughout an
    assembly. A position contested only by fragments, with nothing longer over it to settle the
    question, is not evidence of a mixed sample. One long HSP among the two is enough, because a
    bubble branch is short and still worth hearing when something spanning disagrees with it.

    The base each position *takes* is unaffected. The floor decides only whether a disagreement is
    recorded, never which symbol wins.

    Indels are read from the aligned qseq/sseq: a run of '-' in qseq is an insertion in the
    sample (its bases have no reference position, so they are left out of `seq`), and a run of
    '-' in sseq is a deletion (written into `seq` as '-'). An indel is kept only when its own HSP
    supplied the positions it affects: both flanks of an insertion, every position of a
    deletion. Otherwise another HSP's bases are what the row shows there.
    """
    # Four parallel lists, not one list of tuples. `chars` is the row itself, so an
    # ungapped HSP landing on untouched ground is one slice assignment from a string rather
    # than a tuple allocated per reference position; `own` is the HSP that supplied each
    # position, from which the rest is looked up. Depth and identity are not stored at all:
    # depth is `min(p - first, last - p)` for the owning HSP and identity is `pid_of[k]`.
    chars = [None] * ref_len
    own = [-1] * ref_len
    # The best real base at each position, ignoring gaps. A deletion writes '-' into the row
    # like any other symbol, so without this an edge deletion hides a base another HSP
    # supplied.
    bchar = [None] * ref_len
    bown = [-1] * ref_len

    spans = {}                # HSP index -> (first query position, one past the last)
    pid_of = {}               # HSP index -> its identity to the reference
    firsts, lasts = {}, {}    # and its first and last reference position, for the depth rule
    raw = []                  # (HSP index, kind, pos, length, bases)
    conflicts = set()
    single_ge_min_cov = False
    single_full = False

    # Every qualifying HSP's reference span, before any base is placed: the depth rule has to
    # know where an HSP ends while it is still writing its first position. A reference
    # position is consumed by each qseq character that is not a gap, which is what the walk
    # below also does, so this cannot drift from it.
    qualifying = []
    for k, h in enumerate(hsps):
        pident, qcov, qstart, qend, qseq, sseq = h
        if pident < min_pident:
            continue
        if qcov >= (min_coverage * 100.0):
            single_ge_min_cov = True
        if qstart == 1 and qend == ref_len:
            single_full = True
        qualifying.append((k, h))
        spans[k] = (qstart - 1, qstart - 1 + len(qseq) - qseq.count('-'))
        pid_of[k] = pident
        firsts[k] = spans[k][0]
        lasts[k] = spans[k][1] - 1

    # Which HSPs are long enough to contest a position, as a flat list indexed by HSP number:
    # `own[p]` holds an HSP index, and this is read once per contested position in the hot loop
    # below, so it cannot be a set lookup.
    conflict_floor = min(min_conflict_hsp_bp, ref_len // CONFLICT_MIN_HSP_REF_FRACTION)
    long_enough = [False] * len(hsps)
    for k, (a, b) in spans.items():
        long_enough[k] = (b - a) >= conflict_floor

    def depth_at(p, k):
        first, end = spans[k]
        a, b = p - first, end - 1 - p
        return a if a < b else b

    # Positions already carrying a base, as a few disjoint half-open intervals. An HSP that
    # lands where nothing has been placed has no contest to run and no conflict to detect,
    # so those stretches are written by slice assignment. On a shattered gene most positions
    # are like that: the HSPs tile it, and only their k-1 overlaps are contested.
    filled = []

    def untouched(a, b):
        """The sub-intervals of [a, b) that no HSP has written yet."""
        out, cur = [], a
        for fa, fb in filled:
            if fb <= cur:
                continue
            if fa >= b:
                break
            if fa > cur:
                out.append((cur, fa))
            cur = fb if fb > cur else cur
            if cur >= b:
                return out
        if cur < b:
            out.append((cur, b))
        return out

    def mark(a, b):
        filled.append((a, b))
        filled.sort()
        merged = []
        for fa, fb in filled:
            if merged and fa <= merged[-1][1]:
                if fb > merged[-1][1]:
                    merged[-1][1] = fb
            else:
                merged.append([fa, fb])
        filled[:] = [(x, y) for x, y in merged]

    def place(p, ch, k):
        """One contested position: deepest wins, then identity, then input order."""
        ex = own[p]
        if ex < 0:
            chars[p], own[p] = ch, k
        else:
            d_new, d_old = depth_at(p, k), depth_at(p, ex)
            if (chars[p] != ch and (long_enough[k] or long_enough[ex])
                    and d_new >= INDEL_EDGE_WINDOW
                    and d_old >= INDEL_EDGE_WINDOW):
                conflicts.add(p)
            if (d_new, pid_of[k]) > (d_old, pid_of[ex]):
                chars[p], own[p] = ch, k
        if ch != '-':
            b = bown[p]
            if b < 0 or (depth_at(p, k), pid_of[k]) > (depth_at(p, b), pid_of[b]):
                bchar[p], bown[p] = ch, k

    for k, (pident, qcov, qstart, qend, qseq, sseq) in qualifying:
        first, end = spans[k]
        last = end - 1
        k_long = long_enough[k]

        q_pos = qstart - 1
        if '-' not in qseq and '-' not in sseq:
            # Most HSPs carry no gap at all; skip the indel bookkeeping for them.
            lo, hi = (0 if first < 0 else first), (ref_len if end > ref_len else end)
            if lo >= hi:
                continue
            upper = sseq.upper()
            free = untouched(lo, hi)
            cur = lo
            contested = []
            for a, b in free:
                if a > cur:
                    contested.append((cur, a))
                chars[a:b] = upper[a - first:b - first]
                own[a:b] = [k] * (b - a)
                bchar[a:b] = upper[a - first:b - first]
                bown[a:b] = [k] * (b - a)
                cur = b
            if cur < hi:
                contested.append((cur, hi))
            # Inlined rather than calling `place`: `place` and the `depth_at` inside it were
            # 60% of this function's time, at four calls per contested position.
            for a, b in contested:
                for p in range(a, b):
                    ch = upper[p - first]
                    ex = own[p]
                    d_new = p - first
                    if last - p < d_new:
                        d_new = last - p
                    if ex < 0:
                        chars[p], own[p] = ch, k
                    else:
                        d_old = p - firsts[ex]
                        if lasts[ex] - p < d_old:
                            d_old = lasts[ex] - p
                        if (chars[p] != ch and (k_long or long_enough[ex])
                                and d_new >= INDEL_EDGE_WINDOW
                                and d_old >= INDEL_EDGE_WINDOW):
                            conflicts.add(p)
                        if (d_new, pident) > (d_old, pid_of[ex]):
                            chars[p], own[p] = ch, k
                    if ch != '-':
                        bk = bown[p]
                        if bk < 0:
                            bchar[p], bown[p] = ch, k
                        else:
                            d_b = p - firsts[bk]
                            if lasts[bk] - p < d_b:
                                d_b = lasts[bk] - p
                            if (d_new, pident) > (d_b, pid_of[bk]):
                                bchar[p], bown[p] = ch, k
            mark(lo, hi)
            continue

        ins, dele, del_start = [], [], 0
        for qc, sc in zip(qseq, sseq):
            if qc == '-':
                if dele:
                    raw.append((k, "del", del_start, len(dele), "".join(dele).upper()))
                    dele = []
                ins.append(sc)
                continue
            if ins:
                raw.append((k, "ins", q_pos, len(ins), "".join(ins).upper()))
                ins = []
            if sc == '-':
                if not dele:
                    del_start = q_pos
                dele.append(qc)
            elif dele:
                raw.append((k, "del", del_start, len(dele), "".join(dele).upper()))
                dele = []
            if 0 <= q_pos < ref_len:
                place(q_pos, sc.upper(), k)
            q_pos += 1
        if dele:
            raw.append((k, "del", del_start, len(dele), "".join(dele).upper()))
        if ins:
            raw.append((k, "ins", q_pos, len(ins), "".join(ins).upper()))
        lo, hi = (0 if first < 0 else first), (ref_len if end > ref_len else end)
        if lo < hi:
            mark(lo, hi)

    # An insertion occupies no reference position, so two HSPs disagreeing about one cannot
    # show up as two symbols at a position. If another HSP spans the junction well inside
    # itself and carries no matching insertion there, they disagree.
    ins_at = {(k, pos): bases for k, kind, pos, length, bases in raw if kind == "ins"}
    for (k, pos), bases in ins_at.items():
        first, end = spans[k]
        if min(pos - first, end - pos) < INDEL_EDGE_WINDOW:
            continue
        for j, (jf, je) in spans.items():
            if j == k or pos - 1 < jf or pos >= je:
                continue
            if not (long_enough[k] or long_enough[j]):
                continue
            if min(pos - 1 - jf, je - pos) < INDEL_EDGE_WINDOW:
                continue
            if ins_at.get((j, pos)) != bases:
                conflicts.add(pos)
                break

    def supplied_by(p, k):
        return 0 <= p < ref_len and own[p] == k

    def margin_of(k, kind, pos, length):
        first, end = spans[k]
        return (min(pos - first, end - pos) if kind == "ins"
                else min(pos - first, end - (pos + length)))

    # An indel is an end artefact only if EVERY HSP carrying it sees it near one of its own
    # ends. Reading the margin from the owning HSP alone marks an indel 10 bases from one
    # contig's end but deep inside the next contig's overlap as `(edge)`, and which contig
    # owns it depends on input order: a frameshifted gene enters the clean cohort for
    # no better reason than the order of the rows.
    #
    # Measured note: with ownership decided by depth this set cannot change a verdict,
    # because the HSP that owns every position of an indel is necessarily the one with the
    # largest margin. Checked over 20,000 randomized multi-HSP tilings with indels: zero
    # differences against reading the margin from the owner alone. It is kept because it
    # states the rule itself rather than a consequence of how ownership happens to be
    # decided, so a future change to that rule cannot quietly bring the defect back.
    inner = {(kind, pos, length, bases)
             for k, kind, pos, length, bases in raw
             if margin_of(k, kind, pos, length) >= INDEL_EDGE_WINDOW}

    indels = []
    for k, kind, pos, length, bases in raw:
        if kind == "ins":
            owned = supplied_by(pos - 1, k) and supplied_by(pos, k)
        else:
            owned = all(supplied_by(p, k) for p in range(pos, pos + length))
        if owned:
            indels.append(Indel(kind, pos, length, bases,
                                (kind, pos, length, bases) not in inner))
    indels.sort(key=lambda e: e.pos)

    # An edge deletion is an alignment artefact, not a claim about the sample, so its '-'
    # must not hide a base another HSP placed at the same position. Letting the gap into the
    # per-position contest leaves it there even after the deletion is ruled an artefact, so
    # the row shows '-', `covered` is short by its length, and the accession becomes a
    # haplotype of its own.
    # Refilled only when every position of the deletion has a base to put back: a partial
    # refill would leave an indel record that no longer describes what the row shows.
    if any(e.kind == "del" and e.edge for e in indels):
        kept = []
        for e in indels:
            span = range(e.pos, e.pos + e.length)
            if e.kind == "del" and e.edge and all(bown[q] >= 0 for q in span):
                for q in span:
                    chars[q], own[q] = bchar[q], bown[q]
                continue
            kept.append(e)
        indels = kept

    # Mean identity over the positions that actually carry a subject base. Summing pident over
    # gap positions too (subject deletions, where the row holds '-') inflates the mean and can
    # push it above 100%.
    merged_seq = "".join(c if c is not None else '-' for c in chars)
    n_based = ref_len - merged_seq.count('-')
    # `sum` over a list, not a running total: CPython's `sum` compensates for rounding when
    # it is given floats, and a plain `+=` loop differed from it in the last digit on 3 of
    # 30,000 randomized tilings. Reported identity is compared between runs, so it has to be
    # the same number however it was reached.
    mean_pid = ((sum([pid_of[k] for c, k in zip(chars, own) if c is not None and c != '-'])
                 / n_based) if n_based else 0.0)
    # A trusted deletion is a stretch of the reference the sample aligned across and lacks; it
    # counts as covered, so --require-full-length does not reject an accession for carrying one.
    covered = n_based + sum(e.length for e in indels if e.kind == "del" and not e.edge)

    return Tiling(merged_seq, covered, mean_pid, single_ge_min_cov, single_full, indels,
                  tuple(sorted(conflicts)))



def split_by_frame(indels: list) -> tuple:
    """
    Splits an accession's trusted indels by their effect on the reading frame.

    Walks them in reference order, tracking the frame offset (net indel length mod 3). Returns
    (in_frame, shifted): `in_frame` are indels that leave the frame intact, alone or as part of a
    cluster that restores it within FRAMESHIFT_MIN_SPAN reference bases; `shifted` are the indels
    of every stretch that stays out of frame for longer, or to the end of the gene. An accession
    with anything in `shifted` has a frameshift. Edge indels are in neither list.
    """
    in_frame, shifted, cluster = [], [], []
    offset = start = skipped = 0
    for e in sorted((e for e in indels if not e.edge), key=lambda e: e.pos):
        new = (offset + (e.length if e.kind == "ins" else -e.length)) % 3
        if offset == 0 and new == 0:
            in_frame.append(e)
            continue
        opening = offset == 0
        if opening:
            # The frame shifts only AFTER a deletion's own bases, which the sample
            # does not have. Opening the cluster at `e.pos` counted them as read out of
            # frame, so `31_64del,70_71del` (5 bases genuinely out of frame) was
            # FLAGGED_FRAMESHIFT while `31_34del,40_41del`, the same 5 bases, was clean.
            start = e.pos + (e.length if e.kind == "del" else 0)
            cluster, skipped = [], 0
        cluster.append(e)
        offset = new
        if offset == 0:
            span = e.pos - start - skipped
            (in_frame if span < FRAMESHIFT_MIN_SPAN else shifted).extend(cluster)
            cluster = []
        elif not opening and e.kind == "del":
            skipped += e.length
    shifted.extend(cluster)
    return in_frame, shifted


def apply_indels(seq: str, ref_len: int, applied: list):
    """The sample's own sequence, with `applied` indels spliced into the reference-grid row.

    Returns `(sample, starts, pieces)`: `starts[j]` is where piece `j` begins in the sample
    and `pieces[j]` is `(reference position it begins at, is it an insertion)`, which is what
    maps a position in the sample back to a reference codon.

    Split out of `find_internal_stops_with_indels` so the haplotype table can use it too.
    Calling amino-acid changes on the reference-GRID row instead loses them: an insertion is
    absent there, deleted bases are `-`, and `translate_nt` returns `-` for any codon holding
    a gap, so every codon an indel touches is skipped, inserted residues never appear, and a
    haplotype whose only changes are indels prints
    `num_aa_muts 0` and `aa_mutations NA` while `nt_mutations` and the `InFrameIndel` class
    said plainly that something had changed.
    """
    parts, starts, pieces = [], [], []
    sample_len, cur = 0, 0
    for e in applied:
        if e.pos > cur:
            parts.append(seq[cur:e.pos])
            starts.append(sample_len)
            pieces.append((cur, False))
            sample_len += e.pos - cur
            cur = e.pos
        if e.kind == "ins":
            parts.append(e.bases)
            starts.append(sample_len)
            pieces.append((e.pos, True))
            sample_len += e.length
        else:
            cur = max(cur, e.pos + e.length)
    parts.append(seq[cur:ref_len])
    starts.append(sample_len)
    pieces.append((cur, False))
    # 'N' never translates to a stop, which is how an uncovered '-' was always treated.
    return "".join(parts).replace('-', 'N'), starts, pieces


def _aa_same(a: str, b: str) -> bool:
    """Residues that are not evidence of a difference. `X` is an uncovered or ambiguous
    codon and `-` a codon the grid never filled; neither says the sample differs."""
    return a == b or a in "X-" or b in "X-"


def indel_protein_changes(ref_seq: str, row_seq: str, ref_aa: str, indels: list,
                          codon_table: dict = CODON_TABLE_11):
    """Protein changes for trusted in-frame indels, and the reference codons they cover.

    Returns `(changes, covered_codons)`. `covered_codons` are 0-based reference codon indices
    inside a reported window, so the caller does not also report a grid-based substitution
    there and count the same difference twice.

    Each cluster of nearby indels is spliced on its own and compared with the reference
    protein over a local window, so two indels in different parts of the gene are two
    changes rather than one enormous `delins`. The notation is HGVS-like: `P11del`,
    `K10_P11insA`, `R101_E105delinsKRVQQ`. HGVS 3' shifting in a repeat is not applied, so a
    call inside a homopolymer can differ from what another tool prints.
    """
    applied, shifted = split_by_frame(indels)
    if shifted or not applied:
        return [], set()

    # Cluster indels that are close enough for one window to describe both. The same
    # FRAMESHIFT_MIN_SPAN the frame rule uses, so "near" means the same thing in both places.
    clusters = []
    for e in applied:
        if clusters and e.pos - clusters[-1][-1].pos <= FRAMESHIFT_MIN_SPAN:
            clusters[-1].append(e)
        else:
            clusters.append([e])

    changes, covered = [], set()
    n_ref = len(ref_aa)
    for cluster in clusters:
        sample, _starts, _pieces = apply_indels(row_seq, len(ref_seq), cluster)
        h_aa = translate_nt(sample, codon_table)
        n_h = len(h_aa)

        # Bound the comparison to the cluster's own neighbourhood first, so a substitution
        # elsewhere in the row cannot drag the window open from either side, then trim it
        # from both ends to the shortest window that still describes the change. After the
        # cluster the two proteins differ only by `delta` codons, which is how far the right
        # bound has to move to line them up again.
        last = cluster[-1].pos + (cluster[-1].length if cluster[-1].kind == "del" else 0)
        delta = sum(e.length if e.kind == "ins" else -e.length for e in cluster) // 3
        lo = max(0, cluster[0].pos // 3 - 1)
        hi_r = min(n_ref, last // 3 + 2)
        hi_h = max(lo, min(n_h, hi_r + delta))

        while lo < hi_r and lo < hi_h and _aa_same(ref_aa[lo], h_aa[lo]):
            lo += 1
        while hi_r > lo and hi_h > lo and _aa_same(ref_aa[hi_r - 1], h_aa[hi_h - 1]):
            hi_r -= 1
            hi_h -= 1

        ref_win, sam_win = ref_aa[lo:hi_r], h_aa[lo:hi_h]
        if ref_win == sam_win:
            continue
        covered.update(range(lo, hi_r))
        if not ref_win:
            left = f"{ref_aa[lo - 1]}{lo}" if lo else "0"
            changes.append(f"{left}_{ref_aa[lo]}{lo + 1}ins{sam_win}" if lo < n_ref
                           else f"{left}ins{sam_win}")
        elif not sam_win:
            changes.append(f"{ref_aa[lo]}{lo + 1}del" if len(ref_win) == 1
                           else f"{ref_aa[lo]}{lo + 1}_{ref_aa[hi_r - 1]}{hi_r}del")
        elif len(ref_win) == 1 and len(sam_win) == 1:
            changes.append(f"{ref_win}{lo + 1}{sam_win}")
        else:
            span = (f"{ref_aa[lo]}{lo + 1}" if len(ref_win) == 1
                    else f"{ref_aa[lo]}{lo + 1}_{ref_aa[hi_r - 1]}{hi_r}")
            changes.append(f"{span}delins{sam_win}")
    return changes, covered


def find_internal_stops_with_indels(seq: str, ref_len: int, indels: list,
                                    codon_table: dict = CODON_TABLE_11,
                                    skip_last: bool = True) -> list:
    """
    `find_internal_stops` on the sample's own coding sequence: trusted in-frame deletions are
    removed and trusted in-frame insertions spliced back in, so a stop codon formed across an
    indel is found, including one inside a short stretch that a cluster of indels shifts out of
    frame and back. Positions are still reported as 1-based reference codon numbers. Untrusted
    indels and those of a lasting frameshift are left as the row shows them.
    """
    applied, _shifted = split_by_frame(indels)
    if not applied:
        return find_internal_stops(seq, ref_len, codon_table, skip_last)
    sample, starts, pieces = apply_indels(seq, ref_len, applied)

    stops = {c for c, aa in codon_table.items() if aa == '*'}
    stop_positions = []
    for i in range(0, (len(sample) - 3) if skip_last else (len(sample) - 2), 3):
        if sample[i:i + 3] in stops:
            j = bisect.bisect_right(starts, i) - 1
            ref_start, is_insertion = pieces[j]
            ref_pos = ref_start if is_insertion else ref_start + (i - starts[j])
            stop_positions.append(ref_pos // 3 + 1)
    return stop_positions

def stream_lexicmap_to_msa(
    ref_path: str,
    dataset_id: str = None,
    file_path: str = None,
    base_url: str = "https://usegalaxy.org",
    api_key: str = None,
    api_key_file: str = None,
    min_coverage: float = 0.80,
    min_coverage_partial: float = 0.50,
    min_pident: float = 60.0,
    max_internal_stops: int = 0,
    allow_frameshifts: bool = False,
    allow_mixed: bool = False,
    min_conflict_hsp_bp: int = CONFLICT_MIN_HSP_BP,
    genetic_code: int = 11,
    non_coding: bool = False,
    full_length_only: bool = False,
    max_records: int = None,
    sample_cap: int = 10,
    write_ledger: bool = False,
    out_prefix: str = None,
    write_a2m: bool = True
) -> dict:
    """
    Main streaming pipeline:
    - Streams rows, groups by accession, merges multi-HSP alignments
    - Partitions cohort into:
        1. Clean Full-Length (cov >= min_coverage, 0 internal stops, no frameshift,
           no HSPs of the accession disagreeing about a base)
        2. Clean Partial (min_coverage_partial <= cov < min_coverage, same rules)
        3-6. Flagged (internal stops, low coverage, a frameshift, or mixed HSPs)
    - Writes clean MSAs, flagged audit FASTA/TSV, and complete cohort ledger.
    """
    # Which of the three inputs this run read, for `inputs.source`. A summary that does not
    # say what it consumed cannot be checked against anything.
    if dataset_id:
        src_kind, src_id = "galaxy", str(dataset_id)
    elif file_path and file_path != "-":
        src_kind, src_id = "file", str(file_path)
    else:
        src_kind, src_id = "stdin", "-"

    # --require-full-length REPLACES --min-cov as the tier 1 gate, so everything that
    # describes or applies that gate has to read the same value: the banner, the tier 2
    # label, the argument check, and the value handed to `tile_accession_hsps`, which is
    # what `multi_hsp_rescued_accessions` is judged against in a
    # run whose gate was 100%.
    full_cov = 1.0 if full_length_only else min_coverage

    ref_name, ref_header, ref_seq = load_reference(ref_path)
    ref_len = len(ref_seq)
    codon_table = codon_table_for(genetic_code)

    # `ref_len % 3 == 0` is not enough on its own: one non-coding target in three has a
    # length that is a multiple of 3 (E. coli 16S is 1542 nt), and its own in-frame stops
    # would be counted against every sample, perfect matches included, putting the whole
    # cohort in FLAGGED_PREMATURE_STOPS and exiting 1 on an empty clean cohort: the right
    # exit code for the wrong reason, and nothing pointed at the reference.
    ref_is_cds = (not non_coding) and ref_len % 3 == 0
    if not non_coding and ref_len % 3 != 0:
        print(f"[warn] {ref_name} is {ref_len} nt, not a multiple of 3, so it is treated as "
              f"non-coding: the stop and frame rules are off. Pass --non-coding to say so "
              f"deliberately, or give an in-frame CDS.", file=sys.stderr)
    if ref_is_cds:
        ref_stops = find_internal_stops(ref_seq.upper(), ref_len, codon_table)
        if ref_stops:
            raise ValueError(
                f"the reference {ref_name} has {len(ref_stops)} in-frame stop codon(s) under "
                f"genetic code {genetic_code} (codon(s) "
                f"{','.join(str(s) for s in ref_stops[:10])}"
                f"{',...' if len(ref_stops) > 10 else ''}). Every sample would be flagged for "
                f"the reference's own stops. Give an in-frame CDS, select the right code with "
                f"--genetic-code, or pass --non-coding."
            )
    ref_aa = translate_nt(ref_seq, codon_table) if ref_is_cds else ""
    # Only skip the final codon when the reference actually ends in a stop.
    ref_has_term_stop = ref_aa.endswith("*")

    if dataset_id:
        check_dataset_ready(dataset_id,
                            resolve_dataset_state(dataset_id, base_url, api_key, api_key_file))

    print("=" * 80)
    print(" LexicMapStreamer: Multi-HSP Coordinate Tiler & Cohort Partitioner")
    print("=" * 80)
    print(f" Reference Sequence : {ref_name} ({ref_len} bp, {len(ref_aa)} codons)")
    print(f" Reference Header   : {ref_header}")
    print(f" Target QC Settings : full_length_cov={full_cov*100:.1f}%, partial_cov={min_coverage_partial*100:.1f}%, "
          f"min_pident={min_pident:.1f}%, max_stops={max_internal_stops if ref_is_cds else 'N/A'}")

    if dataset_id:
        print(f" Input Source       : Streaming Galaxy dataset: {dataset_id} ({base_url})")
        raw = get_galaxy_stream(dataset_id, base_url, api_key, api_key_file)
    elif file_path == "-":
        print(" Input Source       : Reading from STDIN")
        raw = sys.stdin.buffer
    elif file_path:
        print(f" Input Source       : Local file: {file_path}")
        raw = open(file_path, "rb")
    else:
        raise ValueError("Must specify either --dataset or --file.")
    try:
        stream = text_stream(raw)
    except Exception:
        if raw is not sys.stdin.buffer:
            raw.close()
        raise

    stem = output_stem(ref_name, ref_path)
    if stem != ref_name:
        print(f"[*] Output files are named '{stem}', from the record id {ref_name!r}.")
    out_prefix = resolve_out_prefix(out_prefix, stem)

    # Read header line (skipping comment lines if any)
    while True:
        header_line = stream.readline()
        if not header_line:
            if raw is not sys.stdin.buffer:
                stream.close()
            raise ValueError("Input stream was empty or did not contain a header line.")
        header_line = header_line.strip()
        if header_line and not (header_line.startswith("#") and "\t" not in header_line):
            break

    if header_line.startswith("#"):
        header_line = header_line[1:].strip()

    header_cols = header_line.split("\t")
    col_idx = {c.strip(): i for i, c in enumerate(header_cols)}

    def find_col(candidates):
        for c in candidates:
            if c in col_idx:
                return col_idx[c]
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
    # `sstr` is not read: LexicMap writes a minus-strand hit's sseq already in the query's
    # orientation. On 107,963 HSPs from eight C. auris tables, the 53,193 on the minus strand
    # lined up with qseq as given as often as the plus-strand ones, and never reversed.
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
        # The stream is open and the main try/finally has not been entered yet, so close it
        # here rather than leaving it to the collector.
        if raw is not sys.stdin.buffer:
            stream.close()
        # A KeyError here gave a 15-line traceback for the likeliest user mistake of
        # all, which is running `lexicmap search` without `-a/--all`: its default output
        # carries no qseq/sseq at all. `run` promises one line naming the cause.
        hint = ("" if not ({"qseq", "sseq"} & {m for m in missing}) else
                " The alignment columns come from `lexicmap search -a/--all`; its default "
                "output does not carry them, so re-run the search with that flag.")
        raise ValueError(f"the table is missing the required column(s) "
                         f"{', '.join(missing)}. Header: {header_cols}.{hint}")

    if idx_query is None:
        print(f"[warn] no query-name column in the header ({header_cols}); rows cannot be "
              f"filtered to '{ref_name}'. If this file holds hits for more than one target, "
              f"their HSPs will be merged onto this one reference.", file=sys.stderr)
    if idx_qlen is None:
        print(f"[*] no qlen column; skipping the reference length agreement check.")

    staged = StagedOutputs(out_prefix)
    f_clean_fa = staged.open(".clean.accessions.fasta")
    f_clean_exp_fa = staged.open(".clean_expanded.accessions.fasta")
    f_flagged_fa = staged.open(".flagged.accessions.fasta")
    # `cohort_ledger.tsv` is a tabular rendering of the summary's `records[]`, which carries
    # the same row for every accession. It is therefore **off by default** from 0.2.0 and
    # written only under `--cohort-ledger`, for Galaxy and for spreadsheets. `flagged.tsv` is
    # gone outright: it was a strict subset of this ledger, the rows whose status is not
    # clean. Filter `records[]` (or the ledger) on `status` instead of either file.
    f_cohort_tsv = None
    if write_ledger:
        f_cohort_tsv = staged.open(".cohort_ledger.tsv")
        f_cohort_tsv.write("accession\tcoverage\tmean_pident\tstatus\ttier\tn_stops"
                           "\tstop_codons\thsps_merged\tindels\tconflict_sites\n")
    # One record per accession, appended as the stream runs so the summary costs no more
    # memory than the stage already uses. Spliced into records[] at the end.
    records = summary_io.Records(out_prefix)
    emit_record = records.add

    # Tracking & State. A clean haplotype is keyed by its row AND its trusted indels: an
    # insertion has no column in the row, so two accessions that differ only by one would
    # otherwise collapse into a single haplotype (usually the wild type).
    clean_full_haps = Counter()
    clean_full_pident_sum = defaultdict(float)
    clean_full_samples = defaultdict(list)
    hap_indels = {}   # haplotype key -> its trusted Indel list

    # Counters, not lists: one Python float per accession is about 177 B, so lists would be
    # 384 MB at 2M accessions and break the memory bound this module documents. Coverage is
    # `covered / ref_len` with `covered` an integer in [0, ref_len], so a Counter over it is
    # exact and holds at most ref_len + 1 keys however many accessions there are; identities
    # go in a C array of doubles rather than a list of boxed floats.
    cov_counts = Counter()
    all_pids = array.array('d')
    cov_ge_50_count = 0
    cov_ge_50_stops = 0

    total_rows = 0
    limit_hit = False
    n_malformed = 0
    n_short = 0
    n_inconsistent = 0
    n_other_query = 0
    n_qlen_mismatch = 0
    total_accessions = 0
    count_clean_full = 0
    count_clean_partial = 0
    count_flagged_stops = 0
    count_flagged_low_cov = 0
    count_flagged_frameshift = 0
    count_flagged_mixed = 0
    count_conflict_sites = 0
    count_clean_full_with_indels = 0
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
        nonlocal count_flagged_stops, count_flagged_low_cov, count_flagged_frameshift
        nonlocal count_flagged_mixed, count_conflict_sites
        nonlocal count_clean_full_with_indels
        nonlocal total_rescued_multihsp, total_single_hsp_full
        nonlocal cov_ge_50_count, cov_ge_50_stops

        if not acc or not hsps:
            return

        if acc in seen_accessions:
            raise ValueError(
                f"Accession '{acc}' reappears after its block was already processed, so the "
                f"input is not grouped by accession. All HSPs of an accession must be "
                f"contiguous. Sort the data rows, keeping the header where it is:\n"
                f"  (head -n1 in.tsv; tail -n +2 in.tsv | "
                f"sort -t$'\\t' -k{idx_sgenome + 1},{idx_sgenome + 1} -s) > sorted.tsv\n"
                f"then re-run on sorted.tsv. A bare `sort` will not do: under LC_ALL=C it "
                f"puts the header last, so the first data row is read as the header and the "
                f"re-run dies in the schema check."
            )
        seen_accessions.add(acc)

        total_accessions += 1

        tiling = tile_accession_hsps(hsps, ref_len, min_pident, full_cov,
                                     min_conflict_hsp_bp)
        merged_seq, mean_pid = tiling.seq, tiling.mean_pident
        single_ge_cov, single_full = tiling.single_ge_min_cov, tiling.single_full

        cov = tiling.covered / ref_len
        cov_counts[tiling.covered] += 1
        if cov > 0:
            all_pids.append(mean_pid)

        trusted = [e for e in tiling.indels if not e.edge]
        _in_frame, frameshifts = split_by_frame(tiling.indels)
        indel_str = ",".join(e.notation() for e in tiling.indels) if tiling.indels else "None"
        indel_tag = f" indels={indel_str}" if tiling.indels else ""

        # Check premature stop codons in coding body
        stop_codons = (find_internal_stops_with_indels(merged_seq, ref_len, tiling.indels,
                                                       codon_table, ref_has_term_stop)
                       if ref_is_cds else [])
        n_stops = len(stop_codons)
        stop_str = ",".join(str(p) for p in stop_codons) if stop_codons else "None"

        if cov >= 0.50:
            cov_ge_50_count += 1
            if n_stops > 0:
                cov_ge_50_stops += 1

        # Stratify accession. With --require-full-length, the clean full-length tier requires
        # complete coverage of the reference (every position tiled), not merely min_coverage.
        # A frameshift is judged on the indels themselves: the row stays on the reference grid,
        # so it never shows one as a run of stop codons. Only a shift lasting
        # FRAMESHIFT_MIN_SPAN bases counts. --allow-frameshift-alts waives that tier only;
        # premature stops still flag.
        # Positions where this accession's own HSPs placed different bases. The
        # tiler picked one by depth, but which allele the sample carries is not something the
        # tiler can know, so a conflicted row is reported rather than published as clean.
        conflicts = tiling.conflicts
        conflict_str = (",".join(str(q + 1) for q in conflicts[:CONFLICT_REPORT_LIMIT])
                        + (",..." if len(conflicts) > CONFLICT_REPORT_LIMIT else "")) if conflicts else "None"
        count_conflict_sites += len(conflicts)

        stops_ok = n_stops <= max_internal_stops
        # The frame rule applies to a coding reference only. For an rRNA, ITS, UTR
        # or intergenic marker the tool already switches the stop scan off and has
        # NonCoding_* haplotype classes, but any accession with a 1 or 2 nt indel lasting 30
        # bases was still filed FLAGGED_FRAMESHIFT at exit 0, so the indel allele could never
        # reach the clean cohort and NonCoding_Variant could never describe an indel.
        frame_ok = (not ref_is_cds) or not frameshifts or allow_frameshifts
        mixed_ok = not conflicts or allow_mixed
        if tiling.covered and cov >= full_cov and stops_ok and frame_ok and mixed_ok:
            status = "CLEAN_FULL_LENGTH"
            tier = 1
            count_clean_full += 1
            if single_full:
                total_single_hsp_full += 1
            if not single_ge_cov:
                total_rescued_multihsp += 1
            if trusted:
                count_clean_full_with_indels += 1

            key = (merged_seq, tuple(e.notation() for e in trusted))
            clean_full_haps[key] += 1
            clean_full_pident_sum[key] += mean_pid
            hap_indels.setdefault(key, trusted)
            samples = clean_full_samples[key]
            if len(samples) < sample_cap and acc not in samples:
                samples.append(acc)

            f_clean_fa.write(f">{acc} cov={cov:.4f} pid={mean_pid:.2f} status=CLEAN_FULL_LENGTH{indel_tag}\n{merged_seq}\n")
            f_clean_exp_fa.write(f">{acc} cov={cov:.4f} pid={mean_pid:.2f} status=CLEAN_FULL_LENGTH{indel_tag}\n{merged_seq}\n")

        elif tiling.covered and cov >= min_coverage_partial and stops_ok and frame_ok and mixed_ok:
            status = "CLEAN_PARTIAL"
            tier = 2
            count_clean_partial += 1
            f_clean_exp_fa.write(f">{acc} cov={cov:.4f} pid={mean_pid:.2f} status=CLEAN_PARTIAL{indel_tag}\n{merged_seq}\n")

        # Coverage decides before the frame and stop rules do. Below the partial threshold
        # there is too little of the gene on the grid to tell a real frameshift or premature
        # stop from a gap in the tiling: judging frame first would put 61,897 of 64,034
        # SARS-CoV-2 frameshift calls at a median coverage of 0.18. In this order the flagged
        # tiers mean what they say: covered well enough to judge, and broken.
        #
        # Nothing tiled is a fragment whatever the thresholds say: `--min-cov-partial 0`
        # would otherwise admit an accession whose HSPs all failed --min-pident, and its
        # all-gap row would be counted as a `WT_partial` haplotype.
        else:
            # One writer for the four flagged tiers; they differed only in the label and in
            # what each adds after it, and the same line was written out four times.
            def flag(label, extra=""):
                f_flagged_fa.write(f">{acc} cov={cov:.4f} pid={mean_pid:.2f} flag={label}"
                                   f"{extra}{indel_tag}\n{merged_seq}\n")

            if not tiling.covered or cov < min_coverage_partial:
                status, tier = "LOW_COVERAGE_FRAGMENT", 4
                count_flagged_low_cov += 1
                flag("LOW_COVERAGE")
            elif not frame_ok:
                status, tier = "FLAGGED_FRAMESHIFT", 5
                count_flagged_frameshift += 1
                flag("FRAMESHIFT")
            elif not stops_ok:
                status, tier = "FLAGGED_PREMATURE_STOPS", 3
                count_flagged_stops += 1
                flag("PREMATURE_STOPS", f" stops={n_stops} codons={stop_str}")
            else:  # covered, in frame, no stops: mixed_ok is false by elimination
                status, tier = "FLAGGED_MIXED_HSPS", 6
                count_flagged_mixed += 1
                flag("MIXED_HSPS", f" conflicts={len(conflicts)} sites={conflict_str}")

        if f_cohort_tsv is not None:
            f_cohort_tsv.write(f"{acc}\t{cov:.4f}\t{mean_pid:.2f}\t{status}\t{tier}"
                               f"\t{n_stops}\t{stop_str}\t{len(hsps)}\t{indel_str}"
                               f"\t{conflict_str}\n")
        # The same row, into the summary's records[]. This is the per-accession layer the
        # walker has always had in its own summary and this pathway kept in a side file, so
        # a consumer needed two parsers to ask one question (docs/summaries.md).
        emit_record({
            "acc": acc,
            "status": status,
            "tier": tier,
            "coverage": round(cov, 4),
            "mean_pident": round(mean_pid, 2),
            "hsps_merged": len(hsps),
            "n_stops": n_stops,
            "stop_codons": stop_str,
            "indels": indel_str,
            "n_conflicts": len(conflicts),
            "conflict_sites": conflict_str,
            "seq_sha256": summary_io.seq_sha256(merged_seq),
            "length": len(merged_seq.replace("-", "")),
        })

    def check_download_complete():
        """Raise if a Galaxy body stopped short of its Content-Length.

        Called from two places. A cut lands either between rows, where only the byte count
        reveals it, or inside one, where the short-row check notices first. When it is a
        download, the byte count is the better message: it says what went wrong rather than
        describing the shape of the damage, so it is consulted before that check raises.
        """
        owed = getattr(raw, "length", None) if dataset_id else None
        if owed:
            raise RuntimeError(
                f"the download of Galaxy dataset {dataset_id} ended {owed:,} bytes short of "
                f"its Content-Length (the connection closed early), so the table is truncated "
                f"and every output below it would be incomplete. Rerun it."
            )

    # Hoisted: this was rebuilt for every row.
    max_needed = max([idx_sgenome, idx_pident, idx_qcov, idx_qstart, idx_qend, idx_qseq,
                      idx_sseq] + [i for i in (idx_query, idx_qlen) if i is not None])
    try:
        while True:
            line = stream.readline()
            if not line:
                break
            # `strip()` also removed empty trailing fields, so a row whose last column
            # is empty looked short. Only the line ending comes off.
            unterminated = not line.endswith("\n")
            line = line.rstrip("\r\n")
            if not line.strip():
                continue

            total_rows += 1
            parts = line.split("\t")
            if len(parts) <= max_needed:
                # This used to be a silent `continue`, with no counter and no message,
                # unlike every other skip. A plain table cut inside its last row therefore
                # demoted its final accession at exit 0, and a table whose rows were all short
                # exited 1 reporting an empty clean cohort instead of a schema mismatch.
                if unterminated:
                    check_download_complete()
                    raise ValueError(
                        f"the input looks truncated: its last row has {len(parts)} column(s), "
                        f"fewer than the {max_needed + 1} the header declares, and the file "
                        f"does not end in a newline. Re-fetch it."
                    )
                n_short += 1
                if n_short <= MALFORMED_REPORT_LIMIT:
                    print(f"[warn] row {total_rows}: {len(parts)} column(s), fewer than the "
                          f"{max_needed + 1} the header declares, skipping: {line[:120]}",
                          file=sys.stderr)
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

            # Nothing used to compare `qseq` with the reference slice it is projected
            # onto, so a reference carrying the search's query NAME and LENGTH but another
            # sequence tiled every accession against the wrong sequence at exit 0, and a
            # hand-edited or truncated row tiled partially without a word.
            if (len(qseq) != len(sseq) or qstart < 1 or qend > ref_len
                    or qseq.replace('-', '').upper() != ref_seq[qstart - 1:qend]):
                n_inconsistent += 1
                if n_inconsistent <= MALFORMED_REPORT_LIMIT:
                    print(f"[warn] row {total_rows}: qseq is not the reference over "
                          f"{qstart}-{qend}, skipping: {line[:120]}", file=sys.stderr)
                continue

            hsp = (pident, qcov, qstart, qend, qseq, sseq)

            if acc != curr_acc:
                if limit_hit:
                    total_rows -= 1   # the row that ends the preview is not one we read
                    break
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

            # Finish the open block, then stop before opening another. Cutting through an
            # accession's HSP block judges it from a fraction of its HSPs, so a clean
            # multi-HSP accession comes out LOW_COVERAGE_FRAGMENT, which is the wrong thing
            # for a preview to show.
            if max_records and total_rows >= max_records and not limit_hit:
                print(f"\n[*] Reached row limit (--max-records {max_records:,}); "
                      f"finishing the HSP block of '{curr_acc}'.")
                limit_hit = True

        # A body cut short of its Content-Length reads as a plain end of file:
        # `HTTPResponse.read1` returns b"" rather than raising `IncompleteRead`, deliberately,
        # and `TextIOWrapper` passes that on as EOF, so without this check a half-delivered
        # table is tiled as the whole table and written with `"status": "ok"`. A preview is
        # the one case where bytes are owed on purpose.
        if not limit_hit:
            check_download_complete()

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

        if n_short:
            print(f"[!] {n_short:,} of {total_rows:,} rows had fewer columns than the header "
                  f"declares and were skipped.", file=sys.stderr)
        if n_inconsistent:
            print(f"[!] {n_inconsistent:,} of {total_rows:,} rows carried a qseq that is not "
                  f"the reference over the span they claim, and were skipped. A reference "
                  f"that is not the one the search used would do this to every row.",
                  file=sys.stderr)
            if n_inconsistent > total_rows / 2:
                raise ValueError(
                    f"{n_inconsistent:,} of {total_rows:,} rows disagree with "
                    f"{ref_path} over the reference span they name. This is the wrong "
                    f"reference for this search, not dirty data."
                )
        if n_malformed:
            print(f"[!] {n_malformed:,} of {total_rows:,} rows had an unparseable numeric "
                  f"field and were skipped"
                  f"{' (only the first %d were named)' % MALFORMED_REPORT_LIMIT if n_malformed > MALFORMED_REPORT_LIMIT else ''}.",
                  file=sys.stderr)
            # Every row unusable is a schema mismatch, not dirty data, and the cohort below
            # would be empty rather than incomplete. Fail rather than write a plausible-looking
            # empty output directory.
            if n_malformed + n_short == total_rows:
                raise ValueError(
                    f"all {total_rows:,} data rows had unparseable numeric fields; "
                    f"check the column layout of the input against the resolved schema above"
                )

    except BaseException:
        # Whatever went wrong, this run's half-written outputs do not survive it.
        for fh in (f_clean_fa, f_clean_exp_fa, f_flagged_fa, f_cohort_tsv):
            if fh is not None:
                fh.close()
        staged.abort()
        raise
    finally:
        # For `-f -` this wrapped `sys.stdin.buffer` and skipped closing the wrapper,
        # but never detached it either, so `IOBase.__del__` closed the process's stdin when
        # the wrapper was collected and a library caller reading stdin afterwards got
        # "I/O operation on closed file". The CLI never noticed because fd 0 stays open.
        if raw is sys.stdin.buffer:
            stream.detach()
        else:
            stream.close()
        f_clean_fa.close()
        f_clean_exp_fa.close()
        f_flagged_fa.close()
        if f_cohort_tsv is not None:
            f_cohort_tsv.close()
        records.close()

    total_time = time.time() - start_time

    # Distributions & Diagnostics
    # One quantile helper for the whole suite, so `p50` means `p50` on both pathways. The
    # printed report reads every key, so an empty cohort takes the all-zero row rather than
    # the None that `summary.quantiles` returns; `cohort` in the summary keeps the None,
    # which is the honest answer when nothing was measured.
    cov_stats = summary_io.quantiles_from_counts(cov_counts, 1.0 / ref_len)
    pid_stats = summary_io.quantiles(all_pids)
    med_cov = cov_stats["p50"] if cov_stats else None
    med_pid = pid_stats["p50"] if pid_stats else None
    stop_rate_ge_50 = (cov_ge_50_stops / cov_ge_50_count * 100.0) if cov_ge_50_count > 0 else 0.0

    # An empty list of quantiles must not come back as a row of zeros, or the regime test
    # fires on the placeholders and a header-only table is labelled `divergent` with "median
    # identity is 0.0%". "No data" and "data that is far from the reference" are different
    # answers, and only one of them is a
    # measurement.
    if total_accessions == 0:
        regime, is_divergent = "no_data", False
    else:
        is_divergent = (med_pid is None or med_pid < 85.0
                        or med_cov is None or med_cov < 0.60)
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
    print(f"    2. Clean Partials ({min_coverage_partial:.0%}-{full_cov:.0%} cov): {count_clean_partial:,} accessions")
    print(f"       --> Combined Clean Cohort  : {count_clean_full + count_clean_partial:,} accessions ({pct_of_total(count_clean_full + count_clean_partial):.2f}%)")
    print(f"    3. Flagged for Premature Stops: {count_flagged_stops:,} accessions")
    print(f"    4. Low-Coverage Fragments (<{min_coverage_partial:.0%}): {count_flagged_low_cov:,} accessions")
    print(f"    5. Flagged for Frameshifts    : {count_flagged_frameshift:,} accessions")
    print(f"    6. Flagged for Mixed HSPs     : {count_flagged_mixed:,} accessions"
          f" ({count_conflict_sites:,} conflicting position(s) in all)"
          + ("  [--allow-mixed is on, so these are in the clean tiers]" if allow_mixed else ""))
    print(f"    Clean full-length with in-frame indels: {count_clean_full_with_indels:,} accessions")
    def q(stats, key, fmt):
        return format(stats[key], fmt) if stats else "n/a"

    print(f"    Coverage Distribution (All)   : p10={q(cov_stats,'p10','.2f')}, "
          f"median={q(cov_stats,'p50','.2f')}, p90={q(cov_stats,'p90','.2f')}")
    print(f"    Identity Distribution (All)   : p10={q(pid_stats,'p10','.1f')}%, "
          f"median={q(pid_stats,'p50','.1f')}%, p90={q(pid_stats,'p90','.1f')}%")

    warning_msg = None
    if is_divergent:
        # Name the condition that actually fired. "X and/or Y" used to be printed whichever
        # one it was, so a cohort tripped by coverage alone read as "Median identity is 99.6%
        # (< 85%)", which is simply false and sent people looking at the wrong number.
        fired = []
        if med_pid is None:
            fired.append(f"no HSP reached --min-pident {min_pident:.1f}%")
        elif med_pid < 85.0:
            fired.append(f"median identity is {med_pid:.1f}% (< 85%)")
        if med_cov is None:
            fired.append("nothing was tiled")
        elif med_cov < 0.60:
            fired.append(f"median coverage is {med_cov:.2f} (< 0.60)")
        warning_msg = (
            f"DIVERGENT REGIME DETECTED: {' and '.join(fired)}, with {stop_rate_ge_50:.1f}% premature stops in >=50% covered accessions. "
            f"Reference projection is fragmentary. Clean full-length sequences ({count_clean_full:,}) and clean partials "
            f"({count_clean_partial:,}) have been partitioned from the flagged tiers: {count_flagged_stops:,} with premature stops "
            f"and {count_flagged_frameshift:,} with frameshifts. "
            f"For unbiased whole-cohort haplotype assembly on this target, use disassembler (logan-walker)."
        )
        print("\n" + "!" * 80)
        print(" [!] DIVERGENT TARGET REGIME DETECTED")
        print("!" * 80)
        print(warning_msg)
        print("!" * 80 + "\n")

    sorted_haplotypes = clean_full_haps.most_common()

    # Shannon entropy over the haplotype keys. Every tier 1 accession contributes exactly one
    # key, so these counts sum to count_clean_full and this is the same number the inline loop
    # produced before this helper was shared.
    shannon_entropy = shannon_bits([c for _key, c in sorted_haplotypes])

    # The same cohort measured without the coverage artefact. Additive; the fields
    # above keep their meaning so runs already done stay comparable.
    coverage_independent = coverage_independent_haplotypes(
        clean_full_haps, hap_indels, ref_seq, ref_len)

    msa_clean_fasta = staged.path(".clean.msa.fasta")
    msa_compat_fasta = staged.path(".msa.fasta")
    counts_tsv_path = staged.path(".clean.haplotypes.tsv")
    counts_compat_tsv = staged.path(".haplotypes.tsv")
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
        for rank, ((seq, indels), count) in enumerate(sorted_haplotypes, 1):
            freq = count / count_clean_full if count_clean_full > 0 else 0
            is_wt = (seq == ref_seq and not indels)
            tag = "WT" if is_wt else f"variant_rank_{rank}"
            # An insertion has no column, so the row alone cannot show it; the header does.
            indel_tag = f" indels={','.join(indels)}" if indels else ""
            _w(f">haplotype_{rank:05d} count={count} freq={freq:.5f} tag={tag}{indel_tag}\n{seq}\n")

    if write_a2m:
        with staged.open(".clean.a2m") as f_a2m, staged.open(".a2m") as f_a2m_compat:
            for f in [f_a2m, f_a2m_compat]:
                f.write(f">{ref_name}/1-{ref_len}\n{ref_seq}\n")
                for rank, ((seq, _indels), count) in enumerate(sorted_haplotypes, 1):
                    f.write(f">hap_{rank}_count_{count}/1-{ref_len}\n{seq}\n")

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

        for rank, (key, count) in enumerate(sorted_haplotypes, 1):
            seq, indel_strs = key
            indels = hap_indels[key]
            freq = count / count_clean_full if count_clean_full > 0 else 0
            mean_pid = clean_full_pident_sum[key] / count if count > 0 else 0
            deleted = sum(e.length for e in indels if e.kind == "del")
            cov_pct = (sum(1 for c in seq if c != '-') + deleted) / ref_len * 100.0

            # An ambiguity code is missing information, not a change. `N` used to be
            # reported as a substitution and, through `codon_table.get(codon, 'X')`, as an
            # `X` residue classed Missense.
            nt_muts = []
            for i, (r_char, h_char) in enumerate(zip(ref_seq, seq)):
                if h_char != '-' and h_char != r_char and h_char in "ACGT" and r_char in "ACGT":
                    nt_muts.append(f"{i+1}{r_char}>{h_char}")
            nt_muts.extend(indel_strs)

            aa_muts = []
            if seq == ref_seq and not indels:
                mut_class = "WT"
            elif ref_is_cds:
                # The indels' own protein effect, from the sample's sequence rather
                # than from the reference grid, plus the codons it already accounts for.
                indel_aa, covered_codons = indel_protein_changes(
                    ref_seq, seq, ref_aa, indels, codon_table)
                h_aa = translate_nt(seq, codon_table)
                for r_idx, (r_a, h_a) in enumerate(zip(ref_aa, h_aa)):
                    if r_idx in covered_codons:
                        continue
                    if h_a != '-' and h_a != r_a and h_a != 'X' and r_a != 'X':
                        aa_muts.append(f"{r_a}{r_idx+1}{h_a}")
                aa_muts.extend(indel_aa)
                aa_muts.sort(key=lambda s: int("".join(c for c in s if c.isdigit())[:9] or 0))

                if indels:
                    # The class reports the most severe change. A frameshift only reaches the
                    # clean tiers under --allow-frameshift-alts.
                    if any(a.endswith("*") for a in aa_muts):
                        mut_class = "Nonsense"
                    elif split_by_frame(indels)[1]:
                        mut_class = "Frameshift"
                    else:
                        mut_class = "InFrameIndel"
                elif not aa_muts and not nt_muts:
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
            samples_str = ",".join(clean_full_samples[key])

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

    # The shared contract (docs/summaries.md). Without the provenance blocks, two summaries
    # written at different thresholds are indistinguishable on disk.
    n_flagged = (count_flagged_stops + count_flagged_low_cov + count_flagged_frameshift
                 + count_flagged_mixed)
    n_clean = count_clean_full + count_clean_partial
    summary_data = {
        "schema": summary_io.SCHEMA,
        "tool": summary_io.tool("tiler", __version__),
        "inputs": {
            "references": [
                {"id": ref_name, "length": ref_len, "sha256": summary_io.seq_sha256(ref_seq)}
            ],
            "source": summary_io.source(src_kind, src_id, total_accessions),
        },
        "settings": {
            "min_coverage": min_coverage,
            "min_coverage_partial": min_coverage_partial,
            "min_pident": min_pident,
            "max_internal_stops": max_internal_stops,
            "allow_frameshifts": allow_frameshifts,
            "allow_mixed": allow_mixed,
            "min_conflict_hsp_bp": min_conflict_hsp_bp,
            "genetic_code": genetic_code,
            "non_coding": non_coding,
            "reference_is_cds": ref_is_cds,
            "reference_terminal_stop": ref_has_term_stop,
            "require_full_length": full_length_only,
            "max_records": max_records,
            "samples_per_haplotype": sample_cap,
            "is_cds": ref_is_cds,
        },
        "run": summary_io.run_block(
            status=summary_io.status_of(total_accessions, n_clean),
            started=start_time,
            wall_seconds=total_time,
            n_records=total_accessions,
            n_ok=n_clean,
        ),
        # Every row and every accession the run threw away. All three row counters existed
        # before 0.2.0 and were printed but never recorded, so a run that discarded most of
        # its input and a clean one produced the same summary.
        "discards": {
            "rows": {
                "malformed": n_malformed,
                "other_query": n_other_query,
                "qlen_mismatch": n_qlen_mismatch,
                "too_few_columns": n_short,
                "inconsistent_qseq": n_inconsistent,
            },
            "accessions": {
                "flagged_premature_stops": count_flagged_stops,
                "flagged_low_coverage": count_flagged_low_cov,
                "flagged_mixed_hsps": count_flagged_mixed,
                "flagged_frameshift": count_flagged_frameshift,
            },
        },
        "cohort": {
            "regime": regime,
            "regime_warning": warning_msg,
            "total_streamed_rows": total_rows,
            "clean_full_length": count_clean_full,
            "clean_partial": count_clean_partial,
            "combined_clean_cohort": n_clean,
            "flagged": n_flagged,
            "multi_hsp_rescued_accessions": total_rescued_multihsp,
            "single_hsp_full_length_accessions": total_single_hsp_full,
            "clean_full_length_with_indels": count_clean_full_with_indels,
            "premature_stop_rate_pct_at_50cov": round(stop_rate_ge_50, 2),
            "clean_unique_haplotypes": len(sorted_haplotypes),
            "shannon_entropy_bits": round(shannon_entropy, 4),
            "wildtype_count": clean_full_haps.get((ref_seq, ()), 0),
            "wildtype_fraction": (
                round(clean_full_haps.get((ref_seq, ()), 0) / count_clean_full, 6)
                if count_clean_full > 0 else 0
            ),
            **coverage_independent,
            "cov": summary_io.quantiles_from_counts(cov_counts, 1.0 / ref_len),
            "pident": summary_io.quantiles(all_pids),
            "top_haplotypes": top_preview,
        },
    }
    records.close()
    # Everything written; swap the whole set in at once, then record it.
    staged.commit()
    summary_io.write(summary_json_path, records_file=records.path,
                     **{k: v for k, v in summary_data.items() if k != "schema"})
    records.cleanup()


    print(f"     [✓] Clean MSA FASTA       : {msa_clean_fasta} ({len(sorted_haplotypes) + 1} sequences)")
    print(f"     [✓] Clean Full Accessions : {out_prefix}.clean.accessions.fasta ({count_clean_full:,} accessions)")
    print(f"     [✓] Clean Expanded (>={min_coverage_partial:.0%}): {out_prefix}.clean_expanded.accessions.fasta ({count_clean_full + count_clean_partial:,} accessions)")
    print(f"     [✓] Flagged Accessions    : {out_prefix}.flagged.accessions.fasta ({n_flagged:,} accessions)")
    if write_ledger:
        print(f"     [✓] Complete Cohort Ledger: {out_prefix}.cohort_ledger.tsv "
              f"({total_accessions:,} rows)")
    print(f"     [✓] Summary JSON          : {summary_json_path}")

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

def add_parser(subparsers):
    parser = subparsers.add_parser(
        "tile",
        help="Stream LexicMap tabular hits into a tiled cohort + collapsed haplotypes.",
        description="Stream LexicMap search results, tile multi-HSP alignments, stratify cohorts, and collapse haplotypes.",
    )
    parser.add_argument("-q", "-r", "--query", "--ref", dest="ref", type=str, required=True,
                        help="Reference in-frame CDS FASTA (e.g. queries/query_hiv1_pol.fasta). "
                             "Only the first record is read; its length defines the coordinate "
                             "grid every HSP is projected onto.")
    parser.add_argument("-d", "--dataset", type=str, default=None,
                        help="Galaxy dataset ID to stream from over HTTP. Mutually exclusive "
                             "with -f/-i/--file; one of the two is required.")
    parser.add_argument("-f", "-i", "--file", "--input", dest="file", type=str, default=None,
                        help="Local LexicMap tabular file, or '-' to read from stdin. "
                             "Mutually exclusive with -d/--dataset.")
    parser.add_argument("--base-url", type=str, default="https://usegalaxy.org",
                        help="Galaxy instance base URL used with -d/--dataset "
                             "(default: https://usegalaxy.org).")
    parser.add_argument("--api-key", type=str, default=None,
                        help="Galaxy API key as a literal string. If unset, the key is read "
                             "from --api-key-file, then $GALAXY_API_KEY, then /tmp/gxy.txt.")
    parser.add_argument("--api-key-file", type=str, default=None,
                        help="Path to a file whose contents are the Galaxy API key. Used "
                             "before $GALAXY_API_KEY; with neither, /tmp/gxy.txt is read. "
                             "Ignored if --api-key is set.")
    parser.add_argument("--min-cov", type=fraction, default=0.80,
                        help="Clean full-length cohort gate: fraction of the reference an "
                             "accession's merged HSPs must cover (default: 0.80 = 80%%). "
                             "Below this, an accession falls through to the partial or "
                             "low-coverage tier. Distinct from --msa-min-cov, which gates a "
                             "later MSA stage in the composites.")
    parser.add_argument("--min-cov-partial", type=fraction, default=0.50,
                        help="Clean partial cohort gate: minimum coverage fraction for the "
                             "expanded clean cohort (default: 0.50 = 50%%). Must be "
                             "<= --min-cov. Anything below this is filed as a low-coverage "
                             "fragment.")
    parser.add_argument("--min-pident", type=percentage, default=60.0,
                        help="Minimum percent identity for an HSP to contribute to the tiled "
                             "reconstruction (default: 60.0%%). HSPs below this are dropped "
                             "before the coverage gates are evaluated.")
    parser.add_argument("--max-internal-stops", type=non_negative_int, default=0,
                        help="Maximum premature stop codons tolerated in the CDS body before "
                             "an accession is flagged out of the clean cohort (default: 0).")
    parser.add_argument("--allow-frameshift-alts", action="store_true",
                        help="Let accessions whose tiled sequence carries a frameshifting "
                             "indel into the clean tiers instead of FLAGGED_FRAMESHIFT (a "
                             "reading-frame shift lasting 30 or more bases). "
                             "Premature stops still flag. Off by default.")
    parser.add_argument("--genetic-code", type=positive_int, default=11, help="NCBI genetic code id for translation and the stop scan (default: 11, bacterial and plant plastid). The shipped Candida auris ERG11 marker is a CTG-clade yeast and needs 12; TGA is Trp in code 4 and the mitochondrial codes.")
    parser.add_argument("--non-coding", action="store_true", help="Treat the reference as non-coding: no stop scan, no frame rule, and haplotypes are classed NonCoding_*. A length that is not a multiple of 3 already implies this; pass it for a non-coding target whose length happens to be one, which is one target in three.")
    parser.add_argument("--allow-mixed", action="store_true",
                        help="Let accessions whose own HSPs disagree about a base into the clean tiers instead of FLAGGED_MIXED_HSPS. Logan leaves an unpopped SNP bubble as two contigs, both of which align; the tiler picks the one each position lies deeper inside and records the disagreement, because which allele the sample carries is not something it can know. Off by default.")
    parser.add_argument("--min-conflict-hsp-bp", type=non_negative_int,
                        default=CONFLICT_MIN_HSP_BP,
                        help="Record a conflict only when at least one of the two disagreeing "
                             "HSPs spans this many reference bases "
                             f"(default: {CONFLICT_MIN_HSP_BP}; capped at a tenth of the "
                             "reference, so a reference under 1 kb uses less). A conserved "
                             "gene's termini match short fragments throughout a Logan assembly, "
                             "and two such fragments disagreeing say nothing about the sample's "
                             "allele. A short Logan bubble branch is still heard, because the "
                             "contig spanning the gene is the other voter. 0 lets every pair "
                             "count, which is the pre-0.2.0 behaviour. This changes only which "
                             "disagreements are reported, never which base is tiled.")
    parser.add_argument("--require-full-length", action="store_true",
                        help="Require 100%% reference coverage (every position tiled) for the "
                             "clean full-length tier, instead of the --min-cov threshold. "
                             "Off by default.")
    parser.add_argument("--max-records", type=positive_int, default=None,
                        help="Stop after processing N input rows. Useful for a quick preview "
                             "on a large tabular stream.")
    parser.add_argument("--samples-per-haplotype", type=positive_int, default=10,
                        dest="samples_per_haplotype",
                        help="Maximum representative accessions listed per collapsed haplotype "
                             "in <prefix>.clean.haplotypes.tsv (default: 10). Does not affect "
                             "the .accessions.fasta outputs.")
    out = parser.add_mutually_exclusive_group()
    out.add_argument("-o", "--out-prefix", dest="out_prefix", type=str, default=None,
                     help="Output prefix: all tiler outputs land at <prefix>.*. An existing "
                          "directory, or a path ending in '/', is taken as --out-dir "
                          "(default: msa/<ref_name>).")
    # `--out-dir` was once a third spelling of --out-prefix, so a directory that did not
    # exist yet became a prefix: `--out-dir out/new` wrote out/new.clean.msa.fasta and the rest.
    out.add_argument("--out-dir", dest="out_dir", type=str, default=None,
                     help="Output directory, created if needed: outputs land at "
                          "<dir>/<ref_name>.*.")
    parser.add_argument("--cohort-ledger", action="store_true",
                        help="Also write <prefix>.cohort_ledger.tsv, one row per evaluated "
                             "accession. It is a tabular rendering of the summary's "
                             "records[], which carries the same rows, so it is off "
                             "by default; pass this for Galaxy or a spreadsheet.")
    parser.add_argument("--no-a2m", action="store_true",
                        help="Skip writing the .a2m formatted MSA files "
                             "(<prefix>.clean.a2m and the .a2m compat alias).")
    parser.add_argument("--allow-empty", action="store_true",
                        help="Exit 0 even if the clean cohort is empty. Off by default: "
                             "an empty cohort exits non-zero so a wrapping workflow does "
                             "not mistake a whole-target miss for a good run. Mirrors "
                             "`disassembler walk --allow-empty`.")
    parser.set_defaults(func=run)


EX_TEMPFAIL = 75  # sysexits.h: temporary failure, try again later


def _write_failed_summary(out_prefix: str, args, started: float, error: str) -> None:
    """A summary even when the run failed (`docs/summaries.md` rule 4).

    Until 0.2.0 an input problem raised past the writer, so the only signal was the file's
    absence: a caller could not tell a crash from a run that had not started, and nothing on
    disk said what had failed. The blocks it can still fill are filled.
    """
    if args.dataset:
        src_kind, src_id = "galaxy", str(args.dataset)
    elif args.file and args.file != "-":
        src_kind, src_id = "file", str(args.file)
    else:
        src_kind, src_id = "stdin", "-"
    # The same prefix the successful run would have used. Reading the reference may itself
    # be what failed, in which case there is no record name to build it from.
    try:
        out_prefix = resolve_out_prefix(
            out_prefix, output_stem(load_reference(args.ref)[0], args.ref))
    except (OSError, ValueError, IndexError):
        out_prefix = resolve_out_prefix(out_prefix, "tiler")
    # Nothing half-written outlives a failed run.
    StagedOutputs.clear(out_prefix)
    summary_io.discard_records(out_prefix)
    try:
        summary_io.write(
            f"{out_prefix}.summary.json",
            tool=summary_io.tool("tiler", __version__),
            inputs={"references": None, "source": summary_io.source(src_kind, src_id, 0)},
            settings={
                "min_coverage": args.min_cov,
                "min_coverage_partial": args.min_cov_partial,
                "min_pident": args.min_pident,
                "max_internal_stops": args.max_internal_stops,
                "allow_frameshifts": args.allow_frameshift_alts,
                "allow_mixed": args.allow_mixed,
                "min_conflict_hsp_bp": args.min_conflict_hsp_bp,
                "genetic_code": args.genetic_code,
                "non_coding": args.non_coding,
                "require_full_length": args.require_full_length,
                "max_records": args.max_records,
                "samples_per_haplotype": args.samples_per_haplotype,
            },
            run=summary_io.run_block("failed", started, time.time() - started, 0, 0,
                                     error=error),
            discards={"rows": {}, "accessions": {}},
            records=[],
        )
    except OSError:
        # The output directory itself may be the thing that failed. The caller already
        # printed the real error; losing the summary too must not change the exit code.
        pass


def validate_args(args, prog: str = "tile"):
    """The argument checks both front ends make. `None` when they pass, else an exit code.

    `from-lexicmap` had its own copy, which had already drifted in wording, and a
    check added to one did not reach the other. It still has to run before `--out-dir` is
    created, so it stays a plain function rather than something the stage runner calls.
    """
    if not args.dataset and not args.file:
        print(f"{prog}: must specify either -d/--dataset <GALAXY_ID> or -f/-i/--file <PATH> "
              f"(or '-' for stdin).", file=sys.stderr)
        return 2

    # --require-full-length replaces --min-cov as the tier 1 gate, so a run with
    # `--require-full-length --min-cov-partial 0.9` is well formed and used to be refused
    # with a reason that was not true of it.
    gate = 1.0 if args.require_full_length else args.min_cov
    if args.min_cov_partial > gate:
        print(f"{prog}: --min-cov-partial ({args.min_cov_partial}) must be <= the clean "
              f"full-length gate ({gate}, from "
              f"{'--require-full-length' if args.require_full_length else '--min-cov'}); "
              f"otherwise the partial tier is empty by construction.", file=sys.stderr)
        return 2

    # Fail on a genetic code that does not exist before anything is opened.
    try:
        codon_table_for(args.genetic_code)
    except ValueError as e:
        print(f"{prog}: {e}", file=sys.stderr)
        return 2

    # Galaxy API key file: only validate when streaming from a dataset and the file is the
    # key that will be used (only --api-key wins over it, as in get_galaxy_api_key).
    if args.dataset and not args.api_key:
        key_err = validate_api_key_file(args.api_key_file)
        if key_err:
            print(f"{prog}: {key_err}", file=sys.stderr)
            return 1
    return None


def run(args) -> int:
    bad = validate_args(args)
    if bad is not None:
        return bad

    out_prefix = (os.path.join(args.out_dir, "") if args.out_dir else args.out_prefix)
    started = time.time()
    try:
        summary = stream_lexicmap_to_msa(
            ref_path=args.ref,
            dataset_id=args.dataset,
            file_path=args.file,
            base_url=args.base_url,
            api_key=args.api_key,
            api_key_file=args.api_key_file,
            min_coverage=args.min_cov,
            min_coverage_partial=args.min_cov_partial,
            min_pident=args.min_pident,
            max_internal_stops=args.max_internal_stops,
            allow_frameshifts=args.allow_frameshift_alts,
            allow_mixed=args.allow_mixed,
            min_conflict_hsp_bp=args.min_conflict_hsp_bp,
            genetic_code=args.genetic_code,
            non_coding=args.non_coding,
            full_length_only=args.require_full_length,
            max_records=args.max_records,
            sample_cap=args.samples_per_haplotype,
            out_prefix=out_prefix,
            write_a2m=not args.no_a2m,
            write_ledger=args.cohort_ledger,
        )
    except DatasetNotReady as e:
        # Deliberately no summary here. EX_TEMPFAIL says the input is not ready yet and the
        # caller should retry; no run was attempted, so there is nothing to summarise, and
        # creating the output directory would make an unstarted run look like a failed one.
        # This is the one exception to rule 4 in docs/summaries.md.
        print(f"tile: {e}", file=sys.stderr)
        return EX_TEMPFAIL
    # Input problems (a missing or empty file, a bad table, an unreadable Galaxy dataset, a
    # missing API key) end in one line naming the cause, not a traceback.
    # A cut chunked body raises http.client.IncompleteRead and a cut gzip body raises
    # EOFError. Both were leaving as tracebacks. An HTTPException carries its detail in
    # `repr`, not `str`, so an IncompleteRead would otherwise print as an empty message.
    except (ValueError, OSError, RuntimeError, EOFError, http.client.HTTPException) as e:
        detail = repr(e) if isinstance(e, http.client.HTTPException) else str(e)
        print(f"tile: {detail}", file=sys.stderr)
        _write_failed_summary(out_prefix, args, started, detail)
        return 1

    # Empty-cohort hard-fail: a wrapping workflow should treat "no clean accessions
    # streamed" as failure, not as a good run with an empty output. Mirrors the
    # walker's --allow-empty escape hatch.
    if isinstance(summary, dict) and summary.get("run", {}).get("status") != "failed":
        n_clean = int(summary.get("cohort", {}).get("combined_clean_cohort", 0))
        if n_clean == 0 and not args.allow_empty:
            print(f"tile: clean cohort is empty (0 full-length + 0 partial). "
                  f"Exiting non-zero. Pass --allow-empty to accept this as a good "
                  f"run.", file=sys.stderr)
            return 1
    return 0
