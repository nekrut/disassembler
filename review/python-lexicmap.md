# Python: LexicMap route (cmd/lexicmap.py)

`disassembler lexicmap` (and the `from-lexicmap` composite) streams a LexicMap `search --all` table, tiles each accession's HSPs (high-scoring segment pairs, one per contig alignment) onto the reference grid, sorts accessions into five tiers and collapses tier 1 into haplotypes. The streaming design, the query-name and `qlen` guards, the Galaxy dataset-state check and the tier ordering (coverage before frame and stop rules) are sound and well tested. The weak point is the tiler: when HSPs of one accession overlap, a single whole-HSP identity decides every base, so a reference-matching contig hides variants, stops and frameshifts (PL-01), and indel ownership and edge status are decided from one HSP only (PL-06, PL-07). The reference itself is barely validated: coding status is inferred from `length % 3`, only genetic code 11 exists, whitespace and the last codon are mishandled (PL-04, PL-05, PL-08, PL-09). Summary statistics mix assembly fragmentation and placeholder values into diversity and regime calls (PL-02, PL-15, PL-16), and several input edges (a cut Galaxy download, short rows, `--max-records`) are accepted at exit 0. The common theme of the fixes is to stop accepting silently: count, flag or refuse where the code now proceeds as if nothing happened.

| ID | Severity | Where | Title | Fix effort |
|---|---|---|---|---|
| PL-01 | high | `src/disassembler/cmd/lexicmap.py:363` | Overlapping reference-matching HSPs silently override variants, stops and frameshifts | medium |
| PL-02 | high | `src/disassembler/cmd/lexicmap.py:761` | Diversity summary counts coverage differences as distinct haplotypes | small |
| PL-03 | high | `src/disassembler/cmd/lexicmap.py:811` | `-d` reads a Galaxy download cut short of Content-Length as the whole table | small |
| PL-04 | medium | `src/disassembler/cmd/lexicmap.py:537` | Coding status from length alone; genetic code 11 only; reference never checked for stops | small |
| PL-05 | medium | `src/disassembler/cmd/lexicmap.py:164` | Last codon is never stop-scanned, even when the reference has no terminal stop | trivial |
| PL-06 | medium | `src/disassembler/cmd/lexicmap.py:410` | Indel inside an HSP overlap is marked an edge artefact, so a frameshift can pass as clean | trivial |
| PL-07 | medium | `src/disassembler/cmd/lexicmap.py:389` | Untrusted (edge) deletion's `-` displaces an overlapping HSP's real base | small |
| PL-08 | medium | `src/disassembler/cmd/lexicmap.py:221` | Whitespace inside reference sequence lines becomes reference positions | trivial |
| PL-09 | medium | `src/disassembler/cmd/lexicmap.py:749` | Frameshift tier applied to references classified as non-coding | trivial |
| PL-10 | medium | `src/disassembler/cmd/lexicmap.py:578` | Raw FASTA record id used as a file name: `/` ids crash, `../` or absolute ids escape the output directory | small |
| PL-11 | medium | `src/disassembler/cmd/lexicmap.py:56` | `clean.a2m` is a nucleotide alignment, but docs offer it to protein language models | trivial |
| PL-12 | medium | `src/disassembler/cmd/lexicmap.py:710` | Ungrouped input aborts mid-stream leaving misclassified files; the advised `sort` breaks the header | small |
| PL-13 | low | `src/disassembler/cmd/lexicmap.py:890` | `--max-records` cuts through an accession's HSP block and mis-tiers it | trivial |
| PL-14 | low | `src/disassembler/cmd/lexicmap.py:821` | Rows with too few columns are dropped silently | small |
| PL-15 | low | `src/disassembler/cmd/lexicmap.py:178` | Quantiles pick a rounded index; even-cohort medians are wrong and the regime call flips | trivial |
| PL-16 | low | `src/disassembler/cmd/lexicmap.py:941` | Empty or untiled cohort is called "divergent" with invented 0.0 medians | small |
| PL-17 | low | `src/disassembler/cmd/lexicmap.py:452` | Frameshift span counts deleted reference bases | trivial |
| PL-18 | low | `src/disassembler/cmd/lexicmap.py:785` | `--min-cov-partial 0` files zero-coverage accessions as clean all-gap rows | trivial |
| PL-19 | low | `src/disassembler/cmd/lexicmap.py:1284` | Under `--require-full-length`, `--min-cov` still drives a check, labels and the rescue count | small |
| PL-20 | low | `src/disassembler/cmd/lexicmap.py:868` | A row's `qseq` is never checked against the reference slice it is projected onto | small |
| PL-21 | low | `src/disassembler/cmd/lexicmap.py:636` | Missing `qseq`/`sseq` and truncated or corrupt gzip end in tracebacks | small |
| PL-22 | low | `src/disassembler/cmd/lexicmap.py:646` | Outputs written in place; failed or `--no-a2m` reruns leave new and stale files mixed | small |
| PL-23 | low | `src/disassembler/cmd/lexicmap.py:252` | Galaxy `deferred` and `paused` datasets exit 75 "try again later" | trivial |
| PL-24 | low | `src/disassembler/cmd/lexicmap.py:693` | Memory grows about 180 B per accession, contrary to the documented model | small |
| PL-25 | low | `src/disassembler/cmd/lexicmap.py:343` | Tiling builds a tuple per reference position; slice assignment is about 3.5x faster | small |
| PL-26 | low | `src/disassembler/cmd/lexicmap.py:924` | `stream_lexicmap_to_msa(file_path="-")` closes the process's `sys.stdin` | trivial |
| PL-27 | low | `src/disassembler/cmd/lexicmap.py:128` | Galaxy key lookup and HTTP access duplicated; tool-not-found messages and exit codes differ | small |
| PL-28 | low | `src/disassembler/cmd/lexicmap.py:426` | Test-only wrapper, a counter never read, triplicated flagged-row writes | small |
| PL-29 | low | `src/disassembler/cmd/lexicmap.py:1068` | Haplotype table never reports the protein change of an in-frame indel | small |

### PL-01 · high · Overlapping reference-matching HSPs silently override variants, stops and frameshifts

`src/disassembler/cmd/lexicmap.py:363` · science · confirmed (reproduced) · fix effort: medium · sources: SC2-01, P1-01

**Problem.** `tile_accession_hsps` resolves overlapping HSPs of one accession per position with `if ex is None or pident > ex[1]: grid[q_pos] = (sc.upper(), pident, k)` (`lexicmap.py:363` for ungapped HSPs, `:388-389` for gapped ones). `pident` is the HSP's identity to the reference, so wherever two contigs disagree the one closer to the reference wins every base. A deletion carried by the losing HSP is then dropped as "not owned" (`:404-409`), so premature stops and frameshifts disappear as well. Nothing records that the HSPs disagreed; the accession is tiered and collapsed as if its row were one haplotype. `tests/test_lexicmap.py:49` pins the rule.

**Impact.** Logan assemblies leave unpopped SNP bubbles as separate contigs: both 61-nt branches (2k-1) pass LexicMap's default `-l 50 -i 70`, and the reference branch (100% identity) always beats the alternative one (98.4%) regardless of abundance. Mixed samples become `CLEAN_FULL_LENGTH` wild type, which biases `wildtype_count` and haplotype frequencies towards the reference; a stop or 1-nt deletion next to a reference-matching contig never reaches tiers 3 or 5; a homolog overlapping a truncated ortholog yields a clean Missense chimera. No output column reveals any of it. An emulation of LexicMap HSPs (≥50 nt, ≥70%) on the real Logan contigs of the 9 Zika accessions used in the docs gave 37-46 wrong bases per accession, mostly spurious 1-3% alleles from bubble branches, with 7 of 9 accessions flagged and the 2 clean rows carrying 43-44 mutations against 6 expected (real LexicMap seeding may report fewer short contigs).

**Evidence.** A synthetic full-column table on `queries/markers/zika_ns5.fasta`, run with `disassembler lexicmap -q queries/markers/zika_ns5.fasta -f hits.tsv -o out/z`: accessions carrying a SNP, a stop or a 1-nt deletion alone are tiered correctly (`CLEAN_FULL_LENGTH D334H`, `FLAGGED_PREMATURE_STOPS 341`, `FLAGGED_FRAMESHIFT 1501del`); the same three changes with a 61-nt reference-matching bubble contig, or a full contig plus a 61-nt reference contig, all come out `CLEAN_FULL_LENGTH`, indels `None`, and six of them collapse into the wild-type haplotype (count 7).

**Proposed fix.** Make the resolution position-aware and record disagreement instead of hiding it. Two prototypes were built; the choice between them needs a test on real Logan HSPs:
- Depth rule plus conflict flag: the owner of a position is the HSP it lies deepest inside (distance to the nearer HSP end), then `pident`, then input order. Where two HSPs place different symbols at a position at least `INDEL_EDGE_WINDOW` inside both, or an interior insertion is spanned by another HSP without it, record a conflict. Return `conflicts` in `Tiling`, add a `conflict_sites` column to the ledger and `flagged.tsv`, a tier 6 `FLAGGED_MIXED_HSPS` after the frame and stop rules, an `--allow-mixed` option (also on `from-lexicmap`) and a count in `summary.json`. On the fixture above, the bubble cases become `FLAGGED_MIXED_HSPS`, the masked stop and deletion return to tiers 3 and 5, and the haplotype table holds WT x1, 1000G>C x1; all 156 existing Python tests pass.
- Per-position majority vote across the HSPs covering a position, writing `N` on a tie and counting ties. On the emulated Zika data this gave 0-1 wrong bases and 9/9 accessions clean; two existing tests that pin the old rule fail and must be updated.

Ranking whole HSPs by length or bitscore instead of `pident` does not help: bubble branches tie on length, bitscore favours the reference, and on the emulated data it gave 45-302 wrong bases.

**Tests and docs.** Add tests for a bubble, stop and deletion masking, the override option, and an end mismatch that must not count as a conflict. Update `docs/lexicmap.md` step 1, the tier, output and option tables, the module docstring (`lexicmap.py:9-11`, `:21-26`), the "Flagged Accessions" count line (`:1174`) and the tools-iuc #8800 wrapper help. Expect more accessions to leave tier 1 on shattered viral cohorts, which can empty the clean cohort (exit 1 without `--allow-empty`). Related: PL-06, PL-07.

### PL-02 · high · Diversity summary counts coverage differences as distinct haplotypes

`src/disassembler/cmd/lexicmap.py:761` · science · confirmed (reproduced) · fix effort: small · sources: P1-04, SC2-06

**Problem.** Tier 1 rows are collapsed on the whole reference-length row, which carries `-` at every untiled position: `key = (merged_seq, tuple(e.notation() for e in trusted))` (`lexicmap.py:761`). With the default `--min-cov 0.80` a tier 1 row can be up to 20% `-` (`docs/lexicmap.md:61`). Two accessions with identical bases wherever covered but different contig ends become separate haplotypes. `clean_unique_haplotypes` (`:1162`), `shannon_entropy_bits` (`:989-995`, computed over those keys) and `wildtype_fraction` (`:1165`, exact `(ref_seq, ())` key only) therefore measure assembly fragmentation as well as sequence diversity. The table labels such rows `WT_partial` (`:1087`), but nothing in `summary.json` or the docs warns that the cohort diversity numbers include coverage patterns.

**Impact.** Default options, whenever tier 1 accessions differ in coverage between 0.80 and 1.0, which is common because Logan contigs often stop short of gene ends. Entropy and haplotype counts go up and the wild-type fraction goes down, so reported diversity tracks sequencing depth and cannot be compared across targets or cohorts. Each coverage pattern is also stored as its own key, row and sample list. `--require-full-length` avoids the split but drops those accessions from tier 1.

**Evidence.** A 600-nt CDS with 5 full-length wild-type accessions, 3 wild-type accessions tiled over 1-570, 31-600 and 1-540, and 1 full-length accession with 300G>A, run through `disassembler lexicmap -q ref.fa -f hits.tsv -o o/def`: `clean_unique_haplotypes 5, shannon_entropy_bits 1.88, wildtype_fraction 0.556`; haplotypes 2-4 are single wild-type accessions at 95/95/90% coverage. Expected: 2 haplotypes, 0.503 bits, wild-type fraction 0.889.

**Proposed fix.** Keep the existing fields for compatibility and add coverage-independent ones in `stream_lexicmap_to_msa`, after `sorted_haplotypes`:
1. Full-coverage subset: keys whose non-`-` bases plus trusted deletions reach `ref_len`. Report `full_coverage_rows`, `clean_unique_haplotypes_full_coverage`, `shannon_entropy_bits_full_coverage`.
2. Fold: each partial key goes to the single full-coverage key with identical trusted indels that agrees at every covered base. A partial with no such key, or several, is counted in `partial_rows_unresolved` and never guessed. Report `partial_rows_folded`, `shannon_entropy_bits_folded`, `wildtype_compatible_count`, `wildtype_fraction_folded`.

A prototype (about 40 lines) gives, on the fixture, 2 full-coverage haplotypes at 0.65 bits and folded results of 3 rows folded, 0 unresolved, 0.503 bits and a fraction of 0.889; the 55 existing lexicmap tests pass. The naive fold is O(partials x fulls x L); for large cohorts index full keys by the substring over each partial's covered span. Do not collapse on variant signatures alone: a partial row would be called wild type whenever its variant falls in the uncovered region.

**Tests and docs.** Add tests: four wild-type accessions tiled over 1-60, 1-57, 4-60 and 2-60 give 1 full-coverage haplotype and folded entropy 0; a partial row compatible with both WT and an SNV haplotype stays unresolved. In `docs/lexicmap.md` step 5 (`:72-74`) state that the plain counts depend on coverage without `--require-full-length`, and list the new fields at `:139`. Decide whether to redefine `shannon_entropy_bits` in 0.2 (changelog note). Related: PL-07 (an edge deletion's `-` also splits haplotypes).

### PL-03 · high · `-d` reads a Galaxy download cut short of Content-Length as the whole table

`src/disassembler/cmd/lexicmap.py:811` · correctness · confirmed (reproduced) · fix effort: small · sources: PC-01, P1-06

**Problem.** `get_galaxy_stream` returns the raw `http.client.HTTPResponse` (`lexicmap.py:246`). It has `peek`, so `textio.text_stream` (`src/disassembler/textio.py:26-32`) wraps it directly and `TextIOWrapper` reads it through `HTTPResponse.read1`. When the connection closes before Content-Length bytes arrive, `read1` returns `b''` (CPython does this deliberately for compatibility) and leaves `raw.length` at the bytes still owed. The read loop (`lexicmap.py:810-812`: `line = stream.readline()` / `if not line:` / `break`) takes that as the end of the table. Nothing compares received bytes with Content-Length, so the partial cohort is tiled, `summary.json` says `"status": "ok"` and the exit status is 0. A cut chunked body raises `http.client.IncompleteRead` and a cut gzip body raises `EOFError`; `run()` catches only `(ValueError, OSError, RuntimeError)` (`:1323`), so those end in tracebacks.

**Impact.** `disassembler lexicmap|from-lexicmap -d <ID>` on a large dataset when the server or a proxy closes cleanly mid-transfer (nginx "upstream prematurely closed", a Galaxy worker restart, a proxy time limit; over HTTPS a close without close_notify reads the same way). The dataset named at `lexicmap.py:71` is served by usegalaxy.org with `content-length: 9824523808`, which keeps one connection open for 15 minutes or more. The user gets the first N accessions in table order, a biased subset; the cut row can be mis-tiered; `from-lexicmap` builds the MSA and tree from it. Nothing in any output shows it. A connection reset or a stall over 120 s raises `OSError` and is reported. The tools-iuc wrapper uses `--file`, so Galaxy tool users are not affected.

**Evidence.** A local stand-in server declared the full Content-Length of a 1,667,717-byte, 600-accession `--all` table, sent 668,496 bytes and closed. `lexicmap -d` exited 0 with empty stderr, `status=ok`, 241 accessions (full run: 600), and the last ledger row `SRR9000240 0.5000 CLEAN_PARTIAL` (full run: `1.0000 CLEAN_FULL_LENGTH`). `from-lexicmap` exited 0 and wrote a 241-record MSA and a tree. The same over HTTPS with a ragged close: exit 0, 241 accessions. An independent reproduction with a different server gave the same result (370 of 1000 accessions, exit 0).

**Proposed fix.** After the read loop, fail if the response still owes bytes, except when `--max-records` stopped the loop:
```python
def check_download_complete(raw, dataset_id):
    owed = getattr(raw, "length", None)        # bytes still owed by Content-Length
    if owed:
        raise RuntimeError(f"the download of Galaxy dataset {dataset_id} stopped "
                           f"{owed:,} bytes short of its Content-Length; rerun it")

# in the loop, set stopped_early = True just before the --max-records `break` (:890-892)
# after the loop:
if dataset_id and not stopped_early:
    check_download_complete(raw, dataset_id)
# run(), :1323
except (ValueError, OSError, RuntimeError, EOFError, http.client.HTTPException) as e:
```
`length` is `None` for `BytesIO`, chunked and close-delimited bodies, so the check fires only on a short Content-Length body. With the fix the cut run exits 1 with "stopped 999,221 bytes short of its Content-Length" over HTTP and HTTPS; the chunked and gzip cuts give one-line errors; complete downloads and `--max-records 100` still exit 0. Later hardening: compare a byte count with `file_size` from `/api/datasets/{id}` (already fetched by `resolve_dataset_state`) for close-delimited bodies, and resume with `Range` plus `If-Range: <ETag>`.

**Tests and docs.** Add `tests/test_galaxy_truncation.py` with a real local HTTP server: complete download, body short of Content-Length, chunked body cut short, and `--max-records`; two of these fail on 0.1.1. Existing `TestGalaxyDataset` stubs use `BytesIO` and are unaffected. Add to `docs/lexicmap.md:87-91`: "A download that ends before the size the server declared exits 1; rerun it." Related: PL-13 (changes the same loop), PL-21 (uncaught `EOFError`), PL-22 (a failed rerun leaves stale files).

### PL-04 · medium · Coding status from length alone; genetic code 11 only; reference never checked for stops

`src/disassembler/cmd/lexicmap.py:537` · correctness · confirmed (reproduced) · fix effort: small · sources: P1-02, SC2-05

**Problem.** The only test of the reference is `ref_is_cds = (ref_len % 3 == 0)` (`lexicmap.py:537`); the reference is never translated and checked. `CODON_TABLE_11` (`:109`) is the only codon table: `translate_nt`, `find_internal_stops` and `find_internal_stops_with_indels` default to it and no option selects another (`docs/lexicmap.md:40` says "genetic code 11"). Non-ACGT bases are also called as changes: `if h_char != '-' and h_char != r_char:` (`:1060`) and `codon_table.get(codon, 'X')` (`:155`) turn an `N` into an `X` residue classed Missense.

**Impact.** A non-coding target whose length is a multiple of 3 (one in three of them, e.g. *E. coli* 16S at 1542 nt), or a CDS that uses another code (TGA = Trp in code 4 and mitochondrial codes; selenoproteins), has its own stops counted against every sample: all accessions, perfect matches included, land in `FLAGGED_PREMATURE_STOPS` and the run exits 1 with "clean cohort is empty", which names the wrong cause. The shipped `queries/markers/candida_auris_erg11.fasta` is from a CTG-clade yeast (code 12, CTG = Ser) and has CTG at codon 74, so that residue is labelled L instead of S and TTG to CTG changes are reported `Synonymous`. The other shipped references are code-11 CDSs and are unaffected.

**Evidence.** A 63-nt non-coding reference with two 100% hits: exit 1, both `FLAGGED_PREMATURE_STOPS 3 5 2,4,6,14,20`. A reference with an in-frame TGA and five identical hits: exit 1, all `FLAGGED_PREMATURE_STOPS ... 22`. A CTG reference: `61T>C WT Synonymous`, and an `N` gives `91T>N C31X Missense`.

**Proposed fix.**
```python
codon_table = codon_table_for(genetic_code)   # from Bio.Data.CodonTable.unambiguous_dna_by_id[N]
ref_is_cds = (not non_coding) and ref_len % 3 == 0
if not non_coding and not ref_is_cds:
    warn("reference length not a multiple of 3: treating it as non-coding")
if ref_is_cds and (st := find_internal_stops(ref_seq.upper(), ref_len, codon_table)):
    raise ValueError(f"reference has {len(st)} in-frame stop(s) under code {genetic_code} "
                     f"(codons {st}); give an in-frame CDS, --genetic-code or --non-coding")
```
Pass `codon_table` to `translate_nt` (`:538`, `:1068`) and `find_internal_stops_with_indels` (`:731`); add `--genetic-code N` (default 11, validated against `unambiguous_dna_by_id`) and `--non-coding` to `lexicmap` and `from-lexicmap` (forward at `pipe_from_lexicmap.py:281`); record `genetic_code` in `summary.json`. Call nucleotide changes only on `ACGT` and skip `X` in amino-acid calls. Prototypes of both parts passed the Python suite: the non-coding reference fails fast naming its stops, `--non-coding` gives 2/2 clean, the TGA CDS is clean under `--genetic-code 4`.

**Tests and docs.** Tests: reference with internal stops raises; `--non-coding` is clean; code 4 with TGA is clean; TTG to CTG is `Missense L21S` under code 12; an `N` is not a mutation. Docs: `docs/lexicmap.md:40`, `:70` (how coding is decided) and the options table at `:178`; note that the C. auris marker needs `--genetic-code 12`. Related: PL-09 (frame rule on non-coding references), PL-05.

### PL-05 · medium · Last codon is never stop-scanned, even when the reference has no terminal stop

`src/disassembler/cmd/lexicmap.py:164` · correctness · confirmed (reproduced) · fix effort: trivial · sources: SWP1-02

**Problem.** Both stop scans stop before the last codon unconditionally: `for i in range(0, ref_len - 3, 3):` (`lexicmap.py:164`) and `for i in range(0, len(sample) - 3, 3):` (`:500`), assuming the reference ends in its natural stop. Six shipped references end in a sense codon because they are polyprotein segments: `queries/markers/dengue_1_ns5.fasta`, `dengue_2_ns5`, `dengue_3_ns5` (TGG), `dengue_4_ns5`, `westnile_ns5`, `zika_ns5` (CTG).

**Impact.** A nonsense change in the last codon of those references passes the tier rule, while the haplotype table translates the full row and classes the same haplotype `Nonsense` (`:1091`), so the outputs contradict each other. Dengue NS5 ends TGG, one G>A away from TAG or TGA.

**Evidence.** On `zika_ns5.fasta`, an accession with TAA at codon 903 is `CLEAN_FULL_LENGTH 1 0 None`, yet the top clean haplotype reads `Nonsense L903*`; the same change at codon 11 is `FLAGGED_PREMATURE_STOPS 3 1 11`.

**Proposed fix.** Decide once from the reference: `ref_has_term_stop = ref_aa.endswith("*")` after `:538`. Give `find_internal_stops` and `find_internal_stops_with_indels` a `skip_last=True` parameter and scan `range(0, n - 3 if skip_last else n - 2, 3)`; pass `skip_last=ref_has_term_stop` at `:731` and through the early return at `:473`. With it, the codon-903 case is `FLAGGED_PREMATURE_STOPS 3 1 903`. Optionally record `reference_terminal_stop` in `summary.json`.

**Tests and docs.** `find_internal_stops("ATGGCTGAATAA", 12, skip_last=False) == [4]`, and an end-to-end test with reference `ATGGCTGAACTG` and a sample ending TAA expecting tier 3 at codon 4. Update the docstring at `:161` and `docs/lexicmap.md:41` ("everything before the final stop codon, when the reference ends in one") and `:178`.

### PL-06 · medium · Indel inside an HSP overlap is marked an edge artefact, so a frameshift can pass as clean

`src/disassembler/cmd/lexicmap.py:410` · correctness · confirmed (reproduced) · fix effort: trivial · sources: P1-05

**Problem.** `tile_accession_hsps` keeps only the copy of an indel carried by the HSP that owns the grid positions (higher `pident`, first on a tie) and judges "edge" from that HSP's span alone: `margin = min(pos - first, end - (pos + length))` (`lexicmap.py:408`), then `indels.append(Indel(kind, pos, length, bases, margin < INDEL_EDGE_WINDOW))` (`:410`). The copy in an overlapping HSP, where the indel may sit far from any end, is discarded. A real indel is then marked `(edge)`, ignored by `split_by_frame` and the stop splice, and an edge deletion does not count as covered (`:420`). The verdict depends on HSP input order.

**Impact.** A frameshifted gene can enter the clean cohort and MSA when its indel lies within 10 nt of one contig's end but inside the overlap with the next contig (adjacent Logan contigs overlap by k-1 = 30 nt). An in-frame 3-nt deletion there loses 3 positions of coverage instead, which matters under `--require-full-length`.

**Evidence.** A 120-nt CDS with a 1-nt deletion at 71, as one HSP: `FLAGGED_FRAMESHIFT 71del`. As two equal-`pident` HSPs 1-75 and 46-120 with 1-75 listed first: `0.9917 CLEAN_FULL_LENGTH 71del(edge)`, haplotype `WT_partial`. Swapping the order, or raising the second HSP's `pident`, gives `FLAGGED_FRAMESHIFT` again.

**Proposed fix.** An indel is an end artefact only if every HSP carrying it sees it near its own end:
```python
def margin_of(k, kind, pos, length):
    first, end = spans[k]
    return min(pos - first, end - pos) if kind == "ins" else min(pos - first, end - (pos + length))
inner = {(kind, pos, length, bases) for k, kind, pos, length, bases in raw
         if margin_of(k, kind, pos, length) >= INDEL_EDGE_WINDOW}
...
if owned:
    indels.append(Indel(kind, pos, length, bases, (kind, pos, length, bases) not in inner))
```
All four orderings then give `FLAGGED_FRAMESHIFT 71del`; the 55 lexicmap tests pass. An indel that two HSPs left-align differently (homopolymer) keeps today's behaviour.

**Tests and docs.** Test the two overlapping HSPs in both orders and with either `pident` higher (all `71del`, frameshift), and a lone indel 4 nt from its only HSP's end (stays edge). Update the `tile_accession_hsps` docstring and the module text on `INDEL_EDGE_WINDOW` ("of the end of every HSP that carries it"). Apply before PL-07's refill step. Related: PL-01.

### PL-07 · medium · Untrusted (edge) deletion's `-` displaces an overlapping HSP's real base

`src/disassembler/cmd/lexicmap.py:389` · correctness · confirmed (reproduced) · fix effort: small · sources: SWP1-01

**Problem.** A deletion's gap enters the per-position contest like a base: `if ex is None or pident > ex[1]: grid[q_pos] = (sc.upper(), pident, k)` (`lexicmap.py:388-389`) also runs for `sc == '-'`. When the deletion is later ruled an edge artefact, nothing removes the `-` from `grid`, so a lower-`pident` overlapping HSP that has the base there is ignored; `covered` (`:420`) excludes the position and the row (`:421`) shows `-`. This contradicts the comment at `:288-292` ("never counts as covered, and never splits a haplotype").

**Impact.** Two overlapping contigs where the higher-identity one closes with a forced 1-nt gap inside the overlap (the C. auris case the comment describes): coverage 59/60 instead of 60/60, a separate clean haplotype (wild-type count too low, entropy too high), and demotion to `CLEAN_PARTIAL` under `--require-full-length`.

**Evidence.** A 60-nt CDS; HSP A (`pident` 99, 1-60, deletion at 55) and HSP B (`pident` 98, 31-60, base present): `indels=['55del(edge)'] covered=59/60`, row `GAAC-ACTAA`; default run `0.9833 CLEAN_FULL_LENGTH`, 2 unique haplotypes, wild-type count 1; with `--require-full-length` `CLEAN_PARTIAL`.

**Proposed fix.** In `tile_accession_hsps`, remember the best non-gap base at every position a `-` holds (`best_base[p] = (base, pident, k)`), in both the gapped loop and the ungapped fast path (`:360-366`). After the ownership pass, for each owned deletion with `edge=True`, put `best_base[p]` back into `grid` and drop the indel when every position is refilled. Trusted deletions keep winning. A prototype (about 35 lines) fixes the case in both HSP orders and passes the 59 lexicmap tests. A lone HSP with an edge deletion still splits from WT until PL-02's key change; until then reword `:291-292`.

**Tests and docs.** `test_edge_deletion_does_not_hide_other_hsps_base` (A+B in both orders give the CDS, covered 60, no indel; A alone gives `55del(edge)`, covered 59) and `test_trusted_deletion_still_wins_over_lower_pident_base`. `docs/lexicmap.md:31-33` can say an edge deletion never hides a base another HSP supplied. Related: PL-01, PL-06.

### PL-08 · medium · Whitespace inside reference sequence lines becomes reference positions

`src/disassembler/cmd/lexicmap.py:221` · correctness · confirmed (reproduced) · fix effort: trivial · sources: PS-01, SWP1-06

**Problem.** `load_reference` does `line = line.strip()` (`lexicmap.py:212`) and then `seq_parts.append(line)` (`:221`), so a space or tab inside a sequence line becomes a reference position. The other readers drop it: `validators.fasta_record_lengths` (`len("".join(line.split()))`, `validators.py:94`) and Biopython in `msa`. One file therefore has two lengths in one `from-lexicmap` run. `load_reference` also skips text before the first `>` silently, while `msa` rejects it.

**Impact.** Without a `qlen` column (allowed, `:643`) the run is silently wrong at exit 0: a 12-nt CDS written `ATGGCT GAATAA` loads as 13 nt, `is_cds` is false so stop and frame checks are off, an identical sample gets coverage 0.9231 and spurious `NonCoding_Variant` changes, a premature TAG enters the clean cohort, and the reference row with the space goes into `clean.msa.fasta`. With `qlen` the run stops with "reference length disagreement ... Use the same reference FASTA the search used", which is wrong advice. A `#` line before `>` lets `from-lexicmap` write 12 files and then fail at stage 3.

**Evidence.** `printf '>geneA\nATGGCT GAATAA\n' > sp.fa`; `disassembler lexicmap -q sp.fa -f noq.tsv --out-dir B` gives exit 0, "geneA (13 bp, 0 codons)", `SRR1 0.9231 CLEAN_FULL_LENGTH`, haplotype `7 >G,8G>A,10A>T,11T>A NonCoding_Variant`; with a `qlen` column, exit 1 "13 vs 12".

**Proposed fix.**
```diff
+            if n_records == 0:
+                raise ValueError(f"{path}: text before the first '>' header ({line[:40]!r})")
             if n_records == 1:
-                seq_parts.append(line)
+                seq_parts.append("".join(line.split()))
```
With it the reference is 12 nt, `is_cds` is true, the stop is flagged and the `qlen` run passes. Later, one `textio.read_fasta()` shared by lexicmap, validators, sieve and msa would keep the rules from drifting.

**Tests and docs.** In the `load_reference` tests (`tests/test_lexicmap.py` near `:260`): internal space and tab give the stripped length; a leading `#` line raises. A stream test: spaced reference plus a stop-carrying hit gives `flagged_premature_stops == 1`. `docs/lexicmap.md:102`: "whitespace inside sequence lines is ignored; nothing may precede the first `>`".

### PL-09 · medium · Frameshift tier applied to references classified as non-coding

`src/disassembler/cmd/lexicmap.py:749` · correctness · confirmed (reproduced) · fix effort: trivial · sources: P1-03, P1-02

**Problem.** For a reference the tool itself reports as `"is_cds": false` (length not a multiple of 3, `:537`), the stop scan is switched off (`:731-732 ... if ref_is_cds else []`) and `NonCoding_*` classes exist (`:1099-1101`, `docs/lexicmap.md:70`), but the frame rule is not gated: `frame_ok = not frameshifts or allow_frameshifts` (`lexicmap.py:749`).

**Impact.** For rRNA, ITS, UTR or intergenic markers, any accession with a 1- or 2-nt indel lasting 30 or more bases is filed `FLAGGED_FRAMESHIFT` at exit 0; the indel allele never reaches the clean cohort, and `NonCoding_Variant` can never describe an indel.

**Evidence.** A 100-nt reference with a perfect hit and a hit carrying a 1-nt deletion at 41: `DEL1 1.0000 99.00 FLAGGED_FRAMESHIFT 5 0 None 1 41del`, `"is_cds": false`. After the fix: `CLEAN_FULL_LENGTH`, haplotype `41del NonCoding_Variant`.

**Proposed fix.**
```diff
-        frame_ok = not frameshifts or allow_frameshifts
+        frame_ok = (not ref_is_cds) or not frameshifts or allow_frameshifts
```
A truncated CDS whose length is not 3n loses frameshift flagging, as it already loses stop flagging; say so in the docs. The 3n non-coding case is covered by PL-04's `--non-coding`.

**Tests and docs.** Next to `test_one_nt_deletion_is_a_frameshift` (`tests/test_lexicmap.py:479`), a 100-nt reference with a 1-nt deletion expects `CLEAN_FULL_LENGTH`, `NonCoding_Variant` and `count_flagged_frameshift == 0`. `docs/lexicmap.md:48-51` and step 3: both rules apply only to a coding reference.

### PL-10 · medium · Raw FASTA record id used as a file name: `/` ids crash, `../` or absolute ids escape the output directory

`src/disassembler/cmd/lexicmap.py:578` · security · confirmed (reproduced) · fix effort: small · sources: PX-01, P1-09

**Problem.** `ref_name = header.split()[0] if header else path.stem` (`lexicmap.py:237`) is not sanitised and becomes a path in directory mode: `out_prefix = f"msa/{ref_name}"` (`:573`, no `-o`) and `out_prefix = str(out_path / ref_name)` (`:578`, for `--out-dir` via `:1314`, `-o DIR/` and `-o EXISTING_DIR`). A `/` splits it into path components, `../` climbs out, and an absolute id replaces the directory. `harvest.safe_target_name` (`harvest.py:56`) already fixed this class for kmindex target names.

**Impact.** Real ids such as `hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124|...` or `A/California/07/2009(H1N1)` exit 1 with ENOENT before any row is read (clear message; `-o PREFIX` works around it); ids over about 220 bytes fail with Errno 36. A crafted id (`../../../d/x`, `/abs/d/x`) writes all 12 outputs into an existing directory elsewhere and overwrites files with those suffixes; because the five cohort files are opened at `:646-650` before the table is checked, even a failing run truncates them. Escape needs a crafted reference header, and overwrites are limited to 12 fixed suffixes, hence medium. `from-lexicmap` (fixed prefix, `pipe_from_lexicmap.py:241`) and the Galaxy wrapper (`--out-prefix out/lexicmap`) are not affected.

**Evidence.** Reference id `../../../victim/lexicmap` with `--out-dir sandbox/out`: exit 0, 0 files in `sandbox/out`, 12 in `victim/`. Absolute id with `-o existing_dir`: all 12 outside. GISAID id: exit 1 `[Errno 2] No such file or directory: 'out/hCoV-19/Wuhan/...'`.

**Proposed fix.** Keep `ref_name` for the query filter, headers and `summary.json`; derive a file stem for `:573` and `:578` only, changing ids as little as `safe_target_name` does so names that work today (`lcl|NC_045512.2_cds_...`) are kept:
```python
def file_safe_name(name, limit=200):
    name = name.replace("/", "_").replace("\\", "_").replace("\0", "_").strip(". ")
    return name.encode()[:limit].decode(errors="ignore").rstrip(". ")

file_stem = file_safe_name(ref_name) or file_safe_name(Path(ref_path).stem) or "lexicmap"
```
Print a notice when `file_stem != ref_name`; put the helper in `validators.py` and share it with harvest. A `resolve()` containment check is unnecessary once separators are gone and breaks symlinked directories.

**Tests and docs.** Run `--out-dir` with `../../outside/x`, an absolute id, the GISAID and influenza ids, `lcl|...` and a 300-char id: 12 files inside `out/`, nothing outside, expected stem. `--help` at `lexicmap.py:1258` and `:1263`, and `docs/lexicmap.md:122`, `:181`: "path separators in the record id become `_`".

### PL-11 · medium · `clean.a2m` is a nucleotide alignment, but docs offer it to protein language models

`src/disassembler/cmd/lexicmap.py:56` · docs · confirmed (reproduced) · fix effort: trivial · sources: SC2-07

**Problem.** The A2M writer emits the nucleotide reference and tier 1 haplotype rows untranslated (`f.write(f">{ref_name}/1-{ref_len}\n{ref_seq}\n")`, `lexicmap.py:1027-1029`), yet the module docstring (`:56`) says "A2M formatted MSA for protein language models / EVcouplings", as do `docs/lexicmap.md:135` and the tools-iuc #8800 wrapper (`emit_a2m` help and tool help). The A2M rows also lack the `indels=` tag, so insertion haplotypes are indistinguishable there.

**Impact.** A, C, G and T are valid amino-acid letters, so ESM-MSA, MSA Transformer, EVE or EVcouplings in its default protein mode accept the file silently and return scores for a 4-letter "protein" three times the real length.

**Evidence.** `disassembler lexicmap ... -o out/v`; `head -2 out/v.clean.a2m` shows a 600-column row over the alphabet A, C, G, T.

**Proposed fix.** Now: describe the file truthfully in the docstring, `docs/lexicmap.md:135` and the wrapper:
```
<out-prefix>.clean.a2m   (nucleotide A2M of the reference + collapsed tier 1 haplotypes;
                          translate before using with protein MSA models / EVcouplings)
```
Later (medium): add `<prefix>.clean.protein.a2m`, translating each tier 1 haplotype with its in-frame indels applied on the codon grid, collapsing on protein identity with summed counts. Do not rename the nucleotide file without a deprecation step; Galaxy's `from_work_dir` depends on it.

**Tests and docs.** Docs only for step 1. For step 2: a synthetic cohort with one synonymous, one missense and one in-frame 3-nt deletion; assert row length `ref_len/3` and that the synonymous haplotype collapses into WT.

### PL-12 · medium · Ungrouped input aborts mid-stream leaving misclassified files; the advised `sort` breaks the header

`src/disassembler/cmd/lexicmap.py:710` · robustness · plausible · fix effort: small · sources: T1-06

**Problem.** `process_accession_block` deliberately raises when an accession reappears (`lexicmap.py:705-711`, documented at `docs/lexicmap.md:99`). The remedy it prints, `` `sort -k<sgenome_col> -s input.tsv` `` (`:710`), sorts the header in with the data: under `LC_ALL=C` the header goes last, the first data row is read as the header, and the run dies with an uncaught `KeyError` traceback. The five cohort files opened at `:646-650` are kept, and the split accession has already been written there, classified from its first block only.

**Impact.** Exit 1 with stray `cohort_ledger.tsv` and clean FASTAs but no `summary.json`; following the advice literally then fails with a traceback. Whether real single-gene LexicMap tables are ungrouped is not settled: unmerged commit 4ad9374 (branch `origin/fix/lexicmap-streamer-accession-grouping`) reports so for a 214k-accession table, but the query-column filter added later (a03a802, `:823-831`) removes the multi-query cause.

**Evidence.** A 3-row table with SRR001 split by SRR002: exit 1, ledger line `SRR001 0.5000 100.00 CLEAN_PARTIAL` left on disk (correctly sorted: `1.0000 CLEAN_FULL_LENGTH`). `LC_ALL=C sort -k3 -s ung.tsv` then gives `KeyError: "Missing required columns [...] in LexicMap header: ['g12', '12', 'SRR001', ...]"`.

**Proposed fix.**
1. Print a header-safe command naming the real column: `(head -n1 in.tsv; tail -n +2 in.tsv | sort -t$'\t' -k{c},{c} -s) > sorted.tsv` with `c = idx_sgenome + 1`. Followed verbatim, it gives exit 0 and the correct ledger.
2. Stage outputs and remove them on failure (PL-22).
3. Raise `ValueError` for missing columns (PL-21).
4. If ungrouped tables prove real: `--unsorted-input`, spooling to a temp file and externally sorting on `sgenome` (sorted chunks plus `heapq.merge`). Do not port 4ad9374 as is: it keeps every row's `qseq`/`sseq` in memory.

**Tests and docs.** Extend `test_ungrouped_input_is_rejected` to assert no cohort file remains; test that the advised command yields a valid table. Update `docs/lexicmap.md:99`; close the stale branch with a recorded decision.

### PL-13 · low · `--max-records` cuts through an accession's HSP block and mis-tiers it

`src/disassembler/cmd/lexicmap.py:890` · correctness · confirmed (reproduced) · fix effort: trivial · sources: P1-11

**Problem.** The row limit is tested after a row joins the open accession block (`if max_records and total_rows >= max_records: ... break`, `lexicmap.py:890-892`), and the partial block is then tiered (`:894-895`). The last accession of a preview is judged from a fraction of its HSPs; a clean multi-HSP accession becomes `LOW_COVERAGE_FRAGMENT` or `CLEAN_PARTIAL`. tools-iuc #8800 exposes this option as a Galaxy "preview" knob.

**Evidence.** SRR2 with HSPs 1-20, 21-40, 41-60: full run `1.0000 CLEAN_FULL_LENGTH hsps_merged 3`; `--max-records 2` gives `0.3333 LOW_COVERAGE_FRAGMENT hsps_merged 1`.

**Proposed fix.** Finish the open block, then stop before opening a new one:
```diff
+            if acc != curr_acc and limit_hit:
+                break
 ...
-            if max_records and total_rows >= max_records:
-                break
+            if max_records and total_rows >= max_records and not limit_hit:
+                print(f"[*] Reached row limit; finishing the HSP block of '{curr_acc}'.")
+                limit_hit = True
```
Do not count the row that triggers the `break` in `total_rows`.

**Tests and docs.** A 3-HSP accession straddling `max_records=2` stays `CLEAN_FULL_LENGTH`, `hsps_merged == 3`, and the next accession is absent. `docs/lexicmap.md:192`, `--help` and the wrapper: "the accession in progress is completed". PL-03's early-stop flag must cover this path.

### PL-14 · low · Rows with too few columns are dropped silently

`src/disassembler/cmd/lexicmap.py:821` · reliability · confirmed (reproduced) · fix effort: small · sources: P1-07, PS-06

**Problem.** `if len(parts) <= max(needed): continue` (`lexicmap.py:821-822`) skips a short row with no counter or message, unlike other-query rows (`:830`) and malformed rows (`:863`, `:917`); `total_rows` still counts it. A plain table cut inside its last row demotes the final accession (e.g. `CLEAN_FULL_LENGTH` to `CLEAN_PARTIAL`) at exit 0; a table with every row short exits 1 with "clean cohort is empty" instead of a schema error. `line.strip()` (`:812`) also removes empty trailing fields.

**Proposed fix.** Read `unterminated = not line.endswith("\n")` and use `line.rstrip("\r\n")`. A short unterminated row raises `ValueError` ("the input looks truncated"); other short rows are counted in `n_short` and warned about for the first `MALFORMED_REPORT_LIMIT`. Raise when `n_malformed + n_short == total_rows`. Add `"skipped_rows": {"too_few_columns", "unparseable", "other_query"}` to `summary.json`.

**Tests and docs.** A short middle row is counted and warned (exit 0); an unterminated short last row raises; all rows short raises the schema message. Document `skipped_rows` in the `docs/lexicmap.md` input section.

### PL-15 · low · Quantiles pick a rounded index; even-cohort medians are wrong and the regime call flips

`src/disassembler/cmd/lexicmap.py:178` · correctness · confirmed (reproduced) · fix effort: trivial · sources: P1-14

**Problem.** `compute_quantiles` uses `idx = int(round(p * (n - 1)))` (`lexicmap.py:178`) with banker's rounding, so an even cohort's median is one middle value, sometimes the lower and sometimes the upper. The p50 values drive `is_divergent = (med_pid < 85.0 or med_cov < 0.60)` (`:941`), which sets `regime` in `summary.json` and prints walker advice. Coverages 0.55 and 0.65 (true median 0.60) report median 0.55 and `divergent`; 0.50/0.55/0.65/0.70 report 0.65 and `conserved`. Advisory only; tiers are unaffected.

**Proposed fix.** Linear interpolation (numpy "linear", `statistics.quantiles(method="inclusive")`):
```python
pos = p * (n - 1); lo = int(pos); hi = min(lo + 1, n - 1)
return round(s[lo] + (s[hi] - s[lo]) * (pos - lo), 4)
```

**Tests and docs.** `TestComputeQuantiles` for n = 0, 1, 2, 4, 6, and a CLI test that the 0.55/0.65 cohort is `conserved`. `docs/lexicmap.md:75-76`: "median (linear interpolation)". Golden summaries in `bench/` may need regenerating; PL-24's count-based quantiles must use the same rule.

### PL-16 · low · Empty or untiled cohort is called "divergent" with invented 0.0 medians

`src/disassembler/cmd/lexicmap.py:941` · usability · confirmed (reproduced) · fix effort: small · sources: SWP1-07

**Problem.** `compute_quantiles` returns 0.0 for every quantile of an empty list (`lexicmap.py:173-174`), and `all_pids` only receives accessions with coverage above 0 (`:722-723`). The regime test (`:941`) fires on those placeholders: a header-only table (a search with no hits) is labelled `divergent` with "median identity is 0.0%" and walker advice, and with `--allow-empty` (as a fan-out over many targets would use) it exits 0 with that record. A cohort whose HSPs all fall below `--min-pident` is rightly divergent but described with a 0.0% median no HSP had.

**Proposed fix.** Return `None` quantiles for an empty list and print `n/a`. If `total_accessions == 0`, set `regime = "no_data"` with no warning; otherwise `is_divergent = med_pid is None or med_pid < 85.0 or med_cov < 0.60`, and when `med_pid is None` say "no HSP reached --min-pident {min_pident}%".

**Tests and docs.** Header-only table with `--allow-empty` gives `regime == "no_data"` and null median; a single 50% HSP gives the `--min-pident` wording. `docs/lexicmap.md:75-78` and `:139`; changelog note for consumers expecting numbers.

### PL-17 · low · Frameshift span counts deleted reference bases

`src/disassembler/cmd/lexicmap.py:452` · correctness · confirmed (reproduced) · fix effort: trivial · sources: SWP1-03

**Problem.** `split_by_frame` opens a cluster at `e.pos` (`lexicmap.py:452`) and closes it with `e.pos - start < FRAMESHIFT_MIN_SPAN` (`:456`), so a deletion's own bases count as out of frame although the frame shifts only after them. A full-length HSP with `31_64del,70_71del` (out of frame for 5 bases) is `FLAGGED_FRAMESHIFT`, while `31_34del,40_41del` (same 5-base stretch) is `CLEAN_FULL_LENGTH`. Rare: needs a long deletion with a nearby compensating indel.

**Proposed fix.** Count only bases read out of frame:
```python
opening = offset == 0
if opening:
    start, cluster, skipped = e.pos + (e.length if e.kind == "del" else 0), [], 0
cluster.append(e); offset = new
if offset == 0:
    (in_frame if e.pos - start - skipped < FRAMESHIFT_MIN_SPAN else shifted).extend(cluster)
    cluster = []
elif not opening and e.kind == "del":
    skipped += e.length
```

**Tests and docs.** Extend `test_split_by_frame`: del34@100 + del2@139 in frame; del34@100 + del2 at `134 + FRAMESHIFT_MIN_SPAN` shifted; ins@100 + del33@105 + del1@145 in frame. `docs/lexicmap.md:35-36`: "(deleted bases do not count)".

### PL-18 · low · `--min-cov-partial 0` files zero-coverage accessions as clean all-gap rows

`src/disassembler/cmd/lexicmap.py:785` · correctness · confirmed (reproduced) · fix effort: trivial · sources: SWP1-04

**Problem.** The tiers test only the coverage fraction (`:750`, `:772`, `elif cov < min_coverage_partial:` at `lexicmap.py:785`), and `validators.fraction` accepts 0.0. With `--min-cov-partial 0`, an accession whose HSPs all fail `--min-pident` (coverage 0, row all `-`) becomes `CLEAN_PARTIAL` and enters `clean_expanded.accessions.fasta`; with `--min-cov 0` as well it becomes a `WT_partial` haplotype of `------------` counted in frequencies and entropy. The tools-iuc #8800 wrapper allows `min="0"` on both gates.

**Proposed fix.**
```diff
-        if cov >= full_cov and stops_ok and frame_ok:
+        if tiling.covered and cov >= full_cov and stops_ok and frame_ok:
-        elif cov >= min_coverage_partial and stops_ok and frame_ok:
+        elif tiling.covered and cov >= min_coverage_partial and stops_ok and frame_ok:
-        elif cov < min_coverage_partial:
+        elif not tiling.covered or cov < min_coverage_partial:
```

**Tests and docs.** With both gates at 0, an accession with only sub-`--min-pident` HSPs lands in tier 4 and no clean row is all `-`. `docs/lexicmap.md:50`: "or nothing tiled".

### PL-19 · low · Under `--require-full-length`, `--min-cov` still drives a check, labels and the rescue count

`src/disassembler/cmd/lexicmap.py:1284` · usability · confirmed (reproduced) · fix effort: small · sources: SWP1-05

**Problem.** `--require-full-length` replaces `--min-cov` as the tier 1 gate (`:747`), but `if args.min_cov_partial > args.min_cov:` (`lexicmap.py:1284`, twin at `pipe_from_lexicmap.py:211`) refuses `--min-cov-partial 0.9` with a false reason (exit 2); the banner (`:549`) prints `full_length_cov=80.0%`; the tier 2 label (`:955`) reads "50%-80% cov"; and `tile_accession_hsps` gets `min_coverage` (`:716`), so `multi_hsp_rescued_accessions` is judged against 80%. Tiering itself is correct.

**Proposed fix.** Compute `full_cov = 1.0 if full_length_only else min_coverage` once and use it at `:716`, `:549` and `:955`; in both `run()` functions check `args.min_cov_partial > (1.0 if args.require_full_length else args.min_cov)` and name the active gate. Also replace the hard-coded "(>=50%)" at `:1173` with `{min_coverage_partial:.0%}`.

**Tests and docs.** `require_full_length=True, min_cov=0.8, min_cov_partial=0.9` runs and tiers a 0.95 accession as 2; a 90%+10% accession counts as rescued. `docs/lexicmap.md:186`; `--min-cov` help "ignored with --require-full-length". Changelog: the rescue metric changes.

### PL-20 · low · A row's `qseq` is never checked against the reference slice it is projected onto

`src/disassembler/cmd/lexicmap.py:868` · input-validation · confirmed (reproduced) · fix effort: small · sources: P1-12

**Problem.** `qseq = parts[idx_qseq]; sseq = parts[idx_sseq]` (`lexicmap.py:868-869`) are used as-is; nothing compares `qseq` with `ref_seq[qstart-1:qend]` or checks `len(qseq) == len(sseq)`. A reference with the search's query name and length but another sequence tiles every accession against the wrong sequence at exit 0 (a 60-nt test gave a clean Missense haplotype with 38 changes); short or hand-edited rows tile partially without warning. Ignoring `sstr` is correct (LexicMap reverse-complements the subject).

**Proposed fix.**
```python
qu = qseq.replace('-', '').upper()
if (len(qseq) != len(sseq) or qstart < 1 or qend > ref_len
        or qu != ref_seq[qstart - 1:qend]):
    n_inconsistent += 1   # warn for the first MALFORMED_REPORT_LIMIT, skip the row
    continue
```
Raise `ValueError` when more than half of the query's rows are inconsistent. Compare upper-cased; consider treating non-ACGT reference bases as wildcards.

**Tests and docs.** Fix `TestMutationClassification` (`tests/test_lexicmap.py:394-400`), which uses the mutant as `qseq`; add a wrong-reference test, a truncated-row test and a fixture from a real `lexicmap search -a` result with `sstr '-'` rows.

### PL-21 · low · Missing `qseq`/`sseq` and truncated or corrupt gzip end in tracebacks

`src/disassembler/cmd/lexicmap.py:636` · usability · confirmed (reproduced) · fix effort: small · sources: P1-08, D2-06, PA-08, PB-06, PC-10

**Problem.** `run()` promises one-line errors but catches only `(ValueError, OSError, RuntimeError)` (`lexicmap.py:1323`). The schema check raises `KeyError` (`:636`), so a table written without LexicMap's `-a/--all` (its default, and the likeliest user mistake) gives a 15-line traceback. `textio.text_stream` returns a bare `gzip.GzipFile` (`textio.py:31`), whose `EOFError` (truncated) and `zlib.error` (corrupt) escape every `except (OSError, ValueError)` in `lexicmap`, `msa`, `sieve`, `stratify`, `reconstruct` and `from-lexicmap`. Exit status is still 1.

**Proposed fix.** Fix at the source: raise `ValueError` at `:636` with "Re-run `lexicmap search` with -a/--all" when `qseq` or `sseq` is missing; wrap the gzip stream so read errors become `ValueError` naming the path:
```python
class _Gunzip(io.RawIOBase):
    def __init__(self, gz, name): self._gz, self._name = gz, name
    def readable(self): return True
    def readinto(self, b):
        try: return self._gz.readinto(b)
        except (EOFError, zlib.error) as e:
            raise ValueError(f"{self._name}: truncated or corrupt gzip data ({e})") from e
```

**Tests and docs.** A table without `qseq`/`sseq` gives exit 1, `--all` in stderr, no `Traceback`; half and byte-flipped gzips raise `ValueError` naming the path. Related: PL-03, PL-12.

### PL-22 · low · Outputs written in place; failed or `--no-a2m` reruns leave new and stale files mixed

`src/disassembler/cmd/lexicmap.py:646` · reliability · confirmed (reproduced) · fix effort: small · sources: P1-10

**Problem.** Five cohort files are opened under their final names before streaming (`lexicmap.py:646-650`); MSA, A2M, haplotype tables and `summary.json` are written only on success (`:1008-1169`). A rerun on the same prefix that fails (ungrouped input, `qlen` mismatch, no matching query) leaves new cohort files beside the previous run's `summary.json` and MSA; a successful `--no-a2m` rerun (`:1024`) leaves the old `.clean.a2m` (`hap_1_count_5`) next to a new MSA (`count=2`) at exit 0. Galaxy jobs are unaffected.

**Proposed fix.** Open every output as `path + ".part"`; on exception unlink all `.part` files and re-raise; on success unlink any of the 12 known output names this run did not write, then `os.replace` each `.part` into place. Replace `✓` in progress prints with ASCII, or call `sys.stdout.reconfigure(errors="replace")` in `cli.main`, so a console encoding cannot fail a finished run.

**Tests and docs.** Rerun on ungrouped input leaves no `.part` and a consistent set; a2m then `--no-a2m` leaves no `.a2m`. One line in `docs/lexicmap.md` outputs: "outputs appear together only on success". Related: PL-10, PL-12.

### PL-23 · low · Galaxy `deferred` and `paused` datasets exit 75 "try again later"

`src/disassembler/cmd/lexicmap.py:252` · reliability · confirmed (reproduced) · fix effort: trivial · sources: SWP1-08

**Problem.** `GALAXY_PENDING_STATES` (`lexicmap.py:252-253`) includes `paused` and `deferred`, so `check_dataset_ready` (`:281-283`) raises `DatasetNotReady` and `run` exits 75 ("Rerun once it is 'ok'"). Galaxy's `Dataset.non_ready_states` is `(NEW, UPLOAD, QUEUED, RUNNING, SETTING_METADATA)`; `DEFERRED` is terminal (it needs a separate materialize request) and `PAUSED` waits for a user to resume the job. A script retrying on 75, as `docs/lexicmap.md:89-91` invites, loops until its own limit.

**Proposed fix.**
```python
GALAXY_PENDING_STATES = {"new", "upload", "queued", "running", "setting_metadata"}
GALAXY_NEEDS_ACTION = {
    "deferred": "it is a deferred dataset; materialize it in Galaxy and pass the new dataset's id",
    "paused": "its job is paused; resume it in Galaxy and rerun once the dataset is 'ok'",
}
# check_dataset_ready: if state in GALAXY_NEEDS_ACTION: raise ValueError(...)   -> exit 1
```

**Tests and docs.** Change the tempfail loop at `tests/test_lexicmap.py:707` to the five pending states; `deferred` and `paused` return 1 without downloading. Update `docs/lexicmap.md:89-91`.

### PL-24 · low · Memory grows about 180 B per accession, contrary to the documented model

`src/disassembler/cmd/lexicmap.py:693` · efficiency · confirmed (reproduced) · fix effort: small · sources: PF-08

**Problem.** `seen_accessions` (`lexicmap.py:693`, contiguity check) and `all_covs`/`all_pids` (`:663-664`, kept only for five quantiles) grow by one entry per accession, while the docstring (`:37-41`) and `docs/lexicmap.md:66` say memory grows with distinct haplotypes only. Measured with one haplotype: 28.9 MB at 1,000 accessions, 383.9 MB at 2,000,000 (about 177 B each); about 0.9 GB at 5 M. A Galaxy allocation sized from the docs is too small.

**Proposed fix.** Keep a `Counter` of `tiling.covered` (at most `ref_len + 1` keys) and compute coverage quantiles from counts with the same rule as `compute_quantiles` (PL-15's rule if that lands first); store identities in `array.array('d')` (2 M accessions: 384 to 322 MB). Keep `seen_accessions` or store an 8-byte digest per accession. Correct the docs either way: "O(accessions) at roughly 100-180 bytes each, plus one sequence per distinct clean haplotype".

**Tests and docs.** Count-based quantiles equal `compute_quantiles` on random data. Update `lexicmap.py:37-41` and `docs/lexicmap.md:66`.

### PL-25 · low · Tiling builds a tuple per reference position; slice assignment is about 3.5x faster

`src/disassembler/cmd/lexicmap.py:343` · efficiency · confirmed (reproduced) · fix effort: small · sources: PF-04

**Problem.** `tile_accession_hsps` keeps `grid = [None] * ref_len` of `(base, pident, k)` tuples filled one character at a time even for ungapped HSPs (`lexicmap.py:343`, `:364`), then makes three more full passes (`:416-421`); `find_internal_stops` (`:164`) is a per-codon loop. On 10,000 accessions (60,103 rows, 3,822-nt reference) the run takes 13.1 s against 3.6 s with the replacement, outputs byte-identical; about 20 versus 6 minutes per million accessions.

**Proposed fix.** Three parallel lists `chars`, `pid`, `own`; write kept HSPs in ascending `(pident, -k)` order so a plain overwrite reproduces "higher pident wins, first wins a tie"; ungapped HSPs by slice (`chars[lo:hi] = sseq.upper()[lo-a:hi-a]`); gapped HSPs keep the per-character walk; re-sort `raw` by HSP index; `merged = "".join(chars)`. Stop scan via `seq.find(stop, i + 1)` per stop codon, keeping `i % 3 == 0`. Land after PL-01, PL-06 and PL-07, which change the same function.

**Tests and docs.** A differential test against the current implementation on randomized HSP sets (overlaps, ties, indels, lower case, out-of-range `qstart`).

### PL-26 · low · `stream_lexicmap_to_msa(file_path="-")` closes the process's `sys.stdin`

`src/disassembler/cmd/lexicmap.py:924` · reliability · confirmed (reproduced) · fix effort: trivial · sources: P1-13

**Problem.** For `-f -` the function wraps `sys.stdin.buffer` (`lexicmap.py:557`, `:564`) and skips closing it (`if raw is not sys.stdin.buffer: stream.close()`, `:924-925`), but never detaches the `TextIOWrapper`; when it is collected, `IOBase.__del__` closes `sys.stdin.buffer`. Library callers that read stdin afterwards get `ValueError: I/O operation on closed file`. The CLI is unaffected (fd 0 stays open).

**Proposed fix.**
```python
    finally:
        if raw is sys.stdin.buffer:
            stream.detach()      # leave the process's stdin open
        else:
            stream.close()
            raw.close()          # GzipFile does not close the file it wraps
```

**Tests and docs.** Patch `sys.stdin` with a `TextIOWrapper(BytesIO(tsv))`, call with `file_path="-"`, assert `not sys.stdin.buffer.closed`.

### PL-27 · low · Galaxy key lookup and HTTP access duplicated; tool-not-found messages and exit codes differ

`src/disassembler/cmd/lexicmap.py:128` · reuse · confirmed (reproduced) · fix effort: small · sources: PS-11

**Problem.** `get_galaxy_api_key` exists twice (`harvest.py:42`, `lexicmap.py:128`) and has already drifted once; `tests/test_lexicmap.py:739-743` pins both. Harvest uses a `requests` session with a 429/Retry-After loop; lexicmap makes single `urlopen` calls with no retry (`:245`, `:268`). The binary-not-found messages are four hand-written strings: `sieve.py:155` says "put it on PATH" without naming the file looked up (`l1_pipeline_serial`) or `$DISASSEMBLER_SIEVE_BIN`, so a user who installs a file named `kmer-diversity-sieve` still fails; the composites exit 1 where `sieve` and `walker` exit 127.

**Proposed fix.** `resolve.py`: `class ToolNotFound(RuntimeError)` and `not_found_message(tool)` naming file, env var, build command and resolution order, used by all four sites with exit 127. New `disassembler/galaxy.py` with one `api_key()`. Optionally one retrying `galaxy.session()`, retrying in one layer only.

**Tests and docs.** Each command's not-found message names `l1_pipeline_serial` and `DISASSEMBLER_SIEVE_BIN` and exits 127; one key-order test. Changelog: composites now exit 127.

### PL-28 · low · Test-only wrapper, a counter never read, triplicated flagged-row writes

`src/disassembler/cmd/lexicmap.py:426` · simplification · confirmed (reproduced) · fix effort: small · sources: PS-13

**Problem.** `merge_accession_hsps` (`lexicmap.py:426-431`) only wraps `tile_accession_hsps` and is called only by tests (`tests/test_lexicmap.py:53`, `:148`). `n_qlen_mismatch` (`:671`) is raised on at 1 and never read. `needed` (`:819-820`) is rebuilt per row. The three flagged tiers write near-identical lines (`:789-804`). `from-lexicmap run()` repeats lexicmap's argument checks (`pipe_from_lexicmap.py:206-215`, `:231-237`).

**Proposed fix.** Delete the wrapper and move the tests to `tile_accession_hsps`; raise directly on a `qlen` mismatch; hoist `max_needed` (keep the `is not None` filter); add a `flag(label, extra="")` helper writing both flagged files; factor the checks into `lexicmap.validate_args(args, prog)` called by both `run()` functions, so `from-lexicmap` still fails before creating `--out-dir`.

**Tests and docs.** `from-lexicmap` with neither `-d` nor `-f` returns 2 and creates no directory; a golden comparison of `.flagged.*` files.

### PL-29 · low · Haplotype table never reports the protein change of an in-frame indel

`src/disassembler/cmd/lexicmap.py:1068` · science · confirmed (reproduced) · fix effort: small · sources: SC2-10

**Problem.** Amino-acid changes are computed on the reference-grid row: `h_aa = translate_nt(seq)`, then `if h_a != '-' and h_a != r_a:` (`lexicmap.py:1068-1071`). In that row insertions are absent and deleted bases are `-`, and `translate_nt` returns `-` for any codon containing a gap (`:152-153`), so every codon touched by an indel is skipped and inserted residues never appear. For a haplotype whose only changes are indels the table prints `num_aa_muts` 0 and `aa_mutations` `NA` (`:1104-1105`, `:1117`). `find_internal_stops_with_indels` (`:462-506`) already splices trusted in-frame indels into the sample's own sequence, but only for the stop scan.

**Impact.** Every tier 1 haplotype with a trusted in-frame indel: a codon-aligned deletion or insertion, a 3-nt deletion off a codon boundary that fuses two codons, and a +1/-1 indel pair that restores the frame within `FRAMESHIFT_MIN_SPAN` (30 nt), which the code deliberately keeps clean. In the last case the residues inside the shifted stretch differ in the real protein but the grid shows them as reference, so up to 10 substitutions are hidden. `nt_mutations` and the `InFrameIndel` class still flag the haplotype and `NA` reads as "not computed", hence low; but anyone summing `num_aa_muts` undercounts, and `docs/lexicmap.md:67-68` promises amino-acid changes against the reference for each haplotype.

**Evidence.** Two hits on a test CDS, true proteins checked with Biopython: haplotype G (`300_301insA,313del`, protein `LTREFNET` to `LTKRVQQT`) and haplotype H (`302_304del`, `LTREFNET` to `LTKFNETA`). `disassembler lexicmap -q ref.fasta -f hits.tsv -o out/lm` exits 0 with both rows `num_aa_muts 0`, `aa_mutations NA`, class `InFrameIndel`.

**Proposed fix.**
1. Factor the splice loop of `find_internal_stops_with_indels` into `apply_indels(seq, ref_len, applied) -> (sample, starts, pieces)` and reuse it there.
2. In the haplotype loop, when `split_by_frame(indels)` gives applied indels and no shifted ones, translate the spliced sample, take the common protein prefix and suffix with the reference (an `X` from an uncovered codon matches anything), and report the middle window in HGVS form: `del`, `ins`, a single substitution, or `delins`. Keep grid-based substitutions outside that window, and run the prefix/suffix per indel cluster so two distant indels are not merged into one `delins`. A prototype gives `R101_E105delinsKRVQQ` and `R101_E102delinsK` for G and H, and `P11del` and `K10_P11insA` for a codon-aligned deletion and insertion.
3. Until then, write `indel:unresolved` and leave `num_aa_muts` empty rather than 0.

HGVS 3' shifting in repeats can give notation that differs from other tools, and `num_aa_muts` changes meaning for existing users (changelog note).

**Tests and docs.** In `tests/test_lexicmap.py` extend `test_in_frame_deletion_counts_as_covered_under_require_full_length` (`:485`, assert the codon deletion in `aa_mutations`), `test_in_frame_insertion_is_its_own_haplotype` (`:498`, an `ins` call) and `test_indels_that_restore_the_frame_quickly_are_not_a_frameshift` (`:560`, a `delins` over the shifted stretch); add the off-boundary 3-nt deletion. Document indel notation in `aa_mutations` in `docs/lexicmap.md` step 5. Related: PL-11 (a protein A2M needs the same splicing).
