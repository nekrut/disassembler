# Why unitigs and not Logan contigs

Logan publishes both compacted unitigs (`s3://logan-pub/u/`) and assembled contigs
(`s3://logan-pub/c/`). If the walker only ever emitted the primary consensus, a fair question
is why stream the large unitig file at all rather than pull the gene out of the much smaller
contig file. This is the measurement that answers it.

**Short version:** on internally diverse samples the assembler shatters the gene across dozens
of short contigs, so it is not recoverable as a unit; the unitig walk reconstructs it
full-length anyway. On clonal samples the gene assembles into a single clean contig and the two
routes agree, with contigs far cheaper. The advantage of unitigs is therefore concentrated on
diverse samples, and it holds even in primary-only mode.

---

## Method

Nine *Zika* NS5 accessions (the curated `zika_ns5.top10.txt` cohort, minus `SRR6252011` which
the 350 MB size cap refuses at ~1.05 GB compressed) plus two clonal *C. auris* ERG11 controls.
For each:

```bash
# contigs, not unitigs
curl -s "https://logan-pub.s3.amazonaws.com/c/<ACC>/<ACC>.contigs.fa.zst" | zstd -dc > <ACC>.contigs.fa
makeblastdb -in <ACC>.contigs.fa -dbtype nucl -out <ACC>.db
blastn -query <gene>.fasta -db <ACC>.db \
  -outfmt '6 sseqid pident length qstart qend' -evalue 1e-6 -max_target_seqs 2000 > <ACC>.tsv
```

Then measure, from the BLAST hits at &ge; 90% identity: the longest single contiguous
alignment, the largest fraction of the gene any single contig covers, and how many contig
fragments it takes to tile the gene. Compare against the walker primary
(`logan-walker -a <ACC> -q <gene>.fasta --no-emit-alts`) on the same accession.

`blastn`, `makeblastdb`, `GraphAligner` and `zstd` are all in `environment.yml`.

---

## Result: Zika NS5 (deep, internally diverse)

The gene is **not assembled as a contig in any of the nine.** Max contig length in these sets
is 4.5&ndash;6 kb, so the contigs are long enough to hold a 2.7 kb gene &mdash; they just don't.

| accession | longest &ge;90% contig fragment | best single contig covers | walker primary |
|---|---|---|---|
| SRR6251963 | 636 nt | 7.0% | 2,708 nt, 99.8% |
| SRR6251966 | 636 nt | 8.3% | 2,708 nt, 99.8% |
| SRR6251967 | 468 nt | 4.7% | 2,708 nt, 99.8% |
| SRR6251969 | 88 nt | 3.2% | 2,708 nt, 99.8% |
| SRR6251973 | 89 nt | 3.3% | 2,708 nt, 99.8% |
| SRR6251979 | 101 nt | 3.7% | 2,708 nt, 99.8% |
| SRR6251984 | 88 nt | 3.2% | 2,708 nt, 99.4% |
| SRR6251992 | 94 nt | 3.4% | 2,708 nt, 99.8% |
| SRR6251993 | 88 nt | 3.2% | 2,708 nt, 99.8% |

No single contig holds more than 8% of NS5; for six of nine, the longest contiguous piece is
under 100 nt. The fragments do collectively tile the gene (union coverage 100%), but it takes
~50&ndash;63 of them &mdash; which is exactly the naive coordinate stitching that the parent
project measured at 87% chimeric. The walker reconstructed the full 2,708 nt as one connected
sequence for all nine.

## Control: C. auris ERG11 (clonal)

| accession | longest &ge;90% contig fragment | best single contig covers |
|---|---|---|
| ERR899743 | 1,575 nt (whole gene) | 100.0% (pid 99.9%) |
| SRR12363134 | 1,575 nt (whole gene) | 100.0% (pid 100.0%) |

The clonal gene comes out of the contigs as one clean full-length contig. So the fragmentation
above is not a Logan-contig artifact &mdash; contigs hold a full-length gene fine when the
sample is homogeneous at that locus.

---

## What causes it

Intra-host diversity. The Zika samples carry dense low-frequency variation (12&ndash;32 on-path
bubbles per accession), and the assembler breaks at nearly every bubble it cannot cleanly pop.
The unitig graph retains every branch, and `GraphAligner` threads one reference-guided path
straight through the bubble field. So the preservation of variation is what lets the walker
reconstruct the full-length **consensus** at all, on precisely the samples where the assembler
fails &mdash; even though primary-only mode then discards the individual variant calls.

## Costs and when contigs win

Contigs are 15&ndash;33&times; smaller than unitigs (2&ndash;7 MB vs hundreds of MB here). For
a clonal target that assembles cleanly, pulling the contig is cheaper and equivalent, and the
smaller file also dodges the size cap that dropped `SRR6252011`. The unitig walk earns its cost
on diverse samples and by keeping the intra-host variation recoverable
(`--emit-alts`, `--emit-vcf`); on clonal ones a contig pull would do.

See [`issue_intra_accession_variance_lost.md`](issue_intra_accession_variance_lost.md) for the
related finding on how the primary itself is chosen.
