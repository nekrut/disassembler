# Reference panels

Each panel is a FASTA of coding sequences for one gene across several strains, meant for the
walker's multi-anchor mode (`docs/walker.md`). Every record is a complete CDS: length divisible
by 3, an ATG start, one terminal stop codon, no internal stop, and no N or IUPAC ambiguity code.
`tests/test_panels.py` enforces all of that, so a new panel has to meet the same bar.

The ambiguity rule matters more than it looks. Anchoring is exact-k-mer, so an ambiguous base is
not a base that matches several things, it is a hole in the anchor that matches nothing.

## `sarscov2_spike_11strains.fa`

The spike (S) CDS of eleven SARS-CoV-2 strains. Record 1, `Wuhan_Hu1_ref`, is the coordinate
reference used by the case studies.

| record | lineage | source | nt |
|---|---|---|---:|
| `Wuhan_Hu1_ref` | A | NC_045512.2 | 3,822 |
| `Alpha_B117` | B.1.1.7 | PP846148.1 | 3,813 |
| `Beta_B1351` | B.1.351 | MW598419.1 | 3,813 |
| `Gamma_P1` | P.1 | MZ264787.1 | 3,822 |
| `Delta_B16172` | B.1.617.2 | PP515675.1 | 3,816 |
| `Omicron_BA1` | BA.1 | OL672836.1 | 3,813 |
| `Omicron_BA2` | BA.2 | PP846152.1 | 3,813 |
| `Omicron_BA5` | BA.5 | PV800153.1 | 3,813 |
| `Omicron_XBB15` | XBB.1.5 | PP846142.1 | 3,810 |
| `Omicron_EG51` | EG.5.1 | OR310579.1 | 3,822 |
| `Omicron_JN1` | JN.1 | PX943838.1 | 3,810 |

Every record is the `[gene=S]` CDS of its GenBank source, byte for byte:

```bash
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=<ACCESSION>&rettype=fasta_cds_na&retmode=text"
```

Take the record whose header carries `[gene=S]`, and do not alter it. Copying a sequence by hand,
or stripping characters it does not expect, is what broke the panel once already.

### Three records were replaced in 0.1.1

`Delta_B16172` was 3,821 nt, which is not a multiple of 3. Its GenBank source, OK091006.1, has a
clean 3,822 nt S CDS, but it carries one IUPAC `Y` at position 284 and that character was deleted
rather than resolved when the panel was built. The deletion shifted the reading frame from codon
95 onward, so a reconstruction anchored on the record inherited a broken frame. It was the chosen
anchor for 4,071 of the 9,448 primaries in the 20,000-accession walker run. The source record is
also a poor representative: it has six Ns where Delta's 157-158 deletion belongs, so it does not
carry that deletion at all.

`Omicron_BA2` and `Omicron_XBB15` were faithful copies of their sources, but those genomes have
dropouts inside spike: 268 and 477 Ns. XBB.1.5 was the shortest record in the panel and so won
the 0.1.0 `cov` metric often, which means the panel's most-used anchor had 12.5% of its length
unmatched.

The replacements are complete, unambiguous spike CDSs carrying the canonical marker set for their
lineage:

| record | new source | spike markers |
|---|---|---|
| `Delta_B16172` | PP515675.1 | T19R, T95I, G142D, E156G + 157-158del, L452R, T478K, D614G, P681R, D950N (plus G261V) |
| `Omicron_BA2` | PP846152.1 | T19I, 24-26del + A27S, G142D, V213G, G339D, S371F, S373P, S375F, T376A, D405N, R408S, K417N, N440K, S477N, T478K, E484A, Q493R, Q498R, N501Y, Y505H, D614G, H655Y, N679K, P681H, N764K, D796Y, Q954H, N969K |
| `Omicron_XBB15` | PP846142.1 | the BA.2 core plus V83A, Y144del + H146Q, Q183E, V213E, G252V, G339H, R346T, L368I, V445P, G446S, N460K, F486P, F490S |

The other eight records were not touched and are byte-identical to their 0.1.0 versions.

Replacing an anchor changes which reference wins for a given accession, so results that quote an
`anchor=` tag are not comparable across this change. Delta in particular went from 3,821 nt
without its defining deletion to 3,816 nt with it, which makes it fit real Delta samples far
better than it did.
