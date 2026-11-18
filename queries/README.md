# Reference queries

The reference CDSs the suite ships with. Every file here is a query: the `-q` of the walker and
`reconstruct`, or the `--ref` of `lexicmap` and `sieve`. Multi-strain panels for the walker's
multi-anchor mode live in [`../panels/`](../panels/README.md) instead.

Every record is a coding sequence in frame from its first base, so the codon MSA and the
lexicmap route's stop scan can be run against it without preparation. Anchoring is
exact-k-mer, so an ambiguous base is not a base that matches several things, it is a hole in
the anchor that matches nothing; prefer an unambiguous allele when adding one.

| path | records | what it is |
|---|---|---|
| `markers/*.fasta` | 1 each | 26 single-gene marker CDSs, one gene per file, with the gene id as the whole header |
| `master_planetary_targets.fasta` | 142 | a catalogue of 25 genes with their observed alleles, **not a panel** |
| `pfcrt_only.fasta` | 1 | *P. falciparum* PfCRT, the golden gate's `ERR2028820` reference |
| `query_hiv1_pol.fasta` | 1 | HIV-1 *pol*, the golden gate's `SRR8209082` reference |

## One gene per query

The walker's multi-anchor mode keeps one primary per accession, so a query holding two
different genes reports whichever scored higher and leaves the other as an unselected
candidate. The walker warns when it sees this, by grouping references that share at least 1% of
their canonical 21-mers, but it is a warning rather than a refusal, because a deliberately
divergent panel of one gene splits the same way. Use the marker files as they are, and cut one
group at a time out of the catalogue.

## `master_planetary_targets.fasta`

A catalogue of what a planetary-scale search is looking for, not something to run as one query:
its 142 records fall into 31 k-mer-sharing groups. Records follow a `<gene>_reference` plus
`<gene>_varN` convention, so a group is cut out by prefix:

```bash
grep -A1 '^>candida_auris_erg11_' queries/master_planetary_targets.fasta > erg11_panel.fasta
```

The 25 gene prefixes are `bla_ctx_m`, `bla_kpc`, `bla_ndm`, `bla_oxa48`, `bla_vim`,
`candida_auris_erg11`, `dengue_1_ns5` through `dengue_4_ns5`, `ebola_gp`, `gyr_a`, `h5n1_ha`,
`h5n1_pb2`, `measles_h`, `m_oryzae_avrpik`, `mpox_f13l`, `p_infestans_avr3a`, `plasmid_mcr1`,
`rabies_g`, `rsv_a_f`, `rsv_b_f`, `van_a`, `westnile_ns5` and `zika_ns5`. (`markers/mcr_3.fasta`
has no group here; `plasmid_mcr1` carries the mcr-3 to mcr-10 homologs.)

`tests/test_master_panel_homology.py` enforces the convention and checks that each `_varN` is
genuinely a variant of its `_reference`, at 0.25 translated identity. That bar exists because an
exact-matching unrelated gene beats a real allele under the walker's default `matched` metric,
so a group that is not one gene returns the wrong gene with `anchor=` naming it and no warning.
The gap it sits in is wide: the most divergent genuine group, `plasmid_mcr1`, scores 0.338, and
every other group is 0.57 or more.

Two records in the catalogue are known to be ambiguity-heavy: `ebola_gp_var2` is 42% ambiguous
and `zika_ns5_var3` is 24%. The walker subtracts ambiguous bases from the match count before
ranking, so they cannot win a multi-anchor contest they have not earned, but they anchor poorly.
