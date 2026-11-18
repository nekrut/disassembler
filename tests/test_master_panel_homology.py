#!/usr/bin/env python3
"""test_master_panel_homology.py -- every `<gene>_varN` in the master target file must
actually be a variant of its `<gene>_reference`.

The file follows the `<gene>_reference` + `<gene>_varN` convention documented in
`queries/README.md`, and a user cuts a group out of it as a multi-anchor panel
(`grep -A1 '^>candida_auris_erg11_'`). Two groups were not alleles at all:
`candida_auris_erg11_var1..9` were other CDSs of the same contig (an F-box protein, a
tRNA adenylyltransferase, an NADH kinase, Bud22p), and `p_infestans_avr3a_var1..9` were
PITG_14293-14328 while AVR3a is PITG_14371. Under the walker's default `matched` metric an
exact-matching unrelated gene beats a real allele that differs from the reference, so
exactly the samples such a panel exists to find (ERG11 Y132F or K143R, AVR3a-KI) came back
as the wrong gene, with `anchor=` naming it and no warning.

The threshold is 0.25 translated identity, and the gap it sits in is wide. Measured over
this file: the worst genuine group is `plasmid_mcr1` at 0.338, which is the mcr-3..10
homolog set and is meant to be divergent; every other group is 0.57 or more; the two bad
groups were 0.183 and 0.141. A nucleotide k-mer test would not do: divergent homologs
share no 21-mer either.
"""
import re
import unittest
from pathlib import Path

PANEL = Path(__file__).resolve().parent.parent / "queries" / "master_planetary_targets.fasta"
MIN_IDENTITY = 0.25


def _records():
    recs, rid, seq = [], None, []
    for line in PANEL.read_text().splitlines():
        if line.startswith(">"):
            if rid:
                recs.append((rid, "".join(seq)))
            rid, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if rid:
        recs.append((rid, "".join(seq)))
    return recs


class TestMasterPanelHomology(unittest.TestCase):

    def setUp(self):
        try:
            from Bio import Align  # noqa: F401
        except ImportError:
            self.skipTest("biopython not installed")

    def test_every_var_is_a_variant_of_its_reference(self):
        from Bio import Align
        from Bio.Align import substitution_matrices
        from Bio.Seq import Seq

        aligner = Align.PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score, aligner.extend_gap_score = -11, -1
        aligner.mode = "global"

        def protein(nt):
            nt = nt[: len(nt) // 3 * 3]
            return str(Seq(nt).translate()).rstrip("*")

        groups = {}
        for rid, seq in _records():
            m = re.match(r"^(.*)_(reference|var\d+)$", rid)
            if m:
                groups.setdefault(m.group(1), {})[m.group(2)] = seq

        bad = []
        for gene, group in sorted(groups.items()):
            ref = group.get("reference")
            if ref is None:
                continue
            a = protein(ref)
            for name, nt in sorted(group.items()):
                if name == "reference":
                    continue
                b = protein(nt)
                if not a or not b:
                    continue
                aln = aligner.align(a, b)[0]
                ident = sum(1 for x, y in zip(aln[0], aln[1]) if x == y and x != "-")
                frac = ident / min(len(a), len(b))
                if frac < MIN_IDENTITY:
                    bad.append(f"{gene}_{name}: {frac:.3f}")
        self.assertEqual(bad, [], "records that are not variants of their reference")

    def test_the_convention_itself_holds(self):
        """Every `_varN` has the `_reference` it claims, so a group cut out by prefix is
        a panel and not a fragment."""
        ids = [rid for rid, _ in _records()]
        refs = {rid[: -len("_reference")] for rid in ids if rid.endswith("_reference")}
        orphans = [rid for rid in ids
                   if re.match(r"^(.*)_var\d+$", rid)
                   and re.match(r"^(.*)_var\d+$", rid).group(1) not in refs]
        self.assertEqual(orphans, [])
