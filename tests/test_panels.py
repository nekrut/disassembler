#!/usr/bin/env python3
"""test_panels.py -- every reference panel shipped in panels/ must be a usable CDS panel.

The panels are data, not code, so nothing else checks them. `panels/sarscov2_spike_11strains.fa`
shipped in 0.1.0 with three defective records and none of the test suite noticed:

  - `Delta_B16172` was one base short of a codon boundary, because the record was copied from
    GenBank with its one IUPAC ambiguity code deleted rather than resolved. That shifted its
    frame from codon 95 onward, so a reconstruction anchored on it inherited a broken reading
    frame. It was the chosen anchor for 4,071 of the 9,448 primaries in the 20,000-accession
    walker run.
  - `Omicron_BA2` and `Omicron_XBB15` carried 268 and 477 Ns, from dropouts in the source
    genomes. A run of Ns anchors nothing (no k-mer over it can match) and counts against
    identity, which made XBB.1.5 a poor anchor precisely where it was picked most often.

Each check below is one of those failures. The frame checks apply to any CDS panel; the
ambiguity checks apply because an anchor is matched by exact k-mers.
"""
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
PANEL_DIR = REPO_ROOT / "panels"
STOPS = {"TAA", "TAG", "TGA"}


def read_fasta(path):
    """[(record id, description, sequence)], sequence uppercased."""
    recs = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            head = line[1:].split(None, 1)
            recs.append((head[0], head[1] if len(head) > 1 else "", []))
        elif line.strip() and recs:
            recs[-1][2].append(line.strip().upper())
    return [(i, d, "".join(s)) for i, d, s in recs]


def panels():
    return sorted(p for p in PANEL_DIR.glob("*.fa")) + sorted(PANEL_DIR.glob("*.fasta"))


class TestShippedPanels(unittest.TestCase):

    def test_there_is_at_least_one_panel_to_check(self):
        # Guards the glob: a rename of panels/ would otherwise make every test below vacuous.
        self.assertTrue(panels(), f"no panel FASTA found in {PANEL_DIR}")

    def test_every_record_is_a_complete_in_frame_cds(self):
        for path in panels():
            for rec_id, _, seq in read_fasta(path):
                with self.subTest(panel=path.name, record=rec_id):
                    self.assertTrue(seq, "empty record")
                    self.assertEqual(len(seq) % 3, 0,
                                     f"length {len(seq)} is not a multiple of 3")
                    self.assertEqual(seq[:3], "ATG", f"starts {seq[:3]}, not ATG")
                    self.assertIn(seq[-3:], STOPS, f"ends {seq[-3:]}, not a stop codon")

    def test_no_record_has_an_internal_stop(self):
        for path in panels():
            for rec_id, _, seq in read_fasta(path):
                with self.subTest(panel=path.name, record=rec_id):
                    codons = [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]
                    at = [i + 1 for i, c in enumerate(codons[:-1]) if c in STOPS]
                    self.assertEqual(at, [], f"internal stop codon(s) at codon {at[:5]}")

    def test_no_record_contains_n_or_an_ambiguity_code(self):
        # Anchoring is exact-k-mer, so an ambiguous base is a hole in the anchor rather than a
        # base that matches several things.
        for path in panels():
            for rec_id, _, seq in read_fasta(path):
                with self.subTest(panel=path.name, record=rec_id):
                    bad = sorted({c for c in seq} - set("ACGT"))
                    self.assertEqual(bad, [], f"non-ACGT character(s) {bad}")

    def test_record_ids_are_unique_within_a_panel(self):
        # Multi-anchor mode keys each candidate on the record id, and GraphAligner names a
        # result by the first whitespace token of the header, so two records sharing an id
        # cannot be told apart in the GAF.
        for path in panels():
            ids = [rec_id for rec_id, _, _ in read_fasta(path)]
            dupes = sorted({i for i in ids if ids.count(i) > 1})
            self.assertEqual(dupes, [], f"{path.name}: duplicate record id(s) {dupes}")

    def test_no_two_records_hold_the_same_sequence(self):
        # Two identical anchors are two candidates that can never be distinguished, so the
        # primary choice between them is arbitrary.
        for path in panels():
            seen = {}
            for rec_id, _, seq in read_fasta(path):
                with self.subTest(panel=path.name, record=rec_id):
                    self.assertNotIn(seq, seen,
                                     f"identical to {seen.get(seq)}")
                seen.setdefault(seq, rec_id)


if __name__ == "__main__":
    unittest.main()
