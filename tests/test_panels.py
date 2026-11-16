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


# DD-01 and DD-02. Two records of the SARS-CoV-2 spike panel named a lineage they were not.
# `Omicron_EG51` was the spike of a 2020 B.1.2 genome, one base (D614G) from Wuhan, and
# `Omicron_BA5` was a BA.2 nanoluciferase reporter clone. The sequences were faithful copies
# of what their headers cited; the citations were the wrong genomes. Nothing in the suite
# looked at what a record *is*, only at whether it was a usable CDS, so a mislabelled anchor
# published its label into `anchor=`, the VCF `CHROM` and `snvs.tsv` `gene` with no way for a
# user to notice.
#
# The defining amino-acid changes below are the canonical, published sets for each lineage.
# They are deliberately not derived from the panel: a test that reads its expectation out of
# the data it is checking cannot catch the data being wrong.
LINEAGE_MARKERS = {
    "Alpha_B117": {"present": ["del69-70", "N501Y", "P681H", "D614G"], "absent": ["L452R"]},
    "Beta_B1351": {"present": ["K417N", "E484K", "N501Y", "D614G"], "absent": ["P681R"]},
    # K417T is deliberately NOT required. The shipped record carries 11 of P.1's 12
    # canonical spike markers, including the N-terminal set (L18F, T20N, P26S, D138Y,
    # R190S) and T1027I/V1176F, which together are highly specific to P.1, but it has no
    # change at 417 at all. K417T is absent in a minority of real P.1 genomes, so this is a
    # genuine Gamma spike missing one RBD marker rather than a mislabelled record. It is
    # worth knowing if the panel is ever used to call Gamma by RBD changes.
    "Gamma_P1": {"present": ["E484K", "N501Y", "D614G", "L18F", "T1027I", "V1176F"],
                 "absent": ["del69-70", "L452R"]},
    "Delta_B16172": {"present": ["L452R", "T478K", "P681R", "D614G"], "absent": ["N501Y"]},
    "Omicron_BA1": {"present": ["del69-70", "G339D", "S371L", "N501Y", "P681H"],
                    "absent": ["L452R"]},
    "Omicron_BA2": {"present": ["G339D", "S371F", "T376A", "Q493R", "N501Y"],
                    "absent": ["del69-70", "L452R", "F486V"]},
    "Omicron_BA5": {"present": ["del69-70", "L452R", "F486V", "G339D", "N501Y"],
                    "absent": ["Q493R"]},
    "Omicron_XBB15": {"present": ["G339H", "R346T", "F486P", "N501Y", "V445P"],
                      "absent": ["L452R", "F486V", "F456L"]},
    "Omicron_EG51": {"present": ["F456L", "Q52H", "F486P", "G339H", "N501Y"],
                     "absent": ["F486V"]},
    # JN.1 is BA.2.86 + L455S. F456L belongs to its descendants (KP.2 and relatives) and
    # to EG.5.1, not to JN.1 itself, so requiring it here would have been wrong.
    "Omicron_JN1": {"present": ["L455S", "N501Y", "F486P"], "absent": ["F486V", "F456L"]},
}


def _protein(nt: str) -> str:
    from Bio.Seq import Seq
    return str(Seq(nt).translate(table=1)).rstrip("*")


def protein_changes(wuhan_aa: str, query_aa: str) -> list:
    """Amino-acid changes of `query_aa` against `wuhan_aa`, in Wuhan numbering.

    Substitutions as `D614G`, a one-residue deletion as `Y144del`, a run as `del69-70`, an
    insertion as `ins214EPE`. The Wuhan position of each alignment column is computed once
    up front, because deriving it inside the branches double-counted after a deletion and
    shifted every marker past it by one or two.
    """
    from Bio.Align import PairwiseAligner
    al = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                         open_gap_score=-11, extend_gap_score=-1)
    a = al.align(wuhan_aa, query_aa)[0]
    ref_row, qry_row = str(a[0]), str(a[1])
    pos, n = [], 0
    for ch in ref_row:
        if ch != "-":
            n += 1
        pos.append(n)

    out, i = [], 0
    while i < len(ref_row):
        r, q = ref_row[i], qry_row[i]
        if r == q:
            i += 1
            continue
        if q == "-":
            j = i
            while j < len(ref_row) and qry_row[j] == "-":
                j += 1
            lo, hi = pos[i], pos[j - 1]
            out.append(f"{ref_row[i]}{lo}del" if lo == hi else f"del{lo}-{hi}")
            i = j
        elif r == "-":
            j = i
            while j < len(ref_row) and ref_row[j] == "-":
                j += 1
            out.append(f"ins{pos[i]}{qry_row[i:j]}")
            i = j
        else:
            out.append(f"{r}{pos[i]}{q}")
            i += 1
    return out


class TestSpikePanelLineages(unittest.TestCase):
    """What each spike panel record IS, not merely that it is a usable CDS (DD-01, DD-02)."""

    PANEL = PANEL_DIR / "sarscov2_spike_11strains.fa"

    @classmethod
    def setUpClass(cls):
        try:
            import Bio  # noqa: F401
        except ImportError:  # pragma: no cover
            raise unittest.SkipTest("biopython is needed to translate and align")
        if not cls.PANEL.is_file():
            raise unittest.SkipTest(f"{cls.PANEL} is not present")
        cls.recs = {i: s for i, _d, s in read_fasta(cls.PANEL)}
        cls.wuhan = _protein(cls.recs["Wuhan_Hu1_ref"])

    def test_every_record_carries_its_lineage_markers(self):
        for name, want in LINEAGE_MARKERS.items():
            with self.subTest(record=name):
                self.assertIn(name, self.recs, "the panel lost a record")
                got = set(protein_changes(self.wuhan, _protein(self.recs[name])))
                missing = [m for m in want["present"] if m not in got]
                self.assertEqual(missing, [], f"{name} is missing {missing}; it has {sorted(got)}")
                wrong = [m for m in want["absent"] if m in got]
                self.assertEqual(wrong, [], f"{name} carries {wrong}, which its lineage does not")

    def test_the_header_lineage_matches_the_table(self):
        # The header is the label a run publishes into `anchor=`, so it has to agree with
        # the markers the record actually carries.
        for ident, desc, _seq in read_fasta(self.PANEL):
            if ident == "Wuhan_Hu1_ref":
                continue
            with self.subTest(record=ident):
                self.assertIn(ident, LINEAGE_MARKERS,
                              "a record with no entry in LINEAGE_MARKERS is unchecked")
                self.assertIn("lineage=", desc)

    def test_no_two_records_are_near_duplicates(self):
        # DD-01's shipped `Omicron_EG51` differed from `Wuhan_Hu1_ref` by ONE base, so any
        # D614G sample matched it better than Wuhan and was published as a 2023 variant.
        # An anchor panel whose members are one base apart cannot distinguish them.
        K = 21
        kmers = {i: {s[j:j + K] for j in range(len(s) - K + 1)} for i, s in self.recs.items()}
        worst = None
        for a in sorted(kmers):
            for b in sorted(kmers):
                if a >= b:
                    continue
                private = len(kmers[a] - kmers[b])
                if worst is None or private < worst[0]:
                    worst = (private, a, b)
        self.assertGreaterEqual(
            worst[0], 3 * K,
            f"{worst[1]} and {worst[2]} share all but {worst[0]} of their {K}-mers, which is "
            f"too close for them to be separate anchors")

    def test_lengths_are_what_the_lineages_imply(self):
        # A deletion panel record that is Wuhan-length has not got its deletions.
        self.assertEqual(len(self.recs["Wuhan_Hu1_ref"]), 3822)
        for name, nt_len in (("Omicron_BA5", 3807), ("Omicron_EG51", 3810)):
            with self.subTest(record=name):
                self.assertEqual(len(self.recs[name]), nt_len)


if __name__ == "__main__":
    unittest.main()
