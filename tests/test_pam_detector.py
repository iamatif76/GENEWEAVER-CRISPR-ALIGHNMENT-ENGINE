"""
Unit tests for PAM detection on forward and reverse strands,
overlapping motif handling, and degenerate IUPAC codes.
"""

from app.core.pam_detector import PAMDetector, CandidateSite
from app.core.sequence_utils import reverse_complement


def test_pam_detector_forward_strand():
    # Protospacer (20 nt) + PAM (AGG)
    protospacer = "ACGTACGTACGTACGTACGT"
    pam = "AGG"
    seq = "AAAAA" + protospacer + pam + "TTTTT"

    detector = PAMDetector(pam_pattern="NGG", guide_length=20)
    sites = list(detector.find_candidate_sites("chr1", seq))

    # Verify forward strand site found
    fwd_sites = [s for s in sites if s.strand == "+"]
    assert len(fwd_sites) >= 1
    hit = fwd_sites[0]
    assert hit.protospacer == protospacer
    assert hit.pam == pam
    assert hit.start_position == 6  # 1-based start after 5 'A's
    assert hit.end_position == 28


def test_pam_detector_reverse_strand():
    protospacer = "ACGTACGTACGTACGTACGT"
    pam = "TGG"
    # Target on minus strand: 5' to 3' is protospacer + pam
    # On plus reference strand, this appears as reverse_complement(protospacer + pam)
    target_locus_rc = reverse_complement(protospacer + pam)

    seq = "GGGGG" + target_locus_rc + "CCCCC"
    detector = PAMDetector(pam_pattern="NGG", guide_length=20)
    sites = list(detector.find_candidate_sites("chr1", seq))

    rev_sites = [s for s in sites if s.strand == "-"]
    assert len(rev_sites) >= 1
    hit = rev_sites[0]
    # Candidate protospacer on the minus strand should match our original protospacer
    assert hit.protospacer == protospacer
    assert hit.pam == pam


def test_pam_detector_overlapping_pams():
    # Sequence with consecutive Gs: TGGGG has overlapping NGG motifs:
    # index 0: TGG, index 1: GGG, index 2: GGG
    seq = ("A" * 20) + "TGGGG"
    detector = PAMDetector(pam_pattern="NGG", guide_length=20)
    sites = list(detector.find_candidate_sites("chrX", seq))

    fwd_sites = [s for s in sites if s.strand == "+"]
    assert len(fwd_sites) == 3


def test_pam_detector_degenerate_iupac():
    # Test with NAG PAM (common alternative SpCas9 PAM)
    protospacer = "TGCA" * 5
    pam = "TAG"
    seq = protospacer + pam
    detector = PAMDetector(pam_pattern="NAG", guide_length=20)
    sites = list(detector.find_candidate_sites("chr1", seq))

    fwd_sites = [s for s in sites if s.strand == "+"]
    assert len(fwd_sites) == 1
    assert fwd_sites[0].pam == "TAG"

