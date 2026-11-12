"""
Unit tests for CPU alignment, mismatch scoring, seed region weighting,
and risk classification.
"""

from app.core.aligner_cpu import CPUAligner, AlignmentHit
from app.core.pam_detector import CandidateSite
from app.services.scoring_service import calculate_alignment_score, classify_risk


def test_scoring_service_formulas():
    # Exact match: 100 - (0*10) - (0*15) = 100.0
    assert calculate_alignment_score(0, 0) == 100.0
    assert classify_risk(0) == "Exact match"

    # 1 mismatch outside seed: 100 - (1*10) - (0*15) = 90.0
    assert calculate_alignment_score(1, 0) == 90.0
    assert classify_risk(1) == "High-risk off-target"

    # 1 mismatch inside seed: 100 - (1*10) - (1*15) = 75.0
    assert calculate_alignment_score(1, 1) == 75.0

    # 2 mismatches (1 in seed, 1 outside): 100 - (2*10) - (1*15) = 65.0
    assert calculate_alignment_score(2, 1) == 65.0
    assert classify_risk(2) == "Moderate-risk off-target"

    # 3 mismatches:
    assert calculate_alignment_score(3, 2) == 40.0
    assert classify_risk(3) == "Low-risk off-target"

    # Clamp lower bound
    assert calculate_alignment_score(10, 5) == 0.0


def test_cpu_aligner_exact_match():
    guide = "GAGTCCGAGCAGAAGAAGAA"
    aligner = CPUAligner(guide_sequence=guide, max_mismatches=3)

    site = CandidateSite(
        chromosome="chr1",
        start_position=100,
        end_position=123,
        strand="+",
        protospacer=guide,
        pam="TGG",
    )

    hit = aligner.align_single(site)
    assert hit is not None
    assert hit.mismatch_count == 0
    assert hit.score == 100.0
    assert hit.risk_category == "Exact match"
    assert hit.seed_mismatch_count == 0


def test_cpu_aligner_filtering_above_threshold():
    guide = "GAGTCCGAGCAGAAGAAGAA"
    aligner = CPUAligner(guide_sequence=guide, max_mismatches=2)

    # 3 mismatches: should be excluded
    target_3m = "TTTTCCGAGCAGAAGAAGAA"
    site = CandidateSite(
        chromosome="chr1",
        start_position=100,
        end_position=123,
        strand="+",
        protospacer=target_3m,
        pam="CGG",
    )

    hit = aligner.align_single(site)
    assert hit is None


def test_cpu_aligner_batch_matches_single():
    guide = "GAGTCCGAGCAGAAGAAGAA"
    aligner = CPUAligner(guide_sequence=guide, max_mismatches=3)

    sites = [
        CandidateSite("chr1", 10, 33, "+", guide, "TGG"),
        CandidateSite("chr1", 50, 73, "+", "C" + guide[1:], "AGG"),
        CandidateSite("chr2", 100, 123, "-", guide[:18] + "CC", "GGG"),
    ]

    batch_hits = aligner.align_batch(sites)
    assert len(batch_hits) == 3

    # Ranking: exact match first (highest score)
    assert batch_hits[0].score == 100.0
    assert batch_hits[0].mismatch_count == 0

