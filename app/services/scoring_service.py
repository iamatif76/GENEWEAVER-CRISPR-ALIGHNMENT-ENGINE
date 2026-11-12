"""
Scoring and classification service for CRISPR off-target alignments.
Implements transparent scoring formulas and risk categorization.
"""

from typing import Tuple


def calculate_alignment_score(total_mismatches: int, seed_mismatches: int) -> float:
    """
    Assign an off-target risk score using the transparent formula:
        score = 100 - (total_mismatches * 10) - (seed_mismatches * 15)
    Clamped between 0.0 and 100.0.

    Higher scores indicate higher binding affinity / higher risk of off-target cleavage.
    An exact match (0 mismatches) scores 100.0.
    """
    raw_score = 100.0 - (total_mismatches * 10.0) - (seed_mismatches * 15.0)
    clamped_score = max(0.0, min(100.0, raw_score))
    return round(clamped_score, 2)


def classify_risk(mismatch_count: int) -> str:
    """
    Classify candidate alignment into risk categories based on mismatch count:
        - 0 mismatches: Exact match
        - 1 mismatch:   High-risk off-target
        - 2 mismatches: Moderate-risk off-target
        - 3 mismatches: Low-risk off-target
        - >3 mismatches: Minimal-risk off-target
    """
    if mismatch_count == 0:
        return "Exact match"
    elif mismatch_count == 1:
        return "High-risk off-target"
    elif mismatch_count == 2:
        return "Moderate-risk off-target"
    elif mismatch_count == 3:
        return "Low-risk off-target"
    else:
        return "Minimal-risk off-target"


def evaluate_alignment(
    total_mismatches: int,
    seed_mismatches: int
) -> Tuple[float, str]:
    """Helper returning (score, risk_category)."""
    score = calculate_alignment_score(total_mismatches, seed_mismatches)
    category = classify_risk(total_mismatches)
    return score, category
