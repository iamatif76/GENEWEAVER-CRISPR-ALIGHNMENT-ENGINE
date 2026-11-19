"""
Tests for benchmark ground-truth validation.
"""

from types import SimpleNamespace

from benchmark.benchmark import EXPECTED_TARGETS, validate_ground_truth


def make_hit(signature):
    """Create a lightweight alignment-hit object from a target signature."""
    (
        chromosome,
        start_position,
        end_position,
        strand,
        matched_sequence,
        pam_sequence,
        mismatch_count,
        mismatch_positions,
    ) = signature

    return SimpleNamespace(
        chromosome=chromosome,
        start_position=start_position,
        end_position=end_position,
        strand=strand,
        matched_sequence=matched_sequence,
        pam_sequence=pam_sequence,
        mismatch_count=mismatch_count,
        mismatch_positions=mismatch_positions,
    )


def test_all_ground_truth_targets_match():
    hits = [make_hit(target) for target in EXPECTED_TARGETS]

    result = validate_ground_truth(hits)

    assert result["passed"] == 8
    assert result["total"] == 8
    assert result["accuracy_percent"] == 100.0
    assert result["missing_targets"] == []
    assert result["unexpected_hits"] == []


def test_missing_target_reduces_accuracy():
    targets = list(EXPECTED_TARGETS)
    hits = [make_hit(target) for target in targets[:-1]]

    result = validate_ground_truth(hits)

    assert result["passed"] == 7
    assert result["total"] == 8
    assert result["accuracy_percent"] == 87.5
    assert len(result["missing_targets"]) == 1
    assert result["unexpected_hits"] == []


def test_unexpected_hit_is_detected():
    hits = [make_hit(target) for target in EXPECTED_TARGETS]

    unexpected_target = (
        "chrZ",
        200,
        222,
        "+",
        "GAGTCCGAGCAGAAGAAGAA",
        "TGG",
        0,
        (),
    )
    hits.append(make_hit(unexpected_target))

    result = validate_ground_truth(hits)

    assert result["passed"] == 8
    assert result["accuracy_percent"] == 100.0
    assert result["missing_targets"] == []
    assert len(result["unexpected_hits"]) == 1
    assert result["unexpected_hits"][0] == unexpected_target


def test_empty_hits_return_zero_accuracy():
    result = validate_ground_truth([])

    assert result["passed"] == 0
    assert result["total"] == 8
    assert result["accuracy_percent"] == 0.0
    assert len(result["missing_targets"]) == 8
    assert result["unexpected_hits"] == []