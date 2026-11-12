"""
Unit tests for sequence utilities, IUPAC parsing, and mismatch calculation.
"""

import pytest
from app.core.sequence_utils import (
    clean_and_validate_sequence,
    reverse_complement,
    iupac_to_regex,
    calculate_mismatches,
)


def test_clean_and_validate_sequence_standard():
    seq = "agtccgagcagaagaagaa"
    cleaned = clean_and_validate_sequence(seq, min_len=15, max_len=30)
    assert cleaned == "AGTCCGAGCAGAAGAAGAA"


def test_clean_and_validate_sequence_rna_conversion():
    rna = "AGUCCGAGCAGAAGAAGAA"
    cleaned = clean_and_validate_sequence(rna)
    assert cleaned == "AGTCCGAGCAGAAGAAGAA"
    assert "U" not in cleaned


def test_clean_and_validate_sequence_invalid_characters():
    with pytest.raises(ValueError, match="Invalid nucleotide"):
        clean_and_validate_sequence("AGTCCGAGCXGAAGAAGAA")


def test_clean_and_validate_sequence_length_bounds():
    # Too short
    with pytest.raises(ValueError, match="Sequence length"):
        clean_and_validate_sequence("ATCG", min_len=15)

    # Too long
    with pytest.raises(ValueError, match="Sequence length"):
        clean_and_validate_sequence("A" * 35, max_len=30)


def test_reverse_complement_basic():
    assert reverse_complement("ATCG") == "CGAT"
    assert reverse_complement("AAAA") == "TTTT"
    assert reverse_complement("CCGG") == "CCGG"


def test_reverse_complement_iupac():
    # R (A or G) complements to Y (T or C)
    rc = reverse_complement("ANGG")
    assert rc == "CCNT"


def test_iupac_to_regex():
    assert iupac_to_regex("NGG") == "[ACGT]GG"
    assert iupac_to_regex("NRG") == "[ACGT][AG]G"
    assert iupac_to_regex("TTTV") == "TTT[ACG]"

    with pytest.raises(ValueError, match="Unknown IUPAC"):
        iupac_to_regex("NGZ")


def test_calculate_mismatches_exact():
    guide = "GAGTCCGAGCAGAAGAAGAA"
    target = "GAGTCCGAGCAGAAGAAGAA"
    m_count, positions, seed_m = calculate_mismatches(guide, target, seed_length=10)
    assert m_count == 0
    assert positions == []
    assert seed_m == 0


def test_calculate_mismatches_outside_seed():
    guide = "GAGTCCGAGCAGAAGAAGAA"
    # Mutation at position 2 (index 1): A -> C
    target = "GCGTCCGAGCAGAAGAAGAA"
    m_count, positions, seed_m = calculate_mismatches(guide, target, seed_length=10)
    assert m_count == 1
    assert positions == [2]
    assert seed_m == 0  # Position 2 is outside the last 10 nt (positions 11-20)


def test_calculate_mismatches_inside_seed():
    guide = "GAGTCCGAGCAGAAGAAGAA"
    # Mutation at position 18 (index 17): G -> C
    target = "GAGTCCGAGCAGAAGAACAA"
    m_count, positions, seed_m = calculate_mismatches(guide, target, seed_length=10)
    assert m_count == 1
    assert positions == [18]
    assert seed_m == 1  # Position 18 is inside the seed (last 10 nt)
