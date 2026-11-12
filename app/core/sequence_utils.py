"""
Sequence utilities for CRISPR guide RNA manipulation, reverse complementation,
IUPAC nucleotide pattern expansion, and mismatch calculation.
"""

from typing import List, Tuple

# IUPAC nucleotide code mappings
IUPAC_MAP = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "T",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]",
}

COMPLEMENT_MAP = str.maketrans("ACGTUacgtuNnRrYySsWwKkMmBbDdHhVv", "TGCAAtgcaaNnYyRrSsWwMmKkVvHhDdBb")


def clean_and_validate_sequence(seq: str, min_len: int = 15, max_len: int = 30) -> str:
    """
    Clean and validate a DNA or RNA sequence string.
    Converts RNA 'U' to 'T' internally, removes whitespace, and validates nucleotide alphabet.
    """
    if not seq:
        raise ValueError("Sequence cannot be empty.")

    cleaned = seq.strip().upper().replace("U", "T")
    valid_bases = set("ACGT")
    invalid_bases = set(cleaned) - valid_bases
    if invalid_bases:
        raise ValueError(f"Invalid nucleotide characters detected: {sorted(list(invalid_bases))}. Only A, C, G, T, U allowed.")

    if not (min_len <= len(cleaned) <= max_len):
        raise ValueError(
            f"Sequence length must be between {min_len} and {max_len} nucleotides. Got {len(cleaned)}."
        )

    return cleaned


def reverse_complement(seq: str) -> str:
    """
    Generate the reverse complement of a DNA/RNA sequence.
    Handles standard bases and IUPAC ambiguous bases.
    """
    cleaned = seq.strip().upper().replace("U", "T")
    return cleaned.translate(COMPLEMENT_MAP)[::-1]


def iupac_to_regex(pattern: str) -> str:
    """
    Convert an IUPAC nucleotide string (e.g. 'NGG', 'NRG', 'TTTN') into a regular expression.
    """
    upper_pattern = pattern.strip().upper().replace("U", "T")
    regex_parts = []
    for char in upper_pattern:
        if char in IUPAC_MAP:
            regex_parts.append(IUPAC_MAP[char])
        else:
            raise ValueError(f"Unknown IUPAC nucleotide code '{char}' in PAM pattern '{pattern}'")
    return "".join(regex_parts)


def calculate_mismatches(
    guide: str,
    target: str,
    seed_length: int = 10
) -> Tuple[int, List[int], int]:
    """
    Compare guide RNA with candidate genomic target sequence.

    Args:
        guide: 5' to 3' guide sequence (length L).
        target: 5' to 3' genomic protospacer (length L).
        seed_length: Number of nucleotides at the 3' end adjacent to the PAM
                     treated as the high-sensitivity seed region.

    Returns:
        (total_mismatches, mismatch_positions_1based, seed_mismatches)
    """
    if len(guide) != len(target):
        raise ValueError(f"Guide length ({len(guide)}) and target length ({len(target)}) must match.")

    length = len(guide)
    mismatch_positions: List[int] = []
    seed_start_idx = max(0, length - seed_length)
    seed_mismatches = 0

    for i in range(length):
        if guide[i] != target[i]:
            pos_1based = i + 1
            mismatch_positions.append(pos_1based)
            if i >= seed_start_idx:
                seed_mismatches += 1

    total_mismatches = len(mismatch_positions)
    return total_mismatches, mismatch_positions, seed_mismatches
