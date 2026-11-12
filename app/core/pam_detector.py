"""
PAM (Protospacer Adjacent Motif) detection engine.
Identifies candidate target sites on both forward (+) and reverse (-) strands.
Supports IUPAC degenerate nucleotide patterns (e.g., NGG for SpCas9, TTTV for Cas12a/Cpf1).
"""

import re
from dataclasses import dataclass
from typing import Generator, List
from app.core.sequence_utils import iupac_to_regex, reverse_complement


@dataclass
class CandidateSite:
    chromosome: str
    start_position: int  # 1-based genomic start
    end_position: int    # 1-based genomic end
    strand: str          # "+" or "-"
    protospacer: str     # 5' to 3' candidate sequence on targeted strand
    pam: str             # 5' to 3' PAM on targeted strand


class PAMDetector:
    """Detects PAM sites and extracts candidate protospacers across genome contigs."""

    def __init__(self, pam_pattern: str = "NGG", guide_length: int = 20):
        self.pam_pattern = pam_pattern.strip().upper().replace("U", "T")
        self.guide_length = guide_length
        self.pam_length = len(self.pam_pattern)

        # Forward strand regex (e.g., NGG -> [ACGT]GG)
        self.fwd_pam_regex_str = iupac_to_regex(self.pam_pattern)
        # Reverse complement PAM regex (e.g., CCN for minus strand target)
        self.rev_pam_pattern = reverse_complement(self.pam_pattern)
        self.rev_pam_regex_str = iupac_to_regex(self.rev_pam_pattern)

        # Lookahead regexes to capture overlapping PAM motifs
        self._fwd_pattern = re.compile(f"(?=({self.fwd_pam_regex_str}))")
        self._rev_pattern = re.compile(f"(?=({self.rev_pam_regex_str}))")

    def find_candidate_sites(
        self,
        chromosome: str,
        sequence: str
    ) -> Generator[CandidateSite, None, None]:
        """
        Scan a single chromosome sequence and yield all candidate sites with valid PAMs
        on both forward (+) and reverse (-) strands.
        """
        seq_len = len(sequence)
        if seq_len < (self.guide_length + self.pam_length):
            return

        # 1. Forward strand (+) matches: [Protospacer: L][PAM: P]
        for match in self._fwd_pattern.finditer(sequence):
            pam_start = match.start()
            if pam_start >= self.guide_length:
                protospacer_start = pam_start - self.guide_length
                protospacer = sequence[protospacer_start:pam_start]
                pam = sequence[pam_start:pam_start + self.pam_length]

                # Ensure candidate contains only standard bases without uncalled gap/N characters
                if "N" not in protospacer:
                    yield CandidateSite(
                        chromosome=chromosome,
                        start_position=protospacer_start + 1,  # 1-based
                        end_position=pam_start + self.pam_length,
                        strand="+",
                        protospacer=protospacer,
                        pam=pam
                    )

        # 2. Reverse strand (-) matches: on (+) reference strand reads [PAM_rc: P][Protospacer_rc: L]
        for match in self._rev_pattern.finditer(sequence):
            pam_start = match.start()
            protospacer_start = pam_start + self.pam_length
            protospacer_end = protospacer_start + self.guide_length

            if protospacer_end <= seq_len:
                protospacer_rc = sequence[protospacer_start:protospacer_end]
                pam_rc = sequence[pam_start:pam_start + self.pam_length]

                # Convert to 5' -> 3' of the targeted (-) strand
                protospacer = reverse_complement(protospacer_rc)
                pam = reverse_complement(pam_rc)

                if "N" not in protospacer:
                    yield CandidateSite(
                        chromosome=chromosome,
                        start_position=pam_start + 1,  # 1-based
                        end_position=protospacer_end,
                        strand="-",
                        protospacer=protospacer,
                        pam=pam
                    )


def find_all_candidates_in_genome(
    chromosomes: dict,
    pam_pattern: str = "NGG",
    guide_length: int = 20
) -> List[CandidateSite]:
    """Helper to collect all candidate sites across all chromosomes."""
    detector = PAMDetector(pam_pattern=pam_pattern, guide_length=guide_length)
    results = []
    for chrom, seq in chromosomes.items():
        for site in detector.find_candidate_sites(chrom, seq):
            results.append(site)
    return results
