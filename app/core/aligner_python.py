
"""
Pure Python baseline CRISPR alignment engine.

Uses ordinary Python loops for mismatch comparison.
This implementation serves as the reference baseline for
NumPy CPU and Numba CUDA performance comparisons.
"""

from typing import Dict, List, Optional

from app.core.aligner_cpu import AlignmentHit
from app.core.pam_detector import CandidateSite, PAMDetector
from app.services.scoring_service import evaluate_alignment


class PythonAligner:
    """Reference CRISPR alignment engine using pure Python comparisons."""

    def __init__(
        self,
        guide_sequence: str,
        max_mismatches: int = 3,
        seed_length: int = 10,
    ):
        self.guide = guide_sequence.strip().upper().replace("U", "T")
        self.guide_len = len(self.guide)
        self.max_mismatches = max_mismatches
        self.seed_length = seed_length
        self.seed_start_idx = max(0, self.guide_len - seed_length)

    def align_single(self, site: CandidateSite) -> Optional[AlignmentHit]:
        """Compare one candidate sequence using Python loops."""

        target = site.protospacer.upper().replace("U", "T")

        if len(target) != self.guide_len:
            return None

        mismatch_positions = []
        seed_mismatches = 0

        for i in range(self.guide_len):
            if self.guide[i] != target[i]:
                mismatch_positions.append(i + 1)

                if i >= self.seed_start_idx:
                    seed_mismatches += 1

        mismatch_count = len(mismatch_positions)

        if mismatch_count > self.max_mismatches:
            return None

        score, risk_category = evaluate_alignment(
            mismatch_count,
            seed_mismatches,
        )

        return AlignmentHit(
            chromosome=site.chromosome,
            start_position=site.start_position,
            end_position=site.end_position,
            strand=site.strand,
            matched_sequence=target,
            pam_sequence=site.pam,
            mismatch_count=mismatch_count,
            mismatch_positions=mismatch_positions,
            seed_mismatch_count=seed_mismatches,
            score=score,
            risk_category=risk_category,
        )

    def align_batch(
        self,
        sites: List[CandidateSite],
        batch_size: int = 50000,
    ) -> List[AlignmentHit]:
        """Align multiple candidates using pure Python."""

        if batch_size <= 0:
            raise ValueError("batch_size must be greater than zero.")

        hits = []

        for start in range(0, len(sites), batch_size):
            chunk = sites[start:start + batch_size]

            for site in chunk:
                hit = self.align_single(site)

                if hit is not None:
                    hits.append(hit)

        hits.sort(
            key=lambda h: (
                -h.score,
                h.mismatch_count,
                h.chromosome,
                h.start_position,
            )
        )

        return hits

    def align_genome(
        self,
        chromosomes: Dict[str, str],
        pam_pattern: str = "NGG",
        batch_size: int = 50000,
    ) -> List[AlignmentHit]:
        """Scan genome sequences and rank candidate alignments."""

        detector = PAMDetector(
            pam_pattern=pam_pattern,
            guide_length=self.guide_len,
        )

        candidates = []

        for chromosome, sequence in chromosomes.items():
            for site in detector.find_candidate_sites(
                chromosome,
                sequence,
            ):
                candidates.append(site)

        return self.align_batch(candidates, batch_size=batch_size)
