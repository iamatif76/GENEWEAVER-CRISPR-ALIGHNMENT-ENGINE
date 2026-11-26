"""
High-performance CPU alignment engine using vectorized NumPy operations.
Scans PAM candidate sites, computes mismatch vectors, evaluates seed-region mutations,
and calculates off-target risk scores.
"""

from dataclasses import dataclass
from typing import Dict, Generator, List, Optional
import numpy as np

from app.core.pam_detector import CandidateSite, PAMDetector
from app.services.scoring_service import evaluate_alignment


@dataclass
class AlignmentHit:
    chromosome: str
    start_position: int
    end_position: int
    strand: str
    matched_sequence: str
    pam_sequence: str
    mismatch_count: int
    mismatch_positions: List[int]
    seed_mismatch_count: int
    score: float
    risk_category: str

    def to_dict(self) -> dict:
        return {
            "chromosome": self.chromosome,
            "start_position": self.start_position,
            "end_position": self.end_position,
            "strand": self.strand,
            "matched_sequence": self.matched_sequence,
            "pam_sequence": self.pam_sequence,
            "mismatch_count": self.mismatch_count,
            "mismatch_positions": self.mismatch_positions,
            "seed_mismatch_count": self.seed_mismatch_count,
            "score": self.score,
            "risk_category": self.risk_category,
        }


class CPUAligner:
    """NumPy-vectorized CRISPR alignment engine."""

    def __init__(self, guide_sequence: str, max_mismatches: int = 3, seed_length: int = 10):
        self.guide = guide_sequence.strip().upper().replace("U", "T")
        self.guide_len = len(self.guide)
        self.max_mismatches = max_mismatches
        self.seed_length = seed_length
        self.seed_start_idx = max(0, self.guide_len - self.seed_length)

        # ASCII byte array representation for vectorized comparison
        self.guide_bytes = np.frombuffer(self.guide.encode("ascii"), dtype=np.uint8)

    def align_single(self, site: CandidateSite) -> Optional[AlignmentHit]:
        """Align a single candidate site against the guide sequence."""
        target = site.protospacer
        if len(target) != self.guide_len:
            return None

        target_bytes = np.frombuffer(target.encode("ascii"), dtype=np.uint8)
        diff_mask = target_bytes != self.guide_bytes
        mismatch_count = int(np.sum(diff_mask))

        if mismatch_count > self.max_mismatches:
            return None

        # 1-based positions
        mismatch_positions = (np.where(diff_mask)[0] + 1).tolist()
        seed_mismatches = int(np.sum(diff_mask[self.seed_start_idx:]))

        score, risk_cat = evaluate_alignment(mismatch_count, seed_mismatches)

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
            risk_category=risk_cat,
        )

    def align_batch(self, sites: List[CandidateSite], batch_size: int = 50000) -> List[AlignmentHit]:
        """
        Process candidate sites in high-throughput vectorized NumPy batches.
        """
        if batch_size <= 0:
            raise ValueError("batch_size must be greater than zero.")
        if not sites:
            return []

        hits: List[AlignmentHit] = []
        n_sites = len(sites)

        for i in range(0, n_sites, batch_size):
            chunk = sites[i : i + batch_size]
            chunk_len = len(chunk)

            # Filter candidates that strictly match guide length
            valid_chunk = [s for s in chunk if len(s.protospacer) == self.guide_len]
            if not valid_chunk:
                continue

            # Join all protospacers into a continuous byte array of shape (N, guide_len)
            seq_concat = "".join(s.protospacer for s in valid_chunk)
            byte_matrix = np.frombuffer(seq_concat.encode("ascii"), dtype=np.uint8).reshape(
                len(valid_chunk), self.guide_len
            )

            # Vectorized comparison against broadcasted guide array
            diffs = byte_matrix != self.guide_bytes
            mismatch_counts = np.sum(diffs, axis=1)

            # Filter candidates with mismatch count <= max_mismatches
            pass_mask = mismatch_counts <= self.max_mismatches
            passing_indices = np.where(pass_mask)[0]

            for idx in passing_indices:
                site = valid_chunk[idx]
                m_count = int(mismatch_counts[idx])
                site_diffs = diffs[idx]
                m_positions = (np.where(site_diffs)[0] + 1).tolist()
                seed_mismatches = int(np.sum(site_diffs[self.seed_start_idx:]))

                score, risk_cat = evaluate_alignment(m_count, seed_mismatches)

                hits.append(
                    AlignmentHit(
                        chromosome=site.chromosome,
                        start_position=site.start_position,
                        end_position=site.end_position,
                        strand=site.strand,
                        matched_sequence=site.protospacer,
                        pam_sequence=site.pam,
                        mismatch_count=m_count,
                        mismatch_positions=m_positions,
                        seed_mismatch_count=seed_mismatches,
                        score=score,
                        risk_category=risk_cat,
                    )
                )

        # Sort hits: highest score first, then lowest mismatch count
        hits.sort(key=lambda h: (-h.score, h.mismatch_count, h.chromosome, h.start_position))
        return hits

    def align_genome(
        self,
        chromosomes: Dict[str, str],
        pam_pattern: str = "NGG",
        batch_size: int = 50000
    ) -> List[AlignmentHit]:
        """Scan entire genome dictionary and return ranked alignment hits."""
        detector = PAMDetector(pam_pattern=pam_pattern, guide_length=self.guide_len)
        all_candidates: List[CandidateSite] = []
        for chrom, seq in chromosomes.items():
            for site in detector.find_candidate_sites(chrom, seq):
                all_candidates.append(site)

        return self.align_batch(all_candidates, batch_size=batch_size)

