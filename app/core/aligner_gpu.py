"""
Optional GPU-accelerated CRISPR alignment engine using CuPy / CUDA.
Automatically falls back to CPUAligner when CUDA hardware or CuPy is unavailable.
"""

import logging
from typing import Dict, List, Optional
import numpy as np

from app.core.aligner_cpu import AlignmentHit, CPUAligner
from app.core.gpu_utils import is_gpu_ready
from app.core.pam_detector import CandidateSite, PAMDetector
from app.services.scoring_service import evaluate_alignment

logger = logging.getLogger(__name__)


class GPUAligner:
    """
    GPU-accelerated aligner using CuPy tensors with automatic CPU fallback.
    """

    def __init__(
        self,
        guide_sequence: str,
        max_mismatches: int = 3,
        seed_length: int = 10,
        force_cpu: bool = False
    ):
        self.guide = guide_sequence.strip().upper().replace("U", "T")
        self.guide_len = len(self.guide)
        self.max_mismatches = max_mismatches
        self.seed_length = seed_length
        self.seed_start_idx = max(0, self.guide_len - self.seed_length)
        self.force_cpu = force_cpu

        # Check GPU availability
        self.gpu_active = (not force_cpu) and is_gpu_ready()

        if self.gpu_active:
            try:
                import cupy as cp  # type: ignore
                self.cp = cp
                self.guide_gpu = self.cp.asarray(
                    np.frombuffer(self.guide.encode("ascii"), dtype=np.uint8)
                )
                logger.info("GPUAligner initialized on active CUDA device.")
            except Exception as e:
                logger.warning("CuPy initialization failed (%s). Falling back to CPU.", str(e))
                self.gpu_active = False
                self.cpu_fallback = CPUAligner(guide_sequence, max_mismatches, seed_length)
        else:
            self.cpu_fallback = CPUAligner(guide_sequence, max_mismatches, seed_length)

    def align_batch(self, sites: List[CandidateSite], batch_size: int = 100000) -> List[AlignmentHit]:
        """Align a batch of candidate sites using GPU if available, else CPU fallback."""
        if not self.gpu_active:
            return self.cpu_fallback.align_batch(sites, batch_size=batch_size)

        if not sites:
            return []

        hits: List[AlignmentHit] = []
        n_sites = len(sites)

        try:
            for i in range(0, n_sites, batch_size):
                chunk = sites[i : i + batch_size]
                valid_chunk = [s for s in chunk if len(s.protospacer) == self.guide_len]
                if not valid_chunk:
                    continue

                seq_concat = "".join(s.protospacer for s in valid_chunk)
                host_matrix = np.frombuffer(seq_concat.encode("ascii"), dtype=np.uint8).reshape(
                    len(valid_chunk), self.guide_len
                )

                # Transfer to device
                dev_matrix = self.cp.asarray(host_matrix)
                dev_diffs = dev_matrix != self.guide_gpu
                dev_mismatches = self.cp.sum(dev_diffs, axis=1)

                pass_mask = dev_mismatches <= self.max_mismatches
                dev_passing_indices = self.cp.where(pass_mask)[0]

                # Bring passing hits back to host
                passing_indices = self.cp.asnumpy(dev_passing_indices)
                host_diffs = self.cp.asnumpy(dev_diffs[passing_indices])
                host_mismatches = self.cp.asnumpy(dev_mismatches[passing_indices])

                for local_idx, orig_idx in enumerate(passing_indices):
                    site = valid_chunk[orig_idx]
                    m_count = int(host_mismatches[local_idx])
                    diff_row = host_diffs[local_idx]
                    m_positions = (np.where(diff_row)[0] + 1).tolist()
                    seed_mismatches = int(np.sum(diff_row[self.seed_start_idx:]))

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

            hits.sort(key=lambda h: (-h.score, h.mismatch_count, h.chromosome, h.start_position))
            return hits

        except Exception as e:
            logger.error("GPU alignment encountered an error (%s); falling back to CPU.", str(e))
            return self.cpu_fallback.align_batch(sites, batch_size=batch_size)

    def align_genome(
        self,
        chromosomes: Dict[str, str],
        pam_pattern: str = "NGG",
        batch_size: int = 100000
    ) -> List[AlignmentHit]:
        """Align entire genome using GPU or fallback."""
        detector = PAMDetector(pam_pattern=pam_pattern, guide_length=self.guide_len)
        all_candidates: List[CandidateSite] = []
        for chrom, seq in chromosomes.items():
            for site in detector.find_candidate_sites(chrom, seq):
                all_candidates.append(site)

        return self.align_batch(all_candidates, batch_size=batch_size)
