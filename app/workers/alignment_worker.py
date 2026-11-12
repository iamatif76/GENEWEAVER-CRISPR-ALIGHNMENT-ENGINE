"""
Background worker executing CRISPR alignment jobs asynchronously.
Updates job progress, handles errors gracefully, and persists alignment results.
"""

from datetime import datetime
import json
import logging
from pathlib import Path
from sqlalchemy.orm import Session

from app.database import SessionLocal
from app.models import Job, Genome, Alignment
from app.core.genome_loader import load_fasta
from app.core.pam_detector import PAMDetector, CandidateSite
from app.core.aligner_gpu import GPUAligner
from app.config import settings

logger = logging.getLogger(__name__)


def run_alignment_job_worker(job_id: str) -> None:
    """
    Background worker task to execute alignment for a given job ID.
    Designed to be run in FastAPI BackgroundTasks or a separate process.
    """
    db: Session = SessionLocal()
    try:
        job = db.query(Job).filter(Job.id == job_id).first()
        if not job:
            logger.error("Worker could not find job %s", job_id)
            return

        job.status = "running"
        job.progress_percent = 5
        db.commit()

        genome = db.query(Genome).filter(Genome.id == job.genome_id).first()
        if not genome:
            raise ValueError(f"Associated reference genome (ID: {job.genome_id}) not found.")

        genome_path = Path(genome.file_path)
        if not genome_path.is_file():
            raise FileNotFoundError(f"Genome file on disk missing: {genome_path}")

        # 1. Load genome chromosomes
        job.progress_percent = 15
        db.commit()
        chromosomes = load_fasta(genome_path)
        total_chroms = len(chromosomes)
        if total_chroms == 0:
            raise ValueError("Reference genome contains no chromosomes.")

        # 2. Detect PAM candidate sites across all chromosomes
        detector = PAMDetector(
            pam_pattern=job.pam_pattern,
            guide_length=len(job.guide_sequence)
        )

        all_candidates = []
        for idx, (chrom, seq) in enumerate(chromosomes.items()):
            for site in detector.find_candidate_sites(chrom, seq):
                all_candidates.append(site)
            # Update progress between 20% and 55%
            pct = 20 + int(35 * ((idx + 1) / total_chroms))
            job.progress_percent = pct
            db.commit()

        # 3. Perform Alignment (using GPU if available or CPU fallback)
        job.progress_percent = 60
        db.commit()

        aligner = GPUAligner(
            guide_sequence=job.guide_sequence,
            max_mismatches=job.max_mismatches,
            seed_length=settings.DEFAULT_SEED_LENGTH,
            force_cpu=(settings.GPU_ENABLED == "cpu")
        )

        hits = aligner.align_batch(all_candidates)

        job.progress_percent = 85
        db.commit()

        # 4. Bulk insert alignment hits into SQLite
        # Delete any previous alignments for this job
        db.query(Alignment).filter(Alignment.job_id == job.id).delete()

        batch_records = []
        for hit in hits:
            record = Alignment(
                job_id=job.id,
                chromosome=hit.chromosome,
                start_position=hit.start_position,
                end_position=hit.end_position,
                strand=hit.strand,
                matched_sequence=hit.matched_sequence,
                pam_sequence=hit.pam_sequence,
                mismatch_count=hit.mismatch_count,
                mismatch_positions=json.dumps(hit.mismatch_positions),
                seed_mismatch_count=hit.seed_mismatch_count,
                score=hit.score,
                risk_category=hit.risk_category,
            )
            batch_records.append(record)
            if len(batch_records) >= 1000:
                db.bulk_save_objects(batch_records)
                db.commit()
                batch_records = []

        if batch_records:
            db.bulk_save_objects(batch_records)
            db.commit()

        # 5. Mark completed
        job.status = "completed"
        job.progress_percent = 100
        job.completed_at = datetime.utcnow()
        job.error_message = None
        db.commit()
        logger.info("Job %s completed successfully with %d alignment hits.", job_id, len(hits))

    except Exception as e:
        logger.exception("Error executing alignment job %s: %s", job_id, str(e))
        try:
            job = db.query(Job).filter(Job.id == job_id).first()
            if job:
                job.status = "failed"
                job.error_message = str(e)
                job.completed_at = datetime.utcnow()
                db.commit()
        except Exception as db_err:
            logger.error("Failed to update job status to failed: %s", str(db_err))
    finally:
        db.close()
