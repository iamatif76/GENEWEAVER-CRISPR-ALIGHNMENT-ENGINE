"""
Alignment service coordinating job creation, progress tracking,
result querying, and summary calculations.
"""

from typing import List, Optional, Tuple
from sqlalchemy.orm import Session
from sqlalchemy import func

from app.models import Job, Alignment, Genome
from app.schemas import JobSummary, JobResponse, AlignmentResponse


class AlignmentService:
    """Coordinates CRISPR alignment job persistence and querying."""

    @staticmethod
    def create_job(
        db: Session,
        guide_sequence: str,
        genome_id: str,
        max_mismatches: int = 3,
        pam_pattern: str = "NGG"
    ) -> Job:
        """Create and queue a new alignment job."""
        job = Job(
            guide_sequence=guide_sequence,
            genome_id=genome_id,
            max_mismatches=max_mismatches,
            pam_pattern=pam_pattern,
            status="queued",
            progress_percent=0,
        )
        db.add(job)
        db.commit()
        db.refresh(job)
        return job

    @staticmethod
    def get_job(db: Session, job_id: str) -> Optional[Job]:
        """Fetch raw Job model by ID."""
        return db.query(Job).filter(Job.id == job_id).first()

    @staticmethod
    def get_job_summary(db: Session, job_id: str) -> JobSummary:
        """Compute summary statistics for a completed job."""
        total = db.query(func.count(Alignment.id)).filter(Alignment.job_id == job_id).scalar() or 0
        exact = db.query(func.count(Alignment.id)).filter(
            Alignment.job_id == job_id, Alignment.mismatch_count == 0
        ).scalar() or 0
        high = db.query(func.count(Alignment.id)).filter(
            Alignment.job_id == job_id, Alignment.mismatch_count == 1
        ).scalar() or 0
        moderate = db.query(func.count(Alignment.id)).filter(
            Alignment.job_id == job_id, Alignment.mismatch_count == 2
        ).scalar() or 0
        low = db.query(func.count(Alignment.id)).filter(
            Alignment.job_id == job_id, Alignment.mismatch_count >= 3
        ).scalar() or 0

        return JobSummary(
            total_alignments=total,
            exact_matches=exact,
            high_risk_count=high,
            moderate_risk_count=moderate,
            low_risk_count=low,
        )

    @staticmethod
    def get_job_response(db: Session, job_id: str) -> Optional[JobResponse]:
        """Retrieve Job along with computed summary statistics."""
        job = AlignmentService.get_job(db, job_id)
        if not job:
            return None

        summary = None
        if job.status == "completed":
            summary = AlignmentService.get_job_summary(db, job_id)

        return JobResponse(
            id=job.id,
            guide_sequence=job.guide_sequence,
            genome_id=job.genome_id,
            max_mismatches=job.max_mismatches,
            pam_pattern=job.pam_pattern,
            status=job.status,
            progress_percent=job.progress_percent,
            error_message=job.error_message,
            created_at=job.created_at,
            completed_at=job.completed_at,
            summary=summary,
        )

    @staticmethod
    def get_job_alignments(
        db: Session,
        job_id: str,
        limit: Optional[int] = None,
        offset: int = 0,
        chromosome: Optional[str] = None,
        strand: Optional[str] = None
    ) -> List[Alignment]:
        """Retrieve sorted alignments for a given job with optional filters."""
        query = db.query(Alignment).filter(Alignment.job_id == job_id)

        if chromosome:
            query = query.filter(Alignment.chromosome == chromosome)
        if strand:
            query = query.filter(Alignment.strand == strand)

        # Ranked by score DESC, then mismatch_count ASC
        query = query.order_by(
            Alignment.score.desc(),
            Alignment.mismatch_count.asc(),
            Alignment.chromosome.asc(),
            Alignment.start_position.asc()
        )

        if offset > 0:
            query = query.offset(offset)
        if limit is not None:
            query = query.limit(limit)

        return query.all()
