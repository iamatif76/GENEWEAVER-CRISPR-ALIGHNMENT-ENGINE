"""
Background worker executing CRISPR alignment jobs asynchronously.

Processes FASTA genomes in overlapping chunks, detects PAM sites,
aligns candidates, and persists results incrementally.
"""
import json
import logging
from datetime import datetime
from pathlib import Path

from sqlalchemy.orm import Session

from app.database import SessionLocal
from app.models import Job, Genome, Alignment
from app.core.genome_loader import inspect_fasta, iter_fasta_chunks
from app.core.pam_detector import PAMDetector, CandidateSite
from app.core.aligner_gpu import GPUAligner
from app.config import settings


logger = logging.getLogger(__name__)

# Processing limits
FASTA_CHUNK_SIZE = 1_000_000
ALIGNMENT_BATCH_SIZE = 50_000
DB_INSERT_BATCH_SIZE = 1_000


def run_alignment_job_worker(job_id: str) -> None:
    db: Session = SessionLocal()

    try:
        # 1. Find the requested job.
        job = db.query(Job).filter(Job.id == job_id).first()

        if not job:
            logger.error("Worker could not find job %s", job_id)
            return

        job.status = "running"
        job.progress_percent = 5
        job.error_message = None
        db.commit()

        # 2. Find and validate the genome file.
        genome = db.query(Genome).filter(
            Genome.id == job.genome_id
        ).first()

        if not genome:
            raise ValueError(
                f"Genome not found for job {job_id}"
            )

        genome_path = Path(genome.file_path)

        if not genome_path.is_file():
            raise FileNotFoundError(
                f"Genome file does not exist: {genome_path}"
            )

        # 3. Inspect the FASTA file without loading all sequences.
        job.progress_percent = 10
        db.commit()

        total_bases, chromosome_count, chrom_info = inspect_fasta(
            genome_path
        )

        if total_bases <= 0 or chromosome_count == 0:
            raise ValueError("Genome contains no valid sequences.")

        chromosome_lengths = dict(chrom_info)

        # 4. Configure PAM detection.
        detector = PAMDetector(
            pam_pattern=job.pam_pattern,
            guide_length=len(job.guide_sequence)
        )

        overlap = (
            len(job.guide_sequence)
            + len(job.pam_pattern)
            - 1
        )

        # 5. Initialize the existing GPU/CPU aligner.
        aligner = GPUAligner(
            guide_sequence=job.guide_sequence,
            max_mismatches=job.max_mismatches,
            seed_length=settings.DEFAULT_SEED_LENGTH,
            force_cpu=(settings.GPU_ENABLED == "cpu")
        )

        # 6. Remove old results before inserting new results.
        db.query(Alignment).filter(
            Alignment.job_id == job.id
        ).delete(synchronize_session=False)

        db.commit()

        # 7. Process the genome one chunk at a time.
        processed_bases = 0
        processed_chunks = 0
        total_hits = 0
        batch_records = []

        for chromosome, chunk_start, chunk_sequence in iter_fasta_chunks(
            genome_path,
            chunk_size=FASTA_CHUNK_SIZE,
            overlap=overlap
        ):
            processed_chunks += 1

            chromosome_length = chromosome_lengths[chromosome]

            # Number of new genomic bases owned by this chunk.
            # The right-side overlap is not counted twice.
            core_length = min(
                FASTA_CHUNK_SIZE,
                chromosome_length - chunk_start
            )

            if core_length <= 0:
                continue

            # 8. Detect candidate PAM sites in this chunk.
            chunk_candidates = detector.find_candidate_sites(
                chromosome,
                chunk_sequence
            )

            # 9. Keep only candidates whose starting coordinate
            # belongs to this chunk's non-overlapping core.
            global_candidates = []

            for site in chunk_candidates:
                local_start = site.start_position - 1

                if local_start < 0 or local_start >= core_length:
                    continue

                # Convert chunk-local coordinates to chromosome-wide
                # coordinates. chunk_start is zero-based, while
                # CandidateSite positions are one-based.
                global_site = CandidateSite(
                    chromosome=site.chromosome,
                    start_position=site.start_position + chunk_start,
                    end_position=site.end_position + chunk_start,
                    strand=site.strand,
                    protospacer=site.protospacer,
                    pam=site.pam
                )

                global_candidates.append(global_site)

            # 10. Align candidates in manageable batches.
            for batch_start in range(
                0,
                len(global_candidates),
                ALIGNMENT_BATCH_SIZE
            ):
                candidate_batch = global_candidates[
                    batch_start:
                    batch_start + ALIGNMENT_BATCH_SIZE
                ]

                hits = aligner.align_batch(candidate_batch)

                # 11. Convert alignment hits into database records.
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
                        mismatch_positions=json.dumps(
                            hit.mismatch_positions
                        ),
                        seed_mismatch_count=hit.seed_mismatch_count,
                        score=hit.score,
                        risk_category=hit.risk_category
                    )

                    batch_records.append(record)
                    total_hits += 1

                    # Save records incrementally.
                    if len(batch_records) >= DB_INSERT_BATCH_SIZE:
                        db.bulk_save_objects(batch_records)
                        db.commit()
                        batch_records.clear()

            # 12. Update progress using unique genomic bases.
            processed_bases += core_length

            progress = 10 + int(
                80 * processed_bases / total_bases
            )

            job.progress_percent = min(progress, 90)
            db.commit()

            logger.info(
                "Job %s: processed chunk %d, chromosome %s, "
                "%d candidates, %d hits so far",
                job_id,
                processed_chunks,
                chromosome,
                len(global_candidates),
                total_hits
            )

        # 13. Save any remaining database records.
        if batch_records:
            db.bulk_save_objects(batch_records)
            db.commit()
            batch_records.clear()

        if processed_chunks == 0:
            raise ValueError("No FASTA chunks were processed.")

        # 14. Mark the job as completed.
        job.status = "completed"
        job.progress_percent = 100
        job.completed_at = datetime.utcnow()
        job.error_message = None
        db.commit()

        logger.info(
            "Job %s completed successfully. "
            "Processed %d chunks and saved %d alignment hits.",
            job_id,
            processed_chunks,
            total_hits
        )

    except Exception as exc:
        logger.exception(
            "Alignment job %s failed: %s",
            job_id,
            exc
        )

        try:
            db.rollback()

            job = db.query(Job).filter(
                Job.id == job_id
            ).first()

            if job:
                # Remove partial results from a failed run.
                db.query(Alignment).filter(
                    Alignment.job_id == job.id
                ).delete(synchronize_session=False)

                job.status = "failed"
                job.error_message = str(exc)
                job.completed_at = datetime.utcnow()
                db.commit()

        except Exception as db_error:
            db.rollback()
            logger.error(
                "Could not update failed job %s: %s",
                job_id,
                db_error
            )

    finally:
        db.close()
