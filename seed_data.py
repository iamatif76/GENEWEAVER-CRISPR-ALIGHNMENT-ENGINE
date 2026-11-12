"""
Database seeding script for GeneWeaver.
Initializes tables and registers the sample reference genome.
"""

import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from app.database import init_db, SessionLocal
from app.config import settings
from app.services.genome_service import GenomeService
from app.workers.alignment_worker import run_alignment_job_worker
from app.services.alignment_service import AlignmentService


def seed():
    print("Initializing database...")
    init_db()

    db = SessionLocal()
    try:
        sample_path = settings.DATA_DIR / "sample_genome.fasta"
        if not sample_path.is_file():
            print(f"Error: Sample FASTA not found at {sample_path}")
            return

        print(f"Registering sample genome from {sample_path.name}...")
        genome = GenomeService.register_existing_fasta(
            db=db,
            file_path=sample_path,
            name="Synthetic Human Benchmark Genome (Sample)"
        )
        print(f"Genome registered successfully! ID: {genome.id} ({genome.total_bases} bases, {genome.chromosome_count} contigs)")

        # Create a sample demo alignment job for immediate exploration
        guide_seq = "GAGTCCGAGCAGAAGAAGAA"  # EMX1 on-target guide
        print(f"Creating demonstration alignment job for guide '{guide_seq}'...")
        job = AlignmentService.create_job(
            db=db,
            guide_sequence=guide_seq,
            genome_id=genome.id,
            max_mismatches=3,
            pam_pattern="NGG",
        )
        print(f"Demo Job created (ID: {job.id}). Running alignment...")
        run_alignment_job_worker(job.id)

        # Print summary
        summary = AlignmentService.get_job_summary(db, job.id)
        print(f"Alignment complete!")
        print(f"  Total Hits: {summary.total_alignments}")
        print(f"  Exact Matches: {summary.exact_matches}")
        print(f"  High-risk Off-targets: {summary.high_risk_count}")
        print(f"  Moderate-risk Off-targets: {summary.moderate_risk_count}")
        print(f"  Low-risk Off-targets: {summary.low_risk_count}")

    finally:
        db.close()


if __name__ == "__main__":
    seed()
