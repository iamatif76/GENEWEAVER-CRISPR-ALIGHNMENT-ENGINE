"""
Service for managing reference genome uploads, validation, and metadata extraction.
"""

import shutil
import uuid
from pathlib import Path
from typing import List, Optional
from fastapi import UploadFile
from sqlalchemy.orm import Session

from app.config import settings
from app.core.genome_loader import inspect_fasta, GenomeLoadError
from app.models import Genome


class GenomeService:
    """Manages reference genome storage and lifecycle."""

    @staticmethod
    def get_genome_by_id(db: Session, genome_id: str) -> Optional[Genome]:
        """Fetch a reference genome record by UUID."""
        return db.query(Genome).filter(Genome.id == genome_id).first()

    @staticmethod
    def list_genomes(db: Session) -> List[Genome]:
        """List all registered reference genomes ordered by upload timestamp."""
        return db.query(Genome).order_by(Genome.uploaded_at.desc()).all()

    @staticmethod
    async def save_uploaded_genome(
        db: Session,
        file: UploadFile,
        custom_name: Optional[str] = None
    ) -> Genome:
        """
        Validate, save to disk, and register an uploaded FASTA genome.
        """
        filename = file.filename or "genome.fasta"
        if not (filename.endswith(".fasta") or filename.endswith(".fa") or filename.endswith(".fna")):
            raise ValueError("Only FASTA files (.fasta, .fa, .fna) are supported.")

        genome_id = str(uuid.uuid4())
        dest_filename = f"{genome_id}_{filename}"
        dest_path = settings.UPLOAD_DIR / dest_filename

        # Write uploaded file in chunks
        try:
            with open(dest_path, "wb") as buffer:
                shutil.copyfileobj(file.file, buffer)
        except Exception as e:
            if dest_path.exists():
                dest_path.unlink()
            raise IOError(f"Failed to write uploaded file to disk: {str(e)}") from e

        # Validate and extract metrics
        try:
            total_bases, chrom_count, _ = inspect_fasta(dest_path)
        except GenomeLoadError as e:
            if dest_path.exists():
                dest_path.unlink()
            raise ValueError(f"Uploaded FASTA file is invalid: {str(e)}") from e

        name = custom_name.strip() if custom_name and custom_name.strip() else Path(filename).stem

        genome_record = Genome(
            id=genome_id,
            name=name,
            filename=filename,
            file_path=str(dest_path),
            total_bases=total_bases,
            chromosome_count=chrom_count,
        )

        db.add(genome_record)
        db.commit()
        db.refresh(genome_record)
        return genome_record

    @staticmethod
    def register_existing_fasta(
        db: Session,
        file_path: Path,
        name: str
    ) -> Genome:
        """Register an existing FASTA file already on disk (e.g. sample data)."""
        path = Path(file_path).resolve()
        if not path.is_file():
            raise FileNotFoundError(f"FASTA file not found: {path}")

        # Check if already registered
        existing = db.query(Genome).filter(Genome.file_path == str(path)).first()
        if existing:
            return existing

        total_bases, chrom_count, _ = inspect_fasta(path)
        genome_record = Genome(
            id=str(uuid.uuid4()),
            name=name,
            filename=path.name,
            file_path=str(path),
            total_bases=total_bases,
            chromosome_count=chrom_count,
        )
        db.add(genome_record)
        db.commit()
        db.refresh(genome_record)
        return genome_record

