"""
SQLAlchemy database models for GeneWeaver.
Defines schemas for jobs, genomes, and alignment results.
"""

import uuid
from datetime import datetime
from sqlalchemy import (
    Column,
    String,
    Integer,
    Float,
    DateTime,
    ForeignKey,
    Text,
)
from sqlalchemy.orm import relationship

from app.database import Base


def generate_uuid() -> str:
    """Generate a UUID4 string."""
    return str(uuid.uuid4())


class Genome(Base):
    """Stores reference genome metadata and file locations."""
    __tablename__ = "genomes"

    id = Column(String(36), primary_key=True, default=generate_uuid)
    name = Column(String(255), nullable=False)
    filename = Column(String(255), nullable=False)
    file_path = Column(String(1024), nullable=False)
    total_bases = Column(Integer, default=0, nullable=False)
    chromosome_count = Column(Integer, default=0, nullable=False)
    uploaded_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    # Relationships
    jobs = relationship("Job", back_populates="genome", cascade="all, delete-orphan")


class Job(Base):
    """Stores CRISPR alignment job parameters, execution state, and progress."""
    __tablename__ = "jobs"

    id = Column(String(36), primary_key=True, default=generate_uuid)
    guide_sequence = Column(String(128), nullable=False)
    genome_id = Column(String(36), ForeignKey("genomes.id", ondelete="CASCADE"), nullable=False)
    max_mismatches = Column(Integer, default=3, nullable=False)
    pam_pattern = Column(String(32), default="NGG", nullable=False)
    status = Column(String(32), default="queued", nullable=False)  # queued, running, completed, failed
    progress_percent = Column(Integer, default=0, nullable=False)
    error_message = Column(Text, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    completed_at = Column(DateTime, nullable=True)

    # Relationships
    genome = relationship("Genome", back_populates="jobs")
    alignments = relationship("Alignment", back_populates="job", cascade="all, delete-orphan")


class Alignment(Base):
    """Stores individual candidate on-target and off-target alignment hits."""
    __tablename__ = "alignments"

    id = Column(String(36), primary_key=True, default=generate_uuid)
    job_id = Column(String(36), ForeignKey("jobs.id", ondelete="CASCADE"), nullable=False)
    chromosome = Column(String(128), nullable=False)
    start_position = Column(Integer, nullable=False)
    end_position = Column(Integer, nullable=False)
    strand = Column(String(1), nullable=False)  # "+" or "-"
    matched_sequence = Column(String(128), nullable=False)
    pam_sequence = Column(String(32), nullable=False)
    mismatch_count = Column(Integer, nullable=False)
    mismatch_positions = Column(Text, default="[]", nullable=False)  # JSON-serialized list of 1-based indices
    seed_mismatch_count = Column(Integer, nullable=False)
    score = Column(Float, nullable=False)
    risk_category = Column(String(64), nullable=False)  # Exact match, High-risk, Moderate-risk, Low-risk

    # Relationships
    job = relationship("Job", back_populates="alignments")

