"""
Pydantic schemas for request and response validation in GeneWeaver API.
"""

from datetime import datetime
import json
from typing import List, Optional, Union
from pydantic import BaseModel, Field, field_validator, ConfigDict


class AlignmentResponse(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: str
    job_id: str
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

    @field_validator("mismatch_positions", mode="before")
    @classmethod
    def parse_mismatch_positions(cls, v: Union[str, List[int]]) -> List[int]:
        if isinstance(v, str):
            try:
                return json.loads(v)
            except Exception:
                return []
        return v or []


class GenomeResponse(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: str
    name: str
    filename: str
    total_bases: int
    chromosome_count: int
    uploaded_at: datetime


class JobCreateRequest(BaseModel):
    guide_sequence: str = Field(..., description="Guide RNA or DNA sequence (~20 nucleotides)")
    genome_id: Optional[str] = Field(None, description="Existing reference genome ID")
    max_mismatches: int = Field(3, ge=0, le=10, description="Maximum allowed mismatches")
    pam_pattern: str = Field("NGG", description="PAM pattern with standard IUPAC codes")

    @field_validator("guide_sequence")
    @classmethod
    def validate_guide(cls, v: str) -> str:
        cleaned = v.strip().upper().replace("U", "T")
        if not cleaned:
            raise ValueError("Guide sequence cannot be empty")
        valid_bases = set("ACGT")
        if not set(cleaned).issubset(valid_bases):
            invalid = set(cleaned) - valid_bases
            raise ValueError(f"Invalid nucleotide characters in guide sequence: {invalid}. Only A, C, G, T, U allowed.")
        if len(cleaned) < 15 or len(cleaned) > 30:
            raise ValueError(f"Guide length must be between 15 and 30 nucleotides (got {len(cleaned)})")
        return cleaned

    @field_validator("pam_pattern")
    @classmethod
    def validate_pam(cls, v: str) -> str:
        cleaned = v.strip().upper().replace("U", "T")
        if not cleaned:
            raise ValueError("PAM pattern cannot be empty")
        valid_iupac = set("ACGTURYSWKMBDHVN")
        if not set(cleaned).issubset(valid_iupac):
            invalid = set(cleaned) - valid_iupac
            raise ValueError(f"Invalid characters in PAM pattern: {invalid}. IUPAC codes allowed.")
        return cleaned


class JobSummary(BaseModel):
    total_alignments: int = 0
    exact_matches: int = 0
    high_risk_count: int = 0
    moderate_risk_count: int = 0
    low_risk_count: int = 0


class JobResponse(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: str
    guide_sequence: str
    genome_id: str
    max_mismatches: int
    pam_pattern: str
    status: str
    progress_percent: int
    error_message: Optional[str] = None
    created_at: datetime
    completed_at: Optional[datetime] = None
    summary: Optional[JobSummary] = None


class HealthResponse(BaseModel):
    status: str
    database: str
    cpu_mode: str
    gpu_available: bool
    gpu_device_name: Optional[str] = None
    version: str
    timestamp: datetime = Field(default_factory=datetime.utcnow)


class BenchmarkMetrics(BaseModel):
    positions_scanned: int
    candidate_pams: int
    alignments_reported: int
    cpu_time_seconds: float
    gpu_time_seconds: Optional[float] = None
    speedup_factor: Optional[float] = None
    memory_peak_mb: float
    accuracy_percent: float
    hardware_info: str
