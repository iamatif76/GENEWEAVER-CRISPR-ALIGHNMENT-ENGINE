"""
Endpoints for reference genome management and FASTA uploads.
"""

from typing import List, Optional
from fastapi import APIRouter, Depends, File, Form, HTTPException, UploadFile, status
from sqlalchemy.orm import Session

from app.database import get_db
from app.schemas import GenomeResponse
from app.services.genome_service import GenomeService
from app.core.security import verify_api_key

router = APIRouter(prefix="/genomes", tags=["Reference Genomes"])


@router.get("", response_model=List[GenomeResponse])
def list_reference_genomes(db: Session = Depends(get_db)):
    """List all available reference genomes uploaded to the system."""
    return GenomeService.list_genomes(db)


@router.post(
    "/upload",
    response_model=GenomeResponse,
    status_code=status.HTTP_201_CREATED,
    dependencies=[Depends(verify_api_key)],
)
async def upload_reference_genome(
    file: UploadFile = File(..., description="FASTA format genome file (.fasta, .fa, .fna)"),
    name: Optional[str] = Form(None, description="Optional custom human-readable name"),
    db: Session = Depends(get_db),
):
    """
    Upload a multi-chromosome FASTA reference genome.
    Validates FASTA formatting, extracts total bases and chromosome count,
    and stores metadata in the database.
    Requires API key authentication.
    """
    try:
        genome = await GenomeService.save_uploaded_genome(db, file, custom_name=name)
        return genome
    except ValueError as ve:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=str(ve)
        )
    except IOError as ioe:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Disk I/O error during file upload: {str(ioe)}"
        )
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Unexpected error while processing genome upload: {str(e)}"
        )


@router.get("/{genome_id}", response_model=GenomeResponse)
def get_reference_genome(genome_id: str, db: Session = Depends(get_db)):
    """Get metadata for a specific reference genome."""
    genome = GenomeService.get_genome_by_id(db, genome_id)
    if not genome:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Genome with ID '{genome_id}' not found."
        )
    return genome

