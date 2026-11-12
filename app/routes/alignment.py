"""
Endpoints for initiating CRISPR alignment jobs, tracking progress,
retrieving ranked hits, and exporting data to CSV/JSON.
"""

from typing import List, Optional
from fastapi import (
    APIRouter,
    BackgroundTasks,
    Depends,
    File,
    Form,
    HTTPException,
    Query,
    Response,
    UploadFile,
    status,
)
from sqlalchemy.orm import Session

from app.database import get_db
from app.schemas import (
    AlignmentResponse,
    JobCreateRequest,
    JobResponse,
)
from app.services.alignment_service import AlignmentService
from app.services.genome_service import GenomeService
from app.services.export_service import export_to_csv, export_to_json
from app.core.sequence_utils import clean_and_validate_sequence, iupac_to_regex
from app.core.security import verify_api_key
from app.workers.alignment_worker import run_alignment_job_worker

router = APIRouter(tags=["CRISPR Alignment Engine"])


@router.post(
    "/align",
    response_model=JobResponse,
    status_code=status.HTTP_202_ACCEPTED,
    dependencies=[Depends(verify_api_key)],
)
async def submit_alignment_job(
    background_tasks: BackgroundTasks,
    guide_sequence: str = Form(..., description="Guide RNA/DNA sequence (~20 nucleotides)"),
    genome_id: Optional[str] = Form(None, description="Registered reference genome ID"),
    file: Optional[UploadFile] = File(None, description="Optional FASTA file upload"),
    max_mismatches: int = Form(3, ge=0, le=10, description="Maximum mismatches allowed"),
    pam_pattern: str = Form("NGG", description="PAM pattern with IUPAC codes"),
    db: Session = Depends(get_db),
):
    """
    Submit a CRISPR guide RNA sequence for genome-wide alignment and off-target risk assessment.
    Accepts either an existing reference genome ID or an uploaded FASTA file.
    Runs asynchronously as a background task.
    Requires API key authentication.
    """
    # 1. Validate Guide Sequence
    try:
        clean_guide = clean_and_validate_sequence(guide_sequence)
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=f"Invalid guide sequence: {str(e)}"
        )

    # 2. Validate PAM pattern
    try:
        iupac_to_regex(pam_pattern)
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=f"Invalid PAM pattern: {str(e)}"
        )

    # 3. Resolve Reference Genome
    target_genome_id = genome_id
    if file and file.filename:
        # User uploaded a new FASTA file with the request
        try:
            uploaded_genome = await GenomeService.save_uploaded_genome(db, file)
            target_genome_id = uploaded_genome.id
        except Exception as e:
            raise HTTPException(
                status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
                detail=f"Could not process uploaded genome: {str(e)}"
            )

    if not target_genome_id:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Either a valid 'genome_id' or an uploaded FASTA 'file' must be provided."
        )

    genome = GenomeService.get_genome_by_id(db, target_genome_id)
    if not genome:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Reference genome with ID '{target_genome_id}' not found."
        )

    # 4. Create Job in Database
    job = AlignmentService.create_job(
        db=db,
        guide_sequence=clean_guide,
        genome_id=target_genome_id,
        max_mismatches=max_mismatches,
        pam_pattern=pam_pattern.strip().upper(),
    )

    # 5. Dispatch Background Worker
    background_tasks.add_task(run_alignment_job_worker, job.id)

    # 6. Return Accepted Job Response
    response = AlignmentService.get_job_response(db, job.id)
    return response


@router.post(
    "/align/json",
    response_model=JobResponse,
    status_code=status.HTTP_202_ACCEPTED,
    dependencies=[Depends(verify_api_key)],
)
def submit_alignment_job_json(
    job_req: JobCreateRequest,
    background_tasks: BackgroundTasks,
    db: Session = Depends(get_db),
):
    """
    Submit an alignment job using a JSON payload with an existing genome_id.
    Requires API key authentication.
    """
    if not job_req.genome_id:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="'genome_id' is required for JSON job submission."
        )

    genome = GenomeService.get_genome_by_id(db, job_req.genome_id)
    if not genome:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Reference genome with ID '{job_req.genome_id}' not found."
        )

    job = AlignmentService.create_job(
        db=db,
        guide_sequence=job_req.guide_sequence,
        genome_id=job_req.genome_id,
        max_mismatches=job_req.max_mismatches,
        pam_pattern=job_req.pam_pattern,
    )

    background_tasks.add_task(run_alignment_job_worker, job.id)
    return AlignmentService.get_job_response(db, job.id)


@router.get("/jobs/{job_id}", response_model=JobResponse)
def get_job_status(job_id: str, db: Session = Depends(get_db)):
    """
    Get alignment job status, progress percentage, error messages, and summary metrics.
    """
    job_resp = AlignmentService.get_job_response(db, job_id)
    if not job_resp:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Job with ID '{job_id}' not found."
        )
    return job_resp


@router.get("/jobs/{job_id}/results", response_model=List[AlignmentResponse])
def get_job_results(
    job_id: str,
    limit: Optional[int] = Query(None, ge=1, le=10000, description="Max records to return"),
    offset: int = Query(0, ge=0, description="Record offset for pagination"),
    chromosome: Optional[str] = Query(None, description="Filter by chromosome"),
    strand: Optional[str] = Query(None, pattern="^[+-]$", description="Filter by strand (+ or -)"),
    db: Session = Depends(get_db),
):
    """
    Retrieve ranked alignment results for a completed job in JSON format.
    Ranked by off-target score descending, then mismatch count ascending.
    """
    job = AlignmentService.get_job(db, job_id)
    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Job with ID '{job_id}' not found."
        )

    alignments = AlignmentService.get_job_alignments(
        db=db,
        job_id=job_id,
        limit=limit,
        offset=offset,
        chromosome=chromosome,
        strand=strand,
    )
    return alignments


@router.get("/jobs/{job_id}/export")
def export_job_results(
    job_id: str,
    format: str = Query("csv", pattern="^(csv|json)$", description="Export format ('csv' or 'json')"),
    db: Session = Depends(get_db),
):
    """
    Download alignment results in CSV or JSON format.
    """
    job = AlignmentService.get_job(db, job_id)
    if not job:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Job with ID '{job_id}' not found."
        )

    alignments = AlignmentService.get_job_alignments(db=db, job_id=job_id)

    if format.lower() == "csv":
        csv_data = export_to_csv(alignments)
        filename = f"geneweaver_alignment_{job_id}.csv"
        return Response(
            content=csv_data,
            media_type="text/csv",
            headers={"Content-Disposition": f"attachment; filename={filename}"},
        )
    else:
        json_data = export_to_json(alignments)
        filename = f"geneweaver_alignment_{job_id}.json"
        return Response(
            content=json_data,
            media_type="application/json",
            headers={"Content-Disposition": f"attachment; filename={filename}"},
        )
