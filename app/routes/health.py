"""
Health check and hardware status endpoints.
"""

from fastapi import APIRouter, Depends
from sqlalchemy import text
from sqlalchemy.orm import Session

from app.database import get_db
from app.schemas import HealthResponse
from app.core.gpu_utils import check_gpu_availability

router = APIRouter(prefix="/health", tags=["Health & Status"])


@router.get("", response_model=HealthResponse)
def get_health(db: Session = Depends(get_db)):
    """
    Check operational status of GeneWeaver services:
    - API operational state
    - Database connectivity
    - CPU vectorization engine
    - GPU acceleration availability
    """
    # Check DB
    db_status = "healthy"
    try:
        db.execute(text("SELECT 1"))
    except Exception as e:
        db_status = f"unhealthy ({str(e)})"

    # Check GPU
    gpu_avail, dev_name, _ = check_gpu_availability()

    return HealthResponse(
        status="operational",
        database=db_status,
        cpu_mode="Vectorized NumPy SIMD",
        gpu_available=gpu_avail,
        gpu_device_name=dev_name if gpu_avail else None,
        version="1.0.0",
    )
