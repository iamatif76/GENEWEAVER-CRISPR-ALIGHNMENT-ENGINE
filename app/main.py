"""
GeneWeaver: GPU-Accelerated CRISPR Alignment Engine
Main FastAPI application entrypoint.
"""

from contextlib import asynccontextmanager
import logging
from pathlib import Path
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import RedirectResponse

from app.config import settings
from app.database import init_db, SessionLocal
from app.routes.health import router as health_router
from app.routes.genomes import router as genomes_router
from app.routes.alignment import router as alignment_router
from app.core.gpu_utils import get_hardware_status
from app.services.genome_service import GenomeService

# Configure structured logging
logging.basicConfig(
    level=getattr(logging, settings.LOG_LEVEL.upper(), logging.INFO),
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("geneweaver")


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application startup and shutdown lifecycle management."""
    logger.info("Initializing GeneWeaver CRISPR Alignment Engine...")

    # 1. Initialize SQLite Database
    init_db()
    logger.info("Database schemas verified and initialized.")

    # 2. Probe hardware acceleration capabilities
    hw = get_hardware_status()
    logger.info("Hardware execution environment: CPU=%s, GPU_Available=%s", hw["cpu_mode"], hw["gpu_available"])
    if hw["gpu_available"]:
        logger.info("CUDA device: %s", hw["gpu_device_name"])
    else:
        logger.info("Running on CPU SIMD fallback mode.")

    # 3. Seed sample reference genome if present on disk
    sample_fasta = settings.DATA_DIR / "sample_genome.fasta"
    if sample_fasta.is_file():
        db = SessionLocal()
        try:
            GenomeService.register_existing_fasta(
                db=db,
                file_path=sample_fasta,
                name="Synthetic Human Benchmark Genome (Sample)"
            )
            logger.info("Registered default sample genome: %s", sample_fasta.name)
        except Exception as e:
            logger.warning("Could not auto-register sample genome: %s", str(e))
        finally:
            db.close()

    yield
    logger.info("GeneWeaver shutting down gracefully.")


app = FastAPI(
    title="GeneWeaver: GPU-Accelerated CRISPR Alignment Engine",
    description=(
        "High-performance bioinformatics engine for CRISPR gRNA on-target and off-target "
        "search, PAM detection, mismatch scoring, seed region weighting, and SQLite tracking."
    ),
    version="1.0.0",
    docs_url="/docs",
    redoc_url="/redoc",
    lifespan=lifespan,
)

# Enable CORS for Streamlit frontend and local API clients
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount API Routers under versioned prefix /api/v1
app.include_router(health_router, prefix="/api/v1")
app.include_router(genomes_router, prefix="/api/v1")
app.include_router(alignment_router, prefix="/api/v1")


@app.get("/", tags=["Root"])
def root():
    """Root endpoint providing service overview and documentation links."""
    return {
        "service": "GeneWeaver: GPU-Accelerated CRISPR Alignment Engine",
        "version": "1.0.0",
        "status": "online",
        "documentation": "/docs",
        "health_check": "/api/v1/health",
        "endpoints": {
            "align": "POST /api/v1/align",
            "job_status": "GET /api/v1/jobs/{job_id}",
            "job_results": "GET /api/v1/jobs/{job_id}/results",
            "job_export": "GET /api/v1/jobs/{job_id}/export?format=csv",
            "genomes": "GET /api/v1/genomes",
            "genome_upload": "POST /api/v1/genomes/upload",
            "health": "GET /api/v1/health",
        },
    }
