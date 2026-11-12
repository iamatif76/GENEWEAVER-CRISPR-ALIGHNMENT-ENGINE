"""
Integration tests for FastAPI REST API endpoints, security, and job workflows.
"""

import io
import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.config import settings
from app.database import init_db


@pytest.fixture(scope="session", autouse=True)
def setup_database():
    """Ensure database tables are initialized before running test suite."""
    init_db()


@pytest.fixture
def client():
    """Provides a FastAPI test client."""
    return TestClient(app)


@pytest.fixture
def auth_headers():
    """Returns headers with the configured API key if set."""
    if settings.API_KEY:
        return {"X-API-Key": settings.API_KEY}
    return {}


def test_root_endpoint(client):
    response = client.get("/")
    assert response.status_code == 200
    data = response.json()
    assert data["service"] == "GeneWeaver: GPU-Accelerated CRISPR Alignment Engine"
    assert "endpoints" in data


def test_health_endpoint(client):
    response = client.get("/api/v1/health")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "operational"
    assert "database" in data
    assert "cpu_mode" in data
    assert "gpu_available" in data


def test_list_genomes(client):
    response = client.get("/api/v1/genomes")
    assert response.status_code == 200
    assert isinstance(response.json(), list)


def test_api_key_security_unauthorized(client):
    """Verify protected endpoints reject requests when missing or invalid API key is provided."""
    if not settings.API_KEY:
        pytest.skip("API_KEY not configured in environment")

    # Missing API key
    response = client.post("/api/v1/genomes/upload")
    assert response.status_code == 401

    # Invalid API key
    response = client.post(
        "/api/v1/genomes/upload",
        headers={"X-API-Key": "invalid_wrong_token"}
    )
    assert response.status_code == 401


def test_upload_genome_and_run_alignment(client, auth_headers):
    """Test full upload, alignment creation, and result retrieval pipeline."""
    # 1. Upload small test FASTA
    fasta_content = (
        b">chrTest Test contig\n"
        b"ACGTACGTGAGTCCGAGCAGAAGAAGAATGGACGTACGT\n"
        b">chrTest2 Test contig 2\n"
        b"TTTTCCATTCTTCTTCTGCTCGGACTCGGGG\n"
    )
    files = {"file": ("test_genome.fasta", io.BytesIO(fasta_content), "text/plain")}
    data = {"name": "Test Integration Genome"}

    upload_resp = client.post(
        "/api/v1/genomes/upload",
        files=files,
        data=data,
        headers=auth_headers
    )
    assert upload_resp.status_code == 201
    genome_data = upload_resp.json()
    genome_id = genome_data["id"]
    assert genome_data["chromosome_count"] == 2

    # 2. Submit alignment job
    align_resp = client.post(
        "/api/v1/align",
        data={
            "guide_sequence": "GAGTCCGAGCAGAAGAAGAA",
            "genome_id": genome_id,
            "max_mismatches": 3,
            "pam_pattern": "NGG",
        },
        headers=auth_headers
    )
    assert align_resp.status_code == 202
    job_data = align_resp.json()
    job_id = job_data["id"]
    assert job_data["status"] in ["queued", "running", "completed"]

    # 3. Retrieve Job Status
    job_status_resp = client.get(f"/api/v1/jobs/{job_id}")
    assert job_status_resp.status_code == 200
    assert job_status_resp.json()["id"] == job_id


def test_align_invalid_guide_sequence(client, auth_headers):
    # Invalid base 'Z'
    resp = client.post(
        "/api/v1/align",
        data={
            "guide_sequence": "GAGTCCGAGCAGAZGAAGAA",
            "max_mismatches": 3,
            "pam_pattern": "NGG",
        },
        headers=auth_headers
    )
    assert resp.status_code == 422
