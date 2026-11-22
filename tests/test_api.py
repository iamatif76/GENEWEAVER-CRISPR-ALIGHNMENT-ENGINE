
"""
Integration tests for FastAPI REST API endpoints, security, and job workflows.
"""

import io

import pytest
import app.workers.alignment_worker as alignment_worker
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
    """Verify protected endpoints reject requests without a valid API key."""
    if not settings.API_KEY:
        pytest.skip("API_KEY not configured in environment")

    response = client.post("/api/v1/genomes/upload")
    assert response.status_code == 401

    response = client.post(
        "/api/v1/genomes/upload",
        headers={"X-API-Key": "invalid_wrong_token"}
    )
    assert response.status_code == 401


def test_upload_genome_and_run_alignment(client, auth_headers):
    """Test full upload, alignment creation, and result retrieval pipeline."""
    fasta_content = (
        b">chrTest Test contig\n"
        b"ACGTACGTGAGTCCGAGCAGAAGAAGAATGGACGTACGT\n"
        b">chrTest2 Test contig 2\n"
        b"TTTTCCATTCTTCTTCTGCTCGGACTCGGGG\n"
    )

    files = {
        "file": ("test_genome.fasta", io.BytesIO(fasta_content), "text/plain")
    }
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

    job_status_resp = client.get(f"/api/v1/jobs/{job_id}")
    assert job_status_resp.status_code == 200
    assert job_status_resp.json()["id"] == job_id


def test_align_invalid_guide_sequence(client, auth_headers):
    """Verify invalid guide sequences are rejected."""
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


def test_chunked_worker_finds_boundary_candidate(client, auth_headers, monkeypatch):
    """Verify that a guide crossing a chunk boundary is found at its global coordinates."""

    monkeypatch.setattr(alignment_worker, "FASTA_CHUNK_SIZE", 10)

    guide = "GAGTCCGAGCAGAAGAAGAA"
    sequence = "A" * 9 + guide + "AGG" + "A" * 10

    fasta_content = (
        f">chrBoundary Test chromosome\n{sequence}\n"
    ).encode()

    upload_resp = client.post(
        "/api/v1/genomes/upload",
        files={
            "file": (
                "boundary_test.fasta",
                io.BytesIO(fasta_content),
                "text/plain",
            )
        },
        data={"name": "Chunk Boundary Test Genome"},
        headers=auth_headers,
    )

    assert upload_resp.status_code == 201
    genome_id = upload_resp.json()["id"]

    align_resp = client.post(
        "/api/v1/align",
        data={
            "guide_sequence": guide,
            "genome_id": genome_id,
            "max_mismatches": 0,
            "pam_pattern": "NGG",
        },
        headers=auth_headers,
    )

    assert align_resp.status_code == 202
    job_id = align_resp.json()["id"]

    status_resp = client.get(f"/api/v1/jobs/{job_id}")

    assert status_resp.status_code == 200
    assert status_resp.json()["status"] == "completed"

    results_resp = client.get(f"/api/v1/jobs/{job_id}/results")

    assert results_resp.status_code == 200
    results = results_resp.json()

    boundary_hits = [
        hit for hit in results
        if hit["chromosome"] == "chrBoundary"
        and hit["strand"] == "+"
        and hit["start_position"] == 10
        and hit["end_position"] == 32
    ]

    assert len(boundary_hits) == 1
    assert boundary_hits[0]["mismatch_count"] == 0
    assert boundary_hits[0]["matched_sequence"] == guide
    assert boundary_hits[0]["pam_sequence"] == "AGG"

def test_worker_streams_candidates_in_bounded_batches(
    client, auth_headers, monkeypatch
):
    """Verify candidates are aligned in small bounded batches."""

    monkeypatch.setattr(
        alignment_worker, "ALIGNMENT_BATCH_SIZE", 2
    )

    original_align_batch = alignment_worker.GPUAligner.align_batch
    observed_batch_sizes = []

    def spy_align_batch(self, sites, batch_size=100000):
        observed_batch_sizes.append(len(sites))
        return original_align_batch(
            self, sites, batch_size=batch_size
        )

    monkeypatch.setattr(
        alignment_worker.GPUAligner,
        "align_batch",
        spy_align_batch
    )

    guide = "A" * 20
    sequence = "TTTTT".join([guide + "AGG"] * 5)

    upload_resp = client.post(
        "/api/v1/genomes/upload",
        files={
            "file": (
                "batch_test.fasta",
                io.BytesIO(f">chrBatchTest\n{sequence}\n".encode()),
                "text/plain"
            )
        },
        data={"name": "Streaming Batch Test Genome"},
        headers=auth_headers
    )

    assert upload_resp.status_code == 201
    genome_id = upload_resp.json()["id"]

    align_resp = client.post(
        "/api/v1/align",
        data={
            "guide_sequence": guide,
            "genome_id": genome_id,
            "max_mismatches": 0,
            "pam_pattern": "NGG"
        },
        headers=auth_headers
    )

    assert align_resp.status_code == 202
    job_id = align_resp.json()["id"]

    status_resp = client.get(f"/api/v1/jobs/{job_id}")

    assert status_resp.status_code == 200
    assert status_resp.json()["status"] == "completed"

    assert observed_batch_sizes == [2, 2, 1]

    results_resp = client.get(
        f"/api/v1/jobs/{job_id}/results"
    )

    assert results_resp.status_code == 200
    assert len(results_resp.json()) == 5


@pytest.mark.parametrize("candidate_count", [257, 1000, 5000])
def test_worker_handles_large_candidate_stream(
    client, auth_headers, monkeypatch, candidate_count
):
    """Verify large candidate streams use bounded alignment batches."""

    monkeypatch.setattr(
        alignment_worker, "ALIGNMENT_BATCH_SIZE", 64
    )

    original_align_batch = alignment_worker.GPUAligner.align_batch
    observed_batch_sizes = []

    def spy_align_batch(self, sites, batch_size=100000):
        observed_batch_sizes.append(len(sites))
        return original_align_batch(
            self, sites, batch_size=batch_size
        )

    monkeypatch.setattr(
        alignment_worker.GPUAligner,
        "align_batch",
        spy_align_batch
    )

    guide = "A" * 20
    sequence = "TTTTT".join(
        [guide + "AGG"] * candidate_count
    )

    upload_resp = client.post(
        "/api/v1/genomes/upload",
        files={
            "file": (
                "stress_test.fasta",
                io.BytesIO(
                    f">chrStress\n{sequence}\n".encode()
                ),
                "text/plain"
            )
        },
        data={"name": "Streaming Stress Test Genome"},
        headers=auth_headers
    )

    assert upload_resp.status_code == 201
    genome_id = upload_resp.json()["id"]

    align_resp = client.post(
        "/api/v1/align",
        data={
            "guide_sequence": guide,
            "genome_id": genome_id,
            "max_mismatches": 0,
            "pam_pattern": "NGG"
        },
        headers=auth_headers
    )

    assert align_resp.status_code == 202
    job_id = align_resp.json()["id"]

    status_resp = client.get(f"/api/v1/jobs/{job_id}")

    assert status_resp.status_code == 200
    assert status_resp.json()["status"] == "completed"

    expected_batches = (
        [64] * (candidate_count // 64)
        + ([candidate_count % 64] if candidate_count % 64 else [])
    )

    assert observed_batch_sizes == expected_batches

    results_resp = client.get(
        f"/api/v1/jobs/{job_id}/results"
    )

    assert results_resp.status_code == 200

    results = results_resp.json()

    assert len(results) == candidate_count
    assert len({
        result["start_position"] for result in results
    }) == candidate_count

    assert all(
        result["mismatch_count"] == 0 for result in results
    )
