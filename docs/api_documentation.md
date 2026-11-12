# GeneWeaver REST API Documentation

Version: `v1`  
Base URL: `/api/v1`  
Interactive Swagger UI: `http://localhost:8000/docs`  
ReDoc: `http://localhost:8000/redoc`

---

## 1. Authentication

Protected endpoints require an API key passed via either:
1. **HTTP Header**: `X-API-Key: <your_api_key>`
2. **Authorization Header**: `Authorization: Bearer <your_api_key>`

If the API key is invalid or absent on a protected route, the server responds with:
```json
{
  "detail": "Missing API key. Provide via 'X-API-Key' header or 'Authorization: Bearer <key>'."
}
```
HTTP Status: `401 Unauthorized`

---

## 2. API Endpoints

### 2.1 Health & Status Check
**`GET /api/v1/health`**
- **Authentication**: None (Public)
- **Description**: Returns server operational status, database connectivity, CPU SIMD vectorization mode, and CUDA GPU availability.

#### Response (200 OK):
```json
{
  "status": "operational",
  "database": "healthy",
  "cpu_mode": "Vectorized NumPy SIMD",
  "gpu_available": false,
  "gpu_device_name": null,
  "version": "1.0.0",
  "timestamp": "2026-11-12T12:00:00Z"
}
```

---

### 2.2 List Reference Genomes
**`GET /api/v1/genomes`**
- **Authentication**: None (Public)
- **Description**: Lists all available reference genomes indexed in SQLite.

#### Response (200 OK):
```json
[
  {
    "id": "b4c3eaab-72ce-4f69-89ad-4173c612736f",
    "name": "Synthetic Human Benchmark Genome (Sample)",
    "filename": "sample_genome.fasta",
    "total_bases": 5617,
    "chromosome_count": 4,
    "uploaded_at": "2026-11-12T12:00:00Z"
  }
]
```

---

### 2.3 Upload Reference Genome
**`POST /api/v1/genomes/upload`**
- **Authentication**: Required (`X-API-Key`)
- **Content-Type**: `multipart/form-data`
- **Parameters**:
  - `file` (File, required): Multi-chromosome FASTA file (`.fasta`, `.fa`, `.fna`).
  - `name` (Form string, optional): Human-readable descriptive name.

#### Response (201 Created):
```json
{
  "id": "c89f5643-98fe-4e0e-9274-1296bf11e9df",
  "name": "Human Chromosome 1 Subset",
  "filename": "chr1_subset.fasta",
  "total_bases": 250000,
  "chromosome_count": 1,
  "uploaded_at": "2026-11-12T12:05:00Z"
}
```

#### Error Codes:
- `401 Unauthorized`: Missing or invalid API key.
- `422 Unprocessable Entity`: Invalid nucleotide characters or malformed FASTA format.

---

### 2.4 Submit Alignment Job
**`POST /api/v1/align`**
- **Authentication**: Required (`X-API-Key`)
- **Content-Type**: `multipart/form-data`
- **Parameters**:
  - `guide_sequence` (string, required): ~20 nt guide RNA or DNA sequence.
  - `genome_id` (string, optional): ID of existing reference genome.
  - `file` (file, optional): Upload new FASTA file directly.
  - `max_mismatches` (integer, default: 3): Max allowed mismatches (0-10).
  - `pam_pattern` (string, default: "NGG"): PAM motif using IUPAC codes.

#### Response (202 Accepted):
```json
{
  "id": "70c956d9-bb2c-4481-bf58-297572a2c9eb",
  "guide_sequence": "GAGTCCGAGCAGAAGAAGAA",
  "genome_id": "b4c3eaab-72ce-4f69-89ad-4173c612736f",
  "max_mismatches": 3,
  "pam_pattern": "NGG",
  "status": "queued",
  "progress_percent": 0,
  "error_message": null,
  "created_at": "2026-11-12T12:10:00Z",
  "completed_at": null,
  "summary": null
}
```

---

### 2.5 Submit Alignment Job (JSON Payload)
**`POST /api/v1/align/json`**
- **Authentication**: Required (`X-API-Key`)
- **Content-Type**: `application/json`
- **Request Body**:
```json
{
  "guide_sequence": "GAGTCCGAGCAGAAGAAGAA",
  "genome_id": "b4c3eaab-72ce-4f69-89ad-4173c612736f",
  "max_mismatches": 3,
  "pam_pattern": "NGG"
}
```

#### Response (202 Accepted): Same format as `POST /api/v1/align`.

---

### 2.6 Get Job Status & Summary
**`GET /api/v1/jobs/{job_id}`**
- **Authentication**: None (Public)
- **Description**: Returns execution state, completion percentage, and summary statistics when complete.

#### Response (200 OK):
```json
{
  "id": "70c956d9-bb2c-4481-bf58-297572a2c9eb",
  "guide_sequence": "GAGTCCGAGCAGAAGAAGAA",
  "genome_id": "b4c3eaab-72ce-4f69-89ad-4173c612736f",
  "max_mismatches": 3,
  "pam_pattern": "NGG",
  "status": "completed",
  "progress_percent": 100,
  "error_message": null,
  "created_at": "2026-11-12T12:10:00Z",
  "completed_at": "2026-11-12T12:10:01Z",
  "summary": {
    "total_alignments": 8,
    "exact_matches": 4,
    "high_risk_count": 2,
    "moderate_risk_count": 1,
    "low_risk_count": 1
  }
}
```

---

### 2.7 Get Ranked Alignment Results
**`GET /api/v1/jobs/{job_id}/results`**
- **Authentication**: None (Public)
- **Query Parameters**:
  - `limit` (int, optional): Maximum records to retrieve.
  - `offset` (int, default: 0): Pagination offset.
  - `chromosome` (string, optional): Filter by contig/chromosome name.
  - `strand` (string, optional): Filter by strand (`+` or `-`).

#### Response (200 OK):
```json
[
  {
    "id": "e427b3d0-3796-41e9-9130-1b777a9e332d",
    "job_id": "70c956d9-bb2c-4481-bf58-297572a2c9eb",
    "chromosome": "chr1",
    "start_position": 101,
    "end_position": 123,
    "strand": "+",
    "matched_sequence": "GAGTCCGAGCAGAAGAAGAA",
    "pam_sequence": "TGG",
    "mismatch_count": 0,
    "mismatch_positions": [],
    "seed_mismatch_count": 0,
    "score": 100.0,
    "risk_category": "Exact match"
  },
  {
    "id": "5fa23e3e-3f59-4d3f-a3c3-305141151a66",
    "job_id": "70c956d9-bb2c-4481-bf58-297572a2c9eb",
    "chromosome": "chr1",
    "start_position": 324,
    "end_position": 346,
    "strand": "+",
    "matched_sequence": "GCGTCCGAGCAGAAGAAGAA",
    "pam_sequence": "TGG",
    "mismatch_count": 1,
    "mismatch_positions": [2],
    "seed_mismatch_count": 0,
    "score": 90.0,
    "risk_category": "High-risk off-target"
  }
]
```

---

### 2.8 Export Results (CSV or JSON File Download)
**`GET /api/v1/jobs/{job_id}/export?format=csv`**
- **Authentication**: None (Public)
- **Query Parameters**: `format` (`csv` or `json`, default: `csv`)
- **Headers Returned**:
  `Content-Disposition: attachment; filename=geneweaver_alignment_{job_id}.csv`
  `Content-Type: text/csv`

