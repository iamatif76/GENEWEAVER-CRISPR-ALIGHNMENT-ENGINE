# GeneWeaver Architecture Documentation

## 1. System Overview

GeneWeaver is a high-performance bioinformatics engine designed to identify and score CRISPR on-target and off-target binding sites across multi-chromosome reference genomes. It combines NumPy SIMD vectorization on CPU with optional CUDA tensor acceleration (CuPy) and provides a versioned REST API (FastAPI) and an interactive web dashboard (Streamlit).

```
                      +-----------------------------+
                      |      Streamlit Dashboard    |
                      |   (Web UI / Visualization)  |
                      +--------------+--------------+
                                     |
                                     | REST API (HTTP / JSON / CSV)
                                     v
                      +-----------------------------+
                      |       FastAPI Gateway       |
                      |   (Auth, Routing, Validation)|
                      +--------------+--------------+
                                     |
                +--------------------+--------------------+
                |                                         |
                v                                         v
   +-------------------------+               +-------------------------+
   |   SQLite Database       |               | Background Worker Task  |
   | (Jobs, Genomes, Hits)   | <-------------+   (Asynchronous Queue)  |
   +-------------------------+               +------------+------------+
                                                          |
                                                          v
                                             +-------------------------+
                                             |     PAM Detector        |
                                             | (Forward & Reverse +/-) |
                                             +------------+------------+
                                                          |
                                                          v
                                             +-------------------------+
                                             | Alignment Engine Router |
                                             +-----+-------------+-----+
                                                   |             |
                                 CUDA Available    |             | CPU Fallback
                                                   v             v
                                        +-------------+   +-------------+
                                        | GPU Aligner |   | CPU Aligner |
                                        |   (CuPy)    |   |   (NumPy)   |
                                        +-------------+   +-------------+
```

---

## 2. Component Architecture

### 2.1 Backend REST API (`app/`)
- **`app/main.py`**: Application bootstrap, CORS configuration, SQLite database migration lifecycle, sample genome registration, and global exception handlers.
- **`app/config.py`**: Pydantic BaseSettings loading environment configurations from `.env` via `python-dotenv`.
- **`app/database.py`**: SQLAlchemy engine configuration with SQLite thread pooling for thread-safe concurrent background tasks.
- **`app/models.py`**: Relational models for `genomes`, `jobs`, and `alignments`.
- **`app/schemas.py`**: Pydantic V2 models for input sanitization, response formatting, and validation.

### 2.2 Core Bioinformatics Pipeline (`app/core/`)
- **`sequence_utils.py`**:
  - DNA/RNA cleaning and normalization (U $\to$ T).
  - Reverse complement translation supporting both standard and IUPAC ambiguous bases.
  - 1-based mismatch indexing and PAM-proximal seed window mutation calculation.
  - IUPAC nucleotide ambiguity translation to standard regular expressions.
- **`genome_loader.py`**:
  - Biopython `SeqIO.parse` wrapper for multi-chromosome FASTA loading.
  - Nucleotide character validation and metric extraction (total bases, contig count).
- **`pam_detector.py`**:
  - Bi-directional PAM scanning.
  - Overlapping motif extraction via regex positive lookahead assertions: `(?=([ACGT]GG))`.
  - Minus-strand locus translation: reference `5'-CCN-3'` maps to targeted strand `5'-NGG-3'`.
- **`aligner_cpu.py`**:
  - Vectorized NumPy alignment using ASCII byte matrices.
  - Broadcast comparison `matrix != guide_bytes` evaluated simultaneously across tens of thousands of sites.
- **`aligner_gpu.py`**:
  - CuPy-accelerated tensor comparison on CUDA devices.
  - Seamless zero-downtime CPU fallback if CUDA devices or CuPy are unavailable.
- **`gpu_utils.py`**:
  - Probes CUDA driver, runtime device count, and device properties without throwing unhandled exceptions.

### 2.3 Services & Workers
- **`alignment_worker.py`**: Runs long-running genome-wide scans in background threads. Updates job completion percentages (0% $\to$ 100%) in SQLite so the frontend can poll asynchronously.
- **`scoring_service.py`**: Computes off-target penalty scores and risk classifications.
- **`export_service.py`**: Formats alignment hits into downloadable CSV and JSON formats using Pandas.

---

## 3. Data Flow Diagram

```
User / UI                     FastAPI                   Background Worker             SQLite Database
    |                            |                              |                            |
    |--- POST /api/v1/align ---->|                              |                            |
    |    (Guide, Genome, PAM)    |                              |                            |
    |                            |--- Insert Job (queued) ----->|--------------------------->|
    |<-- HTTP 202 (Job ID) ------|                              |                            |
    |                            |--- Spawn Worker ------------>|                            |
    |                            |                              |--- Status: running (15%) ->|
    |                            |                              |--- Load FASTA contigs ---->|
    |                            |                              |--- Scan PAM Sites (50%) -->|
    |                            |                              |--- Matrix Alignment (85%)->|
    |                            |                              |--- Bulk Insert Hits ------>|
    |                            |                              |--- Status: completed ----->|
    |--- GET /api/v1/jobs/{id} ->|                              |                            |
    |<-- Status & Summary -------|--- Query Job Summary ------------------------------------>|
    |                            |                              |                            |
    |--- GET .../results ------->|                              |                            |
    |<-- Ranked Hit List --------|--- Query Alignments ------------------------------------->|
```

---

## 4. Algorithmic Principles & Formulae

### 4.1 Target Orientation
For SpCas9:
- Forward (+) Strand:
  $$\text{Target Locus} = [\,\text{Protospacer (20 nt)}\,]\;[\,\text{PAM (3 nt)}\,]$$
- Reverse (-) Strand:
  $$\text{Reference Sequence} = [\,\text{PAM}^{\text{RC}}\,(3\text{ nt})\,]\;[\,\text{Protospacer}^{\text{RC}}\,(20\text{ nt})\,]$$

### 4.2 Seed Region Weighting
The PAM-proximal seed region represents the 10 nucleotides at the 3' end of the protospacer (positions 11 through 20). Mismatches here destabilize R-loop formation and Cas9 endonuclease activation significantly more than distal 5' mutations.

### 4.3 Off-Target Scoring Formula
$$\text{Score} = \max\Big(0,\; \min\big(100,\; 100 - (\text{Total Mismatches} \times 10) - (\text{Seed Mismatches} \times 15)\big)\Big)$$

- **Exact Match** (0 mismatches): Score = 100.0
- **1 Mismatch (Outside Seed)**: Score = 90.0
- **1 Mismatch (Inside Seed)**: Score = 75.0
- **2 Mismatches (1 Seed, 1 Distal)**: Score = 65.0
- **2 Mismatches (2 Seed)**: Score = 50.0

---

## 5. Security & Credential Isolation

1. **Zero Secret Hardcoding**: Secrets are never stored in code, Git history, or container images.
2. **Environment Variable Configuration**: Loaded exclusively via `python-dotenv` into typed Pydantic models.
3. **Constant-Time Comparison**: `hmac.compare_digest` prevents timing attacks on protected endpoints.
4. **Header Obfuscation**: The API never logs or echoes client keys in response bodies or exception messages.
