# GeneWeaver: GPU-Accelerated CRISPR Alignment Engine 🧬⚡

GeneWeaver is a high-performance bioinformatics engine designed to identify on-target and off-target binding sites for CRISPR guide RNAs (gRNAs) across multi-chromosome reference genomes. It performs bi-directional sequence alignment, PAM motif detection, seed-region mutation weighting, transparent off-target risk scoring, SQLite tracking, and downloadable CSV/JSON exports.

The system is architected with a decoupled, asynchronous **FastAPI** REST backend, a responsive **Streamlit** user dashboard, and a hybrid execution engine featuring **NumPy SIMD CPU vectorization** with optional **CUDA GPU acceleration (CuPy)** and automatic zero-downtime CPU fallback.

---

## 🌟 Key Features

- **Strand-Symmetric Alignment**: Accurately scans both the forward (+) and reverse (-) DNA strands for candidate protospacers and PAM motifs.
- **Degenerate PAM Handling**: Supports standard SpCas9 (`NGG`) as well as IUPAC ambiguous motifs (e.g., `NAG`, `NRG`, `TTTV`).
- **Overlapping Motif Detection**: Uses regular expression positive lookahead assertions (`(?=...)`) to avoid missing overlapping PAM sites (e.g. `TGGGG`).
- **Seed Region Weighting**: Separately assesses the PAM-proximal seed region (positions 11–20) where mismatches heavily impair Cas9 binding.
- **Transparent Scoring Formula**:
  $$\text{Score} = \text{clamp}\Big(100 - (\text{Total Mismatches} \times 10) - (\text{Seed Mismatches} \times 15), \; 0, \; 100\Big)$$
- **Hardware Resilience**: High-throughput CPU SIMD vectorization via NumPy byte-matrices; seamlessly utilizes CUDA GPU tensors via CuPy when available, with automatic CPU fallback.
- **Asynchronous Task Queue**: Background task workers prevent API request blocking and stream live progress percentage (0% to 100%) to the frontend.
- **Relational Persistence**: Saves all uploaded genomes, jobs, and alignment hits in SQLite via SQLAlchemy.
- **Export Formats**: One-click downloadable CSV and JSON exports.
- **Security-First Architecture**: API keys are isolated exclusively in environment variables with constant-time cryptographic verification (`hmac.compare_digest`).

---

## 🏗️ Project Architecture

```
geneweaver/
├── app/
│   ├── main.py                     # FastAPI application bootstrap & lifespan
│   ├── config.py                   # Pydantic Settings & environment variables
│   ├── database.py                 # SQLAlchemy engine, session & init_db
│   ├── models.py                   # SQLite tables (genomes, jobs, alignments)
│   ├── schemas.py                  # Pydantic V2 request & response models
│   ├── routes/
│   │   ├── alignment.py            # /api/v1/align, /jobs/{id}, /results, /export
│   │   ├── genomes.py              # /api/v1/genomes, /genomes/upload
│   │   └── health.py               # /api/v1/health (CPU, GPU, DB status)
│   ├── services/
│   │   ├── alignment_service.py    # Job creation, querying & summaries
│   │   ├── genome_service.py       # FASTA upload, indexing & storage
│   │   ├── scoring_service.py      # Off-target scoring & risk tiering
│   │   └── export_service.py       # Pandas CSV & JSON exporters
│   ├── core/
│   │   ├── sequence_utils.py       # Sequence validation, revcomp, IUPAC
│   │   ├── genome_loader.py        # Biopython SeqIO multi-contig parser
│   │   ├── pam_detector.py         # Forward & reverse lookahead PAM detector
│   │   ├── aligner_cpu.py          # Vectorized NumPy SIMD matrix aligner
│   │   ├── aligner_gpu.py          # CuPy CUDA tensor aligner with CPU fallback
│   │   └── gpu_utils.py            # CUDA hardware capability probe
│   └── workers/
│       └── alignment_worker.py     # Asynchronous background job worker
├── frontend/
│   └── app.py                      # Interactive Streamlit dashboard & charts
├── tests/
│   ├── test_sequence_utils.py      # Unit tests for sequence logic
│   ├── test_pam_detector.py        # Unit tests for PAM scanning
│   ├── test_aligner_cpu.py         # Unit tests for alignment & scoring
│   ├── test_api.py                 # Integration tests for FastAPI endpoints
│   └── test_export.py              # Unit tests for CSV/JSON export
├── benchmark/
│   ├── benchmark.py                # Pipeline benchmark & accuracy validator
│   └── benchmark_results.json      # Output performance metrics
├── data/
│   └── sample_genome.fasta         # Synthetic 4-contig benchmark genome
├── docs/
│   ├── architecture.md             # Detailed system & data flow diagrams
│   ├── api_documentation.md        # Comprehensive REST API reference
│   ├── user_manual.md              # Step-by-step operating manual
│   ├── project_report.md           # Full technical report monograph
│   └── presentation_slides.md      # 12-slide presentation outline
├── seed_data.py                    # DB seed & demo job creation script
├── requirements.txt                # Python dependencies
├── Dockerfile                      # Production container image
├── docker-compose.yml              # Multi-service container orchestration
├── .env.example                    # Environment template
└── .gitignore                      # Git exclusion rules
```

---

## 🚀 Quick Start Guide

### 1. Clone & Set Up Virtual Environment

```bash
cd "gen weaver"
python -m venv venv

# Windows
.\venv\Scripts\activate

# Linux / macOS
source venv/bin/activate
```

### 2. Install Dependencies

```bash
pip install -r requirements.txt
```

### 3. Configure Environment Variables

Create your `.env` file from the example template:

```bash
cp .env.example .env
```

Your `.env` file should contain:
```ini
API_KEY= YOUR_API_KEY
DATABASE_URL=sqlite:///./geneweaver.db
UPLOAD_DIR=./uploads
HOST=0.0.0.0
PORT=8000
GPU_ENABLED=auto
LOG_LEVEL=INFO
```

> **Security Note**: Never commit your `.env` file or expose credentials in source code. `.gitignore` is configured to prevent committing `.env` and `*.db` files.

### 4. Initialize Database & Seed Sample Data

Run the seeding script to initialize the SQLite database, register the sample reference genome, and run a baseline demonstration job:

```bash
python seed_data.py
```

### 5. Launch the Application

#### Option A: Running Locally (FastAPI + Streamlit)

1. **Start the FastAPI Backend Server**:
   ```bash
   uvicorn app.main:app --host 0.0.0.0 --port 8000 --reload
   ```
   - API Documentation (Swagger UI): [http://localhost:8000/docs](http://localhost:8000/docs)
   - Health Status: [http://localhost:8000/api/v1/health](http://localhost:8000/api/v1/health)

2. **Start the Streamlit User Dashboard**:
   In a separate terminal window:
   ```bash
   streamlit run frontend/app.py
   ```
   - Dashboard UI: [http://localhost:8501](http://localhost:8501)

#### Option B: Running with Docker Compose

```bash
docker-compose up --build
```
- Backend will be accessible on `http://localhost:8000`
- Streamlit dashboard will be accessible on `http://localhost:8501`

---

## 🧪 Running Automated Tests

Run the automated Pytest test suite:

```bash
pytest -v
```

All 26 tests across sequence utilities, PAM detection, CPU vectorization, FastAPI routes, and export formatting pass with a 100% success rate.

---

## ⚡ Running the Performance Benchmark

To measure scanning throughput, candidate PAM detection speed, memory consumption, and ground-truth accuracy:

```bash
python benchmark/benchmark.py
```

### Benchmark Results (Sample Genome):
- **Genomic Contigs Scanned**: 4 contigs (5,617 base pairs)
- **Candidate PAM Sites Scanned**: 701 candidate loci
- **CPU Alignment Execution Time**: **1.71 milliseconds**
- **Throughput**: **>410,000 candidate sites / second**
- **Peak Memory Usage**: 0.12 MB
- **Ground-Truth Accuracy**: **100.0%** across planted controls

---

## 🔌 API Reference & Usage Examples

### 1. Health & Hardware Check
```bash
curl -X GET http://localhost:8000/api/v1/health
```

### 2. List Available Genomes
```bash
curl -X GET http://localhost:8000/api/v1/genomes
```

### 3. Upload a Reference Genome
```bash
curl -X POST http://localhost:8000/api/v1/genomes/upload \
  -H "X-API-Key: YOUR_API_KEY" \
  -F "file=@data/sample_genome.fasta" \
  -F "name=Human GRCh38 Subset"
```

### 4. Submit an Alignment Job
```bash
curl -X POST http://localhost:8000/api/v1/align \
  -H "X-API-Key: YOUR_API_KEY" \
  -F "guide_sequence=GAGTCCGAGCAGAAGAAGAA" \
  -F "genome_id=<GENOME_ID>" \
  -F "max_mismatches=3" \
  -F "pam_pattern=NGG"
```

### 5. Check Job Status & Summary
```bash
curl -X GET http://localhost:8000/api/v1/jobs/<JOB_ID>
```

### 6. Retrieve Ranked Alignment Results
```bash
curl -X GET "http://localhost:8000/api/v1/jobs/<JOB_ID>/results?limit=50"
```

### 7. Download Results in CSV Format
```bash
curl -X GET "http://localhost:8000/api/v1/jobs/<JOB_ID>/export?format=csv" -o alignment_results.csv
```

---

## 📊 Streamlit Dashboard Walkthrough

The Streamlit UI provides:
1. **CRISPR Alignment Tab**:
   - One-click presets for verified human guides: **EMX1**, **VEGFA**, and **FANCF**.
   - Input fields for custom gRNA and PAM pattern.
   - Genome selector (existing indexed genomes or direct FASTA upload).
   - Maximum mismatch slider (0 to 6).
   - Real-time animated job progress bar.
   - Summary metric cards for exact matches, high-risk, moderate-risk, and low-risk hits.
   - Dynamic interactive charts (mismatch distribution and score distribution).
   - Interactive data table with 1-based genomic coordinates.
   - One-click CSV and JSON download buttons.
2. **Reference Genomes Tab**: View contig counts, total bases, and upload new genomes.
3. **Performance Benchmark Tab**: Run real-time hardware benchmarks measuring SIMD speed and accuracy.
4. **Documentation & Scoring Tab**: Complete interactive explanation of Cas9 strand symmetry, seed weighting, and scoring formulas.

---

## 🛡️ Security Architecture

- **No Secrets in Code**: Environment credentials are strictly loaded via `python-dotenv`.
- **Constant-Time Verification**: API key authentication employs `hmac.compare_digest` to prevent timing attacks.
- **Redaction**: Secrets are never output in logs, exception traces, or HTTP responses.
- **Git Protection**: Sensitive patterns (`.env`, `*.db`, `uploads/`, `__pycache__/`) are enforced in `.gitignore`.

---

## 📚 Technical Documentation

- [Architecture & Data Flow Diagram](docs/architecture.md)
- [Comprehensive API Documentation](docs/api_documentation.md)
- [User Manual & Operating Guide](docs/user_manual.md)
- [Complete Project Report Monograph](docs/project_report.md)
- [12-Slide Master Presentation Deck](docs/presentation_slides.md)

---

## 📜 License
This project is licensed under the MIT License.
## Team Development and Contribution

GeneWeaver is developed collaboratively by a four-member team.
Each contributor works on dedicated feature branches for dashboard
development, API backend, CRISPR alignment logic, and testing.

The develop branch is used for integrating and validating features
before merging stable updates into the main branch.
## Project Maintenance and Quality Assurance

GeneWeaver uses automated testing and performance benchmarking
to maintain the reliability of its CRISPR alignment engine.
The development team continuously improves the dashboard,
backend services, alignment algorithms, and documentation.
All updates are reviewed before integration into the main branch.
## Continuous Improvement and Development Standards

GeneWeaver follows a structured development process to improve
code quality, maintainability, and system reliability. The team
works on modular components, performs functional testing, and
reviews updates before integrating them into the develop branch.
Performance benchmarks help evaluate alignment engine improvements
across supported execution environments.
## Documentation and Code Review Guidelines

Clear documentation helps developers understand the project structure,
features, and implementation process. Code changes should follow
consistent naming conventions, maintain modularity, and be tested
before submission. Team members should review updates to identify
potential issues and ensure compatibility with existing components.
