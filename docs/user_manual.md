# GeneWeaver User Manual & Operating Guide

Welcome to **GeneWeaver: GPU-Accelerated CRISPR Alignment Engine**. This manual provides comprehensive, step-by-step guidance on setting up, configuring, and operating GeneWeaver through the Streamlit graphical interface and FastAPI REST endpoints.

---

## 1. Quick Start

### 1.1 Prerequisites
- Python 3.10, 3.11, 3.12, or 3.13
- Git
- (Optional) NVIDIA CUDA-capable GPU with drivers installed

### 1.2 Environment Setup
1. Clone the repository and navigate into the project directory:
   ```bash
   cd geneweaver
   ```

2. Create and activate a virtual environment:
   ```bash
   python -m venv venv
   # On Windows:
   .\venv\Scripts\activate
   # On Linux/macOS:
   source venv/bin/activate
   ```

3. Install required packages:
   ```bash
   pip install -r requirements.txt
   ```

4. Configure local environment variables:
   Copy `.env.example` to `.env` and set your preferred settings:
   ```bash
   cp .env.example .env
   ```

5. Seed the database with sample reference genomes and demonstration hits:
   ```bash
   python seed_data.py
   ```

---

## 2. Launching the Application

### 2.1 Starting the FastAPI Backend
Start the high-performance ASGI backend server:
```bash
uvicorn app.main:app --host 0.0.0.0 --port 8000 --reload
```
- Open `http://localhost:8000/docs` in your browser to view the interactive OpenAPI Swagger documentation.
- Health check: `http://localhost:8000/api/v1/health`

### 2.2 Starting the Streamlit Dashboard
In a second terminal window (with the virtual environment activated):
```bash
streamlit run frontend/app.py
```
Your browser will open automatically to `http://localhost:8501`.

---

## 3. Using the Streamlit Dashboard

### 3.1 Top Header & Hardware Mode Indicator
- Look at the top right of the dashboard:
  - **🟢 GPU Active: [Device Name]** indicates CuPy and CUDA GPU tensors are active.
  - **🔵 CPU Mode: NumPy SIMD** indicates optimized multi-core CPU SIMD mode is active.
- If the backend indicator is red, ensure your FastAPI uvicorn server is running on port 8000.

### 3.2 Running an Alignment Job
1. **Enter Guide Sequence**: Type or paste a 20-nt guide RNA or DNA sequence (e.g. `GAGTCCGAGCAGAAGAAGAA`).
   - Use the **Quick Presets** buttons (EMX1, VEGFA, FANCF) to load pre-verified human target guides.
   - RNA bases (`U`) are converted to DNA (`T`) automatically.
2. **Select PAM Pattern**: Standard default is `NGG` (SpCas9). Degenerate IUPAC patterns like `NAG`, `NRG`, or `TTTV` are fully supported.
3. **Choose Maximum Mismatches**: Use the slider to choose between 0 and 6 mismatches (default: 3).
4. **Select Reference Genome**:
   - Choose the preloaded "Synthetic Human Benchmark Genome (Sample)", OR
   - Upload your own multi-chromosome FASTA file (`.fasta`, `.fa`, `.fna`).
5. **Click "Run Alignment"**:
   - A real-time progress bar reflects the background worker pipeline.
   - The status updates dynamically as contigs are scanned.

### 3.3 Interpreting Results
- **Summary Metrics**: View total hits, exact matches, high-risk, moderate-risk, and low-risk counts at a glance.
- **Charts**:
  - **Mismatch Distribution**: Bar chart displaying hit counts grouped by mismatch level.
  - **Score Distribution**: Histogram showing the distribution of binding risk scores.
- **Results Table**:
  - Filter and sort candidate sites by chromosome, strand (`+` or `-`), genomic start/end positions, PAM sequence, mismatch positions, and score.
- **Exporting Results**:
  - Click **Download CSV** for tabular spreadsheet analysis (compatible with Excel, R, Python).
  - Click **Download JSON** for downstream bioinformatic pipelines.

---

## 4. REST API Usage Examples

### 4.1 Check API Health
```bash
curl -X GET http://localhost:8000/api/v1/health
```

### 4.2 Submit Alignment via curl
```bash
curl -X POST http://localhost:8000/api/v1/align \
  -H "X-API-Key: pcsk_5Tu4bi_GjEHqWmozkPJfjBqKVNTqZdp82MZcQUxwsqCeis2dieDj7DYg7SDKaHsFAbGN7s" \
  -F "guide_sequence=GAGTCCGAGCAGAAGAAGAA" \
  -F "genome_id=<GENOME_UUID>" \
  -F "max_mismatches=3" \
  -F "pam_pattern=NGG"
```

### 4.3 Check Job Status
```bash
curl -X GET http://localhost:8000/api/v1/jobs/<JOB_UUID>
```

### 4.4 Download CSV Results
```bash
curl -X GET "http://localhost:8000/api/v1/jobs/<JOB_UUID>/export?format=csv" -o results.csv
```

---

## 5. Performance Benchmarking
To run an automated benchmark measuring scanned positions, candidate PAM counts, throughput, and accuracy:
```bash
python benchmark/benchmark.py
```
Output includes memory metrics (MB) and execution times (milliseconds).

---

## 6. Troubleshooting

| Issue | Cause | Solution |
| :--- | :--- | :--- |
| **Backend API Offline** in UI | FastAPI server not started | Run `uvicorn app.main:app --reload` |
| **401 Unauthorized** error | Missing or invalid API key | Set `API_KEY` in `.env` and pass it via header |
| **Invalid guide sequence** | Character outside A, C, G, T, U | Ensure sequence contains valid nucleotide letters |
| **FASTA Parse Error** | File format corruption | Ensure file begins with `>` header lines |
| **GPU not utilized** | CUDA / CuPy not installed | GeneWeaver automatically falls back to CPU SIMD mode without crashing |
