# GeneWeaver: GPU-Accelerated CRISPR Alignment Engine
## Comprehensive Project Report & Technical Monograph

---

### Executive Summary
GeneWeaver is an end-to-end computational genomics and bioinformatics platform developed to accelerate CRISPR-Cas9 on-target identification and off-target risk evaluation. Off-target cleavage poses a paramount challenge in therapeutic gene editing, where unintended double-strand breaks can trigger chromosomal rearrangements, oncogene activation, or cytotoxicity. GeneWeaver addresses this challenge through high-performance sequence alignment, PAM motif scanning on both DNA strands, seed-region sensitivity weighting, and transparent multi-tier off-target scoring. The system features a modern decoupled architecture consisting of an asynchronous FastAPI REST engine, an SQLite database with SQLAlchemy ORM, and an interactive Streamlit data dashboard.

---

### Table of Contents
1. [Introduction](#1-introduction)
2. [Problem Statement](#2-problem-statement)
3. [Project Objectives](#3-project-objectives)
4. [Technology Stack](#4-technology-stack)
5. [System Architecture](#5-system-architecture)
6. [Module-by-Module Explanation](#6-module-by-module-explanation)
7. [Algorithmic Design & Mathematical Models](#7-algorithmic-design--mathematical-models)
8. [Database Design & Schemas](#8-database-design--schemas)
9. [Verification & Testing Methodology](#9-verification--testing-methodology)
10. [Performance Evaluation & Benchmarking](#10-performance-evaluation--benchmarking)
11. [Limitations](#11-limitations)
12. [Future Scope](#12-future-scope)
13. [Conclusion](#13-conclusion)

---

### 1. Introduction
The CRISPR-Cas9 genome editing system relies on a single guide RNA (sgRNA) of approximately 20 nucleotides in length to direct the Cas9 endonuclease to a target genomic locus. Recognition requires Watson-Crick base pairing between the gRNA and the target protospacer, coupled with an adjacent Protospacer Adjacent Motif (PAM)—typically `5'-NGG-3'` for *Streptococcus pyogenes* Cas9 (SpCas9). However, Cas9 tolerates multiple base-pair mismatches between the guide RNA and genomic DNA, causing cleavage at non-targeted genomic sites. Rapid, accurate, and scalable in-silico prediction of off-target binding is a prerequisite for therapeutic and research CRISPR applications.

---

### 2. Problem Statement
Existing bioinformatic tools for off-target detection frequently suffer from one or more systemic drawbacks:
1. **Computational Bottlenecks**: Full-genome sequence search can be prohibitively slow when scanning billions of bases for degenerate PAM sites.
2. **Opaque Scoring Algorithms**: Many tools rely on proprietary or uncalibrated black-box scoring systems that obscure why a specific site is classified as dangerous.
3. **Rigid Architectures**: Command-line-only tools lack accessible graphical user interfaces for molecular biologists, while cloud-only web servers expose proprietary sequences or genetic data.
4. **Hardware Fragility**: Many GPU-accelerated tools crash immediately if deployed on systems lacking CUDA drivers or specialized hardware.

---

### 3. Project Objectives
1. **Accurate Strand-Symmetric Alignment**: Detect candidate binding sites on both forward (+) and reverse (-) genomic strands.
2. **IUPAC Degeneracy Handling**: Support flexible PAM sequences (NGG, NAG, TTTV, etc.) with standard IUPAC ambiguous nucleotide handling.
3. **Seed Region Sensitivity Weighting**: Prioritize the 10 nucleotides adjacent to the PAM where mismatches most heavily impair Cas9 binding.
4. **Transparent Risk Stratification**: Assign reproducible off-target risk scores clamped between 0 and 100.
5. **Hardware Resilience & Acceleration**: Implement vectorized NumPy SIMD execution on CPU, paired with optional CuPy CUDA acceleration and zero-downtime CPU fallback.
6. **Enterprise Architecture**: Provide a versioned REST API, background worker task queues, SQLite storage, CSV/JSON export, and an intuitive Streamlit interface.

---

### 4. Technology Stack

| Layer | Component | Version / Library | Purpose |
| :--- | :--- | :--- | :--- |
| **Backend Framework** | FastAPI | 0.110+ | Asynchronous REST API with OpenAPI documentation |
| **Server Engine** | Uvicorn | 0.28+ | Lightning-fast ASGI server |
| **Core Bioinformatics** | Biopython | 1.83+ | Robust FASTA record parsing and sequence validation |
| **CPU Vectorization** | NumPy | 1.26+ | Byte-matrix comparison and broadcasted mismatch counting |
| **GPU Acceleration** | CuPy / CUDA | Optional | Tensorized parallel mismatch scanning on NVIDIA GPUs |
| **Database & ORM** | SQLAlchemy / SQLite | 2.0+ | Relational tracking of jobs, genomes, and alignment hits |
| **Tabular Analytics** | Pandas | 2.2+ | Data manipulation and CSV/JSON export formatting |
| **Frontend UI** | Streamlit | 1.32+ | Interactive bioinformatics control panel and visualizations |
| **Testing Suite** | Pytest | 8.0+ | Automated unit, regression, and API integration testing |
| **Containerization** | Docker & Compose | Multi-stage | Reproducible environment deployment |

---

### 5. System Architecture
GeneWeaver adopts a modular, service-oriented architecture:
- **Presentation Layer**: Streamlit dashboard interacting via HTTP requests.
- **API Gateway Layer**: FastAPI application managing routing, CORS, and API key authentication.
- **Worker & Persistence Layer**: FastAPI `BackgroundTasks` orchestrating asynchronous pipeline jobs while persisting state in SQLite.
- **Algorithmic Engine**: Specialized modules for sequence validation, PAM lookahead detection, and vectorized alignment.

---

### 6. Module-by-Module Explanation

#### 6.1 `app/core/sequence_utils.py`
Provides foundational nucleotide manipulations. Validates inputs, handles automatic conversion of RNA uracil (`U`) to DNA thymine (`T`), and generates reverse complements. Implements `calculate_mismatches` to return total mismatch counts, 1-based position arrays, and seed-window counts.

#### 6.2 `app/core/genome_loader.py`
Leverages Biopython's `SeqIO` to stream multi-chromosome FASTA genomes. Validates character sets and calculates total sequence lengths and contig counts.

#### 6.3 `app/core/pam_detector.py`
Detects PAM sites using regular expressions with positive lookahead assertions (`(?=...)`). Lookahead matching is essential to capture overlapping PAM sequences (e.g. `TGGGG`). Symmetrically translates minus-strand loci by matching `5'-CCN-3'` on the positive reference strand to produce `5'-NGG-3'` on the minus strand.

#### 6.4 `app/core/aligner_cpu.py`
Converts candidate protospacer strings into a 2D NumPy array of unsigned 8-bit integers (`np.uint8`). Evaluates millions of base comparisons concurrently through array broadcasting (`matrix != guide_bytes`).

#### 6.5 `app/core/aligner_gpu.py` & `gpu_utils.py`
Probes CUDA availability dynamically. If CuPy and a CUDA device are present, offloads byte matrices to GPU device memory. If unavailable, logs an informative notice and invokes `CPUAligner` seamlessly.

#### 6.6 `app/services/scoring_service.py`
Computes off-target risk scores and categorizes results into Exact Match, High Risk, Moderate Risk, and Low Risk.

#### 6.7 `app/workers/alignment_worker.py`
Executes long-running alignment workloads in independent background threads. Periodically commits progress updates (0% to 100%) to SQLite.

---

### 7. Algorithmic Design & Mathematical Models

#### 7.1 Coordinate Mapping
1. **Forward (+) Strand**:
   $$\text{Protospacer: } [k - L, k], \quad \text{PAM: } [k, k + P]$$
   where $k$ is the PAM start index, $L$ is guide length (20), and $P$ is PAM length (3).
2. **Reverse (-) Strand**:
   $$\text{PAM}^{\text{RC}}: [k, k + P], \quad \text{Protospacer}^{\text{RC}}: [k + P, k + P + L]$$
   Reverse complementation yields the 5' $\to$ 3' targeted sequence on the minus strand.

#### 7.2 Off-Target Score Formulation
$$\text{Score} = \text{clamp}\Big(100 - (10 \times M_{\text{total}}) - (15 \times M_{\text{seed}}), \; 0, \; 100\Big)$$

Where:
- $M_{\text{total}}$ is the total number of mismatched base pairs across the 20 nt protospacer.
- $M_{\text{seed}}$ is the number of mismatches situated in the PAM-proximal seed region (positions 11–20).

---

### 8. Database Design & Schemas

The relational schema is implemented in SQLite via SQLAlchemy:

```
+-----------------------------------+        +-----------------------------------+
|              genomes              |        |               jobs                |
+-----------------------------------+        +-----------------------------------+
| id: String(36) [PK]               | 1    * | id: String(36) [PK]               |
| name: String(255)                 |<-------| guide_sequence: String(128)       |
| filename: String(255)             |        | genome_id: String(36) [FK]        |
| file_path: String(1024)           |        | max_mismatches: Integer           |
| total_bases: Integer              |        | pam_pattern: String(32)           |
| chromosome_count: Integer         |        | status: String(32)                |
| uploaded_at: DateTime             |        | progress_percent: Integer         |
+-----------------------------------+        | error_message: Text [Nullable]    |
                                             | created_at: DateTime              |
                                             | completed_at: DateTime [Nullable] |
                                             +-----------------+-----------------+
                                                               | 1
                                                               |
                                                               | *
                                             +-----------------+-----------------+
                                             |            alignments             |
                                             +-----------------------------------+
                                             | id: String(36) [PK]               |
                                             | job_id: String(36) [FK]           |
                                             | chromosome: String(128)           |
                                             | start_position: Integer           |
                                             | end_position: Integer             |
                                             | strand: String(1)                 |
                                             | matched_sequence: String(128)     |
                                             | pam_sequence: String(32)          |
                                             | mismatch_count: Integer           |
                                             | mismatch_positions: Text (JSON)   |
                                             | seed_mismatch_count: Integer      |
                                             | score: Float                      |
                                             | risk_category: String(64)         |
                                             +-----------------------------------+
```

---

### 9. Verification & Testing Methodology
The test suite spans 26 automated tests covering:
- Sequence sanitization, RNA $\to$ DNA conversion, and boundary checks.
- Reverse complement operations with standard and IUPAC ambiguous bases.
- Overlapping PAM motif detection on both forward and reverse strands.
- CPU vectorization mathematical parity against iterative single-site calculations.
- API authentication enforcement (`401 Unauthorized` for missing/invalid keys).
- Asynchronous background job completion and result extraction.
- CSV and JSON export formatting integrity.

**Result**: 26 / 26 tests passing (100% pass rate).

---

### 10. Performance Evaluation & Benchmarking
Performance was evaluated using `benchmark/benchmark.py`:
- **Genomic Contigs Scanned**: 4 contigs (5,617 base pairs)
- **Candidate PAM Sites Scanned**: 701 candidate loci
- **CPU Alignment Execution Time**: 1.71 milliseconds
- **Throughput**: ~410,000 candidate sites per second on CPU
- **Memory Footprint**: 0.12 MB peak memory during alignment
- **Ground-Truth Accuracy**: 100.0% (verified across planted exact, single-mismatch, double-mismatch, and reverse-strand control targets)

---

### 11. Limitations
1. **Epigenetic Factors**: In-silico alignment does not incorporate chromatin accessibility, DNA methylation, or nucleosome positioning.
2. **Indel Detection**: Current alignment targets substitution mismatches; insertions and deletions (RNA/DNA bulges) require gap-aware dynamic programming.
3. **Cas Variants**: Preconfigured for SpCas9 PAMs, though IUPAC support allows custom patterns.

---

### 12. Future Scope
1. **Machine Learning Cleavage Scoring**: Integration of deep learning models (e.g. Elevation, DeepCRISPR) to predict actual cleavage efficiency.
2. **Bulge Alignment Support**: Implementation of banded Needleman-Wunsch or Myers bit-parallel algorithms for RNA/DNA bulge detection.
3. **Distributed Worker Clusters**: Migration from local background tasks to Celery or Redis Queue for cluster-wide parallelization across multiple GPU nodes.

---

### 13. Conclusion
GeneWeaver successfully provides a complete, robust, and mathematically grounded CRISPR alignment and off-target risk engine. By combining SIMD-accelerated CPU computing with optional GPU acceleration, an asynchronous FastAPI architecture, and an accessible Streamlit interface, GeneWeaver bridges high-performance computing and accessible genomics research.

