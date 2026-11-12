"""
GeneWeaver: GPU-Accelerated CRISPR Alignment Engine
Streamlit Interactive Dashboard and Bioinformatics Workbench.
Calls the FastAPI backend REST API.
"""

import json
import os
import time
from typing import Optional
import requests
import streamlit as st
import pandas as pd
from dotenv import load_dotenv

# Load local environment configuration
load_dotenv()

# App Configuration
st.set_page_config(
    page_title="GeneWeaver — CRISPR Alignment Engine",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for modern bioinformatics styling
st.markdown("""
<style>
    .main-title {
        font-size: 2.2rem;
        font-weight: 700;
        background: linear-gradient(90deg, #1E88E5, #00E676);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        margin-bottom: 0.2rem;
    }
    .sub-title {
        font-size: 1.05rem;
        color: #78909C;
        margin-bottom: 1.5rem;
    }
    .metric-card {
        background-color: #1e2530;
        border-radius: 8px;
        padding: 16px;
        border-left: 4px solid #1E88E5;
        margin-bottom: 12px;
    }
    .risk-badge-exact {
        background-color: #00E676;
        color: #000;
        padding: 3px 8px;
        border-radius: 4px;
        font-weight: 600;
    }
    .risk-badge-high {
        background-color: #FF1744;
        color: #fff;
        padding: 3px 8px;
        border-radius: 4px;
        font-weight: 600;
    }
</style>
""", unsafe_allow_html=True)

# Sidebar Configuration
st.sidebar.markdown("### 🧬 GeneWeaver Control")
api_base_url = st.sidebar.text_input("Backend API Base URL", value="http://localhost:8000")
default_api_key = os.getenv("API_KEY", "")
api_key = st.sidebar.text_input("API Key (Authentication)", value=default_api_key, type="password")

headers = {}
if api_key:
    headers["X-API-Key"] = api_key.strip()


def check_api_health():
    try:
        r = requests.get(f"{api_base_url}/api/v1/health", timeout=3)
        if r.status_code == 200:
            return r.json()
        return None
    except Exception:
        return None


def fetch_genomes():
    try:
        r = requests.get(f"{api_base_url}/api/v1/genomes", headers=headers, timeout=5)
        if r.status_code == 200:
            return r.json()
        return []
    except Exception:
        return []


health_info = check_api_health()

# Top Hardware Mode Banner
col_head1, col_head2 = st.columns([3, 1])
with col_head1:
    st.markdown('<div class="main-title">GeneWeaver — GPU-Accelerated CRISPR Alignment Engine</div>', unsafe_allow_html=True)
    st.markdown('<div class="sub-title">Genome-wide on-target and off-target search, PAM detection, seed region weighting, and scoring</div>', unsafe_allow_html=True)

with col_head2:
    if health_info:
        if health_info.get("gpu_available"):
            st.success(f"🟢 GPU Active: {health_info.get('gpu_device_name', 'CUDA Device')}")
        else:
            st.info(f"🔵 CPU Mode: {health_info.get('cpu_mode', 'NumPy SIMD')}")
    else:
        st.error("🔴 Backend API Offline (Start FastAPI server)")

# Main Navigation Tabs
tab_align, tab_genomes, tab_benchmark, tab_about = st.tabs([
    "🎯 CRISPR Alignment",
    "📚 Reference Genomes",
    "⚡ Performance & GPU Benchmark",
    "📖 Documentation & Scoring"
])

# ==============================================================================
# TAB 1: CRISPR ALIGNMENT
# ==============================================================================
with tab_align:
    st.markdown("#### Configure Alignment Job")

    col_inp1, col_inp2 = st.columns([2, 1])

    with col_inp1:
        # Preset buttons for convenience
        st.caption("Quick Presets:")
        preset_cols = st.columns(3)
        if preset_cols[0].button("EMX1 (Human On-Target)"):
            st.session_state["guide_input"] = "GAGTCCGAGCAGAAGAAGAA"
        if preset_cols[1].button("VEGFA (Site 1 Target)"):
            st.session_state["guide_input"] = "GGGTGGGGGGAGTTTGCTCC"
        if preset_cols[2].button("FANCF Target"):
            st.session_state["guide_input"] = "GGAATCCCTTCTGCAGCACC"

        guide_val = st.session_state.get("guide_input", "GAGTCCGAGCAGAAGAAGAA")
        guide_seq = st.text_input(
            "CRISPR Guide RNA Sequence (approx 20 nt, DNA or RNA)",
            value=guide_val,
            help="Accepts DNA bases (A, C, G, T) or RNA (U is converted to T). Length 15-30 nt."
        )

    with col_inp2:
        pam_pattern = st.text_input(
            "PAM Pattern (IUPAC notation)",
            value="NGG",
            help="Default is NGG for SpCas9. Other examples: NAG, TTTV, NNNNGATT."
        )
        max_mismatches = st.slider(
            "Maximum Mismatches Allowed",
            min_value=0,
            max_value=6,
            value=3,
            help="Alignments with more mismatches than this threshold are filtered out."
        )

    st.markdown("##### Reference Genome Source")
    genome_choice = st.radio(
        "Choose Genome Source:",
        ["Select from uploaded genomes", "Upload new FASTA file"],
        horizontal=True
    )

    selected_genome_id = None
    uploaded_file = None

    available_genomes = fetch_genomes()

    if genome_choice == "Select from uploaded genomes":
        if available_genomes:
            genome_options = {
                f"{g['name']} ({g['chromosome_count']} contigs, {g['total_bases']:,} bp)": g["id"]
                for g in available_genomes
            }
            selected_label = st.selectbox("Available Reference Genomes", list(genome_options.keys()))
            selected_genome_id = genome_options[selected_label]
        else:
            st.warning("No reference genomes registered in the database yet. Please upload a FASTA file.")
    else:
        uploaded_file = st.file_uploader("Upload Genome FASTA (.fasta, .fa, .fna)", type=["fasta", "fa", "fna"])

    # Run Alignment Action
    col_btn, _ = st.columns([1, 4])
    with col_btn:
        run_btn = st.button("🚀 Run Alignment", type="primary", use_container_width=True)

    if run_btn:
        if not guide_seq.strip():
            st.error("Please enter a valid guide RNA sequence.")
        elif genome_choice == "Upload new FASTA file" and not uploaded_file:
            st.error("Please upload a FASTA file or select an existing genome.")
        elif genome_choice == "Select from uploaded genomes" and not selected_genome_id:
            st.error("Please select a reference genome.")
        else:
            with st.spinner("Submitting alignment job to backend..."):
                try:
                    payload = {
                        "guide_sequence": guide_seq.strip(),
                        "max_mismatches": str(max_mismatches),
                        "pam_pattern": pam_pattern.strip(),
                    }
                    if selected_genome_id:
                        payload["genome_id"] = selected_genome_id

                    files_dict = None
                    if uploaded_file:
                        files_dict = {"file": (uploaded_file.name, uploaded_file.getvalue(), "text/plain")}

                    res = requests.post(
                        f"{api_base_url}/api/v1/align",
                        data=payload,
                        files=files_dict,
                        headers=headers,
                        timeout=30
                    )

                    if res.status_code in [201, 202]:
                        job_info = res.json()
                        st.session_state["active_job_id"] = job_info["id"]
                        st.success(f"Alignment job submitted successfully! Job ID: `{job_info['id']}`")
                    else:
                        st.error(f"Error {res.status_code}: {res.text}")
                except Exception as e:
                    st.error(f"Failed to communicate with backend: {str(e)}")

    # Active Job Tracking & Results Display
    active_job_id = st.session_state.get("active_job_id")
    if active_job_id:
        st.divider()
        st.markdown(f"#### Active Job: `{active_job_id}`")

        status_placeholder = st.empty()
        progress_bar = st.progress(0)

        # Poll job status
        max_retries = 40
        poll_count = 0
        job_data = None

        with st.spinner("Processing CRISPR alignment across reference contigs..."):
            while poll_count < max_retries:
                try:
                    jr = requests.get(f"{api_base_url}/api/v1/jobs/{active_job_id}", headers=headers, timeout=5)
                    if jr.status_code == 200:
                        job_data = jr.json()
                        status = job_data.get("status")
                        pct = job_data.get("progress_percent", 0)
                        progress_bar.progress(min(pct, 100))

                        if status == "completed":
                            status_placeholder.success("🎉 Alignment Completed!")
                            break
                        elif status == "failed":
                            status_placeholder.error(f"❌ Job failed: {job_data.get('error_message')}")
                            break
                        else:
                            status_placeholder.info(f"⏳ Status: **{status.upper()}** ({pct}% complete)")
                    time.sleep(0.5)
                    poll_count += 1
                except Exception as e:
                    st.warning(f"Polling error: {str(e)}")
                    break

        if job_data and job_data.get("status") == "completed":
            summary = job_data.get("summary") or {}

            # Metric Summary Cards
            st.markdown("##### Alignment Summary")
            m_col1, m_col2, m_col3, m_col4, m_col5 = st.columns(5)
            m_col1.metric("Total Hits", summary.get("total_alignments", 0))
            m_col2.metric("Exact Matches (0 mm)", summary.get("exact_matches", 0))
            m_col3.metric("High-Risk (1 mm)", summary.get("high_risk_count", 0))
            m_col4.metric("Moderate-Risk (2 mm)", summary.get("moderate_risk_count", 0))
            m_col5.metric("Low-Risk (3 mm)", summary.get("low_risk_count", 0))

            # Fetch detailed results
            results_res = requests.get(
                f"{api_base_url}/api/v1/jobs/{active_job_id}/results",
                headers=headers,
                timeout=10
            )

            if results_res.status_code == 200:
                results_list = results_res.json()
                if results_list:
                    df = pd.DataFrame(results_list)

                    # Charts
                    chart_col1, chart_col2 = st.columns(2)
                    with chart_col1:
                        st.markdown("###### Mismatch Distribution")
                        mismatch_counts = df["mismatch_count"].value_counts().sort_index()
                        st.bar_chart(mismatch_counts)

                    with chart_col2:
                        st.markdown("###### Off-Target Risk Score Distribution")
                        st.bar_chart(df["score"].value_counts().sort_index())

                    # Results Table
                    st.markdown("##### Detailed Alignment Hits (Ranked by Score)")
                    display_cols = [
                        "chromosome", "start_position", "end_position", "strand",
                        "matched_sequence", "pam_sequence", "mismatch_count",
                        "mismatch_positions", "seed_mismatch_count", "score", "risk_category"
                    ]
                    st.dataframe(df[display_cols], use_container_width=True)

                    # Export Downloads
                    st.markdown("##### Export Results")
                    col_dl1, col_dl2, _ = st.columns([1, 1, 3])

                    csv_dl_url = f"{api_base_url}/api/v1/jobs/{active_job_id}/export?format=csv"
                    json_dl_url = f"{api_base_url}/api/v1/jobs/{active_job_id}/export?format=json"

                    csv_data = requests.get(csv_dl_url, headers=headers).text
                    json_data = requests.get(json_dl_url, headers=headers).text

                    col_dl1.download_button(
                        label="📥 Download CSV",
                        data=csv_data,
                        file_name=f"geneweaver_{active_job_id}.csv",
                        mime="text/csv",
                        use_container_width=True
                    )
                    col_dl2.download_button(
                        label="📥 Download JSON",
                        data=json_data,
                        file_name=f"geneweaver_{active_job_id}.json",
                        mime="application/json",
                        use_container_width=True
                    )
                else:
                    st.info("No alignment sites met the mismatch criteria.")

# ==============================================================================
# TAB 2: REFERENCE GENOMES
# ==============================================================================
with tab_genomes:
    st.markdown("#### Registered Reference Genomes")
    genomes = fetch_genomes()
    if genomes:
        gdf = pd.DataFrame(genomes)
        display_gdf = gdf[["id", "name", "filename", "total_bases", "chromosome_count", "uploaded_at"]]
        st.dataframe(display_gdf, use_container_width=True)
    else:
        st.info("No genomes uploaded yet.")

    st.divider()
    st.markdown("#### Upload a New Reference Genome")
    with st.form("genome_upload_form", clear_on_submit=True):
        new_genome_name = st.text_input("Genome Name (e.g., Human GRCh38 Subset, E. coli K-12)")
        new_genome_file = st.file_uploader("Upload FASTA file (.fasta, .fa, .fna)", type=["fasta", "fa", "fna"])
        upload_submit = st.form_submit_button("Upload and Index Genome")

        if upload_submit:
            if not new_genome_file:
                st.error("Please choose a FASTA file to upload.")
            else:
                with st.spinner("Uploading and indexing FASTA contigs..."):
                    try:
                        u_res = requests.post(
                            f"{api_base_url}/api/v1/genomes/upload",
                            data={"name": new_genome_name},
                            files={"file": (new_genome_file.name, new_genome_file.getvalue(), "text/plain")},
                            headers=headers,
                            timeout=60
                        )
                        if u_res.status_code == 201:
                            st.success(f"Successfully uploaded and indexed '{new_genome_file.name}'!")
                            st.rerun()
                        else:
                            st.error(f"Upload failed ({u_res.status_code}): {u_res.text}")
                    except Exception as e:
                        st.error(f"Error during upload: {str(e)}")

# ==============================================================================
# TAB 3: PERFORMANCE BENCHMARK & GPU DIAGNOSTICS
# ==============================================================================
with tab_benchmark:
    st.markdown("#### Hardware Execution & Comparative Benchmark")
    st.write(
        "GeneWeaver evaluates candidate off-target sites across the entire reference genome using "
        "vectorized SIMD arrays on CPU, or massively parallel tensor kernels on CUDA GPUs when present."
    )

    hw = check_api_health()
    if hw:
        b_col1, b_col2, b_col3 = st.columns(3)
        b_col1.metric("CPU Vectorization Engine", hw.get("cpu_mode", "NumPy SIMD"))
        b_col2.metric("CUDA GPU Acceleration", "Active" if hw.get("gpu_available") else "Fallback (CPU)")
        b_col3.metric("Device Name", hw.get("gpu_device_name") or "Host CPU")

    st.markdown("##### Run Live Pipeline Benchmark")
    if st.button("⚡ Run Benchmark (Sample Genome)", type="primary"):
        with st.spinner("Running high-throughput benchmark..."):
            try:
                from benchmark.benchmark import run_benchmark
                res = run_benchmark()
                st.success("Benchmark completed successfully!")

                bm_col1, bm_col2, bm_col3, bm_col4 = st.columns(4)
                bm_col1.metric("Positions Scanned", f"{res['positions_scanned']:,} bp")
                bm_col2.metric("Candidate PAMs", f"{res['candidate_pams']:,}")
                bm_col3.metric("CPU Execution Time", f"{res['cpu_time_seconds']*1000:.2f} ms")
                bm_col4.metric("Ground-Truth Accuracy", f"{res['accuracy_percent']:.1f}%")

                st.json(res)
            except Exception as e:
                st.error(f"Could not execute benchmark: {str(e)}")

# ==============================================================================
# TAB 4: DOCUMENTATION & SCORING
# ==============================================================================
with tab_about:
    st.markdown("### CRISPR Alignment & Off-Target Scoring Methodology")

    st.markdown("""
    #### 1. PAM Detection & Reverse-Strand Symmetry
    CRISPR-Cas9 recognizes a 20-nucleotide protospacer adjacent to a downstream Protospacer Adjacent Motif (PAM),
    typically **5'-NGG-3'** for *Streptococcus pyogenes* Cas9 (SpCas9).

    - **Forward Strand (+)**: Targets appear as `5'-[Protospacer: 20 nt][PAM: NGG]-3'`.
    - **Reverse Strand (-)**: On the reference plus strand, the reverse-complement PAM `5'-CCN-3'` precedes the reverse-complement protospacer:
      `5'-[CCN][Protospacer_rc]-3'`. GeneWeaver extracts both strands symmetrically.

    #### 2. PAM-Proximal Seed Region
    The **seed region** corresponds to the 8–12 nucleotides directly adjacent to the PAM (positions 11–20 from the 5' end for a 20 nt guide).
    Mismatches in the seed region drastically reduce Cas9 binding and cleavage activity, making them critical for off-target risk stratification.

    #### 3. Scoring Formula
    Off-target risk is quantified using a transparent linear penalty formula:
    $$\\text{Score} = \\max(0, \\min(100, 100 - (\\text{Total Mismatches} \\times 10) - (\\text{Seed Mismatches} \\times 15)))$$

    #### 4. Risk Classification
    | Mismatches | Classification | Description |
    | :--- | :--- | :--- |
    | **0** | **Exact match** | Primary on-target cleavage site (Score 100.0) |
    | **1** | **High-risk off-target** | High probability of unintended cleavage |
    | **2** | **Moderate-risk off-target** | Potential off-target, sensitive to seed mutations |
    | **3** | **Low-risk off-target** | Low probability of cleavage |
    | **>3** | **Minimal-risk off-target** | Filtered by default |
    """)
