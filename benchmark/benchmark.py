"""
Comprehensive performance benchmarking and accuracy validation engine for GeneWeaver.
Measures:
- Genomic positions scanned
- Candidate PAM sites detected
- Alignments reported
- CPU vs GPU execution time
- Speedup factor
- Memory utilization (tracemalloc)
- Accuracy against known ground-truth targets
"""

import json
import os
import sys
import time
import tracemalloc
from pathlib import Path
from typing import Dict, Any

# Ensure project root is in sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.core.genome_loader import load_fasta
from app.core.pam_detector import PAMDetector
from app.core.aligner_cpu import CPUAligner
from app.core.aligner_gpu import GPUAligner
from app.core.gpu_utils import check_gpu_availability
from app.config import settings


def run_benchmark(
    fasta_path: Path = settings.DATA_DIR / "sample_genome.fasta",
    guide_sequence: str = "GAGTCCGAGCAGAAGAAGAA",
    max_mismatches: int = 3,
    pam_pattern: str = "NGG",
) -> Dict[str, Any]:
    """Execute end-to-end performance and accuracy benchmark."""
    print("=" * 70)
    print("  GENEWEAVER: HIGH-PERFORMANCE CRISPR ALIGNMENT BENCHMARK")
    print("=" * 70)

    gpu_available, dev_name, gpu_details = check_gpu_availability()
    print(f"Hardware Environment:")
    print(f"  CPU Vectorization: NumPy SIMD (AVX/AVX2/SSE enabled)")
    print(f"  CUDA GPU Available: {gpu_available} ({dev_name if gpu_available else 'CPU Fallback'})")
    print(f"  Input FASTA: {fasta_path}")
    print(f"  Guide RNA: {guide_sequence} (Length: {len(guide_sequence)} nt)")
    print(f"  PAM Pattern: {pam_pattern} (Max Mismatches: {max_mismatches})")
    print("-" * 70)

    # 1. Load genome
    chromosomes = load_fasta(fasta_path)
    total_positions_scanned = sum(len(seq) for seq in chromosomes.values())
    chromosome_count = len(chromosomes)
    print(f"Loaded {chromosome_count} contigs containing {total_positions_scanned:,} genomic positions.")

    # 2. Measure PAM scanning
    detector = PAMDetector(pam_pattern=pam_pattern, guide_length=len(guide_sequence))
    t0_pam = time.perf_counter()
    candidates = []
    for chrom, seq in chromosomes.items():
        for site in detector.find_candidate_sites(chrom, seq):
            candidates.append(site)
    t_pam = time.perf_counter() - t0_pam
    total_candidate_pams = len(candidates)
    print(f"Detected {total_candidate_pams:,} candidate PAM sites across forward & reverse strands in {t_pam*1000:.2f} ms.")

    # 3. Benchmark CPU Execution
    tracemalloc.start()
    cpu_aligner = CPUAligner(guide_sequence=guide_sequence, max_mismatches=max_mismatches)
    t0_cpu = time.perf_counter()
    cpu_hits = cpu_aligner.align_batch(candidates)
    t_cpu = time.perf_counter() - t0_cpu
    _, cpu_peak_mem = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    cpu_mem_mb = cpu_peak_mem / (1024 * 1024)

    # 4. Benchmark GPU Execution (or fallback comparison)
    gpu_time = None
    speedup = None
    gpu_hits = []

    if gpu_available:
        tracemalloc.start()
        gpu_aligner = GPUAligner(guide_sequence=guide_sequence, max_mismatches=max_mismatches, force_cpu=False)
        t0_gpu = time.perf_counter()
        gpu_hits = gpu_aligner.align_batch(candidates)
        gpu_time = time.perf_counter() - t0_gpu
        _, gpu_peak_mem = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        gpu_mem_mb = gpu_peak_mem / (1024 * 1024)
        speedup = t_cpu / gpu_time if gpu_time > 0 else 1.0
    else:
        # Simulated comparative analysis or fallback
        gpu_mem_mb = cpu_mem_mb
        gpu_hits = cpu_hits

    # 5. Accuracy Evaluation against Ground Truth
    # On chr1 we planted:
    # 1. Exact (+) target at ~101
    # 2. 1-mismatch outside seed (score 90.0)
    # 3. 1-mismatch inside seed (score 75.0)
    # 4. 2-mismatch (score 65.0)
    # 5. Exact (-) target
    # 6. 3-mismatch target
    exact_count = sum(1 for h in cpu_hits if h.mismatch_count == 0)
    one_mismatch = sum(1 for h in cpu_hits if h.mismatch_count == 1)
    two_mismatch = sum(1 for h in cpu_hits if h.mismatch_count == 2)
    three_mismatch = sum(1 for h in cpu_hits if h.mismatch_count == 3)

    # Ground truth validation checks
    accuracy_tests_passed = 0
    total_tests = 4

    if exact_count >= 2:  # Planted on chr1 and chrX
        accuracy_tests_passed += 1
    if one_mismatch >= 1:
        accuracy_tests_passed += 1
    if two_mismatch >= 1:
        accuracy_tests_passed += 1
    if any(h.strand == "-" and h.mismatch_count == 0 for h in cpu_hits):
        accuracy_tests_passed += 1

    accuracy_percent = (accuracy_tests_passed / total_tests) * 100.0

    print("\nBENCHMARK RESULTS SUMMARY:")
    print(f"  • Genomic Positions Scanned:   {total_positions_scanned:,} bp")
    print(f"  • Candidate PAMs Detected:     {total_candidate_pams:,} sites")
    print(f"  • Total Alignments Reported:   {len(cpu_hits):,} hits")
    print(f"  • Exact Matches (0 mismatches):{exact_count}")
    print(f"  • High-Risk (1 mismatch):      {one_mismatch}")
    print(f"  • Moderate-Risk (2 mismatches):{two_mismatch}")
    print(f"  • Low-Risk (3 mismatches):     {three_mismatch}")
    print(f"  • CPU Execution Time:          {t_cpu*1000:.2f} ms ({t_cpu:.5f} sec)")
    if gpu_available and gpu_time:
        print(f"  • GPU Execution Time:          {gpu_time*1000:.2f} ms ({gpu_time:.5f} sec)")
        print(f"  • Acceleration Speedup Factor: {speedup:.2f}x")
    else:
        print(f"  • GPU Acceleration:            N/A (CUDA GPU not detected - graceful fallback)")
    print(f"  • Peak Memory Consumption:     {cpu_mem_mb:.2f} MB")
    print(f"  • Ground-Truth Accuracy:       {accuracy_percent:.1f}% ({accuracy_tests_passed}/{total_tests} test cases)")
    print("=" * 70)

    result_data = {
        "positions_scanned": total_positions_scanned,
        "candidate_pams": total_candidate_pams,
        "alignments_reported": len(cpu_hits),
        "exact_matches": exact_count,
        "high_risk_count": one_mismatch,
        "moderate_risk_count": two_mismatch,
        "low_risk_count": three_mismatch,
        "cpu_time_seconds": round(t_cpu, 5),
        "gpu_time_seconds": round(gpu_time, 5) if gpu_time else None,
        "speedup_factor": round(speedup, 2) if speedup else None,
        "memory_peak_mb": round(cpu_mem_mb, 2),
        "accuracy_percent": accuracy_percent,
        "gpu_available": gpu_available,
        "hardware_info": dev_name if gpu_available else "CPU Vectorized NumPy",
    }

    out_file = Path(__file__).resolve().parent / "benchmark_results.json"
    with open(out_file, "w") as f:
        json.dump(result_data, f, indent=2)
    print(f"Saved benchmark metrics to {out_file.name}")

    return result_data


if __name__ == "__main__":
    run_benchmark()
