
from app.core.aligner_python import PythonAligner
from app.core.aligner_cpu import CPUAligner
from app.core.pam_detector import CandidateSite


def test_python_aligner_exact_match():
    guide = "GAGTCCGAGCAGAAGAAGAA"

    aligner = PythonAligner(guide)

    site = CandidateSite(
        "chr1", 100, 123, "+", guide, "TGG"
    )

    hit = aligner.align_single(site)

    assert hit is not None
    assert hit.mismatch_count == 0
    assert hit.score == 100.0


def test_python_aligner_seed_mismatch():
    guide = "GAGTCCGAGCAGAAGAAGAA"
    target = guide[:-1] + "C"

    aligner = PythonAligner(guide)

    site = CandidateSite(
        "chr1", 100, 123, "+", target, "TGG"
    )

    hit = aligner.align_single(site)

    assert hit is not None
    assert hit.mismatch_count == 1
    assert hit.seed_mismatch_count == 1
    assert hit.score == 75.0


def test_python_aligner_filters_excess_mismatches():
    guide = "GAGTCCGAGCAGAAGAAGAA"

    aligner = PythonAligner(guide, max_mismatches=2)

    target = "TTTTCCGAGCAGAAGAAGAA"

    site = CandidateSite(
        "chr1", 100, 123, "+", target, "TGG"
    )

    assert aligner.align_single(site) is None


def test_python_and_numpy_results_match():
    guide = "GAGTCCGAGCAGAAGAAGAA"

    python_aligner = PythonAligner(guide)
    numpy_aligner = CPUAligner(guide)

    sites = [
        CandidateSite("chr1", 10, 33, "+", guide, "TGG"),
        CandidateSite("chr1", 50, 73, "+", "C" + guide[1:], "AGG"),
        CandidateSite("chr2", 100, 123, "-", guide[:18] + "CC", "GGG"),
    ]

    python_hits = python_aligner.align_batch(sites)
    numpy_hits = numpy_aligner.align_batch(sites)

    assert [h.to_dict() for h in python_hits] == [
        h.to_dict() for h in numpy_hits
    ]


def test_python_aligner_empty_batch():
    aligner = PythonAligner("GAGTCCGAGCAGAAGAAGAA")

    assert aligner.align_batch([]) == []
