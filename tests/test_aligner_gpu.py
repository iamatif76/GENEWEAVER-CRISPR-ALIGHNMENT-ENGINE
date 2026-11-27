
import pytest

from app.core.aligner_gpu import GPUAligner
from app.core.pam_detector import CandidateSite


GUIDE = "A" * 20


def make_site(position, protospacer=GUIDE):
    return CandidateSite(
        chromosome="chr1",
        start_position=position,
        end_position=position + 22,
        strand="+",
        protospacer=protospacer,
        pam="TGG",
    )


def make_aligner():
    return GPUAligner(
        guide_sequence=GUIDE,
        force_cpu=True,
    )


def test_align_genome_returns_empty_list_for_empty_genome():
    aligner = make_aligner()

    result = aligner.align_genome({}, batch_size=2)

    assert result == []


@pytest.mark.parametrize("batch_size", [0, -1])
def test_iter_align_genome_rejects_invalid_batch_size(batch_size):
    aligner = make_aligner()

    with pytest.raises(ValueError, match="batch_size"):
        list(
            aligner.iter_align_genome(
                {"chr1": "A" * 100},
                batch_size=batch_size,
            )
        )
def test_iter_align_genome_processes_full_and_partial_batches(monkeypatch):
    aligner = make_aligner()
    sites = [make_site(i * 23) for i in range(5)]
    batch_sizes = []

    def fake_find_candidate_sites(self, chromosome, sequence):
        yield from sites

    def fake_align_batch(batch, batch_size=100000):
        batch_sizes.append(len(batch))
        return list(batch)

    monkeypatch.setattr(
        "app.core.aligner_gpu.PAMDetector.find_candidate_sites",
        fake_find_candidate_sites
    )
    monkeypatch.setattr(aligner, "align_batch", fake_align_batch)

    results = list(
        aligner.iter_align_genome({"chr1": "A"}, batch_size=2)
    )

    assert results == sites
    assert batch_sizes == [2, 2, 1]


def test_iter_align_genome_processes_batches_lazily(monkeypatch):
    aligner = make_aligner()
    sites = [make_site(i * 23) for i in range(4)]
    processed_batches = []

    def fake_find_candidate_sites(self, chromosome, sequence):
        yield from sites

    def fake_align_batch(batch, batch_size=100000):
        processed_batches.append(len(batch))
        return list(batch)

    monkeypatch.setattr(
        "app.core.aligner_gpu.PAMDetector.find_candidate_sites",
        fake_find_candidate_sites
    )
    monkeypatch.setattr(aligner, "align_batch", fake_align_batch)

    results = aligner.iter_align_genome(
        {"chr1": "A"},
        batch_size=2
    )

    # Creating the iterator must not process any candidates.
    assert processed_batches == []

    # The first request processes only the first batch.
    assert next(results) == sites[0]
    assert processed_batches == [2]

    # The second result comes from the already processed batch.
    assert next(results) == sites[1]
    assert processed_batches == [2]

    # Requesting another result triggers the next batch.
    assert next(results) == sites[2]
    assert processed_batches == [2, 2]


def test_align_genome_sorts_hits_across_batches(monkeypatch):
    aligner = make_aligner()

    one_mismatch = "C" + GUIDE[1:]

    sites = [
        make_site(1, one_mismatch),
        make_site(2, one_mismatch),
        make_site(3, GUIDE),
        make_site(4, one_mismatch),
    ]

    def fake_find_candidate_sites(self, chromosome, sequence):
        yield from sites

    monkeypatch.setattr(
        "app.core.aligner_gpu.PAMDetector.find_candidate_sites",
        fake_find_candidate_sites,
    )

    results = aligner.align_genome({"chr1": "A" * 100}, batch_size=2)

    assert len(results) == 4
    assert results[0].mismatch_count == 0
    assert results[0].start_position == 3
    assert results[0].score >= results[-1].score
@pytest.mark.parametrize("batch_size", [0, -1])
def test_gpu_aligner_rejects_invalid_batch_size(batch_size):
    aligner = make_aligner()

    with pytest.raises(ValueError, match="batch_size"):
        aligner.align_batch([], batch_size=batch_size)
