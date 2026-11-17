
import pytest
from app.core.genome_loader import (
    GenomeLoadError,
    load_fasta,
    inspect_fasta,
)
from app.core.genome_loader import iter_fasta_chunks
def test_load_fasta_valid_sequences(tmp_path):
    fasta = tmp_path / "valid.fasta"
    fasta.write_text(
        ">chr1\nACGTACGT\n>chr2\nTTGGCCAA\n"
    )

    result = load_fasta(fasta)

    assert result == {
        "chr1": "ACGTACGT",
        "chr2": "TTGGCCAA",
    }


def test_load_fasta_converts_rna_to_dna(tmp_path):
    fasta = tmp_path / "rna.fasta"
    fasta.write_text(">chr1\nACGUACGU\n")

    result = load_fasta(fasta)

    assert result["chr1"] == "ACGTACGT"


def test_load_fasta_detects_invalid_character_after_1000_bases(tmp_path):
    fasta = tmp_path / "invalid.fasta"

    sequence = "A" * 1500 + "X"
    fasta.write_text(f">chr1\n{sequence}\n")

    with pytest.raises(GenomeLoadError, match="invalid nucleotide"):
        load_fasta(fasta)


def test_load_fasta_rejects_empty_file(tmp_path):
    fasta = tmp_path / "empty.fasta"
    fasta.write_text("")

    with pytest.raises(GenomeLoadError):
        load_fasta(fasta)


def test_load_fasta_rejects_empty_sequence_record(tmp_path):
    fasta = tmp_path / "empty_record.fasta"
    fasta.write_text(">chr1\n")

    with pytest.raises(GenomeLoadError, match="empty"):
        load_fasta(fasta)


def test_load_fasta_rejects_duplicate_ids(tmp_path):
    fasta = tmp_path / "duplicate.fasta"
    fasta.write_text(">chr1\nACGT\n>chr1\nTTAA\n")

    with pytest.raises(GenomeLoadError, match="Duplicate"):
        load_fasta(fasta)


def test_load_fasta_rejects_missing_file(tmp_path):
    fasta = tmp_path / "missing.fasta"

    with pytest.raises(GenomeLoadError, match="does not exist"):
        load_fasta(fasta)


def test_inspect_fasta_returns_correct_statistics(tmp_path):
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(
        ">chr1\nACGTACGT\n>chr2\nTTGGCC\n"
    )

    total_bases, chromosome_count, chrom_info = inspect_fasta(fasta)

    assert total_bases == 14
    assert chromosome_count == 2
    assert chrom_info == [
        ("chr1", 8),
        ("chr2", 6),
    ]


def test_inspect_fasta_detects_invalid_character(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    fasta.write_text(">chr1\nACGTX\n")

    with pytest.raises(GenomeLoadError):
        inspect_fasta(fasta)


def test_inspect_fasta_rejects_empty_file(tmp_path):
    fasta = tmp_path / "empty.fasta"
    fasta.write_text("")

    with pytest.raises(GenomeLoadError):
        inspect_fasta(fasta)

def test_iter_fasta_chunks_preserves_sequence_and_overlap(tmp_path):
    fasta = tmp_path / "chunks.fasta"
    sequence = "ACGTACGTACGTACGTACGTACGT"
    fasta.write_text(f">chr1\n{sequence}\n")

    chunks = list(
        iter_fasta_chunks(fasta, chunk_size=10, overlap=5)
    )

    assert chunks == [
        ("chr1", 0, sequence[:15]),
        ("chr1", 10, sequence[10:]),
    ]


def test_iter_fasta_chunks_handles_multiple_records(tmp_path):
    fasta = tmp_path / "multiple.fasta"
    fasta.write_text(
        ">chr1\nACGTACGT\n"
        ">chr2\nTTGGCCAA\n"
    )

    chunks = list(
        iter_fasta_chunks(fasta, chunk_size=4, overlap=2)
    )

    assert chunks == [
        ("chr1", 0, "ACGTAC"),
        ("chr1", 4, "ACGT"),
        ("chr2", 0, "TTGGCC"),
        ("chr2", 4, "CCAA"),
    ]


def test_iter_fasta_chunks_rejects_invalid_nucleotide(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    fasta.write_text(">chr1\nACGTXYZ\n")

    with pytest.raises(GenomeLoadError):
        list(iter_fasta_chunks(fasta, chunk_size=4, overlap=2))


def test_iter_fasta_chunks_rejects_invalid_chunk_size(tmp_path):
    fasta = tmp_path / "valid.fasta"
    fasta.write_text(">chr1\nACGTACGT\n")

    with pytest.raises(ValueError):
        list(iter_fasta_chunks(fasta, chunk_size=0))


def test_iter_fasta_chunks_rejects_duplicate_ids(tmp_path):
    fasta = tmp_path / "duplicate.fasta"
    fasta.write_text(
        ">chr1\nACGT\n"
        ">chr1\nTTGG\n"
    )

    with pytest.raises(GenomeLoadError):
        list(iter_fasta_chunks(fasta))
