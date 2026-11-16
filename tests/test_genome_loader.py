
import pytest
from app.core.genome_loader import (
    GenomeLoadError,
    load_fasta,
    inspect_fasta,
)
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
