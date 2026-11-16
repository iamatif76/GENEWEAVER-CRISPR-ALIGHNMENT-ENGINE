
"""
Genome loading and validation module using Biopython SeqIO.

Supports multi-chromosome and multi-contig FASTA files.
"""

from pathlib import Path
from typing import Dict, List, Tuple

from Bio import SeqIO



class GenomeLoadError(Exception):
    """Raised when a FASTA file fails validation or parsing."""


# DNA IUPAC nucleotide symbols, including ambiguous bases.
VALID_NUCLEOTIDES = set("ACGTNRYKMSWBDHV")


def _validate_sequence(sequence: str, record_id: str) -> str:
    """
    Normalize and validate a DNA sequence.

    RNA U is converted to T. All other unsupported characters
    cause a GenomeLoadError.
    """
    sequence = str(sequence).upper().replace("U", "T")

    invalid_chars = set(sequence) - VALID_NUCLEOTIDES

    if invalid_chars:
        raise GenomeLoadError(
            f"Sequence '{record_id}' contains invalid nucleotide "
            f"characters: {sorted(invalid_chars)}"
        )

    if not sequence:
        raise GenomeLoadError(
            f"Sequence '{record_id}' is empty."
        )

    return sequence






def load_fasta(file_path: Path) -> Dict[str, str]:
    """
    Load a FASTA file into a chromosome/contig dictionary.
    """
    path = Path(file_path)

    if not path.is_file():
        raise GenomeLoadError(f"Genome file does not exist: {path}")

    chromosomes: Dict[str, str] = {}

    try:
        with path.open("r", encoding="utf-8", errors="strict") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id in chromosomes:
                    raise GenomeLoadError(
                        f"Duplicate FASTA sequence ID: '{record.id}'"
                    )

                chromosomes[record.id] = _validate_sequence(
                    record.seq, record.id
                )

    except GenomeLoadError:
        raise
    except (OSError, UnicodeError, ValueError) as exc:
        raise GenomeLoadError(
            f"Failed to parse FASTA file: {exc}"
        ) from exc

    if not chromosomes:
        raise GenomeLoadError(
            "FASTA file contains no valid sequence records."
        )

    return chromosomes



def inspect_fasta(
    file_path: Path,
) -> Tuple[int, int, List[Tuple[str, int]]]:
    """
    Inspect a FASTA file without retaining all sequences in memory.

    Returns:
        (total_bases, chromosome_count, list of (chromosome_id, length))
    """
    path = Path(file_path)

    if not path.is_file():
        raise GenomeLoadError(
            f"Genome file does not exist: {path}"
        )

    total_bases = 0
    chromosome_count = 0
    chrom_info: List[Tuple[str, int]] = []
    seen_ids = set()

    try:
        with path.open("r", encoding="utf-8", errors="strict") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id in seen_ids:
                    raise GenomeLoadError(
                        f"Duplicate FASTA sequence ID: '{record.id}'"
                    )

                sequence = _validate_sequence(record.seq, record.id)

                seen_ids.add(record.id)
                sequence_length = len(sequence)

                total_bases += sequence_length
                chromosome_count += 1
                chrom_info.append((record.id, sequence_length))

    except GenomeLoadError:
        raise
    except (OSError, UnicodeError, ValueError) as exc:
        raise GenomeLoadError(
            f"Failed to inspect FASTA file: {exc}"
        ) from exc

    if chromosome_count == 0:
        raise GenomeLoadError(
            "FASTA file contains no valid sequence records."
        )

    return total_bases, chromosome_count, chrom_info
