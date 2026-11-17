
"""
Genome loading and validation module using Biopython SeqIO.

Supports multi-chromosome and multi-contig FASTA files.
"""

from pathlib import Path
from Bio import SeqIO
from typing import Dict, Iterator, List, Tuple

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

def iter_fasta_chunks(
    file_path: Path,
    chunk_size: int = 1_000_000,
    overlap: int = 22,
) -> Iterator[Tuple[str, int, str]]:
    """
    Read FASTA sequences incrementally and yield overlapping chunks.

    Yields:
        (chromosome_id, chunk_start, chunk_sequence)

        chromosome_id: FASTA record ID.
        chunk_start: Zero-based genomic starting position.
        chunk_sequence: Chunk including its right-side overlap.

    The overlap allows candidate sequences spanning chunk boundaries
    to be detected.
    """
    path = Path(file_path)

    if not path.is_file():
        raise GenomeLoadError(f"Genome file does not exist: {path}")

    if chunk_size <= 0:
        raise ValueError("chunk_size must be greater than zero.")

    if overlap < 0:
        raise ValueError("overlap cannot be negative.")

    seen_ids = set()

    try:
        with path.open("r", encoding="utf-8", errors="strict") as handle:
            record_id = None
            buffer = ""
            chunk_start = 0
            has_sequence = False

            for raw_line in handle:
                line = raw_line.strip()

                if not line:
                    continue

                # A new FASTA record starts.
                if line.startswith(">"):

                    # Finish the previous record.
                    if record_id is not None:
                        if not has_sequence:
                            raise GenomeLoadError(
                                f"Sequence '{record_id}' is empty."
                            )

                        # Do not emit an overlap-only remainder.
                        if chunk_start == 0 or len(buffer) > overlap:
                            yield record_id, chunk_start, buffer

                    header = line[1:].strip()

                    if not header:
                        raise GenomeLoadError(
                            "FASTA record has an empty header."
                        )

                    record_id = header.split()[0]

                    if record_id in seen_ids:
                        raise GenomeLoadError(
                            f"Duplicate FASTA sequence ID: '{record_id}'"
                        )

                    seen_ids.add(record_id)

                    buffer = ""
                    chunk_start = 0
                    has_sequence = False
                    continue

                if record_id is None:
                    raise GenomeLoadError(
                        "FASTA sequence found before the first header."
                    )

                # Remove whitespace and normalize DNA.
                sequence_line = "".join(line.split()).upper().replace("U", "T")

                invalid_chars = set(sequence_line) - VALID_NUCLEOTIDES

                if invalid_chars:
                    raise GenomeLoadError(
                        f"Sequence '{record_id}' contains invalid "
                        f"nucleotide characters: {sorted(invalid_chars)}"
                    )

                if not sequence_line:
                    continue

                has_sequence = True
                buffer += sequence_line

                # Emit complete chunks with right-side overlap.
                while len(buffer) >= chunk_size + overlap:
                    chunk = buffer[:chunk_size + overlap]

                    yield record_id, chunk_start, chunk

                    # Advance by the core size, not the full chunk size.
                    buffer = buffer[chunk_size:]
                    chunk_start += chunk_size

            # Finish the final FASTA record.
            if record_id is not None:
                if not has_sequence:
                    raise GenomeLoadError(
                        f"Sequence '{record_id}' is empty."
                    )

                if chunk_start == 0 or len(buffer) > overlap:
                    yield record_id, chunk_start, buffer

            if record_id is None:
                raise GenomeLoadError(
                    "FASTA file contains no valid sequence records."
                )

    except GenomeLoadError:
        raise
    except (OSError, UnicodeError, ValueError) as exc:
        raise GenomeLoadError(
            f"Failed to read FASTA file: {exc}"
        ) from exc



   
