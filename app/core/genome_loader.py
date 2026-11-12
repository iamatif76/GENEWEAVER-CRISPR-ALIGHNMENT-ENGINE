"""
Genome loading and validation module using Biopython SeqIO.
Supports multi-chromosome and multi-contig FASTA files.
"""

from pathlib import Path
from typing import Dict, List, Tuple
from Bio import SeqIO


class GenomeLoadError(Exception):
    """Raised when a FASTA file fails validation or parsing."""
    pass


def load_fasta(file_path: Path) -> Dict[str, str]:
    """
    Parse a FASTA file and return a dictionary mapping chromosome/contig ID to uppercase sequence.

    Args:
        file_path: Path to the FASTA file.

    Returns:
        Dict[str, str]: Chromosome identifier -> DNA sequence string.

    Raises:
        GenomeLoadError: If file is empty, missing, or has invalid FASTA format.
    """
    path = Path(file_path)
    if not path.is_file():
        raise GenomeLoadError(f"Genome file does not exist: {file_path}")

    chromosomes: Dict[str, str] = {}
    valid_nucleotides = set("ACGTNacgtnRYKMSWBDHVrykmswbdhv")

    try:
        record_count = 0
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                record_count += 1
                seq_str = str(record.seq).upper().replace("U", "T")
                if not seq_str:
                    continue

                # Verify sequence contains biological nucleotide letters
                # (allowing IUPAC ambiguity like N, R, Y)
                sample = seq_str[:1000]
                if not set(sample).issubset(valid_nucleotides):
                    invalid_chars = set(sample) - valid_nucleotides
                    raise GenomeLoadError(
                        f"Chromosome '{record.id}' contains non-nucleotide characters: {invalid_chars}"
                    )

                chromosomes[record.id] = seq_str

        if record_count == 0 or len(chromosomes) == 0:
            raise GenomeLoadError("FASTA file contains no valid sequence records.")

        return chromosomes

    except Exception as e:
        if isinstance(e, GenomeLoadError):
            raise
        raise GenomeLoadError(f"Failed to parse FASTA file: {str(e)}") from e


def inspect_fasta(file_path: Path) -> Tuple[int, int, List[Tuple[str, int]]]:
    """
    Inspect a FASTA file to obtain summary statistics without keeping entire contents in memory.

    Returns:
        (total_bases, chromosome_count, list of (chrom_name, length))
    """
    chroms = load_fasta(file_path)
    chrom_info = [(cid, len(seq)) for cid, seq in chroms.items()]
    total_bases = sum(length for _, length in chrom_info)
    return total_bases, len(chroms), chrom_info
