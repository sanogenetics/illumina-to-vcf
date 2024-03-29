from typing import Union

from fsspec.core import OpenFile
from pyfaidx import Fasta


class ReferenceGenome:
    reference: Union[str, OpenFile]
    reference_index: Union[str, OpenFile]
    reference_fasta: Fasta

    def __init__(self, reference: Union[str, OpenFile], reference_index: Union[str, OpenFile]) -> None:
        self.reference = reference
        self.reference_index = reference_index
        self.reference_fasta = Fasta(
            reference,
            reference_index,
            read_ahead=1024 * 16,  # 16kb read ahead buffer
        )

    def get_reference_bases(self, chrom: str, start: int, end: int) -> str:
        """
        Get the reference bases from start to end
        Args:
            chrom (string): Chromsome to query
            start (int): Start position to query
            end (int): End position (not inclusive)
        Returns:
            string: The genome sequence
        Raises:
            ValueError - Invalid arguments
        """
        if start >= end:
            msg = f"Start/stop coordinates incorrect for: {chrom!s}:{start!s}-{end!s}"
            raise ValueError(msg)

        # force chr prefix
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"

        if chrom not in self.reference_fasta:
            msg = f"FASTA reference is missing entry for chromosome {chrom!s}"
            raise ValueError(msg)
        return str(self.reference_fasta[chrom][start:end])
