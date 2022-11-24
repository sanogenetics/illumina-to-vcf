from typing import Dict, Generator, Iterable, List, Tuple, Union

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

    def get_reference_bases(self, chrom, start, end):
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
            raise ValueError("Start/stop coordinates incorrect for: " + str(chrom) + ":" + str(start) + "-" + str(end))
        # force chr prefix
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        if chrom not in self.reference_fasta:
            raise ValueError("FASTA reference is missing entry for chromosome " + str(chrom))
        ref_base = str(self.reference_fasta[chrom][start:end])
        return ref_base
