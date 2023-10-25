import logging
from typing import Optional, Tuple

from illumina2vcf.bpm.illuminabeadarrayfiles import RefStrand
from illumina2vcf.bpm.referencegenome import ReferenceGenome

logger = logging.getLogger(__name__)

COMPLEMENT_MAP = dict(zip("ABCDGHKMRTVYNID", "TVGHCDMKYABRNID"))
DEGENERACY_MAP = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ATCG",
}


def reverse(sequence: str) -> str:
    """
    Reverse a string
    Args:
        sequence    Sequence to reverse
    Returns:
        The reversed sequence
    """
    return sequence[::-1]


def complement(sequence: str) -> str:
    """
    Complement a nucleotide sequence. Note that complement of D and I are D and I,
    respectively. This is intended to be called on the "SNP" portion of a source sequence.
    Args:
        sequence    The input sequence
    Returns:
        The complemented sequence
    """
    return "".join(COMPLEMENT_MAP[x] for x in sequence)


def determine_left_shift(five_prime: str, indel: str, three_prime: str) -> Tuple[str, str]:
    """
    Adjust 5' and 3' context of indel such that
    indel is fully shifted to 5'
    Args:
        five_prime      Five prime sequence
        indel           Sequence of indel
        three_prime     Three prime sequence
    Returns:
        Tuple of new 5' and 3' sequences
    """
    while five_prime.endswith(indel):
        five_prime = five_prime[: -len(indel)]
        three_prime = indel + three_prime
    # may have not fully shifted homopolymer
    while len(indel) * five_prime[-1] == indel:
        three_prime = five_prime[-1] + three_prime
        five_prime = five_prime[:-1]
    return (five_prime, three_prime)


def reverse_complement(sequence: str) -> str:
    """
    Reverse complement a sequence
    Args:
        sequence    The input sequence
    Returns:
        The reverse-complement of the input sequence
    """
    return reverse(complement(sequence))


def max_suffix_match(str1: str, str2: str) -> int:
    """
    Determine the maximum length of exact suffix
    match between str1 and str2
    str2 may contain degenerate IUPAC characters
    Args:
        str1    First string
        str2    Second string
    Returns:
        Length of maximum suffix match
    """
    result = 0
    for char1, char2 in zip(str1[::-1], str2[::-1]):
        assert char1 in "ACGT"
        if char1 in DEGENERACY_MAP[char2]:
            result += 1
        else:
            break
    return result


def max_prefix_match(str1: str, str2: str) -> int:
    """
    Determine the maximum length of exact prefix
    match between str1 and str2
    str2 may contain degenerate IUPAC characters
    Args:
        str1    First string
        str2    Second string
    Returns:
        Length of maximum prefix match
    """
    result = 0
    for char1, char2 in zip(str1, str2):
        assert char1 in "ACGT"
        if char1 in DEGENERACY_MAP[char2]:
            result += 1
        else:
            break
    return result


class IndelSourceSequence:
    five_prime: str
    indel: str
    three_prime: str
    """
    Represents the source sequence for an indel
    Attributes:
        five_prime      Sequence 5' of indel (on the design strand)
        indel           Indel sequence (on the design strand)
        three_prime     Sequence 3' of indel (on the design strand)
    """

    def __init__(self, source_sequence):
        (self.five_prime, self.indel, self.three_prime) = self.split_source_sequence(source_sequence.upper())

    def get_split_sequence(self, generate_reverse_complement: bool, left_shift: bool) -> Tuple[str, str, str]:
        """
        Return the components of the indel source sequence
        Args:
            generate_reverse_complement         Return reverse complement of original source sequence
            left_shift                          Left shift position of indel on requested strand
        Returns:
            (five_prime, indel, three_prime)    Tuple of three components of indel source sequence
        """
        if generate_reverse_complement:
            (five_prime, indel, three_prime) = (
                reverse_complement(self.three_prime),
                reverse_complement(self.indel),
                reverse_complement(self.five_prime),
            )
        else:
            (five_prime, indel, three_prime) = (self.five_prime, self.indel, self.three_prime)

        if left_shift:
            (five_prime, three_prime) = determine_left_shift(five_prime, indel, three_prime)
        return (five_prime, indel, three_prime)

    @staticmethod
    def split_source_sequence(source_sequence: str) -> Tuple[str, str, str]:
        """
        Break source sequence into different piecdes
        Args:
            source_sequence             Source sequence string (e.g., ACGT[-/AGA]ATAT)
        Returns:
            (string, string, string)    Tuple with 5' sequence, indel sequence, 3' sequence
        """
        left_position = source_sequence.find("/")
        right_position = source_sequence.find("]")
        assert source_sequence[left_position - 1] == "-"
        return (
            source_sequence[: (left_position - 2)],
            source_sequence[(left_position + 1) : right_position],
            source_sequence[(right_position + 1) :],
        )


class BPMRecord:
    name: str
    address_a: str
    probe_a: str
    chromosome: str
    pos: int
    snp: str
    ref_strand: int
    assay_type: int
    indel_source_sequence: Optional[IndelSourceSequence]
    index_num: int
    is_deletion: Optional[bool]
    """
    Represents entry from a manifest.abs
    Attributes:
        name                    Entry name (unique)
        address_a               Address of probe A
        probe_a (string): Sequence of probe A
        chromosome              Chromsome name
        pos                     Mapping position
        snp                     SNP variation (e.g., [A/C])
        ref_strand              Reference strand of snp
        assay_type              0 for Inf II, 1 for Inf I
        indel_source_sequence   Sequence of indel (on design strand), None for SNV
        index_num               Index in original manifest
        is_deletion             Whether indel record represents deletion, None for SNV
    """

    def __init__(
        self,
        name: str,
        address_a: str,
        probe_a: str,
        chromosome: str,
        pos: int,
        snp: str,
        ref_strand: int,
        assay_type: int,
        indel_source_sequence: Optional[IndelSourceSequence],
        source_strand: str,
        ilmn_strand: str,
        genome_reader: Optional[ReferenceGenome],
        index: int,
    ):
        """
        Create a new BPM record
        Args:
            name                    Name field from manifest
            address_a               AddressA_ID field from manifest
            chromosome              Chr field from manifest
            pos                     MapInfo field from manifest
            ref_strand              RefStrand from manifest
            assay_type              0 for Inf II, 1 for Inf I
            indel_source_sequence   Source sequence for indel, may be None for SNV
            source_strand           SourceStrand field from manifest
            ilmn_strand             IlmnStrand field from manifest
            genome_reader           Allows query of genomic sequence, may be None for SNV
            index                   Index of entry within manifest/GTC files
        """
        self.name = name
        self.address_a = address_a
        self.probe_a = probe_a
        self.chromosome = chromosome
        self.pos = pos
        self.snp = snp
        self.ref_strand = ref_strand
        self.assay_type = assay_type
        self.indel_source_sequence = indel_source_sequence
        self._genome_reader = genome_reader
        self.index_num = index

        self.plus_strand_alleles = self._determine_plus_strand_alleles(snp, ref_strand)

        if self.indel_source_sequence:
            source_strand = source_strand[0].upper()
            ilmn_strand = ilmn_strand[0].upper()

            if source_strand == "U" or ilmn_strand == "U":
                raise ValueError('Unable to process indel with customer or ILMN strand value of "U"')

            if source_strand in ["P", "M"]:
                assert ilmn_strand in ["P", "M"]
            else:
                assert ilmn_strand in ["T", "B"]

            self.is_source_on_design_strand = source_strand == ilmn_strand
            self.is_deletion = self._calculate_is_deletion()
        else:
            self.is_source_on_design_strand = None
            self.is_deletion = None

    def is_indel(self) -> bool:
        """
        Check whether a BPM record represents an indel
        Args:
            None
        Returns:
            bool: True if record is indel
        """
        return "D" in self.snp

    def get_indel_source_sequences(self, ref_strand: int) -> Tuple[str, str, str]:
        assert self.indel_source_sequence
        return self.indel_source_sequence.get_split_sequence(
            self.is_source_on_design_strand != (self.ref_strand == ref_strand), True
        )

    def _calculate_is_deletion(self) -> Optional[bool]:
        if self.chromosome == "0" or self.pos == 0:
            return None

        start_index = self.pos - 1
        chromosome = "X" if self.chromosome == "XY" else self.chromosome

        # get indel sequence on the plus strand
        (five_prime, indel_sequence, three_prime) = self.get_indel_source_sequences(RefStrand.Plus)

        assert self._genome_reader

        genomic_sequence = self._genome_reader.get_reference_bases(
            chromosome, start_index, start_index + len(indel_sequence)
        )
        indel_sequence_match = indel_sequence == genomic_sequence

        genomic_deletion_five_prime = self._genome_reader.get_reference_bases(
            chromosome, start_index - len(five_prime), start_index
        )
        genomic_deletion_three_prime = self._genome_reader.get_reference_bases(
            chromosome, start_index + len(indel_sequence), start_index + len(indel_sequence) + len(three_prime)
        )
        (genomic_deletion_five_prime, genomic_deletion_three_prime) = determine_left_shift(
            genomic_deletion_five_prime, indel_sequence, genomic_deletion_three_prime
        )

        genomic_insertion_five_prime = self._genome_reader.get_reference_bases(
            chromosome, start_index - len(five_prime) + 1, start_index + 1
        )
        genomic_insertion_three_prime = self._genome_reader.get_reference_bases(
            chromosome, start_index + 1, start_index + len(three_prime) + 1
        )
        (genomic_insertion_five_prime, genomic_insertion_three_prime) = determine_left_shift(
            genomic_insertion_five_prime, indel_sequence, genomic_insertion_three_prime
        )

        deletion_context_match_lengths = (
            max_suffix_match(genomic_deletion_five_prime, five_prime),
            max_prefix_match(genomic_deletion_three_prime, three_prime),
        )
        max_deletion_context = (
            min(len(genomic_deletion_five_prime), len(five_prime))
            + min(len(genomic_deletion_three_prime), len(three_prime))
            + len(indel_sequence)
        )
        deletion_context_score = (
            sum(deletion_context_match_lengths) + len(indel_sequence) if indel_sequence_match else 0
        ) / float(max_deletion_context)

        insertion_context_match_lengths = (
            max_suffix_match(genomic_insertion_five_prime, five_prime),
            max_prefix_match(genomic_insertion_three_prime, three_prime),
        )
        max_insertion_context = min(len(genomic_insertion_five_prime), len(five_prime)) + min(
            len(genomic_insertion_three_prime), len(three_prime)
        )
        insertion_context_score = sum(insertion_context_match_lengths) / float(max_insertion_context)

        is_deletion = (
            indel_sequence_match
            and deletion_context_score > insertion_context_score
            and min(deletion_context_match_lengths) >= 1
        )

        is_insertion = insertion_context_score > deletion_context_score and min(insertion_context_match_lengths) >= 1

        if is_deletion == is_insertion:
            raise Exception("Unable to determine reference allele for indel")

        if is_deletion:
            if deletion_context_score < 1.0:
                logger.warning("Incomplete match of source sequence to genome for indel " + self.name)

        if is_insertion:
            if insertion_context_score < 1.0:
                logger.warning("Incomplete match of source sequence to genome for indel " + self.name)

        return is_deletion

    def _determine_plus_strand_alleles(self, snp: str, ref_strand: int) -> Tuple[str, str]:
        """
        Return the nucleotides alleles for the record on the plus strand.
        If record is indel, will return alleles in terms of D/I SNP convention
        Args:
            None
        Returns
            tuple of alleles
        Raises:
            Exception - Record does not contains reference strand information
        """
        nucleotides = (snp[1], snp[-2])
        if ref_strand == RefStrand.Plus:
            return nucleotides
        elif ref_strand == RefStrand.Minus:
            return tuple(COMPLEMENT_MAP[nucleotide] for nucleotide in nucleotides)
        else:
            raise Exception("Manifest must contain reference strand information")
