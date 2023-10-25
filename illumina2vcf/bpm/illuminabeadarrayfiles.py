class SourceStrand:
    Unknown = 0
    Forward = 1
    Reverse = 2

    @staticmethod
    def to_string(source_strand):
        """Get an integer representation of source strand annotation
        Args:
            source_strand (str) : string representation of source strand annotation (e.g., "F")

        Returns:
            int : int representation of source strand annotation (e.g. SourceStrand.Forward)
        """
        if source_strand == SourceStrand.Unknown:
            return "U"
        elif source_strand == SourceStrand.Forward:
            return "F"
        elif source_strand == SourceStrand.Reverse:
            return "R"
        else:
            raise Exception("Unexpected value for source strand " + source_strand)

    @staticmethod
    def from_string(source_strand):
        """Get a string representation of source strand annotation
        Args:
            source_strand (int) : int representation of source strand (e.g., SourceStrand.Forward)

        Returns:
            str : string representation of source strand annotation
        """
        if source_strand == "U" or source_strand == "":
            return SourceStrand.Unknown
        if source_strand == "F":
            return SourceStrand.Forward
        elif source_strand == "R":
            return SourceStrand.Reverse
        else:
            raise Exception("Unexpected value for source strand " + source_strand)


class RefStrand:
    Unknown = 0
    Plus = 1
    Minus = 2

    @staticmethod
    def to_string(ref_strand):
        """Get a string reprensetation of ref strand annotation
        Args:
            ref_strand (int) : int representation of ref strand (e.g., RefStrand.Plus)

        Returns:
            str : string representation of reference strand annotation
        """
        if ref_strand == RefStrand.Unknown:
            return "U"
        elif ref_strand == RefStrand.Plus:
            return "+"
        elif ref_strand == RefStrand.Minus:
            return "-"
        else:
            raise Exception("Unexpected value for reference strand " + ref_strand)

    @staticmethod
    def from_string(ref_strand):
        """Get an integer representation of ref strand annotation
        Args:
            ref_strand (str) : string representation of reference strand annotation (e.g., "+")

        Returns:
            int : int representation of reference strand annotation (e.g. RefStrand.Plus)
        """
        if ref_strand == "U" or ref_strand == "":
            return RefStrand.Unknown
        if ref_strand == "+":
            return RefStrand.Plus
        elif ref_strand == "-":
            return RefStrand.Minus
        else:
            raise Exception("Unexpected value for reference strand " + ref_strand)


class LocusEntry:
    """Helper class representing a locus entry within a bead pool manifest. Current only support version
    6,7, and 8.
    Attributes:
        ilmn_id (string) : IlmnID (probe identifier) of locus
        name (string): Name (variant identifier) of locus
        snp (string) : SNP value for locus (e.g., [A/C])
        chrom (string) : Chromosome for the locus (e.g., XY)
        map_info (int) : Mapping location of locus
        assay_type (int) : Identifies type of assay (0 - Infinium II , 1 - Infinium I (A/T), 2 - Infinium I (G/C)
        address_a (int) : AddressA ID of locus
        address_b (int) : AddressB ID of locus (0 if none)
        ref_strand (int) : See RefStrand class
        source_strand (int) : See SourceStrand class
    """

    def __init__(self, handle):
        """Constructor
        Args:
            handle (file handle):  File handle at start of locus entry record
        Returns:
            LocusEntry
        """
        self.ilmn_id = ""
        self.name = ""
        self.snp = ""
        self.chrom = ""
        self.map_info = -1
        self.assay_type = -1
        self.address_a = -1
        self.address_b = -1
        self.ref_strand = RefStrand.Unknown
        self.source_strand = SourceStrand.Unknown
        self.__parse_file(handle)

    def __parse_file(self, handle):
        """Helper function to initialize this object from a file handle
        Args:
            handle (file handle): File handle at start of locus entry record
        Returns:
            None
        """
        version = read_int(handle)
        if version == 6:
            self.__parse_locus_version_6(handle)
        elif version == 7:
            self.__parse_locus_version_7(handle)
        elif version == 8:
            self.__parse_locus_version_8(handle)
        else:
            raise Exception("Manifest format error: unknown version for locus entry (" + str(version) + ")")

    def __parse_locus_version_6(self, handle):
        """Helper function to parse version 6 locus entry
        Args:
            handle (file handle): File handle at start of locus entry record
        Returns:
            None
        """
        self.ilmn_id = read_string(handle)
        self.source_strand = SourceStrand.from_string(self.ilmn_id.split("_")[-2])
        self.name = read_string(handle)
        for idx in range(3):
            read_string(handle)
        handle.read(4)
        for idx in range(2):
            read_string(handle)
        self.snp = read_string(handle)
        self.chrom = read_string(handle)
        for idx in range(2):
            read_string(handle)
        self.map_info = int(read_string(handle))
        for idx in range(2):
            read_string(handle)
        self.address_a = read_int(handle)
        self.address_b = read_int(handle)
        for idx in range(7):
            read_string(handle)
        handle.read(3)
        self.assay_type = read_byte(handle)
        if self.assay_type not in [0, 1, 2]:
            raise Exception("Format error in reading assay type from locus entry")
        if self.address_b == 0:
            if self.assay_type != 0:
                raise Exception("Manifest format error: Assay type is inconsistent with address B")
        else:
            if self.assay_type == 0:
                raise Exception("Manifest format error: Assay type is inconsistent with address B")

    def __parse_locus_version_7(self, handle):
        """Helper function to parse version 7 locus entry
        Args:
            handle (file handle): File handle at start of locus entry record
        Returns:
            None
        """
        self.__parse_locus_version_6(handle)
        handle.read(4 * 4)

    def __parse_locus_version_8(self, handle):
        """Helper function to parse version 8 locus entry
        Args:
            handle (file handle): File handle at start of locus entry record
        Returns:
            None
        """
        self.__parse_locus_version_7(handle)
        self.ref_strand = RefStrand.from_string(read_string(handle))


complement_map = {"A": "T", "T": "A", "C": "G", "G": "C", "D": "D", "I": "I"}


def complement(nucleotide):
    """Complement a single nucleotide. Complements of D(eletion) and I(nsertion) are D and I, respectively.
    Args:
        nucleotide (string) : Nucleotide, must be A, C, T, G, D, or I
    Returns:
        str : Complemented nucleotide
    """
    if nucleotide in complement_map:
        return complement_map[nucleotide]
    raise ValueError("Nucleotide must be one of A, C, T, G, D, or I")


def read_char(handle):
    """Helper function to parse character from file handle
    Args:
        handle (file handle): File handle
    Returns:
        char value read from handle
    """
    return handle.read(1)


def read_ushort(handle):
    """Helper function to parse ushort from file handle
    Args:
        handle (file handle): File handle
    Returns:
        numpy.int16 value read from handle
    """
    return frombuffer(handle.read(2), dtype=uint16)[0]
