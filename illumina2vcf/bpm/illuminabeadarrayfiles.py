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
            msg = f"Unexpected value for reference strand {ref_strand}"
            raise RuntimeError(msg)

    @staticmethod
    def from_string(ref_strand):
        """Get an integer representation of ref strand annotation
        Args:
            ref_strand (str) : string representation of reference strand annotation (e.g., "+")

        Returns:
            int : int representation of reference strand annotation (e.g. RefStrand.Plus)
        """
        if ref_strand in ["U", ""]:
            return RefStrand.Unknown
        if ref_strand == "+":
            return RefStrand.Plus
        elif ref_strand == "-":
            return RefStrand.Minus
        else:
            msg = f"Unexpected value for reference strand {ref_strand}"
            raise RuntimeError(msg)
