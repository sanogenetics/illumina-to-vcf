import csv
import gzip
import random
from dataclasses import dataclass
from datetime import datetime
from io import StringIO
from typing import Generator, List, Tuple

from illumina2vcf.bpm.BPMRecord import COMPLEMENT_MAP


@dataclass
class Probe:
    name: str
    chrom: str
    pos: int
    a: str
    b: str
    strand: str


class IlluminaBuilder:
    _sano: bool = False
    _rng: random.Random

    def __init__(self):
        self._rng = random.Random(42)
        pass

    def sano(self, sano: bool) -> "IlluminaBuilder":
        self._sano = sano
        return self

    def _generate_probes(self) -> Generator[Probe, None, None]:
        with gzip.open("tests/data/GSA-24v3-0_A2.trim.csv.gz", "rt") as infile:
            # read through the header
            for line in infile:
                if line.strip() == "[Assay]":
                    break
            # switch to CSV
            reader = csv.DictReader(infile)
            for row in reader:
                yield Probe(
                    name=row["Name"],
                    chrom=row["Chr"],
                    pos=int(row["MapInfo"]),
                    a=row["SNP"][1],  # [X/.]
                    b=row["SNP"][-2],  # [./X]
                    strand=row["RefStrand"],
                )

    def _generate_sorted_probes(self) -> List[Probe]:
        def probekey(probe: Probe) -> Tuple[str, int, int]:
            if probe.chrom.isnumeric():
                return ("", int(probe.chrom), probe.pos)
            else:
                return (probe.chrom, 0, probe.pos)

        return sorted(self._generate_probes(), key=probekey)

    def _generate_header_lines(self, num_snps=730059, num_samples=24):
        yield "[Header]"
        yield "GSGT Version	2.0.4"
        # e.g. "Processing Date	2022-06-02 9:25 AM"
        yield f"Processing Date	{datetime.now().strftime('%Y-%m-%d %I:%M %p')}"
        yield "Content		GSAMD-24v3-0-EA_20034606_A2.bpm"
        yield f"Num SNPs	{num_snps}"
        yield f"Total SNPs	{num_snps}"
        yield f"Num Samples	{num_samples}"
        yield f"Total Samples	{num_samples}"
        if self._sano:
            yield "SANO"

    @staticmethod
    def _map_to_line(header, data):
        line_items = []
        for heading in header:
            line_items.append(str(data.get(heading, ".")))
        return "\t".join(line_items)

    def _generate_data_lines(self, samples, probes):
        yield "[Data]"
        header = [
            "Sample ID",
            "RsID",
            "GC Score",
            "SNP Name",
            "SNP Index",
            "Sample Index",
            "Sample Name",
            "Sample Group",
            "SNP Aux",
            "Chr",
            "Position",
            "GT Score",
            "Cluster Sep",
            "SNP",
            "ILMN Strand",
            "Customer Strand",
            "Top Genomic Sequence",
            "Plus/Minus Strand",
            "Allele1 - Plus",
            "Allele2 - Plus",
            "Allele1 - Forward",
            "Allele2 - Forward",
            "Theta",
            "R",
            "X",
            "Y",
            "X Raw",
            "Y Raw",
            "B Allele Freq",
            "Log R Ratio",
            "CNV Value",
            "CNV Confidence",
        ]
        yield "\t".join(header)
        for i, sample in enumerate(samples, 1):
            # TODO mark some as chr/pos 0/0
            for probe in probes:
                data = {}
                data["Sample ID"] = sample
                data["Sample Index"] = i
                data["SNP Name"] = probe.name
                data["Chr"] = probe.chrom
                data["Position"] = probe.pos
                data["SNP"] = f"[{probe.a}/{probe.b}]"
                data["Plus/Minus Strand"] = probe.strand
                alleles = (
                    [COMPLEMENT_MAP[probe.a], COMPLEMENT_MAP[probe.b]] if probe.strand == "-" else [probe.a, probe.b]
                )
                data["Allele1 - Plus"] = self._rng.choice(alleles)
                data["Allele2 - Plus"] = self._rng.choice(alleles)
                yield self._map_to_line(header, data)

    def _generate_lines(self, samples=["sample1"]):
        # first generate the probes
        # make sure they are sorted by chromosome and position
        probes = self._generate_sorted_probes()
        # now we know how many probes and how many samples we have can generate header
        for line in self._generate_header_lines(num_snps=len(probes), num_samples=len(samples)):
            yield line

        # now we can generate data lines for each sample
        for line in self._generate_data_lines(samples, probes):
            yield line

    def build_file(self) -> StringIO:
        return StringIO("\n".join(self._generate_lines()))
