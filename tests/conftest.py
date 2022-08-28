import random
from dataclasses import dataclass
from datetime import datetime
from io import StringIO

import pytest


@dataclass
class Probe:
    name: str
    chrom: str
    pos: int
    a: str
    b: str
    strand: str


class IlluminaBuilder:
    def __init__(self):
        pass

    def _generate_chrom_pos(self, seed=42):
        rng = random.Random(seed)
        # human is about 3 billion over 23 chromosome
        # we can simplify a bit, 8 chromosomes and 80m bases
        chroms = [str(i) for i in range(1, 5)] + ["X", "Y", "MT"]
        for chrom in chroms:
            pos = 10000  # minimum position on a chromosome
            while pos < 10000000:
                # GSA+MD is about 700,000 probes over 3b bases, so ~5kb apart
                pos += rng.randint(1, 10000)
                # simulate PAR
                if chrom in ["X", "Y"] and pos < 2781479:
                    chrom = "XY"
                yield chrom, pos

    def _generate_ref(self, salt=42):
        for chrom, pos in self._generate_chrom_pos():
            rng = random.Random(str((chrom, pos, salt)))
            yield chrom, pos, rng.choice(("A", "T", "C", "G"))

    def _generate_probes(self, salt=42):
        for chrom, pos, ref in self._generate_ref():
            rng = random.Random(str((chrom, pos, ref, salt)))
            a = ref
            b = rng.choice(tuple(set(("A", "T", "C", "G")) - set(ref)))
            # TODO sometimes generate multiple probes at same location
            yield Probe(
                name=f"{chrom}:{pos}:{ref}",
                chrom=chrom,
                pos=pos,
                a=a,
                b=b,
                strand="+",  # TODO sometimes negative strand
            )

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
                yield self._map_to_line(header, data)

    def _generate_lines(self, samples=["sample1"]):
        # first generate the probes
        probes = tuple(self._generate_probes())
        # now we know how many probes and how many samples we have can generate header
        for line in self._generate_header_lines(num_snps=len(probes), num_samples=len(samples)):
            yield line

        # now we can generate data lines for each sample
        for line in self._generate_data_lines(samples, probes):
            yield line

    def build_file(self) -> StringIO:
        return StringIO("\n".join(self._generate_lines()))


@pytest.fixture
def fixture_illumina_lines():
    with IlluminaBuilder().build_file() as illumina:
        yield illumina
