"""Stream-level tests for output safety and explicit outstanding requirements."""
from types import SimpleNamespace
import pytest
from illumina2vcf.illumina import IlluminaRow
from illumina2vcf.vcf import VCFMaker, ConverterError

class Genome:
    def __init__(self):
        self.reference_fasta=SimpleNamespace(faidx=SimpleNamespace(index={'1':SimpleNamespace(rlen=1000)}))
    def get_reference_bases(self,chrom,start,end):
        return {99:'A',199:'C',299:'G'}.get(start,'')

def row(sample='S1',pos=100,snp='[A/G]',a='A',b='G',name='rs1'):
    return IlluminaRow(sample,name,'1',pos,snp,a,b,'+')

def data(lines):
    return [x for x in lines if not (x.comment_key or x.comment_raw)]

def test_successful_loci_are_output_in_order():
    maker=VCFMaker(Genome())
    out=data(maker.generate_lines([[row(pos=100)],[row(pos=200,snp='[C/T]',a='C',b='T')],[row(pos=300,snp='[G/A]',a='G',b='A')]]))
    assert [v.pos for v in out]==[100,200,300]

def test_invalid_locus_is_currently_logged_and_dropped(caplog):
    maker=VCFMaker(Genome())
    out=data(maker.generate_lines([[row(pos=100)],[row(pos=200,snp='[A/T]',a='A', b='T')],[row(pos=300,snp='[G/A]',a='G',b='A')]]))
    assert [v.pos for v in out]==[100,300]
    assert 'Reference' in caplog.text

def test_column_header_emitted_once():
    maker=VCFMaker(Genome())
    out=list(maker.generate_lines([[row()],[row(pos=200,snp='[C/T]',a='C',b='T')]]))
    assert sum(bool(v.comment_raw) for v in out)==1

def test_sample_genotypes_correspond_to_header_order():
    maker=VCFMaker(Genome())
    out=list(maker.generate_lines([[row('S2',a='G',b='G'),row('S1',a='A',b='A')]]))
    header=next(str(v) for v in out if v.comment_raw)
    # Current code uses first-seen order within the initial block.
    ids=header.split('\t')[-2:]
    values=[v.sample[i]['GT'] for v in data(out) for i in range(2)]
    assert dict(zip(ids,values))=={'S1':'0/0','S2':'1/1'}


