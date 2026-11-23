"""Additional input/report parsing regression and safety checks."""
from io import StringIO
import pytest
from illumina2vcf.illumina import IlluminaReader, HeaderError, DateError, DataError

COLS = ('Sample ID','SNP Name','Chr','Position','SNP','Allele1 - Plus','Allele2 - Plus','Plus/Minus Strand')

def report(rows, delimiter='\t'):
    return StringIO(delimiter.join(COLS)+'\n'+'\n'.join(delimiter.join(map(str,r)) for r in rows)+'\n')

def r(sample='S1', probe='rs1', chrom='1', pos=100, snp='[A/G]', a='A', b='G', strand='+'):
    return (sample,probe,chrom,pos,snp,a,b,strand)

def blocks(rows, delimiter='\t', blocklist=''):
    return list(IlluminaReader(delimiter,blocklist).generate_line_blocks(report(rows,delimiter)))

@pytest.mark.parametrize('delimiter', ['\t', ','])
def test_delimiter_support(delimiter):
    assert blocks([r(),r('S2')],delimiter)[0][0].sample_id == 'S1'

@pytest.mark.parametrize('date,expected', [('2022-06-02 09:00','20220602'), ('6/2/2022 09:00','20220602'),('2022/6/2','20220602')])
def test_header_supported_date_formats(date, expected):
    data = ['[Header]', 'Content\t GSAMD-test.bpm',f'Processing Date\t{date}','[Data]']
    value, source = IlluminaReader('\t').parse_header(data)
    assert value == expected
    assert source == 'GSAMD-test'

@pytest.mark.parametrize('header', [
    ['[Header]','Processing Date\t2022-06-02','[Data]'],
    ['[Header]','Content\t test.bpm','[Data]'],
    ['[Header]','Content\t test.bpm','Processing Date\t06-2022','[Data]'],
])
def test_header_missing_data_rejected(header):
    with pytest.raises((HeaderError,DateError,KeyError)):
        IlluminaReader('\t').parse_header(header)

def test_overlong_header_rejected():
    lines=['[Header]']+[f'key{i}\tvalue' for i in range(101)]
    with pytest.raises(HeaderError, match='too long'):
        IlluminaReader('\t').parse_header(lines)

@pytest.mark.parametrize('rows,expected', [
    ([r(pos=100),r(pos=200)], [100,200]),
    ([r(pos=200),r(pos=100)], [100,200]),
    ([r(chrom='2'),r(chrom='1')], ['1','2']),
])
def test_locus_sorting(rows,expected):
    result=blocks(rows)
    observed=([b[0].chrom for b in result] if isinstance(expected[0],str) else [b[0].pos for b in result])
    assert observed==expected

def test_same_locus_multiple_samples_grouped():
    got=blocks([r('S2'),r('S1'),r('S3')])
    assert len(got)==1
    assert [x.sample_id for x in got[0]]==['S1','S2','S3']

def test_empty_input_yields_no_loci():
    assert blocks([])==[]

@pytest.mark.parametrize('name', ['bad','badilmndup'])
def test_blocklisted_and_illumina_duplicate_probes(name,tmp_path):
    p=tmp_path/'exclude.txt'
    p.write_text('bad\n',encoding='utf8')
    got=blocks([r(probe=name),r(probe='good')],blocklist=str(p))
    assert [x.snp_name for b in got for x in b]==['good']

def test_blocklist_is_exact_match(tmp_path):
    p=tmp_path/'exclude.txt';p.write_text('rs1\n',encoding='utf8')
    got=blocks([r(probe='rs1_extra')],blocklist=str(p))
    assert got[0][0].snp_name=='rs1_extra'

@pytest.mark.parametrize('position', ['NaN','', '1.5','-'])
def test_invalid_position_fails(position):
    with pytest.raises((ValueError,TypeError)):
        blocks([r(pos=position)])


