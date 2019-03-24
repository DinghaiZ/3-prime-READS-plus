# Test the code in PASS.py
import sys
sys.path.append('/home/dinghai/dev/3-prime-READS-plus/modules')
import pytest
from pathlib import Path
from subprocess import check_output
import PASS


@pytest.fixture(scope='module')
def file_to_trim():
    return "tests/data/rawfastq/siCtrl_1_small.fastq"

#@pytest.mark.skip()
def test_FastqRecord_get_fastq_record(file_to_trim):
    fastq_file_obj = PASS.FastqFile(file_to_trim, 6, 4)
    for i, fastq_record in enumerate(fastq_file_obj.get_fastq_record()):
        if i == 1:
            assert str(fastq_record) == "@NB500952:186:H35F2BGX5:1:11101:1487"\
                + "5:1080 1:N:0:ATCACG\nTATCTCTTTTATTTTTTTTTTTTGAAGGGCAGATTTA"\
                    + "AAATACACTATTAAAATTATTAAATATTAAAAAACAACA\n+\nAAA/A/EEEE"\
                        + "EEEEEEEEEAEE/AAE///E///E</EE/A/EE</6//6//<66//6//6"\
                            + "//////<//6EE/6/6\n"
            assert fastq_record.name == \
                "NB500952:186:H35F2BGX5:1:11101:14875:1080 1:N:0:ATCACG"
            assert fastq_record.seq == \
                "TATCTCTTTTATTTTTTTTTTTTGAAGGGCAGATTTAAAATACACTATTAAAATTATTAA"\
                    + "ATATTAAAAAACAACA"
            assert fastq_record.qual == \
                "AAA/A/EEEEEEEEEEEEEAEE/AAE///E///E</EE/A/EE</6//6//<66//6//6/"\
                    + "/////<//6EE/6/6"
            print(repr(fastq_record), '\n')
            print(str(fastq_record), '\n')



def FastqRecord_test_trim_5p_Ts(file_to_trim):
    fastq_file_obj = PASS.FastqFile(file_to_trim, 6, 4)
    for i, fastq_record in enumerate(fastq_file_obj.get_fastq_record()):
        fastq_record.trim_5p_Ts(fastq_file_obj.randNT5)
        # No trimming
        if i == 0:
            assert fastq_record.name == \
                "TS0GAGGTC:NB500952:186:H35F2BGX5:1:11101:6050:1064 1:N:0:ATCACG"
            assert fastq_record.seq == \
                "CGCGCGCCCCCCTGTATAGAAATCAACATTCTCTCCCAGAATTCTGTATATCTGT"\
                    + "AACTATAGAAATGTC"
            assert fastq_record.qual == \
                "EEEAAEEEEEEE/E/E/A/E6E/AAE/E/</</AEE/6/////A<//<A////6/"\
                    + "//6/A</EEEAA/<6"
        # Trim twice
        if i == 1:
            assert fastq_record.name == \
                "TS17TATCTC::ATTTTTTTTTTTT:NB500952:186:H35F2BGX5:1:11101:148"\
                    + "75:1080 1:N:0:ATCACG"
            assert fastq_record.seq == \
                "GAAGGGCAGATTTAAAATACACTATTAAAATTATTAAATATTAAAAAACAACA"
            assert fastq_record.qual == \
                "AAE///E///E</EE/A/EE</6//6//<66//6//6//////<//6EE/6/6"
        # Trim once
        if i == 2:
            assert fastq_record.name == \
                "TS16AACCAT:NB500952:186:H35F2BGX5:1:11101:13764:1083 1:N:0:ATCACG"
            assert fastq_record.seq == \
                "AAAAAAAAAATAAAAAATATTTTTAAAATTATAAATTTTTTTTTGTAAAATCAC"
            assert fastq_record.qual == \
                "EEEAEE/E///6AEE////A//<<6//66//6//66//66/6/6////////66"


def test_FastqFile_incorrect_file_format():
    file_to_trim = "tests/data/rawfastq/siCtrl_1_small"
    with pytest.raises(AssertionError):
        PASS.FastqFile(file_to_trim, 6, 4)


def test_FastqFile_repr(file_to_trim):
    fastq_file_obj = PASS.FastqFile(file_to_trim, 6, 4)
    assert repr(fastq_file_obj) == "FastqFile('tests/data/rawfastq/siCtrl_1_small.fastq', 6, 4)"
    print(fastq_file_obj)



def test_FastqFile_get_read_num(file_to_trim):
    fastq_file_obj = PASS.FastqFile(file_to_trim, 6, 4)
    assert fastq_file_obj.get_read_num() == 250


def test_FastqFile_create_trimmed_fastq_file(file_to_trim):
    fastq_file_obj = PASS.FastqFile(file_to_trim, 6, 4)
    fastq_file_obj.create_trimmed_fastq_file()
    assert fastq_file_obj.trimmed5T_num == 142
    assert fastq_file_obj.get_read_num() == 250


@pytest.fixture(scope='module')
def genome_mm9():
    print('\nLoading mm9 genome from fasta files...')
    genome_mm9 = PASS.load_fasta_genome(
        '/home/dinghai/projects/fud/ucsc/genomes/mm9')
    print('mm9 genome loaded')
    print('mm9 genome keys: ', sorted(genome_mm9.keys()))
    yield genome_mm9
    del genome_mm9
    print('\nmm9 genome unloaded')


@pytest.mark.slow
def test_load_fasta_genome(genome_mm9):
    assert isinstance(genome_mm9, dict)
    assert 'chr1' in genome_mm9.keys()


@pytest.mark.slow
@pytest.mark.parametrize('chromosome, strand, start, end, expected_output',
                         [
                             ('chr2', "+", 100000, 100020, 
                             'N'*21),
                             ('chr16', "-", 64855900, 64855927, 
                             'ATCACAATGGTCATTAACTGCCTCGTCA'),
                             ('chrM', "+", 100, 120, 
                             'AGAGGTAAAATTACACATGCA')
                         ]
                        )
def test_get_seq(chromosome, strand, start, end, genome_mm9, expected_output):
    seq = PASS.get_seq(chromosome, strand, start, end, genome_mm9) 
    assert seq == expected_output

