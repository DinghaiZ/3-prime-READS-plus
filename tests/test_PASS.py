# Test the code in PASS.py
import sys
sys.path.append('/home/dinghai/dev/3-prime-READS-plus/modules')
import pytest
from pathlib import Path
from subprocess import check_output
import PASS

def test_fastq_reader():
    infile = "tests/data/rawfastq/siCtrl_1_small.fastq"
    for i, fastq_record in enumerate(PASS.fastq_reader(infile, 6, 4)):
        if i == 1:
            assert fastq_record.randNT5 == 6
            assert fastq_record.randNT3 == 4
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
    
    assert i == fastq_record.read_num - 1

def test_trim_5p_Ts():
    infile = "tests/data/rawfastq/siCtrl_1_small.fastq"
    for i, fastq_record in enumerate(PASS.fastq_reader(infile, 6, 4)):
        fastq_record.trim_5p_Ts()
        # No trimming
        if i == 0:
            assert fastq_record.name == \
                "NB500952:186:H35F2BGX5:1:11101:6050:1064 1:N:0:ATCACG"
            assert fastq_record.seq == \
                "GAGGTCCGCGCGCCCCCCTGTATAGAAATCAACATTCTCTCCCAGAATTCTGTATATCTGT"\
                    + "AACTATAGAAATGTC"
            assert fastq_record.qual == \
                "AAAAAEEEEAAEEEEEEE/E/E/A/E6E/AAE/E/</</AEE/6/////A<//<A////6/"\
                    + "//6/A</EEEAA/<6"
            assert fastq_record.trimmed_num == 0
        # Trim twice
        if i == 1:
            assert fastq_record.name == \
                "TS17TATCTC::ATTTTTTTTTTTT:NB500952:186:H35F2BGX5:1:11101:148"\
                    + "75:1080 1:N:0:ATCACG"
            assert fastq_record.seq == \
                "GAAGGGCAGATTTAAAATACACTATTAAAATTATTAAATATTAAAAAACAACA"
            assert fastq_record.qual == \
                "AAE///E///E</EE/A/EE</6//6//<66//6//6//////<//6EE/6/6"
            assert fastq_record.trimmed_num == 1
        # Trim once
        if i == 2:
            assert fastq_record.name == \
                "TS16AACCAT:NB500952:186:H35F2BGX5:1:11101:13764:1083 1:N:0:ATCACG"
            assert fastq_record.seq == \
                "AAAAAAAAAATAAAAAATATTTTTAAAATTATAAATTTTTTTTTGTAAAATCAC"
            assert fastq_record.qual == \
                "EEEAEE/E///6AEE////A//<<6//66//6//66//66/6/6////////66"
            assert fastq_record.trimmed_num == 2
       
    print("\nTotal trimmed fastq reads: ", fastq_record.trimmed_num)


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


def test_genome_mm9(genome_mm9):
    assert isinstance(genome_mm9, dict)
    assert 'chr1' in genome_mm9.keys()


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


def test_fastq_file_trimmer(infile, randNT5, randNT3):
    outfile = infile.replace('.fastq', '.trimmed.fastq') 
    fastq_file_trimme(infile, randNT5, randNT3)
    assert 0


def test_pick_PASS():
    pass



    

 