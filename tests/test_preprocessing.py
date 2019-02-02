# Test the code used for preprocessing fastq files
import pytest
import sys
from pathlib import Path
from subprocess import check_output

sys.path.append('/home/dinghai/dev/3-prime-READS-plus/modules')
import PASS

def test_fastq_reader():
    infile = "data/rawfastq/siCtrl_1_small.fastq"
    for i, fastq_record in enumerate(PASS.fastq_reader(infile, 6, 4)):
        if i == 1:
            assert fastq_record.randNT5 == 6
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
    infile = "data/rawfastq/siCtrl_1_small.fastq"
    for i, fastq_record in enumerate(PASS.fastq_reader(infile, 6, 4)):
        # No trimming
        if i == 0:
            fastq_record.trim_5p_Ts()
            assert fastq_record.name == \
                "NB500952:186:H35F2BGX5:1:11101:6050:1064 1:N:0:ATCACG"
            assert fastq_record.seq == \
                "GAGGTCCGCGCGCCCCCCTGTATAGAAATCAACATTCTCTCCCAGAATTCTGTATATCTGT"\
                    + "AACTATAGAAATGTC"
            assert fastq_record.qual == \
                "AAAAAEEEEAAEEEEEEE/E/E/A/E6E/AAE/E/</</AEE/6/////A<//<A////6/"\
                    + "//6/A</EEEAA/<6"
        # Trim once
        if i == 2:
            fastq_record.trim_5p_Ts()
            assert fastq_record.name == \
                "TS16AACCAT:NB500952:186:H35F2BGX5:1:11101:13764:1083 1:N:0:ATCACG"
            assert fastq_record.seq == \
                "AAAAAAAAAATAAAAAATATTTTTAAAATTATAAATTTTTTTTTGTAAAATCAC"
            assert fastq_record.qual == \
                "EEEAEE/E///6AEE////A//<<6//66//6//66//66/6/6////////66"
        # Trim twice
        if i == 1:
            fastq_record.trim_5p_Ts()
            assert fastq_record.name == \
                "TS17TATCTC::ATTTTTTTTTTTT:NB500952:186:H35F2BGX5:1:11101:148"\
                    + "75:1080 1:N:0:ATCACG"
            assert fastq_record.seq == \
                "GAAGGGCAGATTTAAAATACACTATTAAAATTATTAAATATTAAAAAACAACA"
            assert fastq_record.qual == \
                "AAE///E///E</EE/A/EE</6//6//<66//6//6//////<//6EE/6/6"

    

 