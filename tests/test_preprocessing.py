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

            assert fastq_record.name == \
                "NB500952:186:H35F2BGX5:1:11101:14875:1080 1:N:0:ATCACG"
            assert fastq_record.seq == \
                "TATCTCTTTTATTTTTTTTTTTTGAAGGGCAGATTTAAAATACACTATTAAAATTATTAAATATTAAAAAACAACA"
            assert fastq_record.qual == \
                "AAA/A/EEEEEEEEEEEEEAEE/AAE///E///E</EE/A/EE</6//6//<66//6//6//////<//6EE/6/6"
    
    assert i == 249

def test_trim_5p_Ts():
    infile = "data/rawfastq/siCtrl_1_small.fastq"
    for i, fastq_record in enumerate(PASS.fastq_reader(infile, 6, 4)):
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
                "TS17TATCTC::ATTTTTTTTTTTT:NB500952:186:H35F2BGX5:1:11101:14875:1080 1:N:0:ATCACG"
            assert fastq_record.seq == \
                "GAAGGGCAGATTTAAAATACACTATTAAAATTATTAAATATTAAAAAACAACA"
            assert fastq_record.qual == \
                "AAE///E///E</EE/A/EE</6//6//<66//6//6//////<//6EE/6/6"

    

 