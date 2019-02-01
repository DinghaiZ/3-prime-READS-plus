# Test the code used for preprocessing fastq files
import pytest
import sys
from pathlib import Path
from subprocess import check_output

sys.path.append('/home/dinghai/dev/3-prime-READS-plus/modules')
import PASS

def test_fastq_reader():
    infile = "data/rawfastq/siCtrl_1_small.fastq"
    for i, fastq_record in enumerate(PASS.fastq_reader(infile)):
        if i == 0:
            assert fastq_record.name == \
                "NB500952:186:H35F2BGX5:1:11101:6050:1064 1:N:0:ATCACG"