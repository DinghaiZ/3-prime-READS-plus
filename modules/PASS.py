# -*- coding: utf-8 -*-
'''This module contains functions for preprocessing fastq reads with long 5'
T-stretches, identification of (unique) PASS (Poly(A) Sequence Supporting) reads 
in sam files, clustering of PASS reads into pA sites, generation of QC reports, 
and visualization of PASS reads on UCSC genome browser.
'''

import os
import re
from pathlib import Path
from subprocess import check_output


class FastqRecord():
    '''Fastq record with a method to trim 5' T-stretches'''
    def __init__(self, name, seq, qual):
        self.name = name
        self.seq = seq
        self.qual = qual
        self.trimmed5T = False
        self.trimmed3N = False

    def __repr__(self):
        return f"FastqRecord('{self.name}', '{self.seq}', '{self.qual}')"

    def __str__(self):
        return '\n'.join(['@' + self.name, self.seq, '+', self.qual + '\n'])

    def __len__(self):
        return len(self.seq)
    
    def trim_5p_Ts(self, randNT5):
        '''Trim random nucleotides and T-stretchs from the 5' of the reads. 

        Uses a two-step trimming strategy to deal with sequencing errors in 
        T-stretches. Information about the trimmed nucleotides are saved in the 
        read name. 
        
        Example:

        @NB500952:186:H35F2BGX5:1:11101:14875:1080 1:N:0:ATCACG
        TATCTCTTTTATTTTTTTTTTTTGAAGGGCAGATTTAAAATACACTATTAAAATTATTAAATATTAAAAAA
        CAACA
        
        becomes:

        @TS17TATCTC::ATTTTTTTTTTTT:NB500952:186:H35F2BGX5:1:11101:14875:1080 1:
        N:0:ATCACG
        GAAGGGCAGATTTAAAATACACTATTAAAATTATTAAATATTAAAAAACAACA
        
        TS17: T-stretch length is 17. 
        TATCTC: Random nucleotides
        ATTTTTTTTTTTT: Potential genomic sequence that will be determined later
        '''
        # The first triming step is enough to process reads like
        # "TATCTCTTTTTTTTTTTTTTTTTGAAGGGCAGATTTAAAATACACTATTAAAATTATTAA"
        match1 = re.match('([ATCGN]{%d})(T*)' % randNT5, self.seq)
        T_length1 = len(match1.groups()[1])
        if T_length1 > 1:
            self.seq = self.seq[match1.end():]
            self.qual = self.qual[match1.end():]
            self.name = ''.join(['TS', 
                                 str(T_length1),
                                 match1.groups()[0], 
                                 ':', 
                                 self.name
								 ])
            self.trimmed5T = True

        # The second trimming step is needed to deal with reads like
        # "TATCTCTTTTATTTTTTTTTTTTGAAGGGCAGATTTAAAATACACTATTAAAATTATTAA"
        if T_length1 > 2:
            match2 = re.match('[ACGN](T+)', self.seq)
            if match2:
                T_length2 = len(match2.groups()[0])
                if T_length2 > 2:
                    self.seq = self.seq[match2.end():]
                    self.qual = self.qual[match2.end():]
                    T_length = T_length1 + 1 + T_length2
                    # Attach the trimmed sequence to the read name
                    # for comparison with genomic sequence later
                    self.name = ''.join(['TS', 
                                         str(T_length),
                                         match1.groups()[0], 
                                         '::',  
                                         match2.group(0),
                                         ':',
                                         self.name[self.name.find(':')+1:]
                                         ])

    def trim_3p_Ns(self, randNT3):
        '''Trim randNT3 NT from 3' end if necessary.

		Example:		

		@HISEQ01:507:HBC1DADXX:1:1101:1222:2105 1:N:0:GTCCGC cutadapt    
		AACTTTTGTTTTTTGACAGTCTCAAGTTTTTATTCAGTGGGTCTCTGTGTC
		
		becomes:
		
		@HISEQ01:507:HBC1DADXX:1:1101:1222:2105 1:N:0:GTCCGC <TGTC>
		AACTTTTGTTTTTTGACAGTCTCAAGTTTTTATTCAGTGGGTCTCTG
		'''
		# Trim randNT3 NT from 3' end if the 5' RNA ligation adapter has been 
		# removed
        if re.search('cutadapt$', self.name):
            self.name = self.name.replace('cutadapt', 
			                              '<' + self.seq[-randNT3:] + '>')
            self.seq = self.seq[:-randNT3]
            self.qual = self.qual[:-randNT3]
            self.trimmed3N = True


class FastqFile():
    '''Fastq file'''
    def __init__(self, filename, randNT5 = None, randNT3 = None):
        assert Path(filename).suffix == '.fastq', 'A *.fastq file is needed.'
        assert isinstance(randNT5, int) or randNT5 is None
        assert isinstance(randNT3, int) or randNT3 is None
        self.name = filename
        self.newname = re.sub('\.fastq$', '.trimmed.fastq', self.name)
        self.randNT5 = randNT5
        self.randNT3 = randNT3
        self.read_num = 0
        self.trimmed5T_num = 0

    def __repr__(self):
        return f"FastqFile('{self.name}', {self.randNT5}, {self.randNT3})"
    
    def __str__(self):
        output = (
            f"FastqFile Object with filename '{self.name}', number of 5' "
            f"random neucleotides: {self.randNT5}, number of 3' random "
            f"nucleotides: {self.randNT3}"
        )
        return output

    def __len__(self):
        return self.get_read_num()

    def get_name(self):
        return self.name

    def get_read_num(self):
        if self.read_num == 0:
            line_num = 0
            for _ in open(self.name): line_num += 1
            self.read_num = line_num//4      
        return self.read_num
    
    def get_fastq_record(self):
        '''Generator of FastqRecord objects from a fastq file'''
        i = 0
        name = None
        seq = None
        qual = None
        self.read_num = 0
        for line in open(self.name):
            i += 1
            curr_line = line.strip()
            if i % 4 == 1:
                name = curr_line[1:]
            elif i % 4 == 2:
                seq = curr_line
            elif i % 4 == 0:
                qual = curr_line
                yield FastqRecord(name, seq, qual)
                self.read_num += 1

    def create_trimmed_fastq_file(self):
        '''Trims 5' Ts of fastq records and write to a new file
        '''
        with open(self.newname, 'w') as fout:
            for fastq_record in self.get_fastq_record():
                fastq_record.trim_5p_Ts(self.randNT5)
                if fastq_record.trimmed5T: self.trimmed5T_num += 1
                fout.write(str(fastq_record))


def count_fastq(filename):
    '''A function for counting fastq record number'''
    return Path(filename).name, FastqFile(filename).get_read_num()


def trim_write_count_fastq(filename, randNT5, randNT3):
    '''A function for trimming and writing fastq records. 
    Also counts fastq records in input and output files.
    '''
    FF = FastqFile(filename, randNT5, randNT3)
    FF.create_trimmed_fastq_file()
    return Path(FF.name).name, FF.read_num, FF.trimmed5T_num


def load_fasta_genome(genome_dir):
    '''Read fasta files containing genomic sequence from genome_dir into a dict.
    '''
    fasta_files = Path(genome_dir).glob('ch*.fa')
    genome = {}
    for fasta_file in fasta_files:
        with open(fasta_file, 'r') as f:
            f.readline()
            line = f.read()
            genome[str(fasta_file.name).replace('.fa', '')] = line.replace('\n', '')
    return genome
    
       
def reverse_complement(seq, mol_type = 'DNA'):
    '''Get reverse complement of seq.'''
	# Genomic sequences contain both lower and upper case letters   
    seq = seq.upper()
    if mol_type == 'DNA':
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    elif mol_type == 'RNA':
        complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A', 'N': 'N'}
    return "".join(complement.get(base) for base in reversed(seq))
    

def get_seq(chromosome, strand, start, end, genome):
    '''Get stranded genomic sequence from 'genome' (a dict).'''
    seq = genome[chromosome][(start-1) : end].upper() 
    if strand == '-':
        return reverse_complement(seq, mol_type = 'DNA')
    else: 
        return seq
    
     
def pick_PASS(sam_file, min_mapq=10,  direction='reverse'):
    '''Split sam_file into PASS, noPASS, and yeast (spike-in) sam files.

    Loop through records in the sam file, if the record is PASS, write to
    sam_pass. Mapped non-PASS reads are saved in sam_nopass. Reads mapped to 
    rRNA genes are skipped. Sequencing direction can only be reverse or forward. 

    Attach matched tail length (ML), unmatched tail length (UL), and last mapped
    position (LM) to the end of readname. The head (such as 'TS11AAC') and the 
    tail (such as 'ML:i:4\tUL:i:9\tLM:i:1234') of the read name is the unique 
    identifier of the read.
    '''

    sam_pass = sam_file.replace('.sam', '.pass')
    sam_nopass = sam_file.replace('.sam', '.nopass')
    sam_ref = sam_file.replace('.sam', '.ref')

    sam_pass_file = open(sam_pass, 'w')
    sam_nopass_file = open(sam_nopass, 'w')
    sam_ref_file = open(sam_ref, 'w')
    lap = 0

    with open(sam_file, 'r') as in_file:
        for line in in_file:
            # Skip header
            if line[0] == '@':
                sam_pass_file.write(line)
                sam_nopass_file.write(line)
                sam_ref_file.write(line)
                continue
            # Process each line
            (readname, flag, chromosome, position, mapq, cigar) = line.split()[:6]
            position = int(position)
            mapq = int(mapq)

            # Discard reads with low mapping scores
            if mapq < min_mapq:
                continue
            # Record reads mapped to spiked-in yeast genome 
            if chromosome[:4] == 'ychr':  # or chromosome[:4] == 'BK000964':
                sam_ref_file.write(line)
                continue
            # Ignore reads mapped to other chromosomes, such as 'BK000964'
            if not chromosome[:3] == 'chr':
                continue

            # Extract the T-stretch length encoded in the readname
            match = re.match('TS(\d+)', line)
            if match:
                t_stretch_len = int(match.groups()[0])
            # Skip records with short T-stretches
            if not match or t_stretch_len < 2:
                sam_nopass_file.write(line)
                continue

            # Get genomic sequence downstream of the LAP (last mapped position).
            # (The length of returned sequence is t_stretch_len.)
            if (direction.lower() == 'reverse' and flag == '16') or \
                    (direction.lower() == 'forward' and flag == '0'):
                strand = '+'
                # Process cigar to determine the LAP.
                # Insertion (I) will not affect the covered distance.
                # if cigar = '38M10S' or '38', nums will be 38
                nums = re.split('[MDN]', re.sub('\d+I', '', cigar))[:-1]
                covered = sum(int(x) for x in nums)
                # get_seq() returns reverse complemented sequence if strand == '-'
                downstream_seq = get_seq(chromosome, strand,
                                         start=position + covered,
                                         end=position + covered + t_stretch_len - 1,
                                         genome=genome)
                lap = position + covered - 1
            elif (direction.lower() == 'reverse' and flag == '0') or \
                    (direction.lower() == 'forward' and flag == '16'):
                strand = '-'
                downstream_seq = get_seq(chromosome, strand,
                                         start=position - t_stretch_len, 
                                         end=position - 1,
                                         genome=genome)
                lap = position
            else:
                sam_nopass_file.write(line)
                continue

            # If the 5' T-stretch was trimmed twice (like 'TTTTTTTGTTTT'),
            # check whether GTTTT (reverse-compliment: AAAAC) comes from the genome
            elements = readname.split('::')
            if len(elements) == 2:
                elements[1] = elements[1].split(':')[0]   
                match = re.match(reverse_complement(elements[1]), downstream_seq)
                if match:
                    downstream_seq = downstream_seq[len(elements[1]):]
                    t_stretch_len -= len(elements[1])
                    if strand == '+':
                        lap += len(elements[1])
                    if strand == '-':
                        lap -= len(elements[1])
                    # If the t_stretch_len is reduced enough, no need to analyze downstream sequence
                    if t_stretch_len == 0:
                        sam_nopass_file.write(line.strip() + '\tML:i:%d\tUL:i:%d\tLM:i:%d\n'
                                              % (0, 0, lap))
                        continue

            # Analyze downstream sequence
            if not downstream_seq[:-1] == 'A' * (t_stretch_len - 1):
                match = re.match('A+', downstream_seq)
                matched_len = len(match.group()) if match else 0
                # ML: Matched Length 
                # UL: Unmatched Length
                # LM: Last Mapped position(UL) 
                sam_pass_file.write(line.strip() + '\tML:i:%d\tUL:i:%d\tLM:i:%d\n'
                                    % (matched_len, t_stretch_len - matched_len, lap))
            else:
                match = re.match('A*', downstream_seq)
                matched_len = len(match.group())  # value: 1 ~ t_stretch_len
                sam_nopass_file.write(line.strip() + '\tML:i:%d\tUL:i:%d\tLM:i:%d\n'
                                      % (matched_len, t_stretch_len - matched_len, lap))

    sam_pass_file.close()
    sam_nopass_file.close()
    sam_ref_file.close()
    os.system('rm ' + sam_file)

# List of sam files as input
sam_files = sorted(str(sam_file) for sam_file in sam_dir.glob('*Aligned.out.sam'))

# Genome loaded into memory is available for all processes
genome = load_fasta_genome(genome_dir)

# Try parallel computing, which requires more disk space
try:
    with mp.Pool(processes = WORKERS) as pool:
        pool.starmap(pick_PASS, [(sam_file, 10, 'reverse') for sam_file in sam_files])
except OSError:
    os.system('rm ' + sam_dir + '/*pass')
    print("No space left on device. Trying again with parallel computing turned off.")
    for sam_file in sam_files:
        print("Processing ", sam_file)
        pick_PASS(sam_file = sam_file)
# Release memory    
del genome      


def count_5Ts_in_folder(sam_dir):
    '''
    Count the number of reads with certain 5'T-stretch lengths for pass and 
    nopass reads in files in the sam_dir. Return a DataFrame.
    '''
    sample_names = []
    results = []  # container

    import re
    import numpy as np
    import pandas as pd

    for sam_file in [filename for filename in os.listdir(sam_dir)
                     if re.search('Aligned.out.(no)?pass$', filename)]:
        file_in = os.path.join(sam_dir, sam_file)
        sample_names.append(sam_file.split(
            '.')[0] + '.' + sam_file.split('.')[-1])

        with open(file_in, 'r') as f:
            string = f.read()
            # search patterns in the file
            TS = re.findall('\nTS(\d+)', string)  # 5' t-strech
            # count number of different lengths of T-stretch
            TS = [TS.count(str(i)) for i in range(51)]  # 50 sequencing cycles
            # save the result
            results.append(TS)

    # create DataFrame
    TS = np.array(results).T
    TS = np.insert(TS, 0, values=np.arange(51), axis=1)  # insert T length
    TS = pd.DataFrame(data=TS, columns=['T_Stretch_Length'] + sample_names)

    return TS

def pick_A_stretch(pass_in, astretch_out, non_astretch_out, min_length = 5):
    '''Pick PASS reads with at least min_length mapped 5' Ts in pass_in and copy 
    them into the astretch_out file. '''
    import re
    fout1 = open(astretch_out, 'w')
    fout2 = open(non_astretch_out, 'w')
    
    with open(pass_in, 'r') as fin:
        # process each line       
        for line in fin:
            # skip header
            if line[0] == '@': 
                fout1.write(line)
                fout2.write(line)
                continue
            # read the length of matched T-stretches
            ML = int(re.search('ML:i:(\d+)\t', line).group(1))
            if ML >= min_length:
                fout1.write(line)
            else:
                fout2.write(line)
                
    fout1.close()
    fout2.close()
    

def count_MAP_in_folder(pass_dir, result_dir, pattern = 'pass\.sam$',
                        outfile = 'map.scores.csv'):
    '''Count the number of reads with certain MAP scores for reads in sam 
    files with pattern in the filename in the pass_dir.'''

    results = {}
    
    from collections import defaultdict 
    import re
    
    for sam_file in [filename for filename in os.listdir(pass_dir) \
    if re.search(pattern, filename)]:
        MAPs = defaultdict(int)
        sam_in = os.path.join(pass_dir, sam_file)   
        print('Count MAP scores in ' + sam_in)
        with open(sam_in, 'r') as fin:
            for line in fin:
                if line[0] == '@' or int(line.split()[4]) < 10:
                    continue
                                
                MAPs[line.split()[4]] += 1
        results[sam_file] = MAPs
        
    import pandas as pd
    df = pd.DataFrame(results)   
    df = df/df.sum()
    df.to_csv(os.path.join(result_dir, outfile), index_label = 'MAP')
    
'''
#test 
count_MAP_in_folder(pass_dir = '../data/batch2', result_dir = '../result/batch3', pattern = 'fastq\.sam$',
                        outfile = 'test.map.scores.csv')
'''                        


def select_unique_fragments_2(input_file, random_NT_len):
    '''
    Use the TS\d+[ATCG]{3} string in read name and chromosome, flag, and LM 
    tags to identify potential PCR duplicates.
    '''
    output_file = input_file.replace('.pass', '.unique.pass')
    import itertools
    # make a set to save unique ids
    unique_ids = set()
    # precompile the search patter
    pattern = re.compile('TS(\d+[ATCGN]{%s})' % random_NT_len)
    fout = open(output_file, 'w')
    # open infile, calculate the id
    with open(input_file, 'r') as fin:
        #print('Identifying uPASS in ' + input_file)
        for line in fin:
            if line[0] == '@':
                fout.write(line)
                continue
            # print(line)

            # calculate the id for the read
            l = line.split()
            m = re.match(pattern, l[0])
            # calculate read length
            cigar = l[5]
            nums = re.split('[MDN]', re.sub('\d+I', '', cigar))[:-1]
            covered = sum(int(x) for x in nums)
            if m:
                l = [m.group(1), l[1:3], l[-1][:-1], str(covered)]
                # if 4N has been trimmed from 3' end
                cutadapt = re.search('4N([ATCG]{4})4N', l[0])
                if cutadapt:
                    l.append(cutadapt.groups[0])
                this_id = ''.join(list(itertools.chain(*l)))
                # only copy the line if this_id has not been seen before
                if not this_id in unique_ids:
                    unique_ids.add(this_id)
                    fout.write(line)
                    
def count_unique_reads_in_folder(sam_dir, outfolder, 
                          outfile = 'unique_pass_num.csv'):
    '''Count number of pass and ref_pass reads in sam files in the pass_dir'''
    from subprocess import check_output
    import glob
    
    sample_names = []
    pass_nums = [] 
    ref_nums = []
    
    # count pass reads 
    for file_name in sorted(glob.glob(os.path.join(sam_dir, '*.pass.unique.sam'))):
        sample_names.append(file_name.split('/')[-1].split('.')[0])
        cmd = 'grep -v @ ' + file_name + ' | wc -l'
        cmd_out = check_output(cmd, shell=True)
        pass_nums.append(cmd_out.strip())
    
    # count reference reads 
    for file_name in sorted(glob.glob(os.path.join(sam_dir, '*.ref.unique.sam'))):
        cmd = 'grep -v @ ' + file_name + ' | wc -l'
        cmd_out = check_output(cmd, shell=True)
        ref_nums.append(cmd_out.strip())
   
    # combine
    records = zip(sample_names, pass_nums, ref_nums)
    with open(os.path.join(outfolder, outfile), 'w') as f:
        f.write('Sample, uPASS, uRef\n')
        for rec in records:
            f.write('%s,%s,%s\n' % rec)
   
  
def find_neighboring_indexes(v, max_distance):
    """ 
    v is like [[100395423, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 56]],...],
    sorted by the position 100395423. 
    Yield the indexes of lists containing neighboring positions
    """
    # a container to hold list of liss
    index_list = []
    for i in range(len(v))[:-1]:
        if v[i + 1][0] - v[i][0] <= max_distance:
            #index_list.extend([i, i+1])
            index_list.append(i)
            continue
        if len(index_list) > 0:
            index_list.append(i)  # append the last neighbor
            yield index_list
            index_list = []


def cluster_neighboring_cleavage_sites(cs_cluster, max_distance):
    """
    cs_cluster is like:[position, [num in sam_files]] (Ordered by position)
    [[155603441, [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]],
     [155603444, [3, 1, 3, 14, 13, 8, 6, 1, 1, 10, 1, 6]],
     [155603445, [8, 1, 4, 13, 20, 8, 6, 0, 5, 19, 5, 19]]], with
    Neighboring positions located within max_distance from each other are merged. 
    Recursively merge the CS in the cluster: devide and conqurer!

    """
    #import itertools
    sample_number = len(cs_cluster[0][1]) // 2  # different from python 2.7
    # check base case
    clustered = True
    for i in range(len(cs_cluster))[:-1]:
        if cs_cluster[i + 1][0] - cs_cluster[i][0] <= max_distance:
            clustered = False
    if clustered == True:
        return True
    # when the clustering is not completed:
    else:
        # get max of normalized read numbers from all samples and the index for max
        # print(cs_cluster)
        (m, i) = max((v, i) for i, v in enumerate((sum(pos_num[1][sample_number:])
                                                   for pos_num in cs_cluster)))
        # merge cs within max_distance from the cs with max read number
        left_index = i  # leftmost position merged
        right_index = i  # rightmost position merged
        for j, (pos, nums) in enumerate(cs_cluster):
            if abs(pos - cs_cluster[i][0]) <= max_distance and \
                    abs(pos - cs_cluster[i][0]) > 0:
                # combine the read numbers for two CS, sample by sample
                # cs_cluster[i] is like
                # [155603471, [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]]
                # nums is like [0, 0, 0, 0, 0, 5, 1, 0, 0, 0, 0, 0]
                cs_cluster[i][1] = [sum(x)
                                    for x in zip(cs_cluster[i][1], nums)]
                cs_cluster[j][1] = 0  # don't delete cs_cluster[j]
                left_index = min(left_index, j)
                right_index = max(right_index, j)

        # devide and conqure
        if left_index > 0:
            cluster_neighboring_cleavage_sites(cs_cluster[:left_index],
                                               max_distance)
        if right_index < len(cs_cluster) - 1:
            cluster_neighboring_cleavage_sites(cs_cluster[right_index + 1:],
                                               max_distance)


def cluster_reads_in_sam_dir(file_pattern='pass.unique.sam',
                             outfile='cluster.numbers.csv',
                             direction='reverse',
                             max_distance=24):
    """
    Analyze a sam file and first generate a table with the following columns: 
    chromosome, strand, pos, num. Then recursively combine reads with LAP within 
    24 nt, from pos with the highest read number to the pos with next highest 
    read number.

    Read sam files with file_pattern in infolder. For each sam file, build a dict
    readcounts like {'chr9:-:100395423':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 56]}.
    The [list of int] saves the number of reads from each sam_file. 
    dict is like sample_count[sample1(str)] = num_from_sample1(int)

    """

    import os
    sam_files = [file_name for file_name in os.listdir(sam_dir)
                 if file_name.endswith(file_pattern)]
    sam_file_num = len(sam_files)

    import re
    # regular expression pattern for finding last mapped position (LM) in read names
    re_pattern = re.compile('LM:i:(\d+)')

    import collections
    # readcounts is a read id counter
    # example: {'chr9:-:100395423':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 56]}
    readcounts = {}
    # go through sam files and get read numbers for each read id from each sam file
    for s, sam_file in enumerate(sam_files):
        # read sam file, calculate read id, and count number of reads
        # s will be used as index
        fin = open(os.path.join(sam_dir, sam_file), 'r')
        for line in fin:
            if line[0] == '@':
                continue
            elements = line.split()
            chromosome = elements[2]
            flag = elements[1]
            # calculate strand
            if (direction == 'reverse' and flag == '16') or \
                    (direction == 'forward' and flag == '0'):
                strand = '+'
            elif (direction == 'reverse' and flag == '0') or \
                    (direction == 'forward' and flag == '16'):
                strand = '-'

            position = re.search(re_pattern, elements[-1]).group(1)

            read_id = ':'.join([chromosome, strand, position])
            # initialize the list for holding number of reads for each file
            readcounts.setdefault(read_id, [0] * sam_file_num)[s] += 1
            """
            # when s is 0 and sam_file[s] is the first one analyzed,
            # readcount.popitem() may return something like 
            # 'chr7:-:97599217', [3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            # After analyzing the second files in the folder,
            # readcount.popitem() may return something like:
            # ('chr11:+:31270274', [6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
            """

        fin.close()

    print('Clustering reads...')

    # calculate the normalized read numbers and attach to the read numbers
    import pandas as pd
    df = pd.DataFrame(readcounts)
    df = df.T
    normalized = df / df.sum(0)
    df = pd.concat([df, normalized], axis=1, ignore_index=True)
    readcounts = df.T.to_dict('list')
    # readcount.popitem() may return something like:
    #  ('chr11:+:31270274', [6.0, 7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    #  0.0, 0.0, 9.1427902698768829e-07,3.0861678694165542e-06, 0.0, 0.0 ...]),
    # with the first 12 values as read counts and the next 12 values as normalized
    # read counts

    # separate positions based on chrmosome & strand combination
    poscounts = collections.defaultdict(list)
    while len(readcounts) > 0:
        k, v = readcounts.popitem()  # saves memory
        k = k.split(':')
        poscounts[':'.join(k[:2])].append(
            [int(k[2]), v])  # don't forget int()!!!
        # poscounts.popitem() is like
        # {'chr9:-':[[100395423, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 56]],...]}
    del readcounts

    # sort the list of lists for each chromosome & strand combination
    for k, v in poscounts.items():
        # sort in place to save memory
        v.sort(key=lambda val: val[0])
        # get the indexes of lists containing neighboring positions
        for indexes in find_neighboring_indexes(v, max_distance):
            # print(indexes)
            # extract the cluster
            cs_cluster = v[indexes[0]:(indexes[-1] + 1)]
            # the original list in the dict will be edited in place:
            cluster_neighboring_cleavage_sites(cs_cluster, max_distance)
        # delete positions with 0 read number after clustering
        poscounts[k] = [posi_num for posi_num in v if not posi_num[1] == 0]

    # write the result to disk in csv format
    print('Writing to file...')
    outfile = os.path.join(result_dir, outfile)
    with open(outfile, 'w') as fout:
        # write header
        sample_string = ','.join([sam_file.split('.')[0]
                                  for sam_file in sam_files])
        fout.write('chromosome,strand,position,%s\n' % sample_string)
        # write read counts for each cluster
        for k in poscounts:
            chromosome, strand = k.split(':')
            for record in poscounts[k]:
                position = str(record[0])
                counts = ','.join([str(int(count)) for count
                                   in record[1][:sam_file_num]])
                fout.write('%s,%s,%s,%s\n' %
                           (chromosome, strand, position, counts))
               
def cluster_CS_in_dirs(infolders, outfolder, 
                         cs_file_name = 'CS.all.reads.csv', 
                         #direction = 'reverse', 
                         max_distance = 24,
                         outfile = 'meta.cluster.numbers.csv'):
    '''
    Read "CS.all.reads.csv" files (containing read numbers for each unclustered 
    cleavage site) in input folders. From each CS.all.reads.csv file, build a dict
    named readcounts like {'chr9:-:100395423':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 56]}.
    The [list of int] saves the number of reads from each CS.all.reads.csv file. 
    dict is like sample_count[sample1(str)] = num_from_sample1(int)
    Output a table with the following columns: 
    chromosome, strand, pos, num. Then recursively combine reads with LAP within 24 nt, 
    from pos with the highest read number to the pos with next highest read 
    number.''' 
    '''
    rootfolder = '../../result'
    direction = 'reverse'
    max_distance = 24
    
    import glob
    for sam_in in sorted(glob.glob(os.path.join(pass_dir, '*pass.sam'))):
        print(sam_in)
    
    '''
    #import glob
#    infolders = ['../../result/batch05', '../../result/batch06', 
#                 '../../result/batch09', '../../result/batch13']
            
    #cs_files = [glob.glob(os.path.join(infolder, cs_file_name)) for infolder in infolders]
    cs_files = [os.path.join(infolder, cs_file_name) for infolder in infolders]
    # calculate the total number of columns and get sample names
    num_column = []
    sample_names = []
    for cs_file in cs_files:
        with open(cs_file, 'r') as fin:
            samples = fin.readline().strip().split(',')[3:]
            num_column.append(len(samples))
            sample_names += samples
    total_num_column = sum(num_column)        
    
    
    import collections 
    # readcounts is a read id counter
    #example: {'chr9:-:100395423':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 56]}
    readcounts = {} 
    # loop through cleavage site (CS) files and get read numbers for each CS id   
    for s, cs_file in enumerate(cs_files): #s will be used as index
        # read sam file, calculate read id, and count number of reads
        print('Reading ' + cs_file)        
        fin = open(cs_file, 'r')
        for line in fin:
            elements = line.strip().split(',')
            if elements[0] == 'chromosome': # header
                continue
            chromosome, strand, position = elements[:3]
            # calculate cleavage site id
            cs_id = ':'.join([chromosome, strand, position]) 
            # setdefault(key[, default])
            # If key is in the dictionary, return its value. 
            # If not, insert key with a value of default and return default. 
            readcounts.setdefault(cs_id, [0]*total_num_column)
            # Copy the read numbers to the right columns
            readcounts[cs_id][sum(num_column[:s]):sum(num_column[:s+1])] = [int(x) for x in elements[3:]]
        fin.close()
    
    # calculate the normalized read numbers and attach to the read numbers  
    print('Calculating RPMs...')
    import pandas as pd
    df = pd.DataFrame(readcounts)    
    df = df.T
    normalized = df/df.sum(0)
    df = pd.concat([df, normalized], axis = 1, ignore_index = True)
    readcounts = df.T.to_dict('list')    
    # readcount.popitem() may return something like:
    #  ('chr11:+:31270274', [6.0, 7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
    #  0.0, 0.0, 9.1427902698768829e-07,3.0861678694165542e-06, 0.0, 0.0 ...]),
    # with the first 12 values as read counts and the next 12 values as normalized
    # read counts
    
    print( 'Clustering cleavage sites...')
    # separate positions based on chrmosome & strand combination
    poscounts = collections.defaultdict(list)
    while len(readcounts) > 0:
        k,v = readcounts.popitem() # saves memory
        k = k.split(':')
        poscounts[':'.join(k[:2])].append([int(k[2]),v]) # don't forget int()!!!
        # poscounts.popitem()        
        #example: {'chr9:-':[[100395423, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 56]],...]}
    del readcounts
    
    # sort the list of lists for each chromosome & strand combination    
    for k,v in poscounts.items(): 
        # sort in place to save memory
        v.sort(key = lambda val: val[0]) 
        # get the indexes of lists containing neighboring positions
        for indexes in find_neighboring_indexes(v, max_distance):
            #print(indexes)            
            # extract the cluster
            cs_cluster = v[indexes[0]:(indexes[-1]+1)]
            # the original list in the dict will be edited in place:
            cluster_neighboring_cleavage_sites(cs_cluster, max_distance)
        # delete positions with 0 read number after clustering
        poscounts[k] = [posi_num for posi_num in v if not posi_num[1] == 0]
    
    # write the result to disk in csv format   
    print( 'Writing to file...')
    with open(os.path.join(outfolder, outfile), 'w') as fout:
        # write header        
        #sample_string = ','.join([sam_file + '.num' for sam_file in sam_files])
        sample_string = ','.join(sample_names)
        fout.write('chromosome,strand,position,%s\n' %sample_string)
        # write read counts for each cluster
        for k in poscounts:
            chromosome, strand = k.split(':')
            for record in poscounts[k]:
                position = str(record[0])
                counts = ','.join([str(int(count)) for count in record[1][:total_num_column]])
                fout.write('%s,%s,%s,%s\n'%(chromosome, strand, position, counts))
                
        
# the sam2bigwid() function cannot be defined within make_url().
# functions are only picklable if they are defined at the top-level of a module.
def sam2bigwig(sam_file):
    prefix = sam_file.split('.')[0]
    sam_file = os.path.join(sam_dir, sam_file)
    prefix = os.path.join(sam_dir, prefix)

    # sam -> bam
    cmd = samtools + ' view -uS ' + sam_file + ' | ' + samtools + \
        ' sort - ' + prefix
    # print(cmd)
    os.system(cmd)

    # bam -> bedGraph
    totalReadNum = count_pass(sam_file)[1]
    cmd = genomeCoverageBed + '-bg -split -ibam ' + prefix + '.bam ' + \
        '-strand + -g ' + genome_size + ' -scale ' + str(-10**6/totalReadNum) + \
        ' > ' + prefix + '.m.bedgraph'
    # print(cmd)
    os.system(cmd)
    cmd = genomeCoverageBed + '-bg -split -ibam ' + prefix + '.bam ' + \
        '-strand - -g ' + genome_size + ' -scale ' + str(10**6/totalReadNum) + \
        ' > ' + prefix + '.p.bedgraph'
    #print(cmd)
    os.system(cmd)

    # bedgraph -> bigWig
    cmd = bedGraphToBigWig + prefix + '.p.bedgraph ' + \
        genome_size + ' ' + prefix + '.p.bw'
    # print(cmd)
    os.system(cmd)
    cmd = bedGraphToBigWig + prefix + '.m.bedgraph ' + \
        genome_size + ' ' + prefix + '.m.bw'
    # print(cmd)
    os.system(cmd)

    # remove intermediate files
    cmd = 'rm ' + prefix + '*.bam'
    os.system(cmd)
    cmd = 'rm ' + prefix + '*bedgraph*'
    os.system(cmd)


# sam -> bam -> bigwig
def make_url(project, experiment, sam_dir, genome_size, genomeCoverageBed,
             bedGraphToBigWig):
    from subprocess import check_output

    # parallel computing
    sam_files = sorted(glob.glob(os.path.join(sam_dir, '*.Aligned.out.pass')))
    with mp.Pool(processes=processes) as pool:
        pool.map(sam2bigwig, sam_files)

    # create UCSC genome browser tracks
    bw_files = [filename for filename in os.listdir(sam_dir)
                if filename.endswith('.bw')]
    bw_samples = sorted(list(set([filename.split('.')[0]
                                  for filename in bw_files])))

    strand2str = {'+': 'p', '-': 'm'}
    f = open(os.path.join(sam_dir, 'bigwigCaller.txt'), 'w')
    for strand in ['+', '-']:
        for (i, sample) in enumerate(bw_samples):
            name = sample
            color = list(colors)[int(sample_color[sample]) - 1]
            track = 'track type=bigWig visibility=2 alwaysZero=on color=' + \
                color + ' graphType=bar maxHeightPixels=30:30:30 itemRgb=On group=' + \
                project + ' name="' + name + strand + '" description="' + name + ' ' + \
                strand + '" bigDataUrl=http://intron.njms.rutgers.edu/zhengdh/bigwig/' + \
                project + '/' + experiment + '/' + sample + '.' + strand2str[strand] + \
                '.bw' + '\n'
            f.write(track)
    f.close()




# sam -> bam -> bigwig
def make_CLIP_url(project, batch, sam_dir, result_dir, genome_size, 
             genomeCoverageBed, norm_bedgraph, bedGraphToBigWig):
    for sam_file in [filename for filename in os.listdir(sam_dir) 
    if re.search('\.sam$', filename)]:
        prefix = sam_file.split('.')[0]    
        sam_file = os.path.join(sam_dir, sam_file)
        #prefix = re.sub('\.sam$', '', sam_file)
        prefix = os.path.join(sam_dir, prefix)
        cmd = 'samtools view -uS ' + sam_file + ' | ' + 'samtools sort - ' + prefix # why need "-"?
        print( cmd)
        os.system(cmd)
        # I got the warning: [bam_header_read] EOF marker is absent. The input is probably truncated.
        # The warning should be ignored. Using samtools view -u gives uncompressed bam, and these do not have the EOF marker. 
    
        
        # bam to bigwig 
        cmd = genomeCoverageBed + '-bg -split -ibam ' + prefix + '.bam ' + \
        '-strand + -g ' + genome_size + ' > ' + prefix + '.p.bedgraph'  # + strand: "p", different from 3'READS+, which is reverse sequenced
        print( cmd)
        os.system(cmd)
        cmd = genomeCoverageBed + '-bg -split -ibam ' + prefix + '.bam ' + \
        '-strand - -g ' + genome_size + ' > ' + prefix + '.m.bedgraph' # - strand: "m", different from 3'READS+, which is reverse sequenced
        print( cmd)
        os.system(cmd)
        
        # normolize bedgraph counts     
        cmd = 'grep -v @ ' + sam_file + ' | wc -l'
        from subprocess import check_output
        cmd_out = check_output(cmd, shell=True)
        totalReadNum = int(cmd_out.split(' ')[0]) 
        cmd = norm_bedgraph + '-t ' + str(totalReadNum) + ' -i "' + prefix + \
        '.p.bedgraph ' + prefix + '.m.bedgraph"' #+ ' -m "' + prefix + '.p.bedgraph"'
        print( cmd)
        os.system(cmd)
                
        # convert to bw
        cmd = bedGraphToBigWig + prefix + '.p.bedgraph.normolized ' + \
        genome_size + ' ' + prefix + '.p.bw'
        print( cmd)
        os.system(cmd)
        cmd = bedGraphToBigWig + prefix + '.m.bedgraph.normolized ' + \
        genome_size + ' ' + prefix + '.m.bw'
        print( cmd)
        os.system(cmd)
        
        # remove intermediate files
        #cmd = 'rm ' + prefix + '*.bam'
        #os.system(cmd)
        cmd = 'rm ' + prefix + '*bedgraph*'
        os.system(cmd)
    # copy the bigwig files to the http-enabled intron server
    #cmd = 'scp ' + os.path.join(pass_dir, '*.bw') + ' zhengdh@intron.njms.rutgers.edu:~/../www/zhengdh/bigwig/' + project 
    #os.system(cmd)
    
    # record sample names for the bw files
    bw_files = [filename for filename in os.listdir(sam_dir) 
    if filename.endswith('.bw')]
    bw_samples = sorted(list(set([filename.split('.')[0] for filename in bw_files])))
    #bw_samples.sort(key = lambda x: x[-1])
    
    #cmd = 'rm ' + os.path.join(pass_dir, '*.bw')
    
    # 
    
    f = open(os.path.join(result_dir, 'bigwigCaller.txt'), 'w')
    for strand in ['+', '-']:
        for (i, sample) in enumerate(bw_samples):
            #name = sample.replace('spikein', '') + '_' + batch
            name = sample
            # default:
            strand2str = {'+': 'p', '-':'m'}
            if 'R_' in name: # AdiPR_HuR, meaning that the library is reverse sequenced
                strand2str = {'+': 'm', '-':'p'}
            #print( i, sample)
            color = '000,000,255'
            #if sample.startswith('AS'):
            if re.search('AS', sample):
                color = '255,000,000'
            elif re.search('HP', sample):
                color = '255, 100,000'
                
            track = 'track type=bigWig visibility=2 alwaysZero=on color=' + \
            color +' graphType=bar maxHeightPixels=30:30:30 itemRgb=On group=' + \
            project + ' name="' + name \
             + strand + '" description=\"CLIP"' + \
             ' bigDataUrl=http://intron.njms.rutgers.edu/zhengdh/bigwig/' + \
             project + '/' + batch + '/' + sample + '.' + strand2str[strand] + '.bw' + '\n'
            f.write(track)
    f.close()     

# sam -> bam -> bigwig
def make_RNAseq_url(project, batch, sam_dir, result_dir, genome_size, 
             genomeCoverageBed, norm_bedgraph, bedGraphToBigWig):
    for sam_file in [filename for filename in os.listdir(sam_dir) 
    if re.search('\.sam$', filename)]:
        prefix = sam_file.split('.')[0]    
        sam_file = os.path.join(sam_dir, sam_file)
        #prefix = re.sub('\.sam$', '', sam_file)
        prefix = os.path.join(sam_dir, prefix)
        cmd = 'samtools view -uS ' + sam_file + ' | ' + 'samtools sort - ' + prefix # why need "-"?
        print( cmd)
        os.system(cmd)
        # I got the warning: [bam_header_read] EOF marker is absent. The input is probably truncated.
        # The warning should be ignored. Using samtools view -u gives uncompressed bam, and these do not have the EOF marker. 
    
        
        # bam to bigwig 
        cmd = genomeCoverageBed + '-bg -split -ibam ' + prefix + '.bam ' + \
        '-strand + -g ' + genome_size + ' > ' + prefix + '.p.bedgraph'  # + strand: "p", different from 3'READS+, which is reverse sequenced
        print( cmd)
        os.system(cmd)
        cmd = genomeCoverageBed + '-bg -split -ibam ' + prefix + '.bam ' + \
        '-strand - -g ' + genome_size + ' > ' + prefix + '.m.bedgraph' # - strand: "m", different from 3'READS+, which is reverse sequenced
        print( cmd)
        os.system(cmd)
        
        # normolize bedgraph counts     
        cmd = 'grep -v @ ' + sam_file + ' | wc -l'
        from subprocess import check_output
        cmd_out = check_output(cmd, shell=True)
        totalReadNum = int(cmd_out.split(' ')[0]) 
        cmd = norm_bedgraph + '-t ' + str(totalReadNum) + ' -i "' + prefix + \
        '.p.bedgraph ' + prefix + '.m.bedgraph"' #+ ' -m "' + prefix + '.p.bedgraph"'
        print( cmd)
        os.system(cmd)
                
        # convert to bw
        cmd = bedGraphToBigWig + prefix + '.p.bedgraph.normolized ' + \
        genome_size + ' ' + prefix + '.p.bw'
        print( cmd)
        os.system(cmd)
        cmd = bedGraphToBigWig + prefix + '.m.bedgraph.normolized ' + \
        genome_size + ' ' + prefix + '.m.bw'
        print( cmd)
        os.system(cmd)
        
        # remove intermediate files
        #cmd = 'rm ' + prefix + '*.bam'
        #os.system(cmd)
        cmd = 'rm ' + prefix + '*bedgraph*'
        os.system(cmd)
    # copy the bigwig files to the http-enabled intron server
    #cmd = 'scp ' + os.path.join(pass_dir, '*.bw') + ' zhengdh@intron.njms.rutgers.edu:~/../www/zhengdh/bigwig/' + project 
    #os.system(cmd)
    
    # record sample names for the bw files
    bw_files = [filename for filename in os.listdir(sam_dir) 
    if filename.endswith('.bw')]
    bw_samples = sorted(list(set([filename.split('.')[0] for filename in bw_files])))
    #bw_samples.sort(key = lambda x: x[-1])
    
    #cmd = 'rm ' + os.path.join(pass_dir, '*.bw')
    
    # 
    
    f = open(os.path.join(result_dir, 'bigwigCaller.txt'), 'w')
    for strand in ['+', '-']:
        for (i, sample) in enumerate(bw_samples):
            #name = sample.replace('spikein', '') + '_' + batch
            name = sample
            # default:
            #strand2str = {'+': 'p', '-':'m'} # forward sequencing
            strand2str = {'+': 'm', '-':'p'} # reverse sequencing
#            if 'R_' in name: # AdiPR_HuR, meaning that the library is reverse sequenced
#                strand2str = {'+': 'm', '-':'p'}
            #print( i, sample)
            color = '000,000,255'
            #if sample.startswith('AS'):
            if re.search('AS', sample):
                color = '255,000,000'
            elif re.search('HP', sample):
                color = '255,100,000'
                
            track = 'track type=bigWig visibility=2 alwaysZero=on color=' + \
            color +' graphType=bar maxHeightPixels=30:30:30 itemRgb=On group=' + \
            project + ' name="' + name \
             + strand + '" description=\"CLIP"' + \
             ' bigDataUrl=http://intron.njms.rutgers.edu/zhengdh/bigwig/' + \
             project + '/' + batch + '/' + sample + '.' + strand2str[strand] + '.bw' + '\n'
            f.write(track)
    f.close()  
    



def prcmd(cmd):
    print ("Command to run:", cmd, "\n")        
    os.system(cmd) 
