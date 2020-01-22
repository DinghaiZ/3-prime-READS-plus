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
import collections
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
import multiprocessing as mp


class FastqRecord():
    '''Fastq record'''
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
        self.seq = self.seq[match1.end():]
        self.qual = self.qual[match1.end():]
        self.name = ''.join(['TS', str(T_length1), match1.groups()[0], ':', self.name])

        if T_length1 > 0: self.trimmed5T = True

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
                    self.name = ''.join(['TS', str(T_length), match1.groups()[0], '::',  
                                         match2.group(0), ':',
                                         self.name[self.name.find(':')+1:]
                                         ])

    def trim_3p_Ns(self, randNT3):
        '''Trim randNT3 NT from 3' end.

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
            self.name = self.name.replace('cutadapt', '<' + self.seq[-randNT3:] + '>')
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
        if not self.read_num:
            line_num = 0
            for _ in open(self.name): line_num += 1
            self.read_num = line_num//4      
        return self.read_num
    
    def get_fastq_record(self):
        '''Generator of FastqRecord objects'''
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
        '''Trims 5' Ts of fastq records and write to a new file'''
        with open(self.newname, 'w') as fout:
            for fastq_record in self.get_fastq_record():
                fastq_record.trim_5p_Ts(self.randNT5)
                if fastq_record.trimmed5T: self.trimmed5T_num += 1
                if len(fastq_record) >= 18: fout.write(str(fastq_record))


def merge_and_rename(selected_fastq_files, output_file, rawfastq_dir):
    '''Unzip (if necessary), merge, and rename a list of selected fastq files. 
    
    Arguments:
    selected_fastq_files: a list of selected fastq files. Must all be in either .fastq or .fastq.gz format.
    output_file: output fastq file name.
    '''    
    if len(selected_fastq_files) > 0:
        joint_file_names = ' '.join(selected_fastq_files)
        if selected_fastq_files[0].endswith('.gz'):
            print('\nMerging, unzipping, and renaming fastq files ....')
            cmd = f'zcat {joint_file_names} > {str(rawfastq_dir)}/{output_file}'
        else:
            print('\nMerging and renaming fastq files ....')
            cmd = f'cat {joint_file_names} >> {str(rawfastq_dir)}/{output_file}'
        print(cmd)
        os.system(cmd)
        # Remove downloaded raw fastq files immediately to save space 
        cmd = 'rm ' + joint_file_names
        os.system(cmd)

def count_fastq(filename):
    '''A function for counting fastq record number'''
    return Path(filename).name, FastqFile(filename).get_read_num()


def trim_write_count_fastq(filename, randNT5, randNT3):
    '''A function for trimming and writing fastq records. Also counts fastq 
    records in input and output files.
    '''
    ff = FastqFile(filename, randNT5, randNT3)
    ff.create_trimmed_fastq_file()
    return Path(ff.name).name, ff.read_num, ff.trimmed5T_num


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
    
     
def split_sam(sam_file, min_mapq = 10, direction = 'reverse', spike_in = None):
    '''Split a sam_file into PASS, nonpass, and spike-in sam files.

    Loop through records in the sam file, if the record is PASS, write to
    pass_file. Mapped non-PASS reads are saved in nonpass_file. Reads mapped to 
    rRNA genes are skipped. Sequencing direction can only be reverse or forward. 

    Attach matched tail length (ML), unmatched tail length (UL), and last mapped
    position (LM) to the end of readname. The head (such as 'TS11AAC') and the 
    tail (such as 'ML:i:4\tUL:i:9\tLM:i:1234') of the read name is the unique 
    identifier of the read.
    
    Arguments:
    sam_file: string. Name of input sam file.
    min_mapq: int. Minimum MAPQ score
    direction: string. Either 'reverse' (default) or 'forward'
    spike_in: a 4-charactor string for identifying spike-in RNA. For example,
    spike_in = 'ychr' for spike-in yeast RNA if yeast chromosomes are named as 'ychr1', 
    'ychrM', etc. 
    
    Outputs:
    Separate sam files containing either PASS reads, nonPASS reads, or spike-ins.
    '''
    pass_file = open(sam_file.replace('.sam', '.pass'), 'w')
    nonpass_file = open(sam_file.replace('.sam', '.nonpass'), 'w')
    if spike_in: spike_in_file = open(sam_file.replace('.sam', '.spike_in'), 'w')
    
    lap = 0

    with open(sam_file, 'r') as in_file:
        for line in in_file:
            # Skip header
            if line[0] == '@':
                pass_file.write(line)
                nonpass_file.write(line)
                if spike_in: spike_in_file.write(line)
                continue
            # Process each line
            (readname, flag, chromosome, position, mapq, cigar) = line.split()[:6]
            position = int(position)
            mapq = int(mapq)

            # Discard reads with low mapping scores
            if mapq < min_mapq:
                continue
            # Record reads mapped to spiked-in yeast genome 
            if spike_in and chromosome[:4] == spike_in:  
                spike_in_file.write(line)
                continue
            # Ignore reads mapped to other chromosomes, such as 'BK000964'
            if not chromosome[:3] == 'chr':
                continue

            # Extract the T-stretch length encoded in the readname
            match = re.match('TS(\d+)', line)
            t_stretch_len = int(match.groups()[0])
            # Skip records with short T-stretches
            if t_stretch_len < 2:
                nonpass_file.write(line)
                continue

            # Get genomic sequence downstream of the LAP (last mapped position).
            # (The length of returned sequence is t_stretch_len.)
            if (direction.lower() == 'reverse' and flag == '16') or \
                    (direction.lower() == 'forward' and flag == '0'):
                strand = '+'
                # Process cigar to determine the LAP.
                # Insertion (I) will not affect the covered distance.
                nums = re.split('[MDN]', re.sub('\d+I', '', cigar))[:-1]
                covered = sum(int(x) for x in nums)
                # get_seq() returns reverse complemented sequence if strand == '-'
                downstream_seq = get_seq(chromosome, strand,
                                         start = position + covered,
                                         end = position + covered + t_stretch_len - 1,
                                         genome = genome)
                lap = position + covered - 1
            elif (direction.lower() == 'reverse' and flag == '0') or \
                    (direction.lower() == 'forward' and flag == '16'):
                strand = '-'
                downstream_seq = get_seq(chromosome, strand,
                                         start = position - t_stretch_len, 
                                         end = position - 1,
                                         genome = genome)
                lap = position
            else:
                continue

            # If the 5' T-stretch was trimmed twice (like 'TTTTTTTGTTTT')
            # and if GTTTT (reverse-compliment: AAAAC) comes from the genome,
            # update the downstream_seq, lap, and t_stretch_len
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
                        nonpass_file.write(line.strip() + '\tML:i:%d\tUL:i:%d\tLM:i:%d\n'
                                              % (0, 0, lap))
                        continue

            # Check if downstream sequence matches trimmed T-stretches 
            if not downstream_seq[:-1] == 'A' * (t_stretch_len - 1):
                match = re.match('A+', downstream_seq)
                matched_len = len(match.group()) if match else 0
                # ML: Matched Length 
                # UL: Unmatched Length
                # LM: Last Mapped position(UL) 
                pass_file.write(line.strip() + '\tML:i:%d\tUL:i:%d\tLM:i:%d\n'
                                    % (matched_len, t_stretch_len - matched_len, lap))
            else:
                match = re.match('A*', downstream_seq)
                matched_len = len(match.group())  # value: t_stretch_len-1 or t_stretch_len
                nonpass_file.write(line.strip() + '\tML:i:%d\tUL:i:%d\tLM:i:%d\n'
                                      % (matched_len, t_stretch_len - matched_len, lap))
    pass_file.close()
    nonpass_file.close()
    if spike_in: spike_in_file.close()
    os.system('rm ' + sam_file)


def pick_unique_pass(pass_file, random_NT_len):
    '''Use the TS\d+[ATCG]{random_NT_len} string in read name and chromosome, flag,
    and LM tags to remove potential PCR duplicates. 
    
    Arguments:
    pass_file: A *.pass file generated by pick_PASS() or split_sam()
    random_NT_len: Length of random nucleotide sequence in the 3' ligation 
    adapter for reverse direction sequencing.

    Output:
    Unique PASS reads will be written to a new file.
    '''
    output_file = pass_file.replace('.pass', '.unique.pass')
    fout = open(output_file, 'w')
    # Make a set to save unique ids
    unique_ids = set()
    # Precompile the search patter
    pattern = re.compile('TS(\d+[ATCGN]{%s})' % random_NT_len)
    with open(pass_file, 'r') as fin:
        for line in fin:
            if line[0] == '@':
                fout.write(line)
                continue
            # Calculate the unique identifier for each PASS read
            l = line.split()
            m = re.match(pattern, l[0])
            # Calculate read length
            cigar = l[5]
            nums = re.split('[MDN]', re.sub('\d+I', '', cigar))[:-1]
            covered = sum(int(x) for x in nums)
            if m:
                l = [m.group(1), l[1:3], l[-1][:-1], str(covered)]
                # if 4N has been trimmed from 3' end
                #cutadapt = re.search('4N([ATCG]{4})4N', l[0])
                #if cutadapt:
                #    l.append(cutadapt.groups[0])
                this_id = ''.join(list(itertools.chain(*l)))
                # Copy the line to new file if this_id has not been seen before
                if not this_id in unique_ids:
                    unique_ids.add(this_id)
                    fout.write(line)


def count_sam(sam_file):
    '''Count number or alignment records in a sam_file.'''
    line_num = 0
    for line in open(sam_file): 
        if not line[0] == '@': line_num += 1
    sample_name = Path(sam_file).name.replace('.Aligned.out.pass.sam', '')
    return sample_name, line_num 


def count_bam(bam_file):
    cmd = f'samtools view -c {bam_file}'
    count = int(check_output(cmd, shell=True).decode().strip())
    sample_name = Path(bam_file).name.split('.')[0]
    return sample_name, count


def count_5T_stretch(sam_file, max_TS = 25):
    with open(sam_file, 'r') as fin:
        fstring = fin.read()
        # Search 5' T-strech patterns in the file
        TS = re.findall('\nTS(\d+)', fstring)  
        # Count number of different lengths (up to 30) of T-stretch
        return [TS.count(str(i)) for i in range(max_TS + 1)]  


def summarize_5T_stretch(sam_files, processes, max_TS = 25):
    '''Calculate the distribution of un-aligned 5'T-stretch lengths for records 
    in sam_files, using multiple threads. 
    
    Arguments:
    
    sam_files: a list of pass file names or nonpass file names.
    max_TS: maximum T-stretch length that you care about
    processes: number of threads for parallel computation
    
    Output:
    A DataFrame with T-stretch length as row index, sample names as column names, 
    and values as number of counts of each T-stretch length for each sample. 
    '''
    sample_names = [re.search('.+\/(.+?)\.Aligned.+pass', sam_file).groups()[0] 
                    for sam_file in sam_files]
    # Use multi-processes to count T-stretch lengthes in multiple files
    with mp.Pool(processes = processes) as pool: 
        TS_counts = pool.starmap(count_5T_stretch, 
                                 zip(sam_files, [max_TS]*len(sam_files)))
    # Create DataFrame
    df = np.array(TS_counts).T
    df = pd.DataFrame(data = df, columns = sample_names)
    df.index.name = 'T_Stretch_Length'
    return df


def summarize_5T_stretch_serial(sam_files, max_TS = 30):
    '''Calculate the distribution of un-aligned 5'T-stretch lengths for records 
    in sam_files, using just one thread. 
    
    Arguments:
    'sam_files' is a list of pass file names or nonpass file names.
    'max_TS': maximum T-stretch length that you care about
    Returns the result as a DataFrame.
    '''
    results = []  
    sample_names = [re.search('.+\/(.+?)\.Aligned.+pass', sam_file).groups()[0] 
                    for sam_file in sam_files]
    for sam_file in sam_files:
        with open(sam_file, 'r') as fin:
            fstring = fin.read()
            # Search 5' T-strech patterns in the file
            TS = re.findall('\nTS(\d+)', fstring)  
            # Count number of different lengths (up to 30) of T-stretch
            TS = [TS.count(str(i)) for i in range(max_TS + 1)]  
            results.append(TS)
    # Create DataFrame
    df = np.array(results).T
    df = pd.DataFrame(data = df, columns = sample_names)
    df.index.name = 'T_Stretch_Length'
    return df


def find_nearby_indexes(v, max_distance):
    """ 
    Yield the indexes of lists containing neighboring positions
    
    Arguments:
    v is like [[100395423, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 56]],...],
    sorted by the position (100395423, etc). 

    Output:
    """
    # A container to hold list of liss
    index_list = []
    for i in range(len(v))[:-1]:
        if v[i + 1][0] - v[i][0] <= max_distance:
            index_list.append(i)
            continue
        if len(index_list) > 0:
            index_list.append(i)  
            yield index_list
            index_list = []


def merge_nearby_clusters(cs_cluster, max_distance):
    """
    Recursively merge neighboring clusters using a devide and conqurer algorithm
    Neighboring positions located within max_distance from the peak cluster with 
    max RPM merged into the peak cluster. 

    Arguments:
    cs_cluster is like:[position, [num in pass_files]] (Ordered by position)
    [[155603441, [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]],
     [155603444, [3, 1, 3, 14, 13, 8, 6, 1, 1, 10, 1, 6]],
     [155603445, [8, 1, 4, 13, 20, 8, 6, 0, 5, 19, 5, 19]]]

    Output:
    """
    sample_number = len(cs_cluster[0][1]) // 2  
    # Check base case
    clustered = True
    for i in range(len(cs_cluster))[:-1]:
        if cs_cluster[i + 1][0] - cs_cluster[i][0] <= max_distance:
            clustered = False

    if clustered:
        return True
    # When the clustering is not completed:
    else:
        # Get max of normalized read numbers from all samples and the index for max
        _, i = max((v, i) for i, v in enumerate((sum(pos_num[1][sample_number:])
                                                   for pos_num in cs_cluster)))
        # Merge cs within max_distance from the cs with max read number
        left_index = i  # leftmost position merged
        right_index = i  # rightmost position merged
        for j, (pos, nums) in enumerate(cs_cluster):
            if abs(pos - cs_cluster[i][0]) <= max_distance and \
                    abs(pos - cs_cluster[i][0]) > 0:
                # Combine the read numbers for two CS, sample by sample
                # cs_cluster[i] is like
                #              [155603471, [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]]
                # nums is like [0, 0, 0, 0, 0, 5, 1, 0, 0, 0, 0, 0]
                cs_cluster[i][1] = [sum(x)
                                    for x in zip(cs_cluster[i][1], nums)]
                cs_cluster[j][1] = 0  
                left_index = min(left_index, j)
                right_index = max(right_index, j)
        # Devide and conqure
        if left_index > 0:
            merge_nearby_clusters(cs_cluster[:left_index],
                                               max_distance)
        if right_index < len(cs_cluster) - 1:
            merge_nearby_clusters(cs_cluster[right_index + 1:],
                                               max_distance)

def cluster_dict(readcounts, max_distance):
    '''Converts a dict of readcount into a dict of clustered readcount.
    
    This function is only meant to be called by cluster_pass_reads(),
    rescue_nonpass_reads(), cluster_cleavage_sites().

    Arguments:
    readcounts, a dict containing something like:
    'chr11:+:31270274': [6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    '''
    print('Clustering reads ...')
    # Calculate the normalized read numbers and attach to the read numbers
    df = pd.DataFrame(readcounts).T
    normalized = df / df.sum(0)
    df = pd.concat([df, normalized], axis = 1, ignore_index = True)
    readcounts = df.T.to_dict('list')
    # Separate positions based on chrmosome & strand combination
    cpcounts = collections.defaultdict(list)
    while len(readcounts) > 0:
        k, v = readcounts.popitem()  
        k = k.split(':')
        cpcounts[':'.join(k[:2])].append([int(k[2]), v]) 
        # cpcounts.popitem() is like
        # {'chr9:-':[[100395423, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 56]],...]}
    del readcounts
    # Sort the list of lists for each chromosome & strand combination
    for k, v in cpcounts.items():
        # Sort in place to save memory
        v.sort(key=lambda val: val[0])
        # Get the indexes of lists containing neighboring positions
        for indexes in find_nearby_indexes(v, max_distance):
            cs_cluster = v[indexes[0]:(indexes[-1] + 1)]
            # The original list in the dict will be edited in place:
            merge_nearby_clusters(cs_cluster, max_distance)
        # Delete positions with 0 read number after clustering
        cpcounts[k] = [posi_num for posi_num in v if not posi_num[1] == 0]
    return cpcounts


def cluster_pass_reads(pass_files,
                       output = 'clusters.csv',
                       direction = 'reverse',
                       max_distance = 24):
    """
    Cluster PASS reads within max_distance.

    First generate a table with the following columns: chromosome, strand, pos, 
    num. Then recursively combine reads with LAP within max_distance nt (from 
    position with the highest read number to the position with next highest 
    read number).

    Build a dict of readcounts like {'chr9:-:100395423':[0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 56]}. The [list of int] saves the number of reads from each 
    pass_file. 
    """
    print('Reading sam files containing PASS reads ...')
    file_num = len(pass_files)
    # Pattern for finding last mapped position (LM) in read names
    re_pattern = re.compile('LM:i:(\d+)')
    # Readcounts is a read id counter
    readcounts = {}
    # Go through sam files and get read numbers for each calculated read id
    for i, pass_file in enumerate(pass_files):
       # Read sam file, calculate read id, and count number of reads
        fin = open(pass_file, 'r')
        for line in fin:
            if line[0] == '@': continue
            elements = line.split()
            chromosome = elements[2]
            flag = elements[1]

            if (direction == 'reverse' and flag == '16') or \
               (direction == 'forward' and flag == '0'): strand = '+'
            elif (direction == 'reverse' and flag == '0') or \
                 (direction == 'forward' and flag == '16'): strand = '-'

            position = re.search(re_pattern, elements[-1]).group(1)

            read_id = ':'.join([chromosome, strand, position])
            # Initialize the list for holding number of reads for each file
            readcounts.setdefault(read_id, [0] * file_num)[i] += 1
            # readcount.popitem() may return something like:
            # ('chr11:+:31270274', [6, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        fin.close()

    # Converts a dict of readcount into a dict of clustered readcount.
    cpcounts = cluster_dict(readcounts, max_distance)

    # Write the result to disk in csv format
    print('Writing to output file ...')
    with open(output, 'w') as fout:
        # Write header
        sample_string = ','.join([Path(pass_file).stem.split('.')[0].split('/')[-1]
                                  for pass_file in pass_files])
        fout.write(f'chromosome,strand,position,{sample_string}\n')
        # Write read counts for each cluster
        for k in cpcounts:
            chromosome, strand = k.split(':')
            for record in cpcounts[k]:
                position = str(record[0])
                counts = ','.join([str(int(count)) for count
                                   in record[1][:file_num]])
                fout.write(f'{chromosome},{strand},{position},{counts}\n')
    print('Done!')


def cluster_cleavage_sites(input_cs_files, output = 'combined.clusters.csv', 
                           max_distance = 24):
    '''Cluster cleavage sites from different cleavage site files (cs_files). 
    
    Arguments:
    input_cs_files: A list of file names. Each input_cs_file should be a csv file 
    with the following columns (in order): chromosome, strand, position, 
        read_num_for_each_sample.
    
    Output:
    A csv file with the following columns: chromosome, strand, position, 
    read_num_for_each_sample(from all cs_files).''' 
           
    # Record number of columns and sample names from each input file
    num_column = []
    sample_names = []
    for cs_file in input_cs_files:
        with open(cs_file, 'r') as fin:
            samples = fin.readline().strip().split(',')[3:]
            num_column.append(len(samples))
            sample_names += samples
    total_num_column = sum(num_column)        
    
    # Loop through cleavage site (CS) files and get read numbers for each CS id   
    readcounts = {}
    for i, cs_file in enumerate(input_cs_files): 
        print('Reading ' + cs_file)        
        fin = open(cs_file, 'r')
        for line in fin:
            elements = line.strip().split(',')
            if elements[0] == 'chromosome': # header
                continue
            chromosome, strand, position = elements[:3]
            # Calculate cleavage site id
            cs_id = ':'.join([chromosome, strand, position]) 
            readcounts.setdefault(cs_id, [0]*total_num_column)
            # Copy the read numbers to corresponding columns
            readcounts[cs_id][sum(num_column[:i]):sum(num_column[:i+1])] = \
                [int(x) for x in elements[3:]]
        fin.close()
    # readcounts example: 
    # {'chr9:-:100395423':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 56]}
   
    # Converts a dict of readcount into a dict of clustered readcount.
    cpcounts = cluster_dict(readcounts, max_distance)
    
    # Write the result to disk in csv format   
    print( 'Writing to file...')
    with open(output, 'w') as fout:
        sample_string = ','.join(sample_names)
        fout.write('chromosome,strand,position,%s\n' %sample_string)
        # Write read counts for each cluster
        for k in cpcounts:
            chromosome, strand = k.split(':')
            for record in cpcounts[k]:
                position = str(record[0])
                counts = ','.join([str(int(count)) for count in 
                                   record[1][:total_num_column]])
                fout.write('%s,%s,%s,%s\n'%(chromosome, strand, position, counts))


def sam2bigwig(sam_file, genome_size, keep_bam = False):
    '''sam -> bam -> bigwig'''
    # The sam2bigwid function cannot be defined within make_url, because
    # functions are only picklable if they are defined at the top-level of 
    # a module.
    p = Path(sam_file)
    prefix = str(p.parent/p.stem.split('.')[0])
    # sam -> bam
    cmd = f'samtools view -uS {sam_file} | samtools sort -o {prefix}.bam'
    os.system(cmd)
    # bam -> bedGraph
    totalReadNum = count_sam(sam_file)[1]
    cmd = (f'genomeCoverageBed -bg -split -ibam {prefix}.bam -strand + -g '
           f'{genome_size} -scale {str(-10**6/totalReadNum)} > '
           f'{prefix}.m.bedgraph'
          )
    os.system(cmd)
    # On some systems (with bad LC_COLLATE setting), sorting is required
    cmd = f'bedSort {prefix}.m.bedgraph {prefix}.m.bedgraph'
    os.system(cmd)
    cmd = (f'genomeCoverageBed -bg -split -ibam {prefix}.bam -strand - -g '
           f'{genome_size} -scale {str(10**6/totalReadNum)} > '
           f'{prefix}.p.bedgraph'
          )
    os.system(cmd)
    # On some systems (with bad LC_COLLATE setting), sorting is required
    cmd = f'bedSort {prefix}.p.bedgraph {prefix}.p.bedgraph'
    os.system(cmd)
    # Format the last column of the bedgraph files
    for file_name in (f'{prefix}.m.bedgraph', f'{prefix}.p.bedgraph'):
        bg = pd.read_csv(file_name, header = None, sep = "\t")
        bg.iloc[:,3] = bg.iloc[:,3].round(1)
        bg.to_csv(file_name, header = False, sep = "\t", index = False)
    # bedgraph -> bigWig
    cmd = f'bedGraphToBigWig {prefix}.p.bedgraph {genome_size} {prefix}.p.bw'
    os.system(cmd)
    cmd = f'bedGraphToBigWig {prefix}.m.bedgraph {genome_size} {prefix}.m.bw'
    os.system(cmd)
    # Remove intermediate files
    cmd = f'rm {prefix}*bedgraph*'
    # os.system(cmd)
    if not keep_bam:
        cmd = f'rm {prefix}*.bam'
        os.system(cmd)


def make_url(project, experiment, sam_dir, sam_files, genome_size, 
             keep_bam, sample_description, processes, bigDataUrl):
    '''Creates a file containing UCSC track records''' 
    # Create bigWig files
    l = len(sam_files)
    # sample_description = sample_description.reset_index()
    if re.search('\.nonpass$', sam_files[0]):
        read_type = 'nonPASS' 
    elif re.search('\.pass$', sam_files[0]):
        read_type = 'PASS' 
    else:
        read_type = 'RNASeq'
    # Convert sam files into bigwig files 
    with mp.Pool(processes = processes) as pool:
        pool.starmap(sam2bigwig, zip(sam_files, 
                                     [genome_size]*l, 
                                     [keep_bam]*l))
    bw_files = sorted([str(bw_file) for bw_file in sam_dir.glob('*.bw')])
    bw_samples = sorted(list(set([Path(filename).stem.split('.')[0].split('/')[-1] 
                                  for filename in bw_files])))
    # Calculate colors
    if 'track_color' in sample_description.columns:
        color_max = sample_description.track_color.max(axis=0)
    else:
        color_max = sample_description.shape[0]
    colors = np.array(sns.color_palette("colorblind", color_max))*255
    # Create UCSC genome browser tracks
    strand2str = {'+': 'p', '-': 'm'}
    f = open(sam_dir/'bigwigCaller.txt', 'w')
    for strand in ['+', '-']:
        for sample in bw_samples:
            color = colors[sample_description.loc[sample].track_color - 1, ]
            color = ','.join([str(int(c)) for c in color])
            track = (f'track type=bigWig visibility=2 alwaysZero=on color={color} '
                     f'graphType=bar maxHeightPixels=30:30:30 itemRgb=On group='
                     f'{project} name="{sample}{strand} {read_type}" '
                     f'description="{sample}{strand} {read_type}" '
                     f'bigDataUrl={bigDataUrl}{project}/{experiment}/'
                     f'{read_type}/{sample}.{strand2str[strand]}.bw\n'
                    )
            f.write(track)
    f.close()


def pick_A_stretch(pass_in, astretch_out, non_astretch_out, min_length=5):
    '''Pick PASS reads with at least min_length mapped 5' Ts in pass_in and copy 
    them into the astretch_out file. '''

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

