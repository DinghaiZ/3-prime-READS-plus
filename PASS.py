# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 22:53:05 2015
11/24/2015: replaced line.split('\t') with line.split() to avoid an error in 
Python 2.7
11/24/2015: replaced re.compile('LM:i:(\d+)\n') with re.compile('LM:i:(\d+)')
#11/28/2015: alow rRNA reads to be saved in PASS and noPASS files

@author: dinghai
"""
import os

import re

################################################################################

class Fastq(object):
    """
    Fastq record
    """
    def __init__(self, name, seq, qual):
        self.name = name
        self.seq = seq
        self.qual = qual
    
    def __str__(self):
        return '\n'.join(['@%s' % str(self.name), self.seq, 
            '+%s' % self.name, self.qual])
    
    def as_simple(self):
        return '\n'.join(['@%s' % str(self.name), self.seq, '+', self.qual])
    
    def get_name(self):
        return self.name
    
    def clean_name(self):
        self.name = self.name.split(" ")[0]
    
    def get_seq(self):
        return self.seq
    
    def get_qual(self):
        return self.qual
    
    def get_ascii(self):
        return map(ord, self.qual)
    
    def get_length(self):
        """Return length"""
        return len(self.seq)

    def get_num_N(self, first=28):
        """Return number of Ns"""
        return self.seq.count('N', 0, first)

    def set_name(self, name):
        """Reset name"""
        self.name = name
    
    def remove_tail_N(self):
        """Remove N at the end"""
        i = len(self.seq) - 1
        while i >= 0:
            if self.seq[i] == 'N':
                i -= 1
            else:
                break
        self.seq = self.seq[:i+1]
        self.qual = self.qual[:i+1]

  
################################################################################

def read_barcode(barcode_file):
    """
    Return oldname:newname dict
    """
    assert barcode_file.endswith('.csv')
    barcode2name = {}
    first_line = True
    for line in open(barcode_file):
        #print(line)
        fields = line.rstrip().split(',')
        assert len(fields) > 1
        
        if first_line: # skip the header
            first_line = False	        
            if 'prefix' in fields or 'sample' in fields:
    	        continue
         
        barcode2name[fields[1].rstrip()] = fields[2].replace('.', '_').rstrip()
        #assert re.search('_\d+$', barcode2name[fields[1].rstrip()])
        
    return barcode2name
    
def sample_colors(barcode_file):
    """
    Return newsamplename:color dict for UCSC genome browser
    In the barcode file, color should be in the '255,0,0' format.
    """
    samplecolors = {}
    first_line = True
    for line in open(barcode_file):
        fields = line.rstrip().split(',')
        if first_line: # skip the header
            first_line = False	        
            if 'color' in fields or 'Color' in fields:
    	        continue
        samplecolors[fields[2].replace('.', '_').rstrip()] =  fields[3]
    return samplecolors
    
################################################################################
def count_fastq(input_file):
    """
    Counts the number of fastq records in fastq file
    """
    from subprocess import check_output
   
    cmd = 'wc -l ' + input_file   
    cmd_out = check_output(cmd, shell=True)
    count = str(cmd_out).split('\'')[1].split(' ')[0]
    count = int(count)//4
    sample_name = input_file.split('/')[-1].split('.')[0]
    return (sample_name, count)
    
def count_pass(input_file):
    """
    Counts the number of pass or nonpass records in sam file
    """
    from subprocess import check_output
   
    cmd = 'grep -v @ ' + input_file + ' | wc -l' 
    cmd_out = check_output(cmd, shell=True)
    #print(cmd_out)
    count = int(str(cmd_out).split('\'')[1].strip('\\n'))
    
    return count
    
################################################################################
    
def reader_fastq(infile):
    """Generator of Fastq object from given file"""
    i = 0
    name = None
    seq = None
    qual = None
    if infile.endswith('.gz'):
        import gzip
        for line in gzip.open(infile, 'rb'):
            i += 1
            curr_line = line.strip()
            if i % 4 == 1:
                name = curr_line[1:]
            elif i % 4 == 2:
                seq = curr_line
            elif i % 4 == 0:
                qual = curr_line
                yield Fastq(name, seq, qual)
    else:
        for line in open(infile):
            i += 1
            curr_line = line.strip()
            if i % 4 == 1:
                name = curr_line[1:]
            elif i % 4 == 2:
                seq = curr_line
            elif i % 4 == 0:
                qual = curr_line
                yield Fastq(name, seq, qual)
        
            
################################################################################
    
def count_fastq_length(infolder, outfolder, pattern = 'trimmed.fastq', 
                      outfile = 'fastq_len.csv'):
    """ Calculate the distribution of read length for fastq files in a folder"""
    import os
    from collections import defaultdict
    fastq_len = {}
    
    fastqs = [filename for filename in os.listdir(infolder) if filename.endswith(pattern)]
    for infile in fastqs:
        i = 0    
        fastq_len[infile] = defaultdict(int)
        with open(os.path.join(infolder, infile), 'r') as fin:
            for line in fin:
                i += 1
                if i % 4 == 2:
                    fastq_len[infile][len(line.strip())] += 1
    #return fastq_len   
    rows = defaultdict(list)                
    for filename in fastq_len.keys():
        rows[filename].append(filename)
        for readlength in sorted(fastq_len[filename].keys()):
            rows[filename].append(fastq_len[filename][readlength])
    rows = rows.values()   
    rows.sort(key = lambda x: x[0])
    rows = zip(*rows)
    with open(os.path.join(outfolder, outfile), 'w') as fout:
        #for k in sorted(row.keys()):
         #   string = ','.join([str(num) for num in row[k]]) + '\n'       
          #  fout.write(string)
        for row in rows:
            string = ','.join([str(num) for num in row]) + '\n' 
            fout.write(string)
            
################################################################################
def fastq_trim_Ns(infile, outfile, random_NT_len=4):
    '''
    If the 3' end of the read has been trimmed by cutadapt to remove the 3' end
    of the 5' adapter, trim the last 4 random nucleotides from the 5' adapter 
    and save them in the read name.  
    
    @HISEQ01:507:HBC1DADXX:1:1101:1222:2105 1:N:0:GTCCGC cutadapt    
    AACTTTTGTTTTTTGACAGTCTCAAGTTTTTATTCAGTGGGTCTCTGTGTC
    becomes:
    @HISEQ01:507:HBC1DADXX:1:1101:1222:2105 1:N:0:GTCCGC TGTC
    AACTTTTGTTTTTTGACAGTCTCAAGTTTTTATTCAGTGGGTCTCTG
    '''
    import re
    #pattern = re.compile('([ATCGN]{%d})$'% random_NT_len)#precompile
    with open(outfile, 'w') as outhandle:
        for fq in reader_fastq(infile):
            name = fq.get_name()
            seq = fq.get_seq()
            qual = fq.get_qual()
            
            # trim 4N from 3' end if the 5' adapter (for reverse sequencing) has been removed
            if re.search('cutadapt$', name):
                name = name.replace('cutadapt', '4N' + seq[-random_NT_len:] + '4N')
                seq = seq[:-random_NT_len]
                qual = qual[:-random_NT_len]
                                                  
            # Only keep the remaining sequence if its length is >= 18 nt
            if len(seq) >= 18:
                outhandle.write('\n'.join([name, seq, '+', qual]) + '\n')
                




################################################################################
def fastq_trim_Ts_2(infile, random_NT_len=3):
    '''
    Trim the first random nucleotides from 3' adapter and save them in the read
    name. Also trim 5'Ts and save the number of 5' Ts in the read name. 
    Use a two-step trimming strategy to deal with sequencing errors in T-stretches.
    TTTTGTTVNNNN will become VNNNN. The GTT sequence will be attached to the end 
    of the read name after ::, which will be compared with genomic sequence 
    downstream of the LAP. For example:

    @HISEQ01:507:HBC1DADXX:1:1101:1222:2105 1:N:0:GTCCGC    
    AACTTTTGTTTTTTGACAGTCTCAAGTTTTTATTCAGTGGGTCTCTGTGTC
    
    becomes:
    
    @TS11AAC:507:HBC1DADXX:1:1101:1243:2103 1:N:0:GTCCGC
    GACAGTCTCAAGTTTTTATTCAGTGGGTCTCTGTGTC
    TS11: T-stretch length is 11; AAC: the three random nucleotide is AAC

    '''
    import re
    pattern = re.compile('([ATCGN]{%d})(T*)' % random_NT_len)  # precompile
    outfile = infile.replace('.fastq', '.trimmed.fastq')
    with open(outfile, 'w') as outhandle:
        for fq in reader_fastq(infile):
            seq = fq.get_seq()
            qual = fq.get_qual()

            match1 = re.match(pattern, seq)
            seq = seq[match1.end():]
            qual = qual[match1.end():]
            T_length1 = len(match1.groups()[1])
            read_name = ''.join(['@TS', str(T_length1),
                                 match1.groups()[0], fq.get_name()[7:]])

            # The second step deal with reads like
            # TTTTTGTTTTTTTCCAGTTGTCAAATGATCCTTTAT
            match2 = re.match('[ACGN](T+)', seq)
            if match2:
                T_length2 = len(match2.groups()[0])
                if T_length2 > 2 and (T_length1 + T_length2) > 5:
                    seq = seq[match2.end():]
                    qual = qual[match2.end():]
                    T_length1 = T_length1 + 1 + T_length2
                    # attach the trimmed sequence to the read name
                    # the sequence will be compared with genomic sequence later
                    read_name = ''.join(['@TS', str(T_length1),
                                         match1.groups()[0], '::',  match2.group(0)])
                    #read_name = read_name + '::' + match2.group(0)

            # Only keep the remaining sequence if its length is >= 18 nt
            if len(seq) >= 18:
                outhandle.write('\n'.join([read_name, seq, '+', qual]) + '\n')
    os.system('rm ' + infile)               
            
################################################################################
def fastqgz_trim_Ts(infile, outfile, random_NT_len=3):
    '''
    Similar to gastq_trim_Ts(), but read and write .gz files
    
    '''
    import re
    import gzip
    pattern = re.compile('([ATCGN]{%d})(T*)'% random_NT_len)#precompile
    with gzip.open(outfile, 'wb') as outhandle:
        for fq in reader_fastq(infile):
            seq = fq.get_seq()
            qual = fq.get_qual()
            
            match1 = re.match(pattern,seq)
            seq = seq[match1.end():]
            qual = qual[match1.end():]
            T_length1 = len(match1.groups()[1])
            #read_name = '@'+fq.get_name()+' TS:'+str(match.end() - random_NT_len)
            read_name = ''.join(['@TS', str(T_length1),
                                match1.groups()[0], fq.get_name()[7:]])
            
            # The second step deal with reads like TTTTTGTTTTTTTCCAGTTGTCAAATGATCCTTTAT
            match2 = re.match('[ACGN](T+)',seq)
            if match2:
                T_length2 = len(match2.groups()[0])
                if T_length2 > 2: #and (T_length1 + T_length2) > 5:
                    seq = seq[match2.end():]
                    qual = qual[match2.end():]
                    T_length1 = T_length1 + 1 + T_length2
                    # attach the trimmed sequence to the read name
                    # the sequence will be compared with genomic sequence later
                    read_name = ''.join(['@TS', str(T_length1),
                                match1.groups()[0], '::',  match2.group(0)])
                    #read_name = read_name + '::' + match2.group(0)
                                      
            # Only keep the remaining sequence if its length is >= 15 nt
            if len(seq) >= 15:
                outhandle.write('\n'.join([read_name, seq, '+', qual]) + '\n')
            

################################################################################            
def fastq_trim_As(infile, outfile, random_NT_len=0):
    '''
    Trim the first random nucleotides from 3' adapter and save them in the read
    name. Also trim 3'As and save the number of 3' As in the sequence name. 
    Use a two-step trimming strategy to deal with sequencing errors in T-stretches.
    TTTTGTTVNNNN will become VNNNN. The GTT sequence will be attached to the end 
    of the read name after ::, which will be compared with genomic sequence 
    downstream of the LAP.
    
    @HISEQ01:507:HBC1DADXX:1:1101:1222:2105 1:N:0:GTCCGC    
    AACTTTTGTTTTTTGACAGTCTCAAGTTTTTATTCAGTGGGTCTCTGTGTC
    becomes:
    @TS11AAC:507:HBC1DADXX:1:1101:1243:2103 1:N:0:GTCCGC
    GACAGTCTCAAGTTTTTATTCAGTGGGTCTCTGTGTC
    TS11: T-stretch length is 11; AAC: the three random nucleotide is AAC
    
    '''
    import re
    pattern = re.compile('([ATCGN]{%d})(A*)$'% random_NT_len)#precompile
    with open(outfile, 'w') as outhandle:
        for fq in reader_fastq(infile):
            seq = fq.get_seq()
            qual = fq.get_qual()
            
            match1 = re.search(pattern,seq)
            seq = seq[:match1.start()]
            qual = qual[:match1.start()]
            T_length1 = len(match1.groups()[1]) 
            read_name = ''.join(['@TS', str(T_length1),
                                match1.groups()[0]]) # , fq.get_name()[7:] deleted
            
            # The second step deal with reads like TTTTTGTTTTTTTCCAGTTGTCAAATGATCCTTTAT
            match2 = re.search('[ACGN](A+)',seq)
            if match2:
                T_length2 = len(match2.groups()[0])
                if T_length2 > 2 and (T_length1 + T_length2) > 5:
                    seq = seq[:match1.start()]
                    qual = qual[:match1.start()]
                    T_length1 = T_length1 + 1 + T_length2
                    # attach the trimmed sequence to the read name
                    # the sequence will be compared with genomic sequence later
                    read_name = ''.join(['@TS', str(T_length1),
                                match1.groups()[0], '::',  match2.group(0)])
                    #read_name = read_name + '::' + match2.group(0)
                                      
            # Only keep the remaining sequence if its length is >= 18 nt
            if len(seq) >= 18:
                outhandle.write('\n'.join([read_name, seq, '+', qual]) + '\n')            
################################################################################
def load_fasta_genome(genome_dir = '/HPCTMP_NOBKUP/wl314/data/ucsc/genomes/mm9'):
    '''
    Read specified genomic sequence saved in genome_dir into a dict.
    '''
    fasta_files = os.listdir(genome_dir)
    fasta_files = [x for x in fasta_files if re.search('^chr.+fa$',x)]
    # !du -ch {genome_dir+'/ch*fa'}|grep total # total size of these files
    genome = {}
    for fasta_file in fasta_files:
        #print('Loading ', fasta_file)
        with open(os.path.join(genome_dir, fasta_file), 'r') as f:
            f.readline()#skip the first line
            line = f.read() #much faster than looping through lines!
            genome[fasta_file.replace('.fa', '')] = line.replace('\n', '')
    return genome
    
################################################################################
def load_fasta_genome_2(genome_dir = '/HPCTMP_NOBKUP/wl314/data/ucsc/genomes/mm9'):
    '''
    Read specified genomic sequence saved in genome_dir into a dict in shared memory.
    '''
    from multiprocessing import Manager
    m = Manager()
    
    genome = m.dict() #compare with genome = {}: not in shared memory
    
    fasta_files = os.listdir(genome_dir)
    fasta_files = [x for x in fasta_files if re.search('^chr.+fa$',x)]
    # !du -ch {genome_dir+'/ch*fa'}|grep total # total size of these files
    
    for fasta_file in fasta_files:
        #print('Loading ', fasta_file)
        with open(os.path.join(genome_dir, fasta_file), 'r') as f:
            f.readline()#skip the first line
            line = f.read() #much faster than looping through lines!
            genome[fasta_file.replace('.fa', '')] = line.replace('\n', '')
    return genome
    
   
        
    
    
################################################################################
def reverse_complement(seq, mol_type = 'DNA'):
    #genomic sequence is a mixture of lower and upper cases    
    seq = seq.upper()
    if mol_type == 'DNA':
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    elif mol_type == 'RNA':
        complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A', 'N': 'N'}
    reverse_complement = "".join(complement.get(base) \
    for base in reversed(seq))
    return reverse_complement
    
################################################################################
def complement(seq, mol_type = 'DNA'):
    #genomic sequence is a mixture of lower and upper cases    
    seq = seq.upper()
    if mol_type == 'DNA':
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    elif mol_type == 'RNA':
        complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A', 'N': 'N'}
    new = "".join(complement.get(base) for base in seq)
    return new
    
################################################################################
def get_seq(chromosome, strand, start, end, genome):
    '''Get the genomic sequence specified. 'genome' is a dict.
    Reverse strand sequence is reversecomplemented.'''
    # get the sequence
    seq = genome[chromosome][(start-1) : end].upper() #genomic sequence is a mixture of lower and upper cases!!!
    
    # reverse complement if strand is '-'
    if strand == '-':
        seq = reverse_complement(seq, mol_type = 'DNA')
    return seq
    
        
################################################################################    
def pick_PASS(sam_in, sam_pass, sam_nopass, sam_ref, genome, 
              min_mapq = 10, direction = 'reverse'):
    '''Loop through records in the sam_in file, if the record is PASS, write to
    sam_pass. Mapped non-PASS reads are saved in sam_nopass. Reads mapped to 
    yeast genome are skipped. "genome" is a dict.  Sequencing direction can only 
    be reverse or forward. sam_ref is the file for spike-in reads. 
    Attach matched tail length (ML), unmatched tail length (UL), and last mapped
    position (LM) to the end of readname. The head (such as 'TS11AAC') and the 
    tail (such as 'ML:i:4\tUL:i:9\tLM:i:1234') of the read name is the unique 
    identifier of the read.'''

    sam_pass_file = open(sam_pass, 'w')
    sam_nopass_file = open(sam_nopass, 'w')
    sam_ref_file = open(sam_ref, 'w')
    lap = 0
            
    with open(sam_in, 'r') as sam_in_file:
        for line in sam_in_file:
            # skip header
            if line[0] == '@': 
                sam_pass_file.write(line)
                sam_nopass_file.write(line)
                sam_ref_file.write(line)
                continue
            #process each line
            #print(line)
            (readname, flag, chromosome, position, mapq, cigar) = \
            line.split()[:6]
            
            position = int(position)
            mapq = int(mapq)
            
            # discard reads with low mapping scores
            if mapq < min_mapq: 
                continue             
            # record reads mapped to yeast genome
            if chromosome[:4] == 'ychr': #or chromosome[:4] == 'BK000964': #rRNA 
                sam_ref_file.write(line)
                continue       
            # ignore reads mapped to other chromosomes, such as 'BK000964'
            if not chromosome[:3] == 'chr': 
                continue  
            
                       
            # extract the T-stretch length encoded in the readname
            t_stretch_len = int(re.match('TS(\d+)', line).groups()[0])
            # skip records with short T-stretches
            if t_stretch_len < 2: 
                sam_nopass_file.write(line)
                continue
            
            # get genomic sequence downstream of the LAP (last mapped position). 
            # (The length of returned sequence is t_stretch_len.)
            if (direction.lower() == 'reverse' and flag == '16') or \
            (direction.lower() == 'forward' and flag == '0'): 
                strand = '+'
                # process cigar to determine the LAP
                # Insertion (I) will not affect the covered distance
                # if cigar = '38M10S' or '38', nums will be 38
                nums = re.split('[MDN]', re.sub('\d+I','', cigar))[:-1]
                covered =  sum(int(x) for x in nums)
                # get_seq() returns reverse complemented sequence if strand == '-'
                downstream_seq = get_seq(chromosome, strand, 
                start = position + covered, 
                end = position + covered + t_stretch_len - 1, 
                genome = genome)
                lap = position + covered -1
            elif (direction.lower() == 'reverse' and flag == '0') or \
            (direction.lower() == 'forward' and flag == '16'):
                strand = '-'
                downstream_seq = get_seq(chromosome, strand, 
                start = position - t_stretch_len, end = position -1, 
                genome = genome)     
                lap = position
            else:
                continue 

            # if the 5' T-stretch was trimmed twice (like 'TTTTTTTGTTTT'), 
            # check whether 'GTTTT' comes from the genome           
            #match = re.match('[ACG](T*)', readname.split('::')[-1])
            elements = readname.split('::')            
            if len(elements) == 2 and re.match(reverse_complement(elements[1]), 
                                              downstream_seq):
                downstream_seq = downstream_seq[len(elements[1]):]
                t_stretch_len -= len(elements[1]) 
                if strand == '+': 
                    lap += len(elements[1])
                if strand == '-':
                    lap -= len(elements[1])
            
            # no need to analyze downstream sequence        
            if t_stretch_len == 0:
                sam_nopass_file.write(line.strip()+'\tML:i:%d\tUL:i:%d\tLM:i:%d\n'
                %(len(elements[1]), 0, lap))
                continue
                        
            # analyze downstream sequence
            if not downstream_seq[:-1] == 'A'*(t_stretch_len-1):                
                match = re.match('A+', downstream_seq)
                if match: 
                    matched_len = len(match.group())
                else: 
                    matched_len = 0
                # attach matched length (ML) and unmatched length (UL) to readname
                # t_stretch_len - matched_len could be 1 for the few reads mapped 
                # to the end of chromosome.
                sam_pass_file.write(line.strip()+'\tML:i:%d\tUL:i:%d\tLM:i:%d\n'
                %(matched_len, t_stretch_len - matched_len, lap))                
            else:                           
                match = re.match('A*', downstream_seq) 
                matched_len = len(match.group())# value: 1 ~ t_stretch_len
                #nopass_aligned_tail_len.append(len(match.group()))
                # attach matched length (ML) and unmatched length (UL) to readname
                sam_nopass_file.write(line.strip()+'\tML:i:%d\tUL:i:%d\tLM:i:%d\n'
                %(matched_len, t_stretch_len - matched_len, lap))                 

    sam_pass_file.close() 
    sam_nopass_file.close() 
    sam_ref_file.close()

################################################################################    
def pick_PASS_2(sam_file, genome, min_mapq = 10,  direction = 'reverse'):
    """
    Loop through records in the sam file, if the record is PASS, write to
    sam_pass. Mapped non-PASS reads are saved in sam_nopass. Reads mapped to 
    yeast genome are skipped. "genome" is a dict.  Sequencing direction can only 
    be reverse or forward. sam_ref is the file for spike-in reads. 
    Attach matched tail length (ML), unmatched tail length (UL), and last mapped
    position (LM) to the end of readname. The head (such as 'TS11AAC') and the 
    tail (such as 'ML:i:4\tUL:i:9\tLM:i:1234') of the read name is the unique 
    identifier of the read.
    """
    
    sam_pass = sam_file.replace('.sam', '.pass')
    sam_nopass = sam_file.replace('.sam', '.nopass')
    sam_ref = sam_file.replace('.sam', '.ref')
    
    sam_pass_file = open(sam_pass, 'w')
    sam_nopass_file = open(sam_nopass, 'w')
    sam_ref_file = open(sam_ref, 'w')
    lap = 0
            
    with open(sam_file, 'r') as in_file:
        for line in in_file:
            # skip header
            if line[0] == '@': 
                sam_pass_file.write(line)
                sam_nopass_file.write(line)
                sam_ref_file.write(line)
                continue
            #process each line
            #print(line)
            # line = 'TS14AAG2:51:HYFCNBGXX:2:11111:21123:6214        16      chr3    126547771       255     59M     *       0       0       CGTAAGGGTAAATGACTGATTGATATATTTACGCGTTAATAAATTTGTGATTTCTGCTG     /A/E///E/EE<A//A//</<///E//AE/E///EA/EE/EEAEE<</<AEE///EEE6     NH:i:1  HI:i:1  AS:i:52 nM:i:3'
            (readname, flag, chromosome, position, mapq, cigar) = \
            line.split()[:6]
            
            position = int(position)
            mapq = int(mapq)
            
            # discard reads with low mapping scores
            if mapq < min_mapq: 
                continue             
            # record reads mapped to yeast genome
            if chromosome[:4] == 'ychr': #or chromosome[:4] == 'BK000964': #rRNA 
                sam_ref_file.write(line)
                continue       
            # ignore reads mapped to other chromosomes, such as 'BK000964'
            if not chromosome[:3] == 'chr': 
                continue  
            
                       
            # extract the T-stretch length encoded in the readname
            t_stretch_len = int(re.match('TS(\d+)', line).groups()[0])
            # skip records with short T-stretches
            if t_stretch_len < 2: 
                sam_nopass_file.write(line)
                continue
            
            # get genomic sequence downstream of the LAP (last mapped position). 
            # (The length of returned sequence is t_stretch_len.)
            if (direction.lower() == 'reverse' and flag == '16') or \
            (direction.lower() == 'forward' and flag == '0'): 
                strand = '+'
                # process cigar to determine the LAP
                # Insertion (I) will not affect the covered distance
                # if cigar = '38M10S' or '38', nums will be 38
                nums = re.split('[MDN]', re.sub('\d+I','', cigar))[:-1]
                covered =  sum(int(x) for x in nums)
                # get_seq() returns reverse complemented sequence if strand == '-'
                downstream_seq = get_seq(chromosome, strand, 
                start = position + covered, 
                end = position + covered + t_stretch_len - 1, 
                genome = genome)
                lap = position + covered -1
            elif (direction.lower() == 'reverse' and flag == '0') or \
            (direction.lower() == 'forward' and flag == '16'):
                strand = '-'
                downstream_seq = get_seq(chromosome, strand, 
                start = position - t_stretch_len, end = position -1, 
                genome = genome)     
                lap = position
            else:
                continue 

            # if the 5' T-stretch was trimmed twice (like 'TTTTTTTGTTTT'), 
            # check whether 'GTTTT' comes from the genome           
            #match = re.match('[ACG](T*)', readname.split('::')[-1])
            elements = readname.split('::')            
            if len(elements) == 2 and re.match(reverse_complement(elements[1]), 
                                              downstream_seq):
                downstream_seq = downstream_seq[len(elements[1]):]
                t_stretch_len -= len(elements[1]) 
                if strand == '+': 
                    lap += len(elements[1])
                if strand == '-':
                    lap -= len(elements[1])
            
            # no need to analyze downstream sequence        
            if t_stretch_len == 0:
                sam_nopass_file.write(line.strip()+'\tML:i:%d\tUL:i:%d\tLM:i:%d\n'
                %(len(elements[1]), 0, lap))
                continue
                        
            # analyze downstream sequence
            if not downstream_seq[:-1] == 'A'*(t_stretch_len-1):                
                match = re.match('A+', downstream_seq)
                if match: 
                    matched_len = len(match.group())
                else: 
                    matched_len = 0
                # attach matched length (ML) and unmatched length (UL) to readname
                # t_stretch_len - matched_len could be 1 for the few reads mapped 
                # to the end of chromosome.
                sam_pass_file.write(line.strip()+'\tML:i:%d\tUL:i:%d\tLM:i:%d\n'
                %(matched_len, t_stretch_len - matched_len, lap))                
            else:                           
                match = re.match('A*', downstream_seq) 
                matched_len = len(match.group())# value: 1 ~ t_stretch_len
                #nopass_aligned_tail_len.append(len(match.group()))
                # attach matched length (ML) and unmatched length (UL) to readname
                sam_nopass_file.write(line.strip()+'\tML:i:%d\tUL:i:%d\tLM:i:%d\n'
                %(matched_len, t_stretch_len - matched_len, lap))                 

    sam_pass_file.close() 
    sam_nopass_file.close() 
    sam_ref_file.close()
    

################################################################################    
def pick_PASS_3(sam_file, genome, min_mapq = 10,  direction = 'reverse'):
    """
    Loop through records in the sam file, if the record is PASS, write to
    sam_pass. Mapped non-PASS reads are saved in sam_nopass. Reads mapped to 
    yeast genome are skipped. "genome" is a dict.  Sequencing direction can only 
    be reverse or forward. sam_ref is the file for spike-in reads. 
    Attach matched tail length (ML), unmatched tail length (UL), and last mapped
    position (LM) to the end of readname. The head (such as 'TS11AAC') and the 
    tail (such as 'ML:i:4\tUL:i:9\tLM:i:1234') of the read name is the unique 
    identifier of the read.
    """
    
    sam_pass = sam_file.replace('.sam', '.pass')
    sam_nopass = sam_file.replace('.sam', '.nopass')
    sam_ref = sam_file.replace('.sam', '.ref')
    
    sam_pass_file = open(sam_pass, 'w')
    sam_nopass_file = open(sam_nopass, 'w')
    sam_ref_file = open(sam_ref, 'w')
    lap = 0
            
    with open(sam_file, 'r') as in_file:
        lines = in_file.read().splitlines()  ### faster then pick_PASS_2()?
        for line in lines:
            # skip header
            if line[0] == '@': 
                sam_pass_file.write(line)
                sam_nopass_file.write(line)
                sam_ref_file.write(line)
                continue
            #process each line
            #print(line)
            # line = 'TS14AAG2:51:HYFCNBGXX:2:11111:21123:6214        16      chr3    126547771       255     59M     *       0       0       CGTAAGGGTAAATGACTGATTGATATATTTACGCGTTAATAAATTTGTGATTTCTGCTG     /A/E///E/EE<A//A//</<///E//AE/E///EA/EE/EEAEE<</<AEE///EEE6     NH:i:1  HI:i:1  AS:i:52 nM:i:3'
            (readname, flag, chromosome, position, mapq, cigar) = \
            line.split()[:6]
            
            position = int(position)
            mapq = int(mapq)
            
            # discard reads with low mapping scores
            if mapq < min_mapq: 
                continue             
            # record reads mapped to yeast genome
            if chromosome[:4] == 'ychr': #or chromosome[:4] == 'BK000964': #rRNA 
                sam_ref_file.write(line)
                continue       
            # ignore reads mapped to other chromosomes, such as 'BK000964'
            if not chromosome[:3] == 'chr': 
                continue  
            
                       
            # extract the T-stretch length encoded in the readname
            t_stretch_len = int(re.match('TS(\d+)', line).groups()[0])
            # skip records with short T-stretches
            if t_stretch_len < 2: 
                sam_nopass_file.write(line)
                continue
            
            # get genomic sequence downstream of the LAP (last mapped position). 
            # (The length of returned sequence is t_stretch_len.)
            if (direction.lower() == 'reverse' and flag == '16') or \
            (direction.lower() == 'forward' and flag == '0'): 
                strand = '+'
                # process cigar to determine the LAP
                # Insertion (I) will not affect the covered distance
                # if cigar = '38M10S' or '38', nums will be 38
                nums = re.split('[MDN]', re.sub('\d+I','', cigar))[:-1]
                covered =  sum(int(x) for x in nums)
                # get_seq() returns reverse complemented sequence if strand == '-'
                downstream_seq = get_seq(chromosome, strand, 
                start = position + covered, 
                end = position + covered + t_stretch_len - 1, 
                genome = genome)
                lap = position + covered -1
            elif (direction.lower() == 'reverse' and flag == '0') or \
            (direction.lower() == 'forward' and flag == '16'):
                strand = '-'
                downstream_seq = get_seq(chromosome, strand, 
                start = position - t_stretch_len, end = position -1, 
                genome = genome)     
                lap = position
            else:
                continue 

            # if the 5' T-stretch was trimmed twice (like 'TTTTTTTGTTTT'), 
            # check whether 'GTTTT' comes from the genome           
            #match = re.match('[ACG](T*)', readname.split('::')[-1])
            elements = readname.split('::')            
            if len(elements) == 2 and re.match(reverse_complement(elements[1]), 
                                              downstream_seq):
                downstream_seq = downstream_seq[len(elements[1]):]
                t_stretch_len -= len(elements[1]) 
                if strand == '+': 
                    lap += len(elements[1])
                if strand == '-':
                    lap -= len(elements[1])
            
            # no need to analyze downstream sequence        
            if t_stretch_len == 0:
                sam_nopass_file.write(line.strip()+'\tML:i:%d\tUL:i:%d\tLM:i:%d\n'
                %(len(elements[1]), 0, lap))
                continue
                        
            # analyze downstream sequence
            if not downstream_seq[:-1] == 'A'*(t_stretch_len-1):                
                match = re.match('A+', downstream_seq)
                if match: 
                    matched_len = len(match.group())
                else: 
                    matched_len = 0
                # attach matched length (ML) and unmatched length (UL) to readname
                # t_stretch_len - matched_len could be 1 for the few reads mapped 
                # to the end of chromosome.
                sam_pass_file.write(line.strip()+'\tML:i:%d\tUL:i:%d\tLM:i:%d\n'
                %(matched_len, t_stretch_len - matched_len, lap))                
            else:                           
                match = re.match('A*', downstream_seq) 
                matched_len = len(match.group())# value: 1 ~ t_stretch_len
                #nopass_aligned_tail_len.append(len(match.group()))
                # attach matched length (ML) and unmatched length (UL) to readname
                sam_nopass_file.write(line.strip()+'\tML:i:%d\tUL:i:%d\tLM:i:%d\n'
                %(matched_len, t_stretch_len - matched_len, lap))                 

    sam_pass_file.close() 
    sam_nopass_file.close() 
    sam_ref_file.close()    
    
    
################################################################################    
def pick_PASS_4(sam_file, min_mapq=10,  direction='reverse'):
    """
    Loop through records in the sam file, if the record is PASS, write to
    sam_pass. Mapped non-PASS reads are saved in sam_nopass. Reads mapped to 
    yeast genome are skipped. "genome" is a dict.  Sequencing direction can only 
    be reverse or forward. sam_ref is the file for spike-in reads. 
    Attach matched tail length (ML), unmatched tail length (UL), and last mapped
    position (LM) to the end of readname. The head (such as 'TS11AAC') and the 
    tail (such as 'ML:i:4\tUL:i:9\tLM:i:1234') of the read name is the unique 
    identifier of the read.
    """

    sam_pass = sam_file.replace('.sam', '.pass')
    sam_nopass = sam_file.replace('.sam', '.nopass')
    sam_ref = sam_file.replace('.sam', '.ref')

    sam_pass_file = open(sam_pass, 'w')
    sam_nopass_file = open(sam_nopass, 'w')
    sam_ref_file = open(sam_ref, 'w')
    lap = 0

    with open(sam_file, 'r') as in_file:
        for line in in_file:
            # skip header
            if line[0] == '@':
                sam_pass_file.write(line)
                sam_nopass_file.write(line)
                sam_ref_file.write(line)
                continue
            # process each line
            # line = 'TS14AAG2:51:HYFCNBGXX:2:11111:21123:6214        16      \
            # chr3    126547771       255     59M     *       0       0       \
            #CGTAAGGGTAAATGACTGATTGATATATTTACGCGTTAATAAATTTGTGATTTCTGCTG     \
            #/A/E///E/EE<A//A//</<///E//AE/E///EA/EE/EEAEE<</<AEE///EEE6     \
            # NH:i:1  HI:i:1  AS:i:52 nM:i:3'
            (readname, flag, chromosome, position, mapq, cigar) = \
                line.split()[:6]

            position = int(position)
            mapq = int(mapq)

            # discard reads with low mapping scores
            if mapq < min_mapq:
                continue
            # record reads mapped to spiked-in yeast genome # or rRNA
            if chromosome[:4] == 'ychr':  # or chromosome[:4] == 'BK000964':
                sam_ref_file.write(line)
                continue
            # ignore reads mapped to other chromosomes, such as 'BK000964'
            if not chromosome[:3] == 'chr':
                continue

            # extract the T-stretch length encoded in the readname
            t_stretch_len = int(re.match('TS(\d+)', line).groups()[0])
            # skip records with short T-stretches
            if t_stretch_len < 2:
                sam_nopass_file.write(line)
                continue

            # get genomic sequence downstream of the LAP (last mapped position).
            # (The length of returned sequence is t_stretch_len.)
            if (direction.lower() == 'reverse' and flag == '16') or \
                    (direction.lower() == 'forward' and flag == '0'):
                strand = '+'
                # process cigar to determine the LAP
                # Insertion (I) will not affect the covered distance
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
                                         start=position - t_stretch_len, end=position - 1,
                                         genome=genome)
                lap = position
            else:
                continue

            # if the 5' T-stretch was trimmed twice (like 'TTTTTTTGTTTT'),
            # check whether 'GTTTT' comes from the genome
            #match = re.match('[ACG](T*)', readname.split('::')[-1])
            elements = readname.split('::')
            if len(elements) == 2 and re.match(reverse_complement(elements[1]),
                                               downstream_seq):
                downstream_seq = downstream_seq[len(elements[1]):]
                t_stretch_len -= len(elements[1])
                if strand == '+':
                    lap += len(elements[1])
                if strand == '-':
                    lap -= len(elements[1])

            # no need to analyze downstream sequence
            if t_stretch_len == 0:
                sam_nopass_file.write(line.strip() + '\tML:i:%d\tUL:i:%d\tLM:i:%d\n'
                                      % (len(elements[1]), 0, lap))
                continue

            # analyze downstream sequence
            if not downstream_seq[:-1] == 'A' * (t_stretch_len - 1):
                match = re.match('A+', downstream_seq)
                if match:
                    matched_len = len(match.group())
                else:
                    matched_len = 0
                # attach matched length (ML) and unmatched length (UL) to readname
                # t_stretch_len - matched_len could be 1 for the few reads mapped
                # to the end of chromosome.
                sam_pass_file.write(line.strip() + '\tML:i:%d\tUL:i:%d\tLM:i:%d\n'
                                    % (matched_len, t_stretch_len - matched_len, lap))
            else:
                match = re.match('A*', downstream_seq)
                matched_len = len(match.group())  # value: 1 ~ t_stretch_len
                # nopass_aligned_tail_len.append(len(match.group()))
                # attach matched length (ML) and unmatched length (UL) to readname
                sam_nopass_file.write(line.strip() + '\tML:i:%d\tUL:i:%d\tLM:i:%d\n'
                                      % (matched_len, t_stretch_len - matched_len, lap))

    sam_pass_file.close()
    sam_nopass_file.close()
    sam_ref_file.close()
    os.system('rm ' + sam_file)    
    

    
################################################################################
def count_pass_in_folder(pass_dir, outfolder,  
                          outfile = 'pass_num.csv'):
    '''Count number of pass and nonpass reads in sam files in the pass_dir'''
    from subprocess import check_output
    import glob
    
    sample_names = []
    pass_nums = [] 
    nonpass_nums = [] 
    ref_nums = []
    
    # count pass reads 
    for file_name in sorted(glob.glob(os.path.join(pass_dir, '*.pass.sam'))):
        sample_names.append(file_name.split('/')[-1].split('.')[0])
        cmd = 'grep -v @ ' + file_name + ' | wc -l'
        cmd_out = check_output(cmd, shell=True)
        pass_nums.append(cmd_out.strip())
    
    # count nonpass reads 
    for file_name in sorted(glob.glob(os.path.join(pass_dir, '*.nopass.sam'))):
        cmd = 'grep -v @ ' + file_name + ' | wc -l'
        cmd_out = check_output(cmd, shell=True)
        nonpass_nums.append(cmd_out.strip())
    
    # count reference reads 
    for file_name in sorted(glob.glob(os.path.join(pass_dir, '*.ref.sam'))):
        cmd = 'grep -v @ ' + file_name + ' | wc -l'
        cmd_out = check_output(cmd, shell=True)
        ref_nums.append(cmd_out.strip())
    
    # combine
    records = zip(sample_names, pass_nums, nonpass_nums, ref_nums)
    with open(outfolder+'/' + outfile, 'w') as f:
        f.write('Sample, PASS, nonPASS, Ref\n')
        for rec in records:
            f.write('%s,%s,%s,%s\n' % rec)
            
            
################################################################################    
def count_5Ts_in_folder_2(pass_dir, result_dir):
    '''Count the number of reads with certain 5'T-stretch lengths for pass and 
    nopass reads in files in the pass_dir.'''
    file_names = []
    results = [] #container                     
    
    import re
    import numpy as np
    import pandas as pd
    
    for sam_file in [filename for filename in os.listdir(pass_dir) \
    if re.search('Aligned.out.(no)?pass$', filename)]:
        sam_in = os.path.join(pass_dir, sam_file)    
        sam_file = sam_file.split('.')[0] + '.' + sam_file.split('.')[-1]
        file_names.append(sam_file)  
        
        with open(sam_in, 'r') as f:
            string = f.read()
            # search patterns in the file        
            TS = re.findall('\nTS(\d+)', string) #5' t-strech
            ML = re.findall('ML:i:(\d+)\t', string) # Mapped length
            UL = re.findall('UL:i:(\d+)\t', string) # Unmapped length
            # count number of patterns
            TS = [TS.count(str(i)) for i in range(51)] #50 sequencing cycles
            ML = [ML.count(str(i)) for i in range(51)]
            UL = [UL.count(str(i)) for i in range(51)]
            # save the result
            # zip returns a generater in python3, so list() is necessary
            results.append(list(zip(TS, ML, UL))) 
    
    # get TS, ML, and UL
    
    TS = []
    for i in range(len(file_names)):
        TS.append([row[0] for row in results[i]])   
    out_file = os.path.join(result_dir, 'T-stretch_len_all.csv')
    header = 'T_Stretch_Length,'+','.join(file_names)
    TS = np.transpose(np.array(TS))
    TS = np.insert(TS, 0, values = np.arange(51), axis = 1) #insert T length
    np.savetxt(out_file, TS, header = header, delimiter=",")
    TS = pd.DataFrame(data=TS, columns=['T_Stretch_Length'] + file_names)
    
    ML = []
    for i in range(len(file_names)):
        ML.append([row[1] for row in results[i]])
    out_file = os.path.join(result_dir, 'T-stretch_len_algined.csv')
    header = 'T_Stretch_Length,'+','.join(file_names)
    ML = np.transpose(np.array(ML))
    ML = np.insert(ML, 0, values = np.arange(51), axis = 1) #insert T length
    np.savetxt(out_file, ML, header = header, delimiter=",")
    TS = pd.DataFrame(data=ML, columns=['T_Stretch_Length'] + file_names)
    
    UL = []
    for i in range(len(file_names)):
        UL.append([row[2] for row in results[i]])
    out_file = os.path.join(result_dir, 'T-stretch_len_unalgined.csv')
    header = 'T_Stretch_Length,'+','.join(file_names)
    UL = np.transpose(np.array(UL))
    UL = np.insert(UL, 0, values = np.arange(51), axis = 1) #insert T length
    np.savetxt(out_file, UL, header = header, delimiter=",")
    TS = pd.DataFrame(data=UL, columns=['T_Stretch_Length'] + file_names)
    
    return(TS, ML, UL)


################################################################################    
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
            ML = int(re.search('ML:i:(\d+)\t', line).group(1)) # int() !!!
            if ML >= min_length:
                fout1.write(line)
            else:
                fout2.write(line)
                
    fout1.close()
    fout2.close()
    

################################################################################  
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
################################################################################
def select_unique_fragments(random_NT_len, infile, outfile):
    '''Use the TS\d+[ATCG]{3} string in read name and chromosome, flag, and LM 
    tags to identify potential PCR duplicates.'''
    import itertools
    # make a set to save unique ids    
    unique_ids = set()
    # precompile the search patter
    pattern = re.compile('TS(\d+[ATCGN]{%s})'% random_NT_len)
    fout = open(outfile, 'w')
    # open infile, calculate the id
    with open(infile, 'r') as fin:
        print('Identifying uPASS in ' + infile)        
        for line in fin:
            if line[0] == '@':
                fout.write(line)                
                continue
            #print(line)
            # calculate the id for the read
            l = line.split()
            m = re.match(pattern, l[0])
            # calculate read length
            cigar = l[5]
            nums = re.split('[MDN]', re.sub('\d+I','', cigar))[:-1]
            covered =  sum(int(x) for x in nums)            
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

################################################################################
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
                    
################################################################################
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
   
  
################################################################################ 
def find_neighboring_indexes(v, max_distance):
    """ 
    v is like [[100395423, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 56]],...],
    sorted by the position 100395423. 
    Yield the indexes of lists containing neighboring positions
    """
    # a container to hold list of liss
    index_list = []
    for i in range(len(v))[:-1]:
        if v[i+1][0] - v[i][0] <= max_distance:
            #index_list.extend([i, i+1])
            index_list.append(i)
            continue
        if len(index_list) > 0:
            index_list.append(i) # append the last neighbor
            yield index_list
            index_list = []
            
################################################################################
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
    sample_number = len(cs_cluster[0][1])//2 ### different from python 2.7
    # check base case
    clustered = True
    for i in range(len(cs_cluster))[:-1]:
        if cs_cluster[i+1][0] - cs_cluster[i][0] <= max_distance:
            clustered = False
    if clustered == True:
        return True
    # when the clustering is not completed:
    else:
        # get max of normalized read numbers from all samples and the index for max       
        #print(cs_cluster)
        (m,i) = max((v,i) for i,v in enumerate((sum(pos_num[1][sample_number:]) for pos_num in cs_cluster)))
        # merge cs within max_distance from the cs with max read number
        left_index = i #leftmost position merged
        right_index = i #rightmost position merged
        for j, (pos, nums) in enumerate(cs_cluster):
            if abs(pos - cs_cluster[i][0]) <= max_distance and \
            abs(pos - cs_cluster[i][0]) > 0:
                # combine the read numbers for two CS, sample by sample
                # cs_cluster[i] is like [155603471, [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]]
                # nums is like [0, 0, 0, 0, 0, 5, 1, 0, 0, 0, 0, 0] 
                cs_cluster[i][1] = [sum(x) for x in zip(cs_cluster[i][1], nums)]
                cs_cluster[j][1] = 0 # don't delete cs_cluster[j]
                left_index = min(left_index, j)
                right_index = max(right_index, j)
        
        # devide and conqure        
        if left_index > 0:        
            cluster_neighboring_cleavage_sites(cs_cluster[:left_index], max_distance)
        if right_index < len(cs_cluster) - 1:
            cluster_neighboring_cleavage_sites(cs_cluster[right_index+1:], max_distance)   
            
################################################################################
def cluster_reads_in_dir(infolder, outfolder, 
                         file_pattern = 'pass.unique.sam', 
                         direction = 'reverse', 
                         max_distance = 24,
                         outfile = 'cluster.numbers.csv'):
    """
    Analyze a sam file and first generate a table with the following columns: 
    chromosome, strand, pos, num. Then recursively combine reads with LAP within 24 nt, 
    from pos with the highest read number to the pos with next highest read 
    number.
    
    Read sam files with file_pattern in infolder. For each sam file, build a dict
    readcounts like {'chr9:-:100395423':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 56]}.
    The [list of int] saves the number of reads from each sam_file. 
    dict is like sample_count[sample1(str)] = num_from_sample1(int)
       
    """
    
    import os    
    sam_files = [file_name for file_name in os.listdir(infolder) if file_name.endswith(file_pattern)]
    sam_file_num = len(sam_files)
        
    import re
    #re_pattern = re.compile('LM:i:(\d+)')
    re_pattern = re.compile('LM:i:(\d+)')
    
    import collections 
    # readcounts is a read id counter
    #example: {'chr9:-:100395423':[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 56]}
    readcounts = {} 
    # loop through sam files and get read numbers for each read id from each sam file  
    for s, sam_file in enumerate(sam_files): #s will be used as index
        # read sam file, calculate read id, and count number of reads
        #print('Reading ' + sam_file)        
        fin = open(os.path.join(infolder, sam_file), 'r')
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
            readcounts.setdefault(read_id, [0]*sam_file_num)[s] += 1
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
    normalized = df/df.sum(0)
    df = pd.concat([df, normalized], axis = 1, ignore_index = True)
    readcounts = df.T.to_dict('list')    
    # readcount.popitem() may return something like:
    #  ('chr11:+:31270274', [6.0, 7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
    #  0.0, 0.0, 9.1427902698768829e-07,3.0861678694165542e-06, 0.0, 0.0 ...]),
    # with the first 12 values as read counts and the next 12 values as normalized
    # read counts
    
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
    print('Writing to file...')
    with open(os.path.join(outfolder, outfile), 'w') as fout:
        # write header        
        #sample_string = ','.join([sam_file + '.num' for sam_file in sam_files])
        sample_string = ','.join([sam_file.split('.')[0] for sam_file in sam_files])
        fout.write('chromosome,strand,position,%s\n' %sample_string)
        # write read counts for each cluster
        for k in poscounts:
            chromosome, strand = k.split(':')
            for record in poscounts[k]:
                position = str(record[0])
                counts = ','.join([str(int(count)) for count in record[1][:sam_file_num]])
                fout.write('%s,%s,%s,%s\n'%(chromosome, strand, position, counts))
                
               
################################################################################
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
                
        
################################################################################
################################################################################
# sam -> bam -> bigwig
def make_url(project, batch, pass_dir, result_dir, genome_size, 
             genomeCoverageBed, norm_bedgraph, bedGraphToBigWig):
    for sam_file in [filename for filename in os.listdir(pass_dir) 
    if re.search('\.pass.sam$', filename)]:
        prefix = sam_file.split('.')[0]    
        sam_file = os.path.join(pass_dir, sam_file)
        #prefix = re.sub('\.sam$', '', sam_file)
        prefix = os.path.join(pass_dir, prefix)
        cmd = 'samtools view -uS ' + sam_file + ' | ' + 'samtools sort - ' + prefix # why need "-"?
        print(cmd)
        os.system(cmd)
        # I got the warning: [bam_header_read] EOF marker is absent. The input is probably truncated.
        # The warning should be ignored. Using samtools view -u gives uncompressed bam, and these do not have the EOF marker. 
    
        
        # bam to bigwig
        cmd = genomeCoverageBed + '-bg -split -ibam ' + prefix + '.bam ' + \
        '-strand + -g ' + genome_size + ' > ' + prefix + '.m.bedgraph' 
        print(cmd)
        os.system(cmd)
        cmd = genomeCoverageBed + '-bg -split -ibam ' + prefix + '.bam ' + \
        '-strand - -g ' + genome_size + ' > ' + prefix + '.p.bedgraph' 
        print( cmd)
        os.system(cmd)
        
        # normolize bedgraph counts
        #cmd = 'wc -l ' + sam_file
        #from subprocess import check_output
        #cmd_out = check_output(cmd, shell=True)
        #totalReadNum = int(cmd_out.split(' ')[0]) - 55
        #cmd = norm_bedgraph + '-t ' + str(totalReadNum) + ' -i "' + prefix + \
        #'.p.bedgraph ' + prefix + '.m.bedgraph"' #+ ' -m "' + prefix + '.p.bedgraph"'
        #print( cmd)
        #os.system(cmd)
        
        cmd = 'grep -v @ ' + sam_file + ' | wc -l'
        from subprocess import check_output
        cmd_out = check_output(cmd, shell=True)
        totalReadNum = int(cmd_out.split(' ')[0]) 
        cmd = norm_bedgraph + '-t ' + str(totalReadNum) + ' -i "' + prefix + \
        '.p.bedgraph ' + prefix + '.m.bedgraph"' #+ ' -m "' + prefix + '.p.bedgraph"'
        print(cmd)
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
        cmd = 'rm ' + prefix + '*.bam'
        os.system(cmd)
        cmd = 'rm ' + prefix + '*bedgraph*'
        os.system(cmd)
    # copy the bigwig files to the http-enabled intron server
    #cmd = 'scp ' + os.path.join(pass_dir, '*.bw') + ' zhengdh@intron.njms.rutgers.edu:~/../www/zhengdh/bigwig/' + project 
    #os.system(cmd)
    
    # record sample names for the bw files
    bw_files = [filename for filename in os.listdir(pass_dir) 
    if filename.endswith('.bw')]
    bw_samples = sorted(list(set([filename.split('.')[0] for filename in bw_files])))
    #bw_samples.sort(key = lambda x: x[-1])
    
    #cmd = 'rm ' + os.path.join(pass_dir, '*.bw')
    
    strand2str = {'+': 'p', '-':'m'}
    f = open(os.path.join(result_dir, 'bigwigCaller.txt'), 'w')
    for strand in ['+', '-']:
        for (i, sample) in enumerate(bw_samples):
            #name = sample.replace('spikein', '') + '_' + batch
            name = sample
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
             + strand + '" description=\"3\'READS+\"' + \
             ' bigDataUrl=http://intron.njms.rutgers.edu/zhengdh/bigwig/' + \
             project + '/' + batch + '/' + sample + '.' + strand2str[strand] + '.bw' + '\n'
            f.write(track)
    f.close()     


################################################################################
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

################################################################################
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
