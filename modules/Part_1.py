## 1.1. Load packages
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
from time import sleep
import sys
# sys.path.append('/home/dinghai/dev/3-prime-READS-plus/modules')
import PASS as ps





## 1.2. Set configurations specific to this project and experiment
GENOME = 'mm9' # Genome name, currently either mm9 or hg19

SAMPLE_DESCRIPTION_FILE = 'sample_description.csv' # Sample discription file.
    
SAVE_SPACE = True # To save disk space, delete intermediate files asap.    

SPIKE_IN = None # None if no spike-in, or a 4-charactor string for identifying spike-in RNA. 

# SPIKE_IN = 'ychr' # yeast RNA spike-in (if yeast chromosomes are named as 'ychr*')

DATA_SOURCE = 'intron.njms.rutgers.edu' # Either a server name or path to a local directionary. 

# DATA_SOURCE = '../../../../tmp' # local directionary

if not Path(DATA_SOURCE).exists(): # DATA_SOURCE = 'intron.njms.rutgers.edu'
    DATA_USER = 'zhengdh'  # Use your own user name



## 1.3. Set general configurations
# Directory containing index files for different genomes
GENOME_INDEX_DIR = '/home/dinghai/projects/fud/star_index'

# Specific genome index for STAR
genome_index = Path(GENOME_INDEX_DIR)/GENOME
assert genome_index.exists(), f"Genome index folder {str(genome_index)} does not exist."

# Directory containing fasta files for different genomes
GENOME_FASTA_DIR = '/home/dinghai/projects/fud/ucsc/genomes' 

# Specific genome sequence in fasta format for identifying PASS reads
genome_dir = Path(GENOME_FASTA_DIR)/GENOME
assert genome_dir.exists(), f"Genome fasta file folder {str(genome_dir)} does not exist."

genome_size = Path(GENOME_FASTA_DIR)/f'{GENOME}.chrom.sizes'
assert genome_size.exists(), f"Genome size file {str(genome_size)} does not exist."

# Configurations for using UCSC genome browser to visualize PASS and nonPASS reads   
MAKE_URL = True # Set to True if you have access to an http server.        
if MAKE_URL:
    HTTP_USER = 'zhengdh' # Use your own user name
    HTTP_SERVER = 'intron.njms.rutgers.edu' # Use your own http server name

# Maximum number of processes for parallel computing
MAX_WORKERS = 8  

NUM_COLUMNS = 4  # Number of columns for fastq QC plots

# Minimum read lenth after removing 5' adapter
CUTADAPT_MINLEN = 25  

# Minimum MAPQ score for alignments in sam files
MIN_MAPQ = 10  

# Sequencing direction relative to mRNA sequence. 
SEQ_DIRECTION = 'reverse' # Either 'reverse' (antisense) or 'forward' (sense)

# Maximum distance between clustered reads in the same cluster
MAX_DISTANCE = 24

# Project directory
exp_dir = Path('..')
PROJECT, EXPERIMENT = exp_dir.cwd().parts[-3:-1]

# Data directories
data_dir = exp_dir/'data'
fastq_dir = data_dir/'fastq'
sam_dir = data_dir/'sam'

# Result directory
result_dir = exp_dir/'results'

# Create directories
data_dir.mkdir(parents=True, exist_ok=True)
fastq_dir.mkdir(parents=True, exist_ok=True)
sam_dir.mkdir(parents=True, exist_ok=True)
result_dir.mkdir(parents=True, exist_ok=True)





## 2. Read sample annotation and choose colors for UCSC genome browser tracks
sample_description = pd.read_csv(SAMPLE_DESCRIPTION_FILE, dtype = 'str')

if not set(['folder_string', 'file_string', 'sample', 'read_type']) <= \
       set(sample_description.columns):
    msg = (f'The sample discription file should at least have the following columns:\n'
           f'"folder_string": unique string in the name of the folder '
           f'"file_string": string in the fastq files.\n'
           f'"sample": user defined sample name.'
           f'"read_type": either "single" or "paired"'
          )
    print(msg)
    
# Remove potential whitespace 
sample_description = sample_description.apply(lambda x: x.str.strip(), axis=1)
                      
if 'track_color' in sample_description.columns:
    sample_description = sample_description.sort_values(by=['folder_string', 'track_color'])
else:
    sample_description = sample_description.sort_values(by=['folder_string'])
    
sample_description = sample_description.set_index('sample')

print("\nSample description:")
display(sample_description)    

print('The following colors will be used for UCSC genome browser tracks:')

if 'track_color' in sample_description.columns:
    sample_description.track_color = sample_description.track_color.astype(int)
    color_max = sample_description.track_color.max(axis=0)
else:
    color_max = sample_description.shape[0]
    
sns.palplot(sns.color_palette("colorblind", color_max))





## 3. Find, download, merge, and rename fastq files
for folder_string in sample_description.folder_string.unique():
    fastq_subdir = fastq_dir/folder_string
    fastq_subdir.mkdir(parents=True, exist_ok=True)
    df = sample_description[sample_description.folder_string == folder_string]
    
    # DATA_SOURCE is something like 'intron.njms.rutgers.edu'
    if not re.search(r'\/', DATA_SOURCE):
        
        # Find folders in the warehouse that may contain the fastq files
        cmd = f'ssh {DATA_USER}@{DATA_SOURCE} find ../database* -type d -name {folder_string}*'
        print(f'\nSearching folders with {folder_string} in their names on {DATA_SOURCE}')
        print(f'\nUsing command "{cmd}"')
        fastq_folders = check_output(cmd, shell=True).decode().split('\n')
        print('\nFastq folders found:')
        print('\n'.join(fastq_folders))

        # Find and download fastq files in the folders
        for file_string in df.file_string:
            sample_name = df[df.file_string == file_string].index[0]
            # Some of the folders may not contain fastq data
            for fastq_folder in fastq_folders:
                cmd = (f'ssh {DATA_USER}@{DATA_SOURCE} find {fastq_folder}'
                       f' -type f -name "*{file_string}*R1*"')
                # print(cmd)
                fastq_files = check_output(cmd, shell=True).decode().replace('\n', " ")
                # Escape special characters
                fastq_files = fastq_files.replace("&", "\&").replace(":", "\:")
                # Skip the current fastq_folder if fastq_files have already been found
                if not fastq_files == '': 
                    break
                    
            # Download fastq files from internal server
            print('\nDownloading fastq files ....')
            cmd = f'scp -r {DATA_USER}@{DATA_SOURCE}:"{fastq_files}" {fastq_subdir}'
            print(cmd)
            os.system(cmd)
            
    # If fastq files are in a local folder:    
    elif Path(DATA_SOURCE).exists():
        cmd = f'find {DATA_SOURCE} -type d -name "{folder_string}*"'
        fastq_subdir = check_output(cmd, shell=True).decode().strip() 
        print('\nLocal fastq folders found:', fastq_subdir)

    # Merge fastq files in the folder
    for file_string in df.file_string:
        sample_name = df[df.file_string == file_string].index[0]

        # Get file names
        R1_fastq_gz = [str(filename) for filename in Path(fastq_subdir).iterdir() 
                       if re.search(file_string + '.*fastq.gz$', str(filename))]
        
        R1_fastq = [str(filename) for filename in Path(fastq_subdir).iterdir() 
                    if re.search(file_string + '.*fastq$', str(filename))]

        # In case that both .fastq and .fastq.gz files exist for the same sample
        R1_fastq_gz = [filename for filename in R1_fastq_gz if not 
                       filename.replace('.gz', '') in R1_fastq]

        # Merge files
        ps.merge_and_rename(R1_fastq_gz, f'{sample_name}.fastq', fastq_dir)
        ps.merge_and_rename(R1_fastq, f'{sample_name}.fastq', fastq_dir)
    
     
    print('\nDone!')

# The following code determines number of workers for parallel computing.
fastq_files = sorted([str(filename.absolute()) for filename in fastq_dir.glob('*.fastq')])

if len(fastq_files) == 0:
    print(f"\nPleae put fastq files in the {str(fastq_dir)} folder first.")    
    
WORKERS = min(MAX_WORKERS, len(fastq_files))

print(f'\n{WORKERS} workers will be used for all (except mapping) parallel computing.')

# Count raw fastq records in fastq files in parallel:
with mp.Pool(processes=WORKERS) as pool:
    raw_fastq_counts = pool.map(ps.count_fastq, fastq_files)
    
raw_fastq_counts = pd.DataFrame(raw_fastq_counts, columns = ["File_Name", "Raw"])
raw_fastq_counts["File_Name"] = raw_fastq_counts["File_Name"].str.replace('.fastq', '') 
raw_fastq_counts = raw_fastq_counts.set_index("File_Name")

print('\nRaw fastq count in each file:')
raw_fastq_counts





## 4. FASTQ QC and calculate length of random nucleotides in 3' ligation adapter
# Currently only implemented in the Jupyter Notebook pipeline.
random_NT_lens 




## 5. Clip adapters
cmds = [] 
for sample_name in sample_description.index:
    fastq_file = str(rawfastq_dir/f'{sample_name}.fastq')
    read_type = sample_description.loc[sample_name, 'read_type']
    if read_type == 'single':
        cmd = (f'cutadapt -a NNNNTGGAATTCTCGGGTGCCAAGG -n 1 -O 10 -m {CUTADAPT_MINLEN} -f fastq -q 10 {fastq_file} '
               f'-o {fastq_file.replace(".fastq", ".cut.fastq").replace("rawfastq", "fastq")} '
               f'&& rm {fastq_file}')
    elif read_type == 'paired':
        cmd = (f'cutadapt -u -75 -a NNNNTGGAATTCTCGGGTGCCAAGG -n 1 -O 10 -m {CUTADAPT_MINLEN} -f fastq '
               f'-q 10 {fastq_file} -o {fastq_file.replace(".fastq", ".cut.fastq").replace("rawfastq", "fastq")} '
               f'&& rm {fastq_file}')
    cmds.append(cmd)

print("\nRemoving ligation adapters from fastq reads in the following files:\n") 
print('\n'.join([fastq_file.split('/')[-1] for fastq_file in fastq_files]))

with mp.Pool(processes = WORKERS) as pool: 
    pool.map(os.system, cmds)
    
print('\nDone!')





## 6. Trim 5' Ts
cutadapt_outputs = sorted([str(filename.absolute()) for filename in fastq_dir.glob('*cut.fastq')])

print("\nTrimming 5' T-stretches from fastq reads in the following files:\n") 
print('\n'.join([fastq_file.split('/')[-1] for fastq_file in cutadapt_outputs]))

with mp.Pool(processes = WORKERS) as pool:
    trimmed_fastq_counts = pool.starmap(ps.trim_write_count_fastq, 
                            list(zip(cutadapt_outputs, list(random_NT_lens), [4]*len(random_NT_lens))))
    
if SAVE_SPACE:
    cmd = 'rm ' + ' '.join(cutadapt_outputs)
    os.system(cmd)    

print('\nDone!')

trimmed_fastq_counts = pd.DataFrame(trimmed_fastq_counts, columns = ["File_Name", "Cut", "Trimmed"])
trimmed_fastq_counts["File_Name"] = trimmed_fastq_counts["File_Name"].str.replace('.cut.fastq', '')  
trimmed_fastq_counts = trimmed_fastq_counts.set_index("File_Name")
fastq_counts = pd.concat([raw_fastq_counts, trimmed_fastq_counts], axis = 1)

print('\nFastq counts in each file:')
fastq_counts





## 7. Map the reads to the genome using STAR
trimmed_fastq_files = sorted([str(filename.absolute()) 
                              for filename in fastq_dir.glob('*trimmed.fastq')])
# Commands for mapping reads
cmds = [(f'STAR --runThreadN {os.cpu_count()} '
         f'--genomeDir {genome_index} --genomeLoad LoadAndKeep '
         f'--readFilesIn {fastq_file} --alignEndsType EndToEnd --outFilterType BySJout '
         f'--outFilterMultimapNmax 10 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 '
         f'--outSAMmultNmax 1 '
         f'--alignIntronMin 20 --alignIntronMax 1000000 --outFileNamePrefix '
         f'{fastq_file.replace("/fastq/", "/sam/").replace("cut.trimmed.fastq", "")}') 
         for fastq_file in trimmed_fastq_files]

# Commands for removing fastq files
cmds_rm = ['rm ' + fastq_file for fastq_file in trimmed_fastq_files]  

print(f'\nMapping reads from the following files to {GENOME}:\n') 
print('\n'.join([fastq_file.split('/')[-1] for fastq_file in trimmed_fastq_files]))

# Map reads and then delete fastq files immediately
for i in range(len(cmds)):
    os.system(cmds[i])    
    if SAVE_SPACE:
        os.system(cmds_rm[i])
    
# Remove genome from memory. Note "--genomeDir" is still needed!
cmd = f'STAR --genomeDir {genome_index} --genomeLoad Remove'
os.system(cmd)

print('\nDone!')





## 8.1. Identfy PASS reads
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

    with open(sam_file, 'r') as fin:
        for line in fin:
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
                downstream_seq = ps.get_seq(chromosome, strand,
                                         start = position + covered,
                                         end = position + covered + t_stretch_len - 1,
                                         genome = genome)
                lap = position + covered - 1
            elif (direction.lower() == 'reverse' and flag == '0') or \
                    (direction.lower() == 'forward' and flag == '16'):
                strand = '-'
                downstream_seq = ps.get_seq(chromosome, strand,
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
                match = re.match(ps.reverse_complement(elements[1]), downstream_seq)
                if match:
                    downstream_seq = downstream_seq[len(elements[1]):]
                    t_stretch_len -= len(elements[1])
                    if strand == '+':
                        lap += len(elements[1])
                    if strand == '-':
                        lap -= len(elements[1])
            
            # If the t_stretch_len is reduced enough, no need to analyze downstream sequence
            #assert t_stretch_len >= 0
            
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
            # If t_stretch_len == 1, the following block will also run.
            else:
                match = re.match('A*', downstream_seq)
                matched_len = len(match.group())  
                nonpass_file.write(line.strip() + '\tML:i:%d\tUL:i:%d\tLM:i:%d\n'
                                      % (matched_len, t_stretch_len - matched_len, lap))
    pass_file.close()
    nonpass_file.close()
    if spike_in: spike_in_file.close()
    
    os.system('rm ' + sam_file)

# List of sam files as input
sam_files = sorted(str(sam_file) for sam_file in sam_dir.glob('*Aligned.out.sam'))
l = len(sam_files)

print('\nPicking PASS reads from the following files:\n') 
print('\n'.join([sam_file.split('/')[-1] for sam_file in sam_files]))

# Genome loaded into memory is available for all processes
genome = ps.load_fasta_genome(genome_dir)

# Try parallel computing, which requires more disk space
try:
    with mp.Pool(processes = WORKERS) as pool:
        pool.starmap(split_sam, zip(sam_files, [MIN_MAPQ]*l, [SEQ_DIRECTION]*l, [SPIKE_IN]*l))
except OSError:
    os.system(f'rm {str(sam_dir)}/*pass')
    if SPIKE_IN: os.system(f'rm {str(sam_dir)}/*.spike_in')
    print("\nNo space left on device. Trying again with parallel computing turned off.")
    for sam_file in sam_files:
        print("Processing ", sam_file)
        split_sam(sam_file = sam_file)

# Release memory    
del genome  

print('\nDone!') 





## 8.2. Pick unique PASS reads
pass_files = sorted(str(pass_file) for pass_file in sam_dir.glob('*Aligned.out.pass'))
print('\nPicking unique PASS reads from the following files:\n') 
print('\n'.join([pass_file.split('/')[-1] for pass_file in pass_files]))

with mp.Pool(processes = WORKERS) as pool:
    pool.starmap(ps.pick_unique_pass, list(zip(pass_files, list(random_NT_lens))))
    
print('\nDone!')  





## 9. Cluster all PASS or unique PASS reads
# Cluster all PASS reads
pass_files = sorted(str(pass_file) for pass_file in sam_dir.glob('*Aligned.out.pass'))
pass_reads_output = result_dir/'clusters.using.all.reads.csv'

ps.cluster_pass_reads(pass_files, output = pass_reads_output,
                      direction = SEQ_DIRECTION, max_distance = MAX_DISTANCE)

# Cluster unique PASS reads
upass_files = sorted(str(pass_file) for pass_file in sam_dir.glob('*Aligned.out.unique.pass'))
upass_reads_output = result_dir/'clusters.using.unique.reads.csv'

ps.cluster_pass_reads(upass_files, output = upass_reads_output,
                      direction = SEQ_DIRECTION, max_distance = MAX_DISTANCE)





## 10. Read number statistics
with mp.Pool(processes = WORKERS) as pool: 
    pass_counts = pool.map(ps.count_sam, pass_files)
    
pass_counts = pd.DataFrame(pass_counts, columns = ['Sample', 'PASS'])
pass_counts.iloc[:,0] = pass_counts.iloc[:,0].str.replace('.Aligned.out.pass', '')
pass_counts = pass_counts.set_index('Sample')   


nonpass_files = sorted(str(nonpass_file) for nonpass_file in sam_dir.glob('*Aligned.out.nonpass'))
with mp.Pool(processes = WORKERS) as pool: 
    nonpass_counts = pool.map(ps.count_sam, nonpass_files)
    
nonpass_counts = pd.DataFrame(nonpass_counts, columns = ['Sample', 'nonPASS'])
nonpass_counts.iloc[:,0] = nonpass_counts.iloc[:,0].str.replace('.Aligned.out.nonpass', '')
nonpass_counts = nonpass_counts.set_index('Sample')  

    
with mp.Pool(processes = WORKERS) as pool: 
    upass_counts = pool.map(ps.count_sam, upass_files)
    
upass_counts = pd.DataFrame(upass_counts, columns = ['Sample', 'uPASS'])
upass_counts.iloc[:,0] = upass_counts.iloc[:,0].str.replace('.Aligned.out.unique.pass', '')
upass_counts = upass_counts.set_index('Sample')    


if SPIKE_IN: 
    spike_in_files = sorted(str(nonpass_file) for nonpass_file in sam_dir.glob('*Aligned.out.spike_in'))
    with mp.Pool(processes = WORKERS) as pool: 
        spike_in_counts = pool.map(ps.count_sam, spike_in_files)
    spike_in_counts = pd.DataFrame(spike_in_counts, columns = ['Sample', 'Spike-in'])
    spike_in_counts.iloc[:,0] = spike_in_counts.iloc[:,0].str.replace('.Aligned.out.Spike-in', '')
    spike_in_counts = spike_in_counts.set_index('Sample')


mapped_counts = pd.DataFrame(pd.concat([pass_counts, nonpass_counts], 
                                       axis=1).sum(axis=1), columns = ['Mapped'])

if SPIKE_IN: 
    count_stats = pd.concat([fastq_counts, mapped_counts, pass_counts, 
                             nonpass_counts, upass_counts, spike_in_counts], axis=1)
else:
    count_stats = pd.concat([fastq_counts, mapped_counts, pass_counts, 
                             nonpass_counts, upass_counts], axis=1)

count_stats['Cut%'] = round((1-count_stats.Cut/count_stats.Raw) * 100, 1)
count_stats['Trimmed%'] = round((1-count_stats.Trimmed/count_stats.Cut) * 100, 1)
count_stats['MAP%'] = round(count_stats.Mapped/count_stats.Cut * 100, 1)
count_stats['PASS%'] = round(count_stats.PASS/count_stats.Mapped * 100, 1)
count_stats['uPASS%'] = round(count_stats.uPASS/count_stats.PASS * 100, 1)

count_stats.to_csv(result_dir/'ReadStats.csv')





## 11. T-stretch statistics
pass_files = sorted([str(sam_file) for sam_file in sam_dir.glob('*.Aligned.out.pass')])
nonpass_files = sorted([str(sam_file) for sam_file in sam_dir.glob('*.Aligned.out.nonpass')])

pass_TS = ps.summarize_5T_stretch(sam_files = pass_files, processes = WORKERS, max_TS = 25)
nonpass_TS = ps.summarize_5T_stretch(sam_files = nonpass_files, processes = WORKERS, max_TS = 25)

pass_TS.columns = pass_TS.columns + '.PASS'
nonpass_TS.columns = nonpass_TS.columns + '.nonPASS'

TS = pd.concat([pass_TS, nonpass_TS], axis = 1)
TS.to_csv(result_dir/'TStats.csv')





## 12. Create UCSC genome browser tracks
pass_files = sorted([str(sam_file) for sam_file in sam_dir.glob('*.Aligned.out.pass')])

# Make UCSC genome browser tracks for PASS reads
ps.make_url(project = PROJECT, experiment = EXPERIMENT, sam_dir= sam_dir, sam_files = pass_files, 
            genome_size = str(genome_size), sample_description = sample_description, processes = WORKERS, 
            keep_bam = False, bigDataUrl = f'http://{HTTP_SERVER}/{HTTP_USER}/bigwig/')

# Prepare for uploading files
pass_bw_dir = sam_dir/PROJECT/EXPERIMENT/'PASS'
pass_bw_dir.mkdir(parents=True, exist_ok=True)

os.system(f'mv {str(sam_dir)}/*.bw {str(pass_bw_dir)}')
os.system(f'mv {str(sam_dir)}/bigwigCaller.txt {str(pass_bw_dir)}')

# Upload the bigwig and track definition files to the http-enabled server
cmd = f'ssh {HTTP_USER}@{HTTP_SERVER} mkdir ~/../www/{HTTP_USER}/bigwig/{PROJECT}'
os.system(cmd)

cmd = f'scp -r {str(sam_dir/PROJECT/EXPERIMENT)} {HTTP_USER}@{HTTP_SERVER}:~/../www/{HTTP_USER}/bigwig/{PROJECT}'
if not os.system(cmd):
    print('bigWig files have been uploaded to http server.')
    os.system(f'rm -r {str(sam_dir/PROJECT)}')
    if SAVE_SPACE: os.system(f'rm -r {str(sam_dir)}')
    
    pass_url = (f'https://genome.ucsc.edu/cgi-bin/hgTracks?db={GENOME}'
                f'&position=chr17%3A35092728-35095537&hgct_customText='
                f'http://{HTTP_SERVER}/{HTTP_USER}/bigwig/{PROJECT}/'
                f'{EXPERIMENT}/PASS/bigwigCaller.txt'
                )    

    print('URL for PASS reads:')
    print(pass_url, '\n')