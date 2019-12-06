result_dir = '.'
sub_dir = file.path(result_dir, "bedgraph")
dir.create(sub_dir)
    
# The dataframe "pas" must contain the following columns:
# chr, strand, pA_pos, *_rpm (RPM for different samples)
pas =  read.csv(file.path(result_dir, "pas.csv"), as.is=T)

## Generate bedGraph files
require(dplyr)
for(rpm in grep("_rpm$", names(pas), value=T)){
    p_file_name = sub("_rpm", ".p.bedgraph", rpm)
    m_file_name = sub("_rpm", ".m.bedgraph", rpm)
    
    pas %>%
      dplyr::filter(strand == "+") %>%
      mutate(pA_pos = pA_pos - 1) %>%
      mutate(start_pos = pA_pos - 1) %>%
      # bedGraph columns: "chrom chromStart chromEnd dataValue"
      dplyr::select(chr, start_pos, pA_pos, eval(rpm)) %>%
      distinct() %>%
      write.table(file = file.path(sub_dir, p_file_name), 
                  row.names=F, col.names=F, quote=F, sep="\t")
    
    pas %>%
      dplyr::filter(strand == "-") %>%
      mutate(pA_pos = pA_pos - 1) %>%
      mutate(end_pos = pA_pos + 1) %>%
      # bedGraph columns: "chrom chromStart chromEnd dataValue"
      dplyr::select(chr, pA_pos, end_pos, eval(rpm)) %>%
      mutate_at(vars(matches("_rpm")), list(~ -1*.)) %>%
      distinct() %>%
      write.table(file = file.path(sub_dir, m_file_name), 
                  row.names=F, col.names=F, quote=F, sep="\t")
}


## Create bigWig files
chrom_sizes = '/home/dinghai/projects/fud/ucsc/genomes/mm9.chrom.sizes'
require(R.utils)
for(bedgraph_file in list.files(getAbsolutePath(sub_dir), full.names = T)){
    # Sort records
    cmd = paste("bedSort", bedgraph_file, bedgraph_file, sep=" ")
    print(cmd)
    system(cmd)
    
    # Create bigWig files
    cmd = paste("bedGraphToBigWig", bedgraph_file, chrom_sizes, 
                sub("bedgraph$", "bw", bedgraph_file), sep=" ")
    print(cmd)
    system(cmd)
    
    # Delete bedgraph files
    cmd = paste("rm", bedgraph_file, sep = " ")
    print(cmd)
    system(cmd)
}
