library(MASS) 
library(ggrepel) 
require(tidyr)
require(dplyr)
require(ggplot2)
theme_set(theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5))) 


findPotentialStartsAndStops = function(sequence){
  library("Biostrings")
  # Define a vector with the sequences of potential start and stop codons
  codons = c("ATG", "TAA", "TAG", "TGA")
  sequence = toupper(sequence)
  # Find the number of occurrences of each type of potential start or stop codon
  for(i in 1:length(codons)){
    codon = codons[i]
    # Find all occurrences of codon "codon" in sequence "sequence"
    occurrences = matchPattern(codon, sequence)
    # Find the start positions of all occurrences of "codon" in sequence "sequence"
    codonpositions = start(occurrences)
    # Find the total number of potential start and stop codons in sequence "sequence"
    numoccurrences = length(codonpositions)
    if (i == 1){
      # Make a copy of vector "codonpositions" called "positions"
      positions = codonpositions
      # Make a vector "types" containing "numoccurrences" copies of "codon"
      types = rep(codon, numoccurrences)
    }else{
      # Add the vector "codonpositions" to the end of vector "positions":
      positions   = append(positions, codonpositions, after=length(positions))
      # Add the vector "rep(codon, numoccurrences)" to the end of vector "types":
      types       = append(types, rep(codon, numoccurrences), after=length(types))
    }
  }
  # Sort the vectors "positions" and "types" in order of position along the input sequence:
  indices = order(positions)
  positions = positions[indices]
  types = types[indices]
  # Return a list variable including vectors "positions" and "types":
  mylist = list(positions,types)
  return(mylist)
}


findORFsinSeq = function(sequence){
  require(Biostrings)
  # Make vectors "positions" and "types" containing information on the positions of ATGs 
  # in the sequence:
  mylist = findPotentialStartsAndStops(sequence)
  positions = mylist[[1]]
  types = mylist[[2]]
  # Make vectors "orfstarts" and "orfstops" to store the predicted start and stop codons of ORFs
  orfstarts = numeric()
  orfstops = numeric()
  # Make a vector "orflengths" to store the lengths of the ORFs
  orflengths = numeric()
  # Print out the positions of ORFs in the sequence:
  # Find the length of vector "positions"
  numpositions = length(positions)
  # There must be at least one start codon and one stop codon to have an ORF.
  if (numpositions >= 2){
    for (i in 1:(numpositions-1)){
      posi = positions[i]
      typei = types[i]
      found = 0
      while (found == 0){
        for (j in (i+1):numpositions){
          posj  = positions[j]
          typej = types[j]
          posdiff = posj - posi
          posdiffmod3 = posdiff %% 3
          # Add in the length of the stop codon
          orflength = posj - posi + 3
          if (typei == "ATG" && (typej == "TAA" || typej == "TAG" || typej == "TGA") && posdiffmod3 == 0){
            # Check if we have already used the stop codon at posj+2 in an ORF
            numorfs = length(orfstops)
            usedstop = -1
            if (numorfs > 0){
              for (k in 1:numorfs){
                orfstopk = orfstops[k]
                if (orfstopk == (posj + 2)) { usedstop = 1 }
              }
            }
            if (usedstop == -1){
              orfstarts = append(orfstarts, posi, after=length(orfstarts))
              orfstops = append(orfstops, posj+2, after=length(orfstops)) # Including the stop codon.
              orflengths = append(orflengths, orflength, after=length(orflengths))
            }
            found = 1
            break
          }
          if (j == numpositions) { found = 1 }
        }
      }
    }
  }
  # Sort the final ORFs by start position:
  indices = order(orfstarts)
  orfstarts = orfstarts[indices]
  orfstops = orfstops[indices]
  # Find the lengths of the ORFs that we have
  orflengths = numeric()
  numorfs = length(orfstarts)
  for (i in 1:numorfs)
  {
    orfstart = orfstarts[i]
    orfstop = orfstops[i]
    orflength = orfstop - orfstart + 1
    orflengths = append(orflengths,orflength,after=length(orflengths))
  }
  mylist = list(orfstarts, orfstops, orflengths)
  names(mylist) = c("starts", "ends", "lengths")
  return(mylist)
}


seeFastq = function(fastq, batchsize=10000, klength=8){ 
  ## Random sample N reads from fastq file (N=batchsize)
  f = FastqSampler(fastq, batchsize)
  fq = yield(f)
  nReads = f$status()[["total"]] # Total number of reads in fastq file
  close(f)
  
  ## If reads are not of constant width then inject them into a matrix pre-populated with 
  ## N/NA values and of dimensions N_rows = number_of_reads and N_columns = length_of_longest_read. 
  if(length(unique(width(fq))) == 1) {
    q = as.matrix(PhredQuality(quality(fq))) 
    s = as.matrix(sread(fq))
  } else {
    mymin = min(width(fq)); mymax = max(width(fq))
    s = matrix("N", length(fq), mymax)
    q = matrix(NA, length(fq), mymax)
    for(i in mymin:mymax) {
      index = width(fq)==i
      if(any(index)) {
        s[index, 1:i] = as.matrix(DNAStringSet(sread(fq)[index], start=1, end=i))
        q[index, 1:i] = as.matrix(PhredQuality(quality(fq))[index]) 
      }
    }
  }
  s[s=="N"] = NA
  row.names(q) = paste("s", 1:length(q[,1]), sep=""); colnames(q) = 1:length(q[1,])
  
  ## (A) Per cycle quality box plot
  ## Generate box plot from precomputed stats	
  bpl = boxplot(q, plot=FALSE)
  astats = data.frame(bpl$names, t(matrix(bpl$stats, dim(bpl$stats))))
  colnames(astats) = c("Cycle", "min", "low", "mid", "top", "max")
  astats[,1] = factor(astats[,1], levels=unique(astats[,1]), ordered=TRUE)  
  
  ## (B) Per cycle base proportion
  bstats = apply(s, 2, function(x) table(factor(x, levels=c("A", "C", "G", "T")))) 
  colnames(bstats) = 1:length(bstats[1,])
  bstats = t(apply(bstats, 1, function(x) x/colSums(bstats)))
  bstats = data.frame(Nuc=rownames(bstats), bstats) 
  convertDF = function(df=df, mycolnames) { 
    myfactor = rep(colnames(df)[-1], each=length(df[,1]))
    mydata = as.vector(as.matrix(df[,-1]))
    df = data.frame(df[,1], mydata, myfactor)
    colnames(df) = mycolnames
    return(df) 
  }
  bstats = convertDF(bstats, mycolnames=c("Base", "Frequency", "Cycle"))
  bstats[,3] = as.numeric(gsub("X", "", bstats[,3]))
  bstats[,3] = factor(bstats[,3], levels=unique(bstats[,3]), ordered=TRUE) 
  
  ## (C) Per cycle average quality of each base type
  A = q; A[s %in% c("T", "G", "C")] = NA; A = colMeans(A, na.rm=TRUE)
  T = q; T[s %in% c("A", "G", "C")] = NA; T = colMeans(T, na.rm=TRUE)
  G = q; G[s %in% c("T", "A", "C")] = NA; G = colMeans(G, na.rm=TRUE)
  C = q; C[s %in% c("T", "G", "A")] = NA; C = colMeans(C, na.rm=TRUE)
  cstats = data.frame(Quality=c(A, C, G, T), Base=rep(c("A", "C", "G", "T"), each=length(A)), Cycle=c(names(A), names(C), names(G), names(T)))
  cstats[,3] = factor(cstats[,3], levels=unique(cstats[,3]), ordered=TRUE) 
  
  ## (D) Relative K-mer Diversity 
  dna = sread(fq)
  loopv = 1:(min(width(dna)) - (klength-1))
  kcount = sapply(loopv, function(x) length(unique(DNAStringSet(start=x, end=x+klength-1, dna))))	
  reldiv = kcount/(5^klength) # 5 instead of 4 because of Ns
  reldiv = c(rep(NA, klength-1), reldiv) # Adds dummy NAs to align with sequencing cycles 
  names(reldiv) = 1:length(reldiv)	
  dstats = data.frame(RelDiv=reldiv, Method=rep(c(1), each=length(reldiv)), Cycle=names(reldiv))
  dstats[,3] = factor(dstats[,3], levels=unique(dstats[,3]), ordered=TRUE) 
  
  ## (E) Number of reads where all Phred scores are above a minimum cutoff
  ev = c("0"=0, "1"=10, "2"=20, "3"=30, "4"=40)
  edf = sapply(ev, function(x) sapply(as.numeric(names(ev)), function(y) sum(rowSums(q >= x, na.rm=TRUE) >= (rowSums(!is.na(q))-y))))
  rownames(edf) = names(ev); colnames(edf) = ev
  edf = edf/max(edf)*100
  edf = data.frame(Percent=paste(">", colnames(edf), sep=""), t(edf), check.names=FALSE)
  estats = convertDF(edf, mycolnames=c("minQuality", "Percent", "Outliers"))
  estats[,1] = factor(estats[,1], levels=unique(estats[,1]), ordered=TRUE)
  estats[,3] = factor(estats[,3], levels=unique(estats[,3]), ordered=TRUE)
  
  ## (F) Distribution of mean quality of reads
  qv = table(round(rowMeans(q)))[as.character(0:max(q, na.rm=TRUE))]
  qv[is.na(qv)] = 0; names(qv) = 0:max(q, na.rm=TRUE)
  fstats = data.frame(Quality=names(qv), Percent=qv)
  fstats[,2] = fstats[,2] / length(q[,1]) * 100
  fstats[,1] = factor(fstats[,1], levels=unique(fstats[,1]), ordered=TRUE)
  
  ## (G) Read length distribution
  l = rep(0, max(width(fq))); names(l) = 1:length(l)
  lv = table(width(fq))
  l[names(lv)] = lv
  gstats = data.frame(Cycle=names(l), Percent=l)
  gstats[,2] = gstats[,2] / sum(gstats[,2]) * 100 
  gstats[,1] = factor(gstats[,1], levels=unique(gstats[,1]), ordered=TRUE)
  
  ## (H) Read occurrence distribution
  qa1 = qa(fq, basename(fastq))
  hstats = qa1[["sequenceDistribution"]][,1:2]
  hstats = data.frame(nOccurrences=hstats[,1], Percent=hstats[,1] * hstats[,2] / batchsize * 100)
  hstats[,1] = factor(hstats[,1], levels=unique(hstats[,1]), ordered=TRUE)
  
  ## Assemble results in list
  return(list(fqstats=c(batchsize=batchsize, nReads=nReads, klength=klength), 
              astats=astats, bstats=bstats, cstats=cstats, dstats=dstats, 
              estats=estats, fstats=fstats, gstats=gstats, hstats=hstats))
}


seeFastqPlot = function(fqlist, arrange=c(1,2,3,4,5,6,7,8), ...){
  require(grid)
  require(ggplot2)
  
  ## Create plotting instances from fqlist
  .fastqPlot = function(x=fqlist) {
    ## (A) Per cycle quality box plot
    astats = x[[1]][["astats"]]
    a = ggplot(astats, aes(x=Cycle, ymin = min, lower = low, middle = mid, upper = top, ymax = max)) + 
      geom_boxplot(stat = "identity", color="#606060", fill="#56B4E9") +
      scale_x_discrete(breaks=c(1, seq(0, length(astats[,1]), by=10)[-1])) + 
      ylab("Quality") + 
      theme(legend.position = "none", plot.title = element_text(size = 12)) +
      ggtitle(names(x))
    
    ## (B) Per cycle base proportion
    bstats = x[[1]][["bstats"]]
    b = ggplot(bstats, aes(Cycle, Frequency, fill=Base), color="black") + 
      scale_x_discrete(breaks=c(1, seq(0, length(unique(bstats$Cycle)), by=10)[-1])) + 
      geom_bar(stat="identity") + 
      theme(legend.position = "none") + 
      ylab("Proportion")  
    
    ## (C) Per cycle average quality of each base type
    cstats = x[[1]][["cstats"]]
    c = ggplot(cstats, aes(Cycle, Quality, group=Base, color=Base)) + 
      geom_line() + ylab("Average Quality") +
      scale_x_discrete(breaks=c(1, seq(0, length(unique(bstats$Cycle)), by=10)[-1])) + 
      theme(legend.position="top", legend.key.size=unit(0.2, "cm")) 
    
    ## (D) Relative K-mer Diversity 
    dstats = x[[1]][["dstats"]]
    d = ggplot(dstats, aes(Cycle, RelDiv, group=Method, color=Method)) + 
      geom_line() + 
      scale_x_discrete(breaks=c(1, seq(0, length(unique(bstats$Cycle)), by=10)[-1])) + 
      ylab(paste(x[[1]][["fqstats"]][["klength"]], "-mer Diversity", sep="")) +
      theme(legend.position = "none")
    
    ## (E) Number of reads where all Phred scores are above a minimum cutoff
    estats = x[[1]][["estats"]]
    e = ggplot(estats, aes(minQuality, Percent, fill = Outliers)) +
      geom_bar(position="dodge", stat="identity") +
      theme(legend.position="top", legend.key.size=unit(0.2, "cm")) + 
      xlab("All Bases Above Min Quality") +
      ylab("% Reads")
    
    ## (F) Distribution of mean quality of reads
    fstats = x[[1]][["fstats"]]
    f = ggplot(fstats, aes(Quality, Percent)) +
      geom_bar(fill="#0072B2", stat="identity") +
      theme(legend.position = "none", plot.title = element_text(size = 9)) +
      ggtitle(paste(formatC(x[[1]][["fqstats"]][["batchsize"]], big.mark = ",", format="f", digits=0), "of", 
                    formatC(x[[1]][["fqstats"]][["nReads"]], big.mark = ",", format="f", digits=0), "Reads")) + 
      scale_x_discrete(breaks=c(0, seq(0, length(unique(fstats$Quality)), by=5)[-1])) +
      xlab("Mean Quality") +
      ylab("% Reads")
    
    ## (G) Read length distribution
    gstats = x[[1]][["gstats"]]
    g = ggplot(gstats, aes(Cycle, Percent)) +
      geom_bar(fill="#0072B2", stat="identity") +
      theme(legend.position = "none") +
      scale_x_discrete(breaks=c(1, seq(0, length(unique(gstats$Cycle)), by=10)[-1])) +
      xlab("Read Length") +
      ylab("% Reads")
    
    ## (H) Read occurrence distribution
    hstats = x[[1]][["hstats"]] 
    myintervals = data.frame(labels=c("1", "2-10", "11-100", "101-1k", "1k-10k", ">10k"), 
                             lower=c(1,2,11,101,1001,10001), 
                             upper=c(2,11,101,1001,10001,Inf))
    iv = sapply(seq(along=myintervals[,1]),
                function(x) sum(hstats[as.numeric(as.vector(hstats$nOccurrences)) >= 
                                myintervals[x,2] & as.numeric(as.vector(hstats$nOccurrences)) 
                                < myintervals[x,3], "Percent"]))
    hstats = data.frame(labels=myintervals[,1], Percent=iv)
    hstats[,1] = factor(hstats[,1], levels=unique(hstats[,1]), ordered=TRUE)
    h = ggplot(hstats, aes(labels, Percent)) +
      geom_bar(fill="#0072B2", stat="identity") +
      theme(legend.position = "none") +
      xlab("Read Occurrence") +
      ylab("% Reads")
    
    ## Assemble results in list
    return(list(a=a, b=b, c=c, d=d, g=g, e=e, f=f, h=h))
  }
  ## Loop to run .fastqPlot and store instances in list 
  fqplot = lapply(names(fqlist), function(z) .fastqPlot(x=fqlist[z]))
  names(fqplot) = names(fqlist)
  
  ## Final plot
  grid.newpage() # Open a new page on grid device
  #	pushViewport(viewport(layout = grid.layout(length(arrange), length(fqplot)))) 
  pushViewport(viewport(layout = grid.layout(length(arrange), n)))
  for(i in seq(along=fqplot)){
    for(j in seq(along=arrange)) {
      suppressWarnings(print(fqplot[[i]][[arrange[j]]], vp = viewport(layout.pos.row = j, layout.pos.col = i)))
    }
  }
}


create_3UTRs_from_pAs = function(pas){
  #` Create 3'UTR GRanges
  required_columns = c("region", "gene_symbol", "pAid", "cds_start", "cds_end")
  if(!all(required_columns %in% names(pas))) {
    stop(cat("The input dataframe must contain the following columns: ", 
             paste(required_columns, collapse=", ")))
  }

  df = subset(pas, region == "3UTR", c("gene_symbol", "pAid", "cds_start", "cds_end"))

  df$chr = sub("(chr.+)[+-]\\d+$", "\\1", df$pAid)
  df$strand = sub("chr.+([+-])\\d+$", "\\1", df$pAid)
  df$pos = as.numeric(sub("chr.+[+-](\\d+)$", "\\1", df$pAid))
  df$start = as.numeric(ifelse(df$strand == "+", df$cds_end + 1, df$pos))
  df$end = as.numeric(ifelse(df$strand == "+", df$pos, df$cds_start - 1))

  df = df[df$end >= df$start, ]
  df[, c("pos", "cds_end", "cds_start")] = NULL
  df = df[complete.cases(df),]
  
  makeGRangesFromDataFrame(df, keep.extra.columns=T)
}


countAllMotif = function(pas, geno = "mm9", search_from = 0, search_len = 50, motif_width = 4){
  # pas must contain "chr", "pA_pos", and "strand" columns
  if(grepl("^mm", geno)){
    require(BSgenome.Mmusculus.UCSC.mm9)
  }else if(grepl("^hg", geno)){
    require(BSgenome.Hsapiens.UCSC.hg19)
  }
  
  require(GenomicRanges)
  require(Biostrings)
  
  upseq.gr = GRanges(seqnames = pas$chr, 
                     ranges = IRanges(start = pas$pA_pos, end = pas$pA_pos),
                     strand = pas$strand)

  start(upseq.gr[strand(upseq.gr) == "+", ]) =
        start(upseq.gr[strand(upseq.gr) == "+", ]) - search_from - search_len + 1
  end(upseq.gr[strand(upseq.gr) == "+", ]) =
        end(upseq.gr[strand(upseq.gr) == "+", ]) - search_from + 1
  end(upseq.gr[strand(upseq.gr) == "-", ]) =
        end(upseq.gr[strand(upseq.gr) == "-", ]) + search_from + search_len - 1
  start(upseq.gr[strand(upseq.gr) == "-", ]) =
        start(upseq.gr[strand(upseq.gr) == "-", ]) + search_from - 1
  
  if(grepl("^mm", geno)){
    upseq = getSeq(Mmusculus, upseq.gr)
  }else if(grepl("^hg", geno)){
    upseq = getSeq(Hsapiens, upseq.gr)
  }
  
  motif_counts = oligonucleotideFrequency(upseq, width = motif_width, step=1,
                                          as.prob=FALSE, as.array=F,
                                          fast.moving.side="right", with.labels=TRUE,
                                          simplify.as="matrix")
  
  motif_counts
}


plot_nucleotide_profile = function(pA, BSgeno = "BSgenome.Hsapiens.UCSC.hg19", window_size = 100){
	#` Plot nucleotide profile near pAs
	require(GenomicRanges)
	require(BSgeno, character.only = T)

	# Get genomic sequences
	seq = getSeq(eval(parse(text=BSgeno)), pA + window_size)

	# Remove sequences with "N"
	seq = seq[!grepl("N", seq)]

	# Computes the consensus matrix of a set of sequences
	m = consensusMatrix(seq)[1:4,]

	# Normalization
	m = scale(m, center = F, scale = colSums(m))

	# Column names indicates position relative to pA
	colnames(m) = -window_size:window_size

	# Write to csv, so that you can plot it in Excel
	# write.csv(m, "nucleotide_profiles.csv")


	# You can also plot using ggplot2:
	require(tidyr)
	require(dplyr)
	require(ggplot2)
	df = as.data.frame(t(m))
	df$Position = as.numeric(rownames(df))
	# png("Nucleotide_Profile.png", 600, 550)
	p = df %>% 
	  gather(key="Base", value = "Fraction", -Position) %>%
	  ggplot(aes(x=Position, y=Fraction, color=Base)) +
	  geom_line() 
	# dev.off()

	p
} 

row.fisher = function(rowdat){
#` Fisher's exact test, accepting a row of 4 integers and returning a p value. 
# Can be used like apply(x, 1, row.fisher), where x is a 4 column matrix
  m = matrix(rowdat, nrow = 2, byrow=T)
  fisher.test(m)$p.value
}

genewise.APA.fisher = function(gene.df){
  #` gene.df is a dataframe with 4 columns containing integers
  # all the cells in the 3rd/4th column contain the sum of the 1st/2nd column
  if(nrow(gene.df) == 1 | any(is.na(gene.df))){
    p_val = rep(1, nrow(gene.df))
  }else{
    gene.df[, 3] = gene.df[, 3] - gene.df[, 1]
    gene.df[, 4] = gene.df[, 4] - gene.df[, 2]
    if(any(gene.df < 0)){
      p_val = rep(1, nrow(gene.df))
    }else{
      p_val = p.adjust(apply(gene.df, 1, row.fisher), "BH")
    }
  }
  p_val
}

APA.fisher = function(dat){
  #` dat must contain "gene_symbol", two "_count" an two corresponding "_gene" columns
  # the "_gene" column contains total read counts for each gene
  # split the data 
  dat.lst = split(dat, dat$gene_symbol)
  # work on one gene each time
  p.vals = sapply(dat.lst, function(x) genewise.APA.fisher(x[, -1])) 
  # recombind the data with results
  p.vals = unlist(p.vals)
}

row.chi = function(rowdat){
  require(MASS)
  m = matrix(rowdat, nrow = 2, byrow=T)
  if(any(is.na(m)) | any(m < 0)){
    p_val = 1
  }else{
    chisq.test(m)$p.value
  }
}

genewise.chi = function(dat){
  #` dat is a dataframe with 4 columns containing integers
  # all the cells in the 3rd/4th column contain the sum of the 1st/2nd column
  dat[, 3] = dat[, 3] - dat[, 1]
  dat[, 4] = dat[, 4] - dat[, 2]
  
  p_val = p.adjust(apply(dat, 1, row.chi), "BH")
}


#### for a gene with n pA isoforms, return n-1 combination of neighboring pA sites
gene2neighorpAPairs = function(gene){
  if(nrow(gene) > 1){
    # set up the columns
  gene_id = rep(gene$gene_id[1],nrow(gene)-1) 
  gene_symbol = rep(gene$gene_symbol[1], nrow(gene)-1) 
  comparedPAi = rep(NA, length(gene_id)) 
  comparedPAj = rep(NA, length(gene_id)) 
  
  # Sort the pAs according to their position. 
  # This ensure that the pAi is always the Prx pA compared to pAj
  if (gene$strand[1] == "-"){
    gene = gene[order(-1*gene$pA_pos), ]
  }else if (gene$strand[1] == "+"){
    gene = gene[order(gene$pA_pos), ]
  }
  
  # create neighboring pA pairs
  k = 1
  for (i in 1:(nrow(gene) -1)){
    comparedPAi[k] = gene$pAid[i]
    comparedPAj[k] = gene$pAid[i+1]
    k = k + 1
  }
  cbind(gene_id, gene_symbol, comparedPAi, comparedPAj)
  }
  
}

#### For a gene with n pA isoforms, return n*(n-1)/2 combination of pA sites
gene2pAPairs = function(gene){ # pAid is needed!
  if(nrow(gene) > 1){
    # set up the columns
  gene_id = rep(gene$gene_id[1], choose(nrow(gene), 2)) 
  gene_symbol = rep(gene$gene_symbol[1], choose(nrow(gene), 2)) 
  comparedPAi = rep(NA, length(gene_id)) 
  comparedPAj = rep(NA, length(gene_id)) 
  
  # Sort the pAs according to their position. 
  # This ensure that the pAi is always the Prx pA compared to pAj
  if(gene$strand[1] == "-"){
    gene = gene[order(-1*gene$pA_pos), ]
  }else if(gene$strand[1] == "+"){
    gene = gene[order(gene$pA_pos), ]
  }
  
  # create pairwise combinations of pAs
  k = 1
  for (i in 1:(nrow(gene) -1)){
    for (j in (i+1):nrow(gene)){ # compare all combinations of the APA isoforms
      comparedPAi[k] = gene$pAid[i]
      comparedPAj[k] = gene$pAid[j]
      k = k + 1
    }
  }
  cbind(gene_id, gene_symbol, comparedPAi, comparedPAj)
  }
  
}

#### pAs from the same gene will get paired, but no further calculations will be done
pairwise.pAs = function(data, neighbor = F, toptwo = F, extra_cols = NULL, match_only=F){
  # pAid is needed!
  # extra_cols are column names other than those defined below that you want to compare between the isoforms
  if(is.null(data$gene_id)){
    data$gene_id = data$gene_symbol
  }
  # cols are column names used later for merging data frames
  if(any(grep("^RAI", names(data)))){
    cols = grep("^RAI", names(data), value = T)
  }
  if(any(grep("APA", names(data)))){
    cols = grep("APA", names(data), value = T)
  }
  if(any(grep("^RPM", names(data)))){
    cols = grep("^RPM", names(data), value = T)
  }
  if(any(grep("^num", names(data)))){
    cols = grep("^num", names(data), value = T)
  }
  if(any(grep("_count$", names(data)))){
    cols = grep("_count$", names(data), value = T)
  }
  cols = c(cols, extra_cols)
  
  # split the data frame according to gene ids
  gene.lst = split(data, data$gene_id)
  # only keep the genes with at least two APA isoforms
  gene.lst = gene.lst[sapply(gene.lst, nrow) >= 2]
  # compare read numbers between combinations of the isoforms in the two fractions
  if(neighbor == F & toptwo == F){ 
    cmp.lst = lapply(gene.lst, function(gene) gene2pAPairs(gene))
  }else if(neighbor == F & toptwo == T){
    cmp.lst = lapply(gene.lst, function(gene) gene2pAPairs(gene)) ## same as toptwo == F !!!!!!!!
  }else if (neighbor == T & toptwo == F){
    cmp.lst = lapply(gene.lst, function(gene) gene2neighorpAPairs(gene))
  }else{
    stop()
  }
  
  # only keep non-empty lists
  #cmp.lst = cmp.lst[sapply(cmp.lst, length) > 0] 
  
  # convert the list of df into one df
  cmp = do.call("rbind", cmp.lst) #### matrix
  cmp = as.data.frame(cmp, row.names = NULL, stringsAsFactors = FALSE)#"row.names = NULL" avoids duplicated row names; Without "stringsAsFactors = F", everything will be factors  
  rownames(cmp) = NULL # I thought the row names are already deleted
  
  # get strand and pAid information from the original data
  cmp = merge(cmp, data[,c("pAid", cols)], by.x = "comparedPAj", by.y = "pAid", sort = F)
  cmp = merge(cmp, data[,c("pAid", cols)], by.x = "comparedPAi", by.y = "pAid", sort = F)
  
  names(cmp) = sub("\\.x$", "_Dis", names(cmp))
  names(cmp) = sub("\\.y$", "_Prx", names(cmp))
  names(cmp) = sub("comparedPAi", "Prx_pA", names(cmp))
  names(cmp) = sub("comparedPAj", "Dis_pA", names(cmp))
  
  if(match_only == F){
    if(any(grep("^RAI", names(data)))){
      # calculate delta RAIs
      for(string in sub("_Dis", "", grep("(RAI.*)_Dis", names(cmp), value = T, perl = T))){
        Dis.col.name = paste0(string, "_Dis")
        Prx.col.name = paste0(string, "_Prx")
        cmp[, paste0("delt_", string)] = cmp[, Dis.col.name] - cmp[, Prx.col.name]
      }
    }
    
    if(any(grep("^RPM", names(data)))){
      # calculate ratio of Dis to Prx pA RPMs
      for(string in sub("_Dis", "", grep("(RPM.*)_Dis", names(cmp), value = T, perl = T))){
        
        Dis.col.name = paste0(string, "_Dis")
        Prx.col.name = paste0(string, "_Prx")
        cmp[, paste0("ratio_", sub("^RPM_", "", string))] = log2(cmp[, Dis.col.name]) - log2(cmp[, Prx.col.name])
      }
    }
    
    if(any(grep("^num", names(data)))){
      # calculate ratio of Dis to Prx pA read numbers
      for(string in sub("_Dis", "", grep("(num.*)_Dis", names(cmp), value = T, perl = T))){
        
        Dis.col.name = paste0(string, "_Dis")
        Prx.col.name = paste0(string, "_Prx")
        cmp[, paste0("ratio_", sub("^num_", "", string))] = log2(cmp[, Dis.col.name]) - log2(cmp[, Prx.col.name])
      }
      #cmp = cmp[, -grep("^num", names(cmp))]
      #names(cmp) = sub("ratio_", "", names(cmp))
    }
    
    if(any(grep("count$", names(data)))){
      # calculate ratio of Dis to Prx pA read numbers
      for(string in sub("_Dis", "", grep("(count.*)_Dis", names(cmp), value = T, perl = T))){
        
        Dis.col.name = paste0(string, "_Dis")
        Prx.col.name = paste0(string, "_Prx")
        cmp[, sub("counts?", "d2p_ratio", string)] = log2(cmp[, Dis.col.name]) - log2(cmp[, Prx.col.name])
      }
    }
    
    if(any(grep("^HL\\.", names(data)))){
      # calculate ratio of Dis to Prx pA halflifes
      for(string in sub("_Dis", "", grep("(^HL\\..*)_Dis", names(cmp), value = T, perl = T))){
        
        Dis.col.name = paste0(string, "_Dis")
        Prx.col.name = paste0(string, "_Prx")
        cmp[, paste0(string, "_ratio")] = log2(cmp[, Dis.col.name]) - log2(cmp[, Prx.col.name])
        cmp[, paste0(string, "_difference")] = cmp[, Dis.col.name] - cmp[, Prx.col.name]
      }
    }
  }
  cmp
}

get.pA.pairs = function(data, cols = "UTR3_size", neighbor = F, toptwo = F, match_only=F){ 
  # gene_id (or gene_symbol) and pAid are needed in input data!
  # when toptwo is set to T, rpm values for each sample should also be provided to select the
  # top two highly expressed isoforms
  # cols are columns that you want to compare between the isoforms
  names(data) =  sub("position", "pA_pos", names(data))
  if(is.null(data$gene_id)){
    data$gene_id = data$gene_symbol
  }
  
  if(any(grep("_count$", names(data)))){
    cols = c(grep("_count$", names(data), value = T), cols)
  }
  
  #### for a gene with n pA isoforms, return 1 pair pA isoforms with highest read counts
  if(toptwo){
    if(any(grepl("rpm$", names(data)))){
      # focus on rpm
      tmp = data[, grep("pAid|gene_id|rpm$", names(data))]
      # calculate total rpm for each pA isoform
      tmp$total_rpm = rowSums(data[, grep("rpm$", names(data))])
      # split the data frame according to gene ids
      gene.lst = split(tmp, tmp$gene_id)
      # only keep the genes with at least two APA isoforms
      gene.lst = gene.lst[sapply(gene.lst, nrow) >= 2]
      # reorder each dataframe that have > 2 pA isoforms
      for(i in 1:length(gene.lst)){
        if(nrow(gene.lst[[i]]) > 2){
          gene.lst[[i]] = gene.lst[[i]][order(gene.lst[[i]]$total_rpm, decreasing = T)[1:2],]
        }
      }
      # recombine the dataframes, and only keep the pAids
      tmp = do.call("rbind", gene.lst)$pAid
      # only keep the pA isoforms in tmp
      data = subset(data, pAid %in% tmp)
      rm(tmp)
      # split the data frame again according to gene ids
      gene.lst = split(data, data$gene_id)
      # compare read numbers between combinations of the isoforms in the two fractions
      cmp.lst = lapply(gene.lst, function(gene) gene2pAPairs(gene))
    }else{
      stop("To select top two isoforms, RPM columns should exist in the input.")
    } 
  }else{
    # split the data frame according to gene ids
    gene.lst = split(data, data$gene_id)
    # only keep the genes with at least two APA isoforms
    gene.lst = gene.lst[sapply(gene.lst, nrow) >= 2]
    # compare read numbers between combinations of the isoforms in the two fractions
    if(neighbor == F){ 
      cmp.lst = lapply(gene.lst, function(gene) gene2pAPairs(gene))
    }else if(neighbor == T){
      cmp.lst = lapply(gene.lst, function(gene) gene2neighorpAPairs(gene))
    }
  }
  
  # convert the list of df into one df
  cmp = do.call("rbind", cmp.lst) #### matrix
  cmp = as.data.frame(cmp, row.names = NULL, stringsAsFactors = F)#"row.names = NULL" avoids duplicated row names; Without "stringsAsFactors = F", everything will be factors  
  rownames(cmp) = NULL # I thought the row names are already deleted
  
  # get strand and pAid information from the original data
  cmp = merge(cmp, subset(data, select = c("pAid", cols)), by.x = "comparedPAj", by.y = "pAid", sort = F)
  cmp = merge(cmp, subset(data, select = c("pAid", cols)), by.x = "comparedPAi", by.y = "pAid", sort = F)
  
  # change column names
  names(cmp) = sub("\\.x$", "_Dis", names(cmp))
  names(cmp) = sub("\\.y$", "_Prx", names(cmp))
  names(cmp) = sub("comparedPAi", "Prx_pA", names(cmp))
  names(cmp) = sub("comparedPAj", "Dis_pA", names(cmp))
  
  # calculate Dis to Prx ratios within each sample
  if(match_only == F){
    if(any(grep("^RAI", names(data)))){
      # calculate delta RAIs
      for(string in sub("_Dis", "", grep("(RAI.*)_Dis", names(cmp), value = T, perl = T))){
        Dis.col.name = paste0(string, "_Dis")
        Prx.col.name = paste0(string, "_Prx")
        cmp[, paste0("delt_", string)] = cmp[, Dis.col.name] - cmp[, Prx.col.name]
      }
    }
    
    if(any(grep("^RPM", names(data)))){
      # calculate ratio of Dis to Prx pA RPMs
      for(string in sub("_Dis", "", grep("(RPM.*)_Dis", names(cmp), value = T, perl = T))){
        
        Dis.col.name = paste0(string, "_Dis")
        Prx.col.name = paste0(string, "_Prx")
        cmp[, paste0("ratio_", sub("^RPM_", "", string))] = log2(cmp[, Dis.col.name]) - log2(cmp[, Prx.col.name])
      }
    }
    
    if(any(grep("^num", names(data)))){
      # calculate ratio of Dis to Prx pA read numbers
      for(string in sub("_Dis", "", grep("(num.*)_Dis", names(cmp), value = T, perl = T))){
        
        Dis.col.name = paste0(string, "_Dis")
        Prx.col.name = paste0(string, "_Prx")
        cmp[, paste0("ratio_", sub("^num_", "", string))] = log2(cmp[, Dis.col.name]) - log2(cmp[, Prx.col.name])
      }
      #cmp = cmp[, -grep("^num", names(cmp))]
      #names(cmp) = sub("ratio_", "", names(cmp))
    }
    
    if(any(grep("count$", names(data)))){
      # calculate ratio of Dis to Prx pA read numbers
      for(string in sub("_Dis", "", grep("(count.*)_Dis", names(cmp), value = T, perl = T))){
        
        Dis.col.name = paste0(string, "_Dis")
        Prx.col.name = paste0(string, "_Prx")
        cmp[, sub("counts?", "d2p_ratio", string)] = log2(cmp[, Dis.col.name]) - log2(cmp[, Prx.col.name])
      }
    }
    
    if(any(grep("^HL\\.", names(data)))){
      # calculate ratio of Dis to Prx pA halflifes
      for(string in sub("_Dis", "", grep("(^HL\\..*)_Dis", names(cmp), value = T, perl = T))){
        
        Dis.col.name = paste0(string, "_Dis")
        Prx.col.name = paste0(string, "_Prx")
        cmp[, paste0(string, "_ratio")] = log2(cmp[, Dis.col.name]) - log2(cmp[, Prx.col.name])
        cmp[, paste0(string, "_difference")] = cmp[, Dis.col.name] - cmp[, Prx.col.name]
      }
    }
  }
  
  # return
  cmp = merge(cmp, unique(data[, c("gene_symbol", "description")]), all.x = T, sort=F)
}

get_red = function(data, cols = "UTR3_size", neighbor = F, toptwo = F, match_only=F){ 
  # gene_id (or gene_symbol) and pAid are needed in input data!
  # when toptwo is set to T, rpm values for each sample should also be provided to select the
  # top two highly expressed isoforms
  # cols are columns that you want to compare between the isoforms
  names(data) =  sub("position", "pA_pos", names(data))
  if(is.null(data$gene_id)){
    data$gene_id = data$gene_symbol
  }
  
  if(any(grep("_count$", names(data)))){
    cols = c(grep("_count$", names(data), value = T), cols)
  }
  
  # For a gene with n pA isoforms, return 1 pair pA isoforms with highest read counts
  if(toptwo){
    if(any(grepl("rpm$", names(data)))){
      # focus on rpm
      tmp = data[, grep("pAid|gene_id|rpm$", names(data))]
      # calculate total rpm for each pA isoform
      tmp$total_rpm = rowSums(data[, grep("rpm$", names(data))])
      # split the data frame according to gene ids
      gene.lst = split(tmp, tmp$gene_id)
      # only keep the genes with at least two APA isoforms
      gene.lst = gene.lst[sapply(gene.lst, nrow) >= 2]
      # reorder each dataframe that have > 2 pA isoforms
      for(i in 1:length(gene.lst)){
        if(nrow(gene.lst[[i]]) > 2){
          gene.lst[[i]] = gene.lst[[i]][order(gene.lst[[i]]$total_rpm, decreasing = T)[1:2],]
        }
      }
      # recombine the dataframes, and only keep the pAids
      tmp = do.call("rbind", gene.lst)$pAid
      # only keep the pA isoforms in tmp
      data = subset(data, pAid %in% tmp)
      rm(tmp)
      # split the data frame again according to gene ids
      gene.lst = split(data, data$gene_id)
      # compare read numbers between combinations of the isoforms in the two fractions
      cmp.lst = lapply(gene.lst, function(gene) gene2pAPairs(gene))
    }else{
      stop("To select top two isoforms, RPM columns should exist in the input.")
    } 
  }else{
    # split the data frame according to gene ids
    gene.lst = split(data, data$gene_id)
    # only keep the genes with at least two APA isoforms
    gene.lst = gene.lst[sapply(gene.lst, nrow) >= 2]
    # compare read numbers between combinations of the isoforms in the two fractions
    if(neighbor == F){ 
      cmp.lst = lapply(gene.lst, function(gene) gene2pAPairs(gene))
    }else if(neighbor == T){
      cmp.lst = lapply(gene.lst, function(gene) gene2neighorpAPairs(gene))
    }
  }
  
  # convert the list of df into one df
  cmp = do.call("rbind", cmp.lst) #### matrix
  cmp = as.data.frame(cmp, row.names = NULL, stringsAsFactors = F)
  rownames(cmp) = NULL # I thought the row names are already deleted
  
  # get strand and pAid information from the original data
  cmp = merge(cmp, subset(data, select = c("pAid", cols)), by.x = "comparedPAj", by.y = "pAid", sort = F)
  cmp = merge(cmp, subset(data, select = c("pAid", cols)), by.x = "comparedPAi", by.y = "pAid", sort = F)
  
  # change column names
  names(cmp) = sub("\\.x$", "_Dis", names(cmp))
  names(cmp) = sub("\\.y$", "_Prx", names(cmp))
  names(cmp) = sub("comparedPAi", "Prx_pA", names(cmp))
  names(cmp) = sub("comparedPAj", "Dis_pA", names(cmp))
  
  # calculate Dis to Prx ratios within each sample
  if(match_only == F){
    if(any(grep("^RPM", names(data)))){
      # calculate ratio of Dis to Prx pA RPMs
      for(string in sub("_Dis", "", grep("(RPM.*)_Dis", names(cmp), value = T, perl = T))){
        Dis.col.name = paste0(string, "_Dis")
        Prx.col.name = paste0(string, "_Prx")
        cmp[, paste0("ratio_", sub("^RPM_", "", string))] = log2(cmp[, Dis.col.name]) - log2(cmp[, Prx.col.name])
      }
    }
    if(any(grep("count$", names(data)))){
      # calculate ratio of Dis to Prx pA read numbers
      for(string in sub("_Dis", "", grep("(count.*)_Dis", names(cmp), value = T, perl = T))){
        Dis.col.name = paste0(string, "_Dis")
        Prx.col.name = paste0(string, "_Prx")
        cmp[, sub("counts?", "d2p_ratio", string)] = log2(cmp[, Dis.col.name]) - log2(cmp[, Prx.col.name])
      }
    }
  }
  
  # return
  cmp = merge(cmp, unique(data[, c("gene_symbol", "description")]), all.x = T, sort=F)
}


#### Convert pAid into a dataframe containing chr, pA_pos, and strand info
pAid2pos = function(pairwise.pas){
  pAid.i = pairwise.pas$Prx_pA
  pAid.j = pairwise.pas$Dis_pA
  pA.pos.df = as.data.frame(do.call(rbind, strsplit(pAid.i, "[+-]")), 
                            stringsAsFactors = F)
  #pos.i.string = "Prx_pA_pos"
  names(pA.pos.df) = c("chr", "Prx_pA_pos")
  pA.pos.df$chr = as.character(pA.pos.df$chr)
  pA.pos.df$Prx_pA_pos = as.integer(pA.pos.df$Prx_pA_pos)
  pA.pos.df$strand = as.vector(sub(".*?([+-])\\d*", "\\1", pAid.i, perl = T))
  
  pA.pos.df$Dis_pA_pos = as.data.frame(do.call(rbind, strsplit(pAid.j, "[+-]")), 
                                          stringsAsFactors = F)[,2]
  pA.pos.df$Dis_pA_pos = as.integer(pA.pos.df$Dis_pA_pos)
  
  cbind(pA.pos.df, pairwise.pas)
}

#### If the input is just pas (with only one pAid column)
split.pAid.in.pas = function(pas){
  pA.pos.df = as.data.frame(do.call(rbind, strsplit(pas$pAid, "[+-]")), 
                            stringsAsFactors = F)
  names(pA.pos.df) = c("chr", "pA_pos")
  pA.pos.df$chr = as.character(pA.pos.df$chr)
  pA.pos.df$pA_pos = as.integer(pA.pos.df$pA_pos)
  pA.pos.df$strand = as.vector(sub(".*?([+-])\\d*", "\\1", pas$pAid, perl = T))
  cbind(pA.pos.df, pas)
}

#### Match aUTRs with slightly different pA sites in two data frames allowing distance of "radius" nt
match.pA.pair = function(df1, df2, radius = 12, keep = "intersect", rm_dup_col = T){
  # "keep" can be one of "union", "intersect", "df1", "df2")
  # "rm_dup_col": should duplicated columns be removed?
  
  if(!keep %in% c("union", "intersect", "df1", "df2")){
    stop('"keep" can only be one of "union", "intersect", "df1", and "df2"')
  }
  
  # calculate chr, strand, and position
  if("Prx_pA" %in% names(df1) & "Dis_pA" %in% names(df1)){
    df1 = pAid2pos(df1)
  }
  if("Prx_pA" %in% names(df2) & "Dis_pA" %in% names(df2)){
    df2 = pAid2pos(df2)
  }
  
  # calculate pA position in case it is not provided (for backward compatibility, when pAid is not supplied)
  if(!"Prx_pA_pos" %in% names(df1)){
    df1$Prx_pA_pos = as.numeric(sub("chr.+[-+]", "", df1$Prx_pA))
  }
  if(!"Dis_pA_pos" %in% names(df1)){
    df1$Dis_pA_pos = as.numeric(sub("chr.+[-+]", "", df1$Dis_pA))
  }
  if(!"Prx_pA_pos" %in% names(df2)){
    df2$Prx_pA_pos = as.numeric(sub("chr.+[-+]", "", df2$Prx_pA))
  }
  if(!"Dis_pA_pos" %in% names(df2)){
    df2$Dis_pA_pos = as.numeric(sub("chr.+[-+]", "", df2$Dis_pA))
  }
  
  # split genes
  if("gene_symbol" %in% names(df1) & "gene_symbol" %in% names(df2)){
    # only keep common genes
    df1_intersect = subset(df1, toupper(gene_symbol) %in% toupper(df2$gene_symbol))
    df2_intersect = subset(df2, toupper(gene_symbol) %in% toupper(df1$gene_symbol))
    # order the dataframe so that the list will be in the same order
    df1_intersect = df1_intersect[order(df1_intersect$gene_symbol),]
    df2_intersect = df2_intersect[order(df2_intersect$gene_symbol),]
    # split into individual genes for comparison
    df1_intersect.lst = split(df1_intersect, df1_intersect$gene_symbol)
    df2_intersect.lst = split(df2_intersect, df2_intersect$gene_symbol)
  }else if("gene_id" %in% names(df1) & "gene_id" %in% names(df2)){
    # only keep common genes
    df1$gene_id = as.numeric(df1$gene_id)
    df2$gene_id = as.numeric(df2$gene_id)
    df1_intersect = subset(df1, gene_id %in% df2$gene_id)
    df2_intersect = subset(df2, gene_id %in% df1$gene_id)
    # order the dataframe so that the list will be in the same order
    df1_intersect = df1_intersect[order(as.numeric(df1_intersect$gene_id)),]
    df2_intersect = df2_intersect[order(as.numeric(df2_intersect$gene_id)),]
    # split into individual genes for comparison
    df1_intersect.lst = split(df1_intersect, df1_intersect$gene_id)
    df2_intersect.lst = split(df2_intersect, df2_intersect$gene_id)
  }
  
  # determine if the pA pairs are the same, gene by gene
  df = data.frame()
  unmatched_genes = vector()
  for(i in 1:length(df1_intersect.lst)){
    unmatched = T
    for(j in 1:nrow(df1_intersect.lst[[i]])){
      for(k in 1:nrow(df2_intersect.lst[[i]])){
        if(abs(df1_intersect.lst[[i]][j,"Prx_pA_pos"] - df2_intersect.lst[[i]][k, "Prx_pA_pos"]) <= radius*2 &
           abs(df1_intersect.lst[[i]][j,"Dis_pA_pos"] - df2_intersect.lst[[i]][k,"Dis_pA_pos"]) <= radius*2){
          row = cbind(df1_intersect.lst[[i]][j,], df2_intersect.lst[[i]][k,])
          df = rbind(df, row)
          
          unmatched = F
          # next
        }
      }
    }
    if(unmatched){
      unmatched_genes = c(unmatched_genes, df1_intersect[i, "gene_symbol"])
    }
  }
  
  if(keep == "union"){
    df1_unique = subset(df1, !toupper(gene_symbol) %in% toupper(df2$gene_symbol))
    if(nrow(df1_unique) > 0){
      df2_empty = matrix(NA, ncol = ncol(df2), nrow = nrow(df1_unique))
      colnames(df2_empty) = names(df2)
      if(length(unmatched_genes)){
        df2_empty = rbind(df2_empty, subset(df2, gene_symbol %in% unmatched_genes))
      }
      colnames(df2_empty) = names(df)[-(1:ncol(df1))]
      df = rbind(df, cbind(df1_unique, df2_empty))
    }
    df2_unique = subset(df2, !toupper(gene_symbol) %in% toupper(df1$gene_symbol))
    if(nrow(df2_unique) > 0){
      df1_empty = matrix(NA, ncol = ncol(df1), nrow = nrow(df2_unique))
      colnames(df1_empty) = names(df1)
      if(length(unmatched_genes)){
        df1_empty = rbind(df1_empty, subset(df1, gene_symbol %in% unmatched_genes))
      }
      colnames(df1_empty) = names(df)[1:ncol(df1)]
      df = rbind(df, cbind(df1_empty, df2_unique))
    }
  }else if(keep == "df1"){
    df1_unique = subset(df1, !toupper(gene_symbol) %in% toupper(df2$gene_symbol))
    if(nrow(df1_unique) > 0){
      df2_empty = matrix(NA, ncol = ncol(df2), nrow = nrow(df1_unique))
      colnames(df2_empty) = names(df2)
      if(length(unmatched_genes)){
        df2_empty = rbind(df2_empty, subset(df2, gene_symbol %in% unmatched_genes))
      }
      colnames(df2_empty) = names(df)[-(1:ncol(df1))]
      df = rbind(df, cbind(df1_unique, df2_empty))
    }
  }else if(keep == "df2"){
    df2_unique = subset(df2, !toupper(gene_symbol) %in% toupper(df1$gene_symbol))
    if(nrow(df2_unique) > 0){
      df1_empty = matrix(NA, ncol = ncol(df1), nrow = nrow(df2_unique))
      colnames(df1_empty) = names(df1)
      if(length(unmatched_genes)){
        df1_empty = rbind(df1_empty, subset(df1, gene_symbol %in% unmatched_genes))
      }
      colnames(df1_empty) = names(df)[1:ncol(df1)]
      df = rbind(df, cbind(df1_empty, df2_unique))
    }
  }
  
  # Remove duplicated columns
  duplicated_columns = intersect(names(df2), names(df1))
  names(df) = make.unique(names(df))
  if(rm_dup_col){
    for(duplicated_column in duplicated_columns){
      df[, duplicated_column] = ifelse(is.na(df[, duplicated_column]), df[, paste0(duplicated_column, ".1")], df[, duplicated_column])
      df[, paste0(duplicated_column, ".1")] = NULL
    }
  }
  df
}
match.pA.pairs = match.pA.pair


#### Match pA sites in two data frames, allowing distance of "radius" nt
match.pA = function(df1, df2, radius = 12, left_join=T, copy_right_columns = T){
  # df1 and df2 are data frames containing the pAid (chr1-196924906) column.
  # pA sites whose genomic positions <= 2*radius bp will be considered as the same pA site
  # all rows and columns of df1 will be preserved
  df1$chr = NULL
  df1$strand = NULL
  df1$pA_pos = NULL
  df1 = split.pAid.in.pas(as.data.frame(df1))
  df2$chr = NULL
  df2$strand = NULL
  df2$pA_pos = NULL
  df2 = split.pAid.in.pas(as.data.frame(df2))
  
  require(GenomicRanges)
  df1.gr = GRanges(seqnames = df1$chr,
                   strand = df1$strand, 
                   ranges = IRanges(start = df1$pA_pos - radius, end = df1$pA_pos + radius),
                   df1.pAid = df1$pAid)
  
  df2.gr = GRanges(seqnames = df2$chr,
                   strand = df2$strand, 
                   ranges = IRanges(start = df2$pA_pos - radius, end = df2$pA_pos + radius),
                   df2.pAid = df2$pAid)
  if(left_join){
    df2.gr = df2.gr[df2.gr %over% df1.gr]
  }
  # make a data frame for matching rows
  olp = findOverlaps(df1.gr, df2.gr)
  # some pAs in one gr can match >= 1 pAs in another gr!!!!!!!
  #   sum(duplicated(queryHits(olp)))
  #   sum(duplicated(subjectHits(olp)))
  
  matched.pA = cbind(df1.pAid = df1.gr$df1.pAid[queryHits(olp)], 
                     df2.pAid = df2.gr$df2.pAid[subjectHits(olp)])
  
  data.cmn = merge(matched.pA, df1, by.y = "pAid", by.x = "df1.pAid",  
                   sort = F, stringAsFactors = F, all=F) ## all = F is important to remove unmatched pAids
  data.cmn = merge(data.cmn, df2, by.x = "df2.pAid", by.y = "pAid", 
                   sort = F, stringAsFactors = F, all=F) ## all = F is important to remove unmatched pAids
  
  
  # clean up the dataframe
  # change colnames so that they can be processed like other common columns
  names(data.cmn)[names(data.cmn) == "df1.pAid"] = "pAid.x"
  names(data.cmn)[names(data.cmn) == "pAid"] = "pAid.y"
  #data.cmn = apply(data.cmn, 2, as.vector) # the "pAd.x" column is factor
  data.cmn[, "pAid.x"] = as.vector(data.cmn[, "pAid.x"])
  # copy info
  if(copy_right_columns == T){
    for (col in colnames(df1)){
      if (col %in% colnames(df2)){
        colx = paste0(col, ".x")
        coly = paste0(col, ".y")
        #data.cmn[, colx][is.na(data.cmn[, colx])] = data.cmn[, coly][is.na(data.cmn[, colx])]
        #print(data.cmn[is.na(data.cmn[, colx]), colx])
        #print(as.vector(data.cmn[is.na(data.cmn[, colx]), coly]))
        data.cmn[is.na(data.cmn[, colx]), colx] = as.vector(data.cmn[is.na(data.cmn[, colx]), coly])
      }
    }
  }
  
  ### Some pAs in df1 match to >=1 pAs in df2. Only keep the one closest to the pA in df1
  data.dup = data.cmn[data.cmn$pAid %in% matched.pA[,1][duplicated(matched.pA[,1])], ]
  data.dup$pA_distance = abs(data.dup$pA_pos.x - data.dup$pA_pos.y)
  data.dup = split(data.dup, data.dup$pAid.x)
  data.dup = lapply(data.dup, function(gene) gene[which.min(gene$pA_distance), ])
  data.dup = as.data.frame(do.call(rbind, data.dup))
  data.dup$pA_distance = NULL
  data.cmn = rbind(data.cmn[!data.cmn$pAid %in% matched.pA[,1][duplicated(matched.pA[,1])], ], data.dup)
  
  # remove .y columns 
  data.cmn = data.cmn[, -grep("\\.y", names(data.cmn))]
  names(data.cmn) = sub("\\.x", "", names(data.cmn))
  
  ### combine data for pAs in df1 that match to >=1 pAs in df2 ??
  ### Not necessary for now
  
  # return
  data.cmn
}
match.pAs = match.pA
match.pAs.2 = match.pA  




