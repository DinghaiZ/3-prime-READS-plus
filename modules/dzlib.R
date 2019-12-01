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


