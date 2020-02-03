library(MASS) 
library(ggrepel) 
require(tidyr)
require(dplyr)
require(ggplot2)
theme_set(theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5))) 

rename_grl = function(grl, old_key = "ACCNUM", new_key = "SYMBOL"){
    # A function to rename GRanges in a GRangeList.
    require(AnnotationDbi)
    require(GenomicRanges)
    
    # Convert a GRangeList to a GRanges for faster calculation
    gr = unlist(grl)

    # Convert the accession numbers of the genes into gene symbols
    names(gr) = mapIds(org.db, keys = names(gr), keytype = old_key, column = new_key)
    
    # Rebuild GRanges by gene symbol 
    grl = split(gr, names(gr))   

    # Collapse features
    grl = reduce(grl)
    
    # Remove GRanges that are on two strands
    grl[sapply(grl, function(x) length(strand(x)@values) == 1)]
}

findPotentialStartsAndStops = function(sequence){
  require(Biostrings)
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
          if (typei == "ATG" && (typej == "TAA" || typej == "TAG" || typej == "TGA")
                        && posdiffmod3 == 0){
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
              # Including the stop codon.
              orfstops = append(orfstops, posj+2, after=length(orfstops)) 
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
  for (i in 1:numorfs){
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
  cstats = data.frame(Quality=c(A, C, G, T), Base=rep(c("A", "C", "G", "T"), 
                      each=length(A)), Cycle=c(names(A), names(C), names(G), names(T)))
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
  edf = sapply(ev, function(x) sapply(as.numeric(names(ev)), 
                   function(y) sum(rowSums(q >= x, na.rm=TRUE) >= (rowSums(!is.na(q))-y))))
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

get_protein_localization = function(df = pas, localization_file){
  #` Add protein localization info to the dataframe df
    
  # df can contain poly(A) sites or pairs of poly(A) sites
  # The gene_symbol column in df is required for joining the dataframes
    
  require(dplyr)
  
  localizations = read.csv(localization_file, as.is = T) %>%
    distinct()

  m = localizations[, -grep("gene_symbol|localizations", names(localizations))]
    
  localizations = localizations %>%
    mutate(Other = as.integer(rowSums(m) == 0)) %>%
    dplyr::select(-localizations)
  
  names(localizations) = paste0("YSU.", names(localizations)) 

  names(localizations) = sub("YSU.gene_symbol", "gene_symbol", names(localizations)) 

  localizations = subset(localizations, !is.na(gene_symbol))
  localizations = subset(localizations, gene_symbol != "")
  localizations = subset(localizations, gene_symbol != "Unknown")
  
  # Genes not found in localizations will be NA
  df$gene_symbol = toupper(df$gene_symbol)
  localizations$gene_symbol = toupper(localizations$gene_symbol)
  df = left_join(df, localizations, by="gene_symbol") 

  df$gene_symbol = Hmisc::capitalize(tolower(df$gene_symbol))
  
  # pAs in introns, UA, etc should not be annotated with protein localization
  if("region" %in% names(df)){ 
    for(loc in grep("gene_symbol", names(localizations), value=T, invert=T)){
      df[df$region != "3UTR", loc] = NA
    }
  }
  df
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

get_hg19_Alu = function(df=pas, txdb_path,
                        repeat.table=file.path(SHARED_DATA_DIR, "hg19_rmsk"), 
                        ignore.intron.Alu = T, 
                        keep.alu.direction = F){
  #` A function to annotate Alu elements in hg19
    
  rmsk = read.table(repeat.table, header = T, sep = "\t", as.is = T, comment.char = "")
  alu = subset(rmsk, repFamily == "Alu", select = names(rmsk)[c(2:11, 14:16)])
  # If "strand" is "+": 
    # The sequence between "genoStart" and "genoEnd" on the positive strand 
    # is similar to the "repName"
  # If "strand" is "-": 
    # The sequence between "genoStart" and "genoEnd" on the minus strand 
    # is similar to the "repName"
  
  # If a 3'UTR is on the "+" strand and it contains an Alu element 
    # defined by ("genoName", "genoStart", "genoEnd", "strand"):
  #    if the Alu "strand" is "+": the Alu in 3'UTR is "S" (sense)
  #    if the Alu "strand" is "-": the Alu in 3'UTR is "AS" (anti-sense)
  # If a 3'UTR is on the "-" strand and it contains an Alu element
    # defined by ("genoName", "genoStart", "genoEnd", "strand"):
  #    if the Alu "strand" is "+": the Alu in 3'UTR is "AS" (anti-sense)
  #    if the Alu "strand" is "-": the Alu in 3'UTR is "S" (sense)
  
  alu$alu.id = 1:nrow(alu)
  
  # Map Alu to 3'UTRs
  require(GenomicRanges)
  pA.gr = create_3UTRs_from_pAs(df)
    
  # Sense
  alu.s = alu[, c("genoName", "genoStart", "genoEnd", "strand", "repName", 
                  "repStart", "repEnd", "alu.id")] 
  names(alu.s) = c("chr", "start", "end", "strand", "repName", 
                   "repStart", "repEnd", "alu.id")
  alu.s.gr = makeGRangesFromDataFrame(alu.s, keep.extra.columns = T) 
    
  # Anti-sense
  alu.a = alu.s
  alu.a$strand = ifelse(alu.a$strand == "+", "-", "+")
  alu.a.gr = makeGRangesFromDataFrame(alu.a, keep.extra.columns = T) 
  
  # Remove alu in introns 
  if(ignore.intron.Alu){
    txdb = loadDb(txdb_path)
    introns = intronsByTranscript(txdb)
    alu.s.gr = subsetByOverlaps(alu.s.gr, introns, type = "within", invert = T)
    alu.a.gr = subsetByOverlaps(alu.a.gr, introns, type = "within", invert = T)
  }
  
  pA.alu.s = mergeByOverlaps(alu.s.gr, pA.gr, type="within")
  pA.alu.s$direction = "s"
  names(pA.alu.s)[1] = "alu.gr"
  
  pA.alu.a = mergeByOverlaps(alu.a.gr, pA.gr, type="within")
  pA.alu.a$direction = "a"
  names(pA.alu.a)[1] = "alu.gr"
  
  pA.alu = rbind(pA.alu.s, pA.alu.a)
  names(pA.alu) = sub("^pA.gr$", "utr3.gr", names(pA.alu))
  pA.alu$alu.gr = as.character(pA.alu$alu.gr)     
  pA.alu$utr3.gr = as.character(pA.alu$utr3.gr)  
  
  require(dplyr)
  require(tidyr)
  alu.df = as.data.frame(pA.alu) %>% 
    mutate(signed.start = as.numeric(sub(".+?:(\\d+?)-\\d+:([+-])", "\\2\\1", alu.gr))) %>% 
    group_by(gene_symbol, pAid, utr3.gr) %>%
    arrange(signed.start) %>%
    summarize(num_Alu = n(), 
              alu.directions = paste(direction, collapse=";"), 
              repNames = paste(repName, collapse=";"), 
              alu.ids = paste(alu.id, collapse=";"),
              alu.grs = paste(alu.gr, collapse=";")) %>%
    dplyr::filter(num_Alu < 100)
  
  
  alu.df$alu.s = !is.na(alu.df$alu.directions) & !grepl("a", alu.df$alu.directions) 
  alu.df$alu.a = !is.na(alu.df$alu.directions) & !grepl("s", alu.df$alu.directions)
  alu.df$alu.a.s = !is.na(alu.df$alu.directions) & grepl("a;s", alu.df$alu.directions)
  alu.df$alu.s.a = !is.na(alu.df$alu.directions) & grepl("s;a", alu.df$alu.directions) 
  
  require(stringr)
  # Number of possible pairs of a/s or s/a is the min of number of s and number of a.
  alu.df$num_IRAlu = pmin(str_count(alu.df$alu.directions, pattern = "s"),
                          str_count(alu.df$alu.directions, pattern = "a"))
  
  alu.df$Alu_type = ifelse(alu.df$alu.a.s | alu.df$alu.s.a, "IRAlu", alu.df$num_Alu)
  alu.df$Alu_type[!alu.df$Alu_type %in% c("1", "2", "IRAlu")] = ">=3"
  
  cols = c("pAid", "Alu_type", "num_Alu", "num_IRAlu", "alu.s", "alu.a", "alu.a.s", "alu.s.a")
  if(keep.alu.direction){
    cols = c(cols, "alu.directions")
  }
  
  alu.df = alu.df[, cols]
  df = merge(df, alu.df, by = "pAid", all.x=T)
  df$Alu_type[is.na(df$Alu_type)] = "0"
  df$num_Alu[is.na(df$num_Alu)] = 0
  df$num_IRAlu[is.na(df$num_IRAlu)] = 0
  
  df
}

countAllMotif = function(pas, geno = "mm9", search_from = 0, search_len = 50, motif_width = 4){
  #` Count number of motifs between search_from and search_len downstream of search_from

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
  # Fisher's Exact Test for APA
  # dat must contain "gene_symbol", two "_count" an two corresponding "_gene" columns
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


gene2neighorpAPairs = function(gene){
  #` For a gene with n pA isoforms, return n-1 combination of neighboring pA sites
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
  
  # Create neighboring pA pairs
  k = 1
  for (i in 1:(nrow(gene) -1)){
    comparedPAi[k] = gene$pAid[i]
    comparedPAj[k] = gene$pAid[i+1]
    k = k + 1
  }
  cbind(gene_id, gene_symbol, comparedPAi, comparedPAj)
  }
}

gene2pAPairs = function(gene){ 
  #` For a gene with n pA isoforms, return n*(n-1)/2 combination of pA sites
  if(nrow(gene) > 1){
  # Set up the columns
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
  
  # Compare all combinations of the APA isoforms
  k = 1
  for (i in 1:(nrow(gene) -1)){
    for (j in (i+1):nrow(gene)){ 
      comparedPAi[k] = gene$pAid[i]
      comparedPAj[k] = gene$pAid[j]
      k = k + 1
    }
  }
  cbind(gene_id, gene_symbol, comparedPAi, comparedPAj)
  }
  
}

get_red = function(data, cols = "UTR3_size", neighbor = F, toptwo = F, match_only=F){ 
  #` Calculate RED
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
  rownames(cmp) = NULL 
  
  # get strand and pAid information from the original data
  cmp = merge(cmp, subset(data, select = c("pAid", cols)), by.x = "comparedPAj", by.y = "pAid", sort = F)
  cmp = merge(cmp, subset(data, select = c("pAid", cols)), by.x = "comparedPAi", by.y = "pAid", sort = F)
  
  # change column names
  names(cmp) = sub("\\.x$", "_Dis", names(cmp))
  names(cmp) = sub("\\.y$", "_Prx", names(cmp))
  names(cmp) = sub("comparedPAi", "Prx_pA", names(cmp))
  names(cmp) = sub("comparedPAj", "Dis_pA", names(cmp))
  
  # calculate Dis to Prx ratios within each sample
  if(!match_only){
     if(any(grep("^n_", names(data)))){
      # calculate ratio of Dis to Prx pA read numbers
      for(string in sub("_Dis", "", grep("^n_.*_Dis", names(cmp), value = T, perl = T))){
        Dis.col.name = paste0(string, "_Dis")
        Prx.col.name = paste0(string, "_Prx")
        cmp[, paste0(string, "_d2p_ratio")] =
           log2(cmp[, Dis.col.name]) - log2(cmp[, Prx.col.name])
      }
    }

    if(any(grep("count$", names(data)))){
      # calculate ratio of Dis to Prx pA read numbers
      for(string in sub("_Dis", "", grep("(count.*)_Dis", names(cmp), value = T, perl = T))){
        Dis.col.name = paste0(string, "_Dis")
        Prx.col.name = paste0(string, "_Prx")
        cmp[, sub("counts?", "d2p_ratio", string)] =
           log2(cmp[, Dis.col.name]) - log2(cmp[, Prx.col.name])
      }
    }
  }
  
  # return
  cmp = merge(cmp, unique(data[, c("gene_symbol", "description")]), all.x = T, sort=F)
}

pAid2pos = function(pap){
  #` Convert pAid into a dataframe containing chr, pA_pos, and strand info
  pAid.i = pap$Prx_pA
  pAid.j = pap$Dis_pA
  
  pA.pos.df = as.data.frame(do.call(rbind, strsplit(pAid.i, "[+-]")), 
                            stringsAsFactors = F)
  names(pA.pos.df) = c("chr", "Prx_pA_pos")

  pA.pos.df$chr = as.character(pA.pos.df$chr)
  pA.pos.df$Prx_pA_pos = as.integer(pA.pos.df$Prx_pA_pos)
  pA.pos.df$strand = as.vector(sub(".*?([+-])\\d*", "\\1", pAid.i, perl = T))
  
  pA.pos.df$Dis_pA_pos = as.data.frame(do.call(rbind, strsplit(pAid.j, "[+-]")), 
                                          stringsAsFactors = F)[,2]
  pA.pos.df$Dis_pA_pos = as.integer(pA.pos.df$Dis_pA_pos)
  
  cbind(pA.pos.df, pap)
}

split.pAid.in.pas = function(pas){
  pA.pos.df = as.data.frame(do.call(rbind, strsplit(pas$pAid, "[+-]")), 
                            stringsAsFactors = F)
  names(pA.pos.df) = c("chr", "pA_pos")
  pA.pos.df$chr = as.character(pA.pos.df$chr)
  pA.pos.df$pA_pos = as.integer(pA.pos.df$pA_pos)
  pA.pos.df$strand = as.vector(sub(".*?([+-])\\d*", "\\1", pas$pAid, perl = T))
  cbind(pA.pos.df, pas)
}

match.pA.pair = function(df1, df2, radius = 12, keep = "intersect", rm_dup_col = T){
  #` Match aUTRs with slightly different pA sites in two data frames allowing distance of "radius" nt
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

match.pA = function(df1, df2, radius = 12, left_join=T, copy_right_columns = T){
  #` Match pA sites in two data frames, allowing distance of "radius" nt
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
                   sort = F, stringAsFactors = F, all=F) 
  data.cmn = merge(data.cmn, df2, by.x = "df2.pAid", by.y = "pAid", 
                   sort = F, stringAsFactors = F, all=F) 
  
  
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
  
  # return
  data.cmn
}
match.pAs = match.pA
match.pAs.2 = match.pA  


get_eCLIP_scores = function(grl = UTR5, 
                            grl_region = "UTR5",
                            bg_grl = NULL, # background grl (such as introns)
                            cell = "HepG2",
                            eCLIP_bed_file=file.path(SHARED_DATA_DIR, "BED_IDR_HepG2.csv")
                            ){
    # A function that sum eCLIP scores for eCLIP peaks mapped to each gr in grl but not in
    # bg_gr
    
    oldw = getOption("warn")
    options(warn = -1)

    require(dplyr)
    require(GenomicRanges)
    # Container for saving the scores
    scores = data.frame(tx_id = names(grl), stringsAsFactors=F)
    
    df = read.csv(eCLIP_bed_file)
    
    rbps = sub(paste0("_", cell, "_IDR"), "", unique(df$dataset_label))
    rbps = toupper(trimws(rbps))
    
    for(rbp in rbps){
        rbp_label = paste0(rbp, "_", cell, "_IDR")
        rbp_scores_df = subset(df, toupper(dataset_label)==toupper(rbp_label))
        # Create RBP GR
        rbp_gr = makeGRangesFromDataFrame(rbp_scores_df, keep.extra.columns=T)
        # Remove overlaps with background (such as introns)
        if(!is.null(bg_grl)){
            rbp_gr = rbp_gr[!rbp_gr %within% bg_grl] 
        }
        # Find overlaps
        olp = findOverlaps(rbp_gr, grl, type="within")
        # Alignment 
        if(length(olp)>0){# There could be no overlap! 
            # Align GRs
            olp_df = data.frame(tx_id = scores[subjectHits(olp), "tx_id"], 
                                rbp_score = rbp_scores_df[queryHits(olp), 
                                            'log2.eCLIP.fold.enrichment.over.size.matched.input.'])
            # Aggregate RBP scores
            olp_df %<>% group_by(tx_id) %>% mutate(rbp_score=sum(rbp_score)) %>% distinct()
        }else{
            olp_df = data.frame(tx_id = scores$tx_id, rbp_score = 0)
        }
        # Rename column 
        names(olp_df) = sub("^rbp_score$", paste0(grl_region, "_", rbp, "_eCLIP"), names(olp_df))
        # Add to scores  
        scores = merge(scores, olp_df, by="tx_id", all.x=T, sort=F)
    }
    options(warn = oldw)
    scores[is.na(scores)] = 0
    scores
}
