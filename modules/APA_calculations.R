require(MASS)
require(grid)
require(ggplot2)
require(gridExtra)


#### Calculate isoform ratios between fractions 
if(length(fractions) > 1){
  for(batch in batches){
    for(treatment in treatments){
      for(i in 1:length(fractions[-1])){
        for(j in (i+1):length(fractions)){
          # only keep desirable comparisons
          if((exists("comparisons") && paste(fractions[j], fractions[i], sep="_") %in% comparisons) | 
             !exists("comparisons")){
            col1 = paste(treatment, fractions[i], batch, "rpm", sep="_")
            col2 = paste(treatment, fractions[j], batch, "rpm", sep="_")
            if(all(c(col1, col2) %in% names(pas))){
              if(!paste(treatment, fractions[j], fractions[i], batch, "RPMR", sep="_") %in% names(pas)){ 
                pas[, paste(treatment, fractions[j], fractions[i], batch, "RPMR", sep = "_")] = 
                  log2(pas[,col2]) - log2(pas[,col1])
              }
            }
          }
        }
      }
    }
  }
  if(length(batches) > 1){
    for(treatment in treatments){
      col_to_be_avaraged = grep(paste0(treatment, "_", fractions[j], "_", fractions[i], ".+_RPMR"), names(pas), value=T)
      col_name_avarage = paste0(treatment, "_", fractions[j], "_", fractions[i], "_mean_RPMR")
      if(length(col_to_be_avaraged) > 1){
        pas[, col_name_avarage] = rowMeans(pas[, col_to_be_avaraged])
      }
    }
  }
}


#### Calculate isoform ratios between treatments 
if(length(treatments) > 1){
  for(batch in batches){
    for(fraction in fractions){
      for(i in 1:length(treatments[-1])){
        for(j in (i+1):length(treatments)){
          # only keep desirable comparisons
          if((exists("comparisons") && paste(treatments[j], treatments[i], sep="_") %in% comparisons) | 
             !exists("comparisons")){
            col1 = paste(treatments[i], fraction, batch, "rpm", sep="_")
            col2 = paste(treatments[j], fraction, batch, "rpm", sep="_")
            if(all(c(col1, col2) %in% names(pas))){
              if(!paste(treatments[j], treatments[i], fraction, batch, "RPMR", sep="_") %in% names(pas)){
                pas[, paste(treatments[j], treatments[i], fraction, batch, "RPMR", sep = "_")] = 
                  log2(pas[,col2]) - log2(pas[,col1])
              }
            }
          }
        }
      }
    }
  }
  if(length(batches) > 1){
    for(fraction in fractions){
      col_to_be_avaraged = grep(paste0(treatments[j], "_", treatments[i], "_", fraction,  ".+_RPMR"), names(pas), value=T)
      col_name_avarage = paste0(treatments[j], "_", treatments[i], "_", fraction, "_mean_RPMR")
      if(length(col_to_be_avaraged) > 1){
        pas[, col_name_avarage] = rowMeans(pas[, col_to_be_avaraged])
      }
    }
  }
}


#### Calculate delta abundance between fractions
if(length(fractions) > 1){
  for(i in 1:(length(fractions)-1)){
    for(j in (i+1):length(fractions)){
      # only keep desirable comparisons
      if((exists("comparisons") && paste(fractions[j], fractions[i], sep="_") %in% comparisons) | 
         !exists("comparisons")){
        for(treatment in treatments){
          for(batch in batches){
            col_name_1 =  paste(treatment, fractions[i], batch, "usage", sep="_")
            col_name_2 =  paste(treatment, fractions[j], batch, "usage", sep="_")
            if(all(c(col_name_1, col_name_2) %in% names(pas))){
              pas[, paste("Delta_abn", treatment, fractions[j], "vs", fractions[i], batch, sep="_")] = pas[, col_name_2]-pas[,col_name_1]
            }
          }
          if(length(batches) > 1){
            col_name_avarage = paste("Delta_abn", treatment, fractions[j], "vs", fractions[i], sep="_")
            if(length(grep(col_name_avarage, names(pas))) > 1){
              pas[, col_name_avarage] = rowMeans(pas[, grep(col_name_avarage, names(pas))])
            }
          }
        }
      }
    }
  }
}


#### Calculate delta abundance between treatments
if(length(treatments) > 1){
  for(i in 1:(length(treatments)-1)){
    for (j in (i+1):length(treatments)){
      # only keep desirable comparisons
      if((exists("comparisons") && paste(treatments[j], treatments[i], sep="_") %in% comparisons) | 
         !exists("comparisons")){
        for (fraction in fractions){
          for (batch in batches){
            col_name_1 =  paste(treatments[i], fraction, batch, "usage", sep="_")
            col_name_2 =  paste(treatments[j], fraction, batch, "usage", sep="_")
            if(all(c(col_name_1, col_name_2) %in% names(pas))){
              pas[, paste("Delta_abn", treatments[j], "vs", treatments[i], fraction, batch, sep="_")] = pas[, col_name_2]-pas[,col_name_1]
            }
          }
          if (length(batches) > 1){
            col_name_avarage = paste("Delta_abn", treatments[j], "vs", treatments[i], fraction, sep="_")
            if(length(grep(col_name_avarage, names(pas))) > 1){
              pas[, col_name_avarage] = rowMeans(pas[, grep(col_name_avarage, names(pas))])
            }
          }
        }
      }
    }
  }
}


#### Calculate APA p_values using fisher's exact test accross fractions
if(length(fractions) > 1){
  for(i in 1:(length(fractions)-1)){
    for (j in (i+1):length(fractions)){
      # only keep desirable comparisons
      if((exists("comparisons") && paste(fractions[j], fractions[i], sep="_") %in% comparisons) | 
         !exists("comparisons")){
        for(treatment in treatments){
          for (batch in batches){
            col_name_1 =  paste(treatment, fractions[i], batch, "count", sep="_")
            col_name_2 =  paste(treatment, fractions[j], batch, "count", sep="_")
            if(all(c(col_name_1, col_name_2) %in% names(pas))){
              col_name_1_gene = sub("count", "gene", col_name_1)
              col_name_2_gene = sub("count", "gene", col_name_2)
              # only keep useful columns
              dat = pas[,c("gene_symbol", col_name_1, col_name_2, col_name_1_gene, col_name_2_gene)]
              # calculate p values (adjusted within gene)
              p.vals = APA.fisher(dat)
              p.val.signs = sign(pas[, paste("Delta_abn", treatment, fractions[j], "vs", fractions[i], batch, sep="_")])
              pas[, paste("APA_SS", treatment, fractions[j], "vs", fractions[i], batch, sep="_")] = -1*p.val.signs * log10(p.vals)
            }
          }
        }
      }
    }
  }
}

#### Calculate APA p_values using fisher's exact test accross fractions accross treatments:
if(length(treatments) > 1){
  for(i in 1:(length(treatments)-1)){
    for(j in (i+1):length(treatments)){
      # only keep desirable comparisons
      if((exists("comparisons") && paste(treatments[j], treatments[i], sep="_") %in% comparisons) | 
         !exists("comparisons")){
        for(fraction in fractions){
          for(batch in batches){
            col_name_1 =  paste(treatments[i], fraction, batch, "count", sep="_")
            col_name_2 =  paste(treatments[j], fraction, batch, "count", sep="_")
            if(all(c(col_name_1, col_name_2) %in% names(pas))){
              col_name_1_gene = sub("count", "gene", col_name_1)
              col_name_2_gene = sub("count", "gene", col_name_2)
              # only keep useful columns
              dat = pas[,c("gene_symbol", col_name_1, col_name_2, col_name_1_gene, col_name_2_gene)]
              # calculate p values (adjusted within gene)
              p.vals = APA.fisher(dat)
              p.val.signs = sign(pas[, paste("Delta_abn", treatments[j], "vs", treatments[i], fraction, batch, sep="_")])
              pas[, paste("APA_SS", treatments[j], "vs", treatments[i], fraction, batch, sep="_")] = -1*p.val.signs * log10(p.vals)
            }
          }
        }
      }
    }
  }
}


#### Calculate genewise foldchange using the RPM accross fractions:
if(length(fractions) > 1){
  for(i in 1:(length(fractions)-1)){
    for(j in (i+1):length(fractions)){
      # only keep desirable comparisons
      if((exists("comparisons") && paste(fractions[j], fractions[i], sep="_") %in% comparisons) | 
         !exists("comparisons")){
        for(treatment in treatments){
          for(batch in batches){
            col_name_1 =  paste("Gene_RPM", treatment, fractions[i], batch, sep="_")
            col_name_2 =  paste("Gene_RPM", treatment, fractions[j], batch, sep="_")
            if(all(c(col_name_1, col_name_2) %in% names(pas))){
              pas[, paste("Gene_FC", treatment, fractions[j], "vs", fractions[i], batch, sep="_")] = log2(pas[, col_name_2])-log2(pas[,col_name_1]) 
            }
          }
          if(length(batches) > 1){
            col_name_avarage = paste("Gene_FC", treatment, fractions[j], "vs", fractions[i], sep="_")
            if(length(grep(col_name_avarage, names(pas)))>1){
              pas[, col_name_avarage] = rowMeans(pas[, grep(col_name_avarage, names(pas))])
            }
          }
        }
      }
    }
  }
}


#### Calculate genewise foldchange using the RPM accross treatments:
if(length(treatments) > 1){
  for(i in 1:(length(treatments)-1)){
    for(j in (i+1):length(treatments)){
      # only keep desirable comparisons
      if((exists("comparisons") && paste(treatments[j], treatments[i], sep="_") %in% comparisons) | 
         !exists("comparisons")){
        for(fraction in fractions){
          for(batch in batches){
            col_name_1 =  paste("Gene_RPM", treatments[i], fraction, batch, sep="_")
            col_name_2 =  paste("Gene_RPM", treatments[j], fraction, batch, sep="_")
            if(all(c(col_name_1, col_name_2) %in% names(pas))){
              pas[, paste("Gene_FC", treatments[j], "vs", treatments[i], fraction, batch, sep="_")] = log2(pas[, col_name_2])-log2(pas[,col_name_1])
            }
          }
          if (length(batches) > 1){
            col_name_avarage = paste("Gene_FC", treatments[j], "vs", treatments[i], fraction, sep="_")
            if(length(grep(col_name_avarage, names(pas)))>1){
              pas[, col_name_avarage] = rowMeans(pas[, grep(col_name_avarage, names(pas))])
            }
          }
        }
      }
    }
  }
}


#### 3'UTR RPM accross treatments:
if(length(treatments) > 1){
  for(i in 1:(length(treatments)-1)){
    for (j in (i+1):length(treatments)){
      # only keep desirable comparisons
      if((exists("comparisons") && paste(treatments[j], treatments[i], sep="_") %in% comparisons) | 
         !exists("comparisons")){
        for(fraction in fractions){
          for (batch in batches){
            col_name_1 =  paste("UTR_RPM", treatments[i], fraction, batch, sep="_")
            col_name_2 =  paste("UTR_RPM", treatments[j], fraction, batch, sep="_")
            if(all(c(col_name_1, col_name_2) %in% names(pas))){
              pas[, paste("UTR_FC", treatments[j], "vs", treatments[i], fraction, batch, sep="_")] = log2(pas[, col_name_2])-log2(pas[,col_name_1])
            }
          }
          if(length(batches) > 1){
            col_name_avarage = paste("UTR_FC", treatments[j], "vs", treatments[i], fraction, sep="_")
            if(length(grep(col_name_avarage, names(pas)))>1){
              pas[, col_name_avarage] = rowMeans(pas[, grep(col_name_avarage, names(pas))])
            }
          }
        }
      }
    }
  }
}


#### Calculate p values for genewise foldchange using chi-squared test (too slow with fisher's exact)
gene.df = unique(pas[,c("gene_symbol", grep("_gene", names(pas), value = T))])
#### Accross fractions:
if(length(fractions) > 1){
  for(i in 1:(length(fractions)-1)){
    for (j in (i+1):length(fractions)){
      # only keep desirable comparisons
      if((exists("comparisons") && paste(fractions[j], fractions[i], sep="_") %in% comparisons) | 
         !exists("comparisons")){
        for(treatment in treatments){
          for (batch in batches){
            col_name_1 =  paste(treatment, fractions[i], batch, "gene", sep="_")
            col_name_2 =  paste(treatment, fractions[j], batch, "gene", sep="_")
            
            if(all(c(col_name_1, col_name_2) %in% names(pas))){
              # prepare input
              dat = cbind(gene.df[,c(col_name_1, col_name_2)], 
                          col_name_1_sum = sum(gene.df[,col_name_1], na.rm = T),
                          col_name_2_sum = sum(gene.df[,col_name_2], na.rm = T))
              
              # calculate p values (adjusted)
              gene.df$p.vals = genewise.chi(dat)
              pas = merge(pas, gene.df[,c("gene_symbol", "p.vals")], by = "gene_symbol", sort=F)
              
              p.val.signs = sign(pas[, paste("Gene_FC", treatment, fractions[j], "vs", fractions[i], batch, sep="_")])
              pas$p.vals = -1*p.val.signs * log10(pas$p.vals)
              names(pas) = sub("^p\\.vals$", paste("Gene_SS", treatment, fractions[j], "vs", fractions[i], batch, sep="_"), names(pas))
            }
          }
        }
      }
    }
  }
}
rm(dat)
#### Accross treatments:
if(length(treatments) > 1){
  for(i in 1:(length(treatments)-1)){
    for(j in (i+1):length(treatments)){
      # only keep desirable comparisons
      if((exists("comparisons") && paste(treatments[j], treatments[i], sep="_") %in% comparisons) | 
         !exists("comparisons")){
        for(fraction in fractions){
          for(batch in batches){
            col_name_1 =  paste(treatments[i], fraction, batch, "gene", sep="_")
            col_name_2 =  paste(treatments[j], fraction, batch, "gene", sep="_")
            
            if(all(c(col_name_1, col_name_2) %in% names(pas))){
              # prepare input
              dat = cbind(gene.df[,c(col_name_1, col_name_2)], 
                          col_name_1_sum = sum(gene.df[,col_name_1], na.rm = T),
                          col_name_2_sum = sum(gene.df[,col_name_2], na.rm = T))
              
              # calculate p values (adjusted)
              gene.df$p.vals = genewise.chi(dat)
              pas = merge(pas, gene.df[,c("gene_symbol", "p.vals")], by = "gene_symbol", sort=F)
              
              p.val.signs = sign(pas[, paste("Gene_FC", treatments[j], "vs", treatments[i], fraction, batch, sep="_")])
              pas$p.vals = -1*p.val.signs * log10(pas$p.vals)
              names(pas) = sub("^p\\.vals$", paste("Gene_SS", treatments[j], "vs", treatments[i], fraction, batch, sep="_"), names(pas))
            }
          }
        }
      }
    }
  }
}
rm(dat)
rm(gene.df)


#### Replace the NaN values in the APA_SS columns with 0
pas[, grep("_SS", names(pas))][is.na(pas[, grep("_SS", names(pas))])] = 0

#### Change the column names containing "usage"
names(pas) = sub("(.+)_usage", "Isoform_abn_\\1", names(pas))


#### Calculate REDs
pA3utr = subset(pas, CDS_size > 0 & region == "3UTR" & !is.na(exonic_3UTR_seq)) 
# Add pseudocounts
if(!exists("USE_PSEUDOCOUNTS") || USE_PSEUDOCOUNTS == T){
  pA3utr[,grep("_count$", names(pA3utr))] = pA3utr[, grep("_count$", names(pA3utr))] + 1
}
colname_pattern = "exonic_3UTR_seq|UTR3_size|_counts$|_sites$"
if(exists("additional_RED_pattern")){
  colname_pattern = paste0(colname_pattern, "|", additional_RED_pattern)
}
pap = get_red(data = pA3utr, cols = grep(colname_pattern, names(pA3utr), value=T), 
              neighbor = neighbor, toptwo = toptwo, match_only=F)

# Calculate aUTR lengths
pap = pAid2pos(pap)
# The distribution of aUTR length
pap$len = abs(pap$Dis_pA_pos - pap$Prx_pA_pos) # default method
index = which(pap$UTR3_size_Dis > 0 & pap$UTR3_size_Prx > 0)
pap[index,]$len = pap[index,]$UTR3_size_Dis - pap[index,]$UTR3_size_Prx
pap = pap[pap$len > 0, ]
pap$len_log10 = log10(pap$len)

ggplot(data = pap, aes(x = pap$len_log10)) + geom_histogram(fill="blue", colour="blue", binwidth = 0.05)

# Calculate relative aUTR length
pap = merge(pap, unique(pas[, c("gene_symbol", "CDS_end", "CDS_size")]), all.x=T, sort=F)
pap = subset(pap, CDS_size > 0)
pap = pap[!is.na(pap$CDS_end),]
pap$cUTRlen = abs(pap$CDS_end - pap$Prx_pA_pos)
pap$relative_len = log2(abs(pap$CDS_end - pap$Dis_pA_pos)) - log2(pap$cUTRlen)
pap$len_bin = cut(pap$len, breaks = quantile(pap$len,probs = seq(0, 1, 0.2), na.rm = T), include.lowest = T)
pap$relative_len_bin = cut(pap$relative_len, breaks = quantile(pap$relative_len,probs = seq(0, 1, 0.2), na.rm = T), include.lowest = T)

Num_pA = as.data.frame(table(pA3utr$gene_symbol))
names(Num_pA) = c("gene_symbol", "Num_pA")
pap = merge(pap, Num_pA, by = "gene_symbol", all.x = T, sort=F)
rm(Num_pA)


#### calculate and plot RED accross fractions:
if(length(fractions) > 1){
  p.list = list()
  k = 1
  for(i in 1:(length(fractions)-1)){
    for(j in (i+1):length(fractions)){
      # only keep desirable comparisons
      if((exists("comparisons") && paste(fractions[j], fractions[i], sep="_") %in% comparisons) | 
         !exists("comparisons")){
        for(treatment in treatments){
          for(batch in batches){
            col_name_1 =  paste(treatment, fractions[i], batch, "DPR", sep="_")
            col_name_2 =  paste(treatment, fractions[j], batch, "DPR", sep="_")
            if(all(c(col_name_1, col_name_2) %in% names(pap))){
              pap[, paste(treatment, fractions[j], "vs", fractions[i], batch, "RED", sep="_")] = pap[, col_name_2]-pap[,col_name_1] 
            }
          }
          if(length(batches) > 1){
            col_to_be_avaraged = grep(paste0(treatment, "_", fractions[j], "_vs_", fractions[i], ".+_RED"), names(pap), value=T)
            col_name_avarage = paste0(treatment, "_", fractions[j], "_vs_", fractions[i], "_mean_RED")
            if(length(col_to_be_avaraged) > 1){
              pap[, col_name_avarage] = rowMeans(pap[, col_to_be_avaraged])
            }
          }
        }
      }
    }
  }
}


#### Calculate APA accross treatments:
if(length(treatments) > 1){
  for(i in 1:(length(treatments)-1)){
    for(j in (i+1):length(treatments)){
      # only keep desirable comparisons
      if((exists("comparisons") && paste(treatments[j], treatments[i], sep="_") %in% comparisons) | 
         !exists("comparisons")){
        for(fraction in fractions){
          for(batch in batches){
            col_name_1 =  paste(treatments[i], fraction, batch, "DPR", sep="_")
            col_name_2 =  paste(treatments[j], fraction, batch, "DPR", sep="_")
            if(all(c(col_name_1, col_name_2) %in% names(pap))){
              pap[, paste(treatments[j], "vs", treatments[i], fraction, batch, "RED", sep="_")] = pap[, col_name_2]-pap[,col_name_1]
            }
          }
          if(length(batches) > 1){
            col_to_be_avaraged = grep(paste0(treatments[j], "_vs_", treatments[i], "_", fraction,  ".+_RED"), names(pap), value=T)
            col_name_avarage = paste0(treatments[j], "_vs_", treatments[i], "_", fraction, "_mean_RED")
            if(length(col_to_be_avaraged) > 1){
              pap[, col_name_avarage] = rowMeans(pap[, col_to_be_avaraged])
            }
          }
        }
      }
    }
  }
}


keeped_columns = c("gene_symbol", "uORF", "UTR5_size", "UTR5_GC", "CDS_size", "CDS_GC", "max_intron_size", "min_intron_size", 
                   "total_intron_size", "num_intron", "description")
pap = merge(pap, unique(pA3utr[, keeped_columns]), all.x = T, sort=F)
pap = pap[, -c(grep("^chr$|pA_pos$|strand|_log10$", names(pap)))]

if(exists("additional_RED_pattern")){
  pap = pap[, c("gene_symbol","gene_id", "Prx_pA", "Dis_pA", "UTR3_size_Prx", "UTR3_size_Dis", "len", 
                        grep(paste0(c(sub("\\$$", "", additional_RED_pattern), "^exonic_", "_RED$", "_count_"), collapse = "|"), names(pap), value=T), 
                        keeped_columns)]
}else{
  pap = pap[, c("gene_symbol","gene_id", "Prx_pA", "Dis_pA", "UTR3_size_Prx", "UTR3_size_Dis", "len", 
                        grep(paste0(c("^exonic_", "_RED$", "_count_"), collapse = "|"), names(pap), value=T), 
                        keeped_columns)]
}

# Add 3'UTR GC
require(Biostrings)
cUTR_seq = DNAStringSet(pap$exonic_3UTR_seq_Prx)
base_fraction = alphabetFrequency(cUTR_seq, baseOnly=T, as.prob=T)
pap$cUTR_GC = base_fraction[, "C"] + base_fraction[, "G"]

fUTR_seq = DNAStringSet(pap$exonic_3UTR_seq_Dis)
base_fraction = alphabetFrequency(fUTR_seq, baseOnly=T, as.prob=T)
pap$fUTR_GC = base_fraction[, "C"] + base_fraction[, "G"]

aUTR_seq = apply(pap, 1, 
                 function(row) substr(row["exonic_3UTR_seq_Dis"], 
                                      start = as.integer(row["UTR3_size_Prx"])+1, 
                                      stop = row["UTR3_size_Dis"]))
aUTR_seq = DNAStringSet(aUTR_seq)
base_fraction = alphabetFrequency(aUTR_seq, baseOnly=T, as.prob=T)
pap$aUTR_GC = base_fraction[, "C"] + base_fraction[, "G"]


names(pap) = sub("^len$", "aUTR_size", names(pap))
write.csv(pap, file.path(result_dir, "pap.csv"), row.names = F)

rm(pA3utr)

