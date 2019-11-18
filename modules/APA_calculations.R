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
            if(all(c(col1, col2) %in% names(pA.df))){
              if(!paste(treatment, fractions[j], fractions[i], batch, "rpm_ratio", sep="_") %in% names(pA.df)){ 
                pA.df[, paste(treatment, fractions[j], fractions[i], batch, "rpm_ratio", sep = "_")] = 
                  log2(pA.df[,col2]) - log2(pA.df[,col1])
              }
            }
          }
        }
      }
    }
  }
  if(length(batches) > 1){
    for(treatment in treatments){
      col_to_be_avaraged = grep(paste0(treatment, "_", fractions[j], "_", fractions[i], ".+_rpm_ratio"), names(pA.df), value=T)
      col_name_avarage = paste0(treatment, "_", fractions[j], "_", fractions[i], "_mean_rpm_ratio")
      if(length(col_to_be_avaraged) > 1){
        pA.df[, col_name_avarage] = rowMeans(pA.df[, col_to_be_avaraged])
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
            if(all(c(col1, col2) %in% names(pA.df))){
              if(!paste(treatments[j], treatments[i], fraction, batch, "rpm_ratio", sep="_") %in% names(pA.df)){
                pA.df[, paste(treatments[j], treatments[i], fraction, batch, "rpm_ratio", sep = "_")] = 
                  log2(pA.df[,col2]) - log2(pA.df[,col1])
              }
            }
          }
        }
      }
    }
  }
  if(length(batches) > 1){
    for(fraction in fractions){
      col_to_be_avaraged = grep(paste0(treatments[j], "_", treatments[i], "_", fraction,  ".+_rpm_ratio"), names(pA.df), value=T)
      col_name_avarage = paste0(treatments[j], "_", treatments[i], "_", fraction, "_mean_rpm_ratio")
      if(length(col_to_be_avaraged) > 1){
        pA.df[, col_name_avarage] = rowMeans(pA.df[, col_to_be_avaraged])
      }
    }
  }
}


#### Calculate delta abundance between accross fractions
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
            if(all(c(col_name_1, col_name_2) %in% names(pA.df))){
              pA.df[, paste("Delta_abn", treatment, fractions[j], "vs", fractions[i], batch, sep="_")] = pA.df[, col_name_2]-pA.df[,col_name_1]
            }
          }
          if(length(batches) > 1){
            col_name_avarage = paste("Delta_abn", treatment, fractions[j], "vs", fractions[i], sep="_")
            if(length(grep(col_name_avarage, names(pA.df))) > 1){
              pA.df[, col_name_avarage] = rowMeans(pA.df[, grep(col_name_avarage, names(pA.df))])
            }
          }
        }
      }
    }
  }
}


#### Calculate delta abundance between accross treatments
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
            if(all(c(col_name_1, col_name_2) %in% names(pA.df))){
              pA.df[, paste("Delta_abn", treatments[j], "vs", treatments[i], fraction, batch, sep="_")] = pA.df[, col_name_2]-pA.df[,col_name_1]
            }
          }
          if (length(batches) > 1){
            col_name_avarage = paste("Delta_abn", treatments[j], "vs", treatments[i], fraction, sep="_")
            if(length(grep(col_name_avarage, names(pA.df))) > 1){
              pA.df[, col_name_avarage] = rowMeans(pA.df[, grep(col_name_avarage, names(pA.df))])
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
            if(all(c(col_name_1, col_name_2) %in% names(pA.df))){
              col_name_1_gene = sub("count", "gene", col_name_1)
              col_name_2_gene = sub("count", "gene", col_name_2)
              # only keep useful columns
              dat = pA.df[,c("gene_symbol", col_name_1, col_name_2, col_name_1_gene, col_name_2_gene)]
              # calculate p values (adjusted within gene)
              p.vals = APA.fisher(dat)
              p.val.signs = sign(pA.df[, paste("Delta_abn", treatment, fractions[j], "vs", fractions[i], batch, sep="_")])
              pA.df[, paste("APA_SS", treatment, fractions[j], "vs", fractions[i], batch, sep="_")] = -1*p.val.signs * log10(p.vals)
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
            if(all(c(col_name_1, col_name_2) %in% names(pA.df))){
              col_name_1_gene = sub("count", "gene", col_name_1)
              col_name_2_gene = sub("count", "gene", col_name_2)
              # only keep useful columns
              dat = pA.df[,c("gene_symbol", col_name_1, col_name_2, col_name_1_gene, col_name_2_gene)]
              # calculate p values (adjusted within gene)
              p.vals = APA.fisher(dat)
              p.val.signs = sign(pA.df[, paste("Delta_abn", treatments[j], "vs", treatments[i], fraction, batch, sep="_")])
              pA.df[, paste("APA_SS", treatments[j], "vs", treatments[i], fraction, batch, sep="_")] = -1*p.val.signs * log10(p.vals)
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
            if(all(c(col_name_1, col_name_2) %in% names(pA.df))){
              pA.df[, paste("Gene_FC", treatment, fractions[j], "vs", fractions[i], batch, sep="_")] = log2(pA.df[, col_name_2])-log2(pA.df[,col_name_1]) 
            }
          }
          if(length(batches) > 1){
            col_name_avarage = paste("Gene_FC", treatment, fractions[j], "vs", fractions[i], sep="_")
            if(length(grep(col_name_avarage, names(pA.df)))>1){
              pA.df[, col_name_avarage] = rowMeans(pA.df[, grep(col_name_avarage, names(pA.df))])
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
            if(all(c(col_name_1, col_name_2) %in% names(pA.df))){
              pA.df[, paste("Gene_FC", treatments[j], "vs", treatments[i], fraction, batch, sep="_")] = log2(pA.df[, col_name_2])-log2(pA.df[,col_name_1])
            }
          }
          if (length(batches) > 1){
            col_name_avarage = paste("Gene_FC", treatments[j], "vs", treatments[i], fraction, sep="_")
            if(length(grep(col_name_avarage, names(pA.df)))>1){
              pA.df[, col_name_avarage] = rowMeans(pA.df[, grep(col_name_avarage, names(pA.df))])
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
            if(all(c(col_name_1, col_name_2) %in% names(pA.df))){
              pA.df[, paste("UTR_FC", treatments[j], "vs", treatments[i], fraction, batch, sep="_")] = log2(pA.df[, col_name_2])-log2(pA.df[,col_name_1])
            }
          }
          if(length(batches) > 1){
            col_name_avarage = paste("UTR_FC", treatments[j], "vs", treatments[i], fraction, sep="_")
            if(length(grep(col_name_avarage, names(pA.df)))>1){
              pA.df[, col_name_avarage] = rowMeans(pA.df[, grep(col_name_avarage, names(pA.df))])
            }
          }
        }
      }
    }
  }
}


#### Calculate p values for genewise foldchange using chi-squared test (too slow with fisher's exact)
gene.df = unique(pA.df[,c("gene_symbol", grep("_gene", names(pA.df), value = T))])
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
            
            if(all(c(col_name_1, col_name_2) %in% names(pA.df))){
              # prepare input
              dat = cbind(gene.df[,c(col_name_1, col_name_2)], 
                          col_name_1_sum = sum(gene.df[,col_name_1], na.rm = T),
                          col_name_2_sum = sum(gene.df[,col_name_2], na.rm = T))
              
              # calculate p values (adjusted)
              gene.df$p.vals = genewise.chi(dat)
              pA.df = merge(pA.df, gene.df[,c("gene_symbol", "p.vals")], by = "gene_symbol", sort=F)
              
              p.val.signs = sign(pA.df[, paste("Gene_FC", treatment, fractions[j], "vs", fractions[i], batch, sep="_")])
              pA.df$p.vals = -1*p.val.signs * log10(pA.df$p.vals)
              names(pA.df) = sub("^p\\.vals$", paste("Gene_SS", treatment, fractions[j], "vs", fractions[i], batch, sep="_"), names(pA.df))
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
            
            if(all(c(col_name_1, col_name_2) %in% names(pA.df))){
              # prepare input
              dat = cbind(gene.df[,c(col_name_1, col_name_2)], 
                          col_name_1_sum = sum(gene.df[,col_name_1], na.rm = T),
                          col_name_2_sum = sum(gene.df[,col_name_2], na.rm = T))
              
              # calculate p values (adjusted)
              gene.df$p.vals = genewise.chi(dat)
              pA.df = merge(pA.df, gene.df[,c("gene_symbol", "p.vals")], by = "gene_symbol", sort=F)
              
              p.val.signs = sign(pA.df[, paste("Gene_FC", treatments[j], "vs", treatments[i], fraction, batch, sep="_")])
              pA.df$p.vals = -1*p.val.signs * log10(pA.df$p.vals)
              names(pA.df) = sub("^p\\.vals$", paste("Gene_SS", treatments[j], "vs", treatments[i], fraction, batch, sep="_"), names(pA.df))
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
pA.df[, grep("_SS", names(pA.df))][is.na(pA.df[, grep("_SS", names(pA.df))])] = 0

#### Change the column names containing "usage"
names(pA.df) = sub("(.+)_usage", "Isoform_abn_\\1", names(pA.df))


#### Calculate REDs
pA3utr = subset(pA.df, CDS_size > 0 & region == "3UTR" & !is.na(exonic_3UTR_seq)) 
# Add pseudocounts
if(!exists("USE_PSEUDOCOUNTS") || USE_PSEUDOCOUNTS == T){
  pA3utr[,grep("_count$", names(pA3utr))] = pA3utr[, grep("_count$", names(pA3utr))] + 1
}
colname_pattern = "exonic_3UTR_seq|UTR_length|_counts$|_sites$"
if(exists("additional_RED_pattern")){
  colname_pattern = paste0(colname_pattern, "|", additional_RED_pattern)
}
pair.pA = get_red(data = pA3utr, cols = grep(colname_pattern, names(pA3utr), value=T), 
                  neighbor = neighbor, toptwo = toptwo, match_only=F)

# Calculate aUTR lengths
pair.pA = pAid2pos(pair.pA)
# The distribution of aUTR length
pair.pA$len = abs(pair.pA$distal_pA_pos - pair.pA$proximal_pA_pos) # default method
index = which(pair.pA$UTR_length_distal > 0 & pair.pA$UTR_length_proximal > 0)
pair.pA[index,]$len = pair.pA[index,]$UTR_length_distal - pair.pA[index,]$UTR_length_proximal
pair.pA = pair.pA[pair.pA$len > 0, ]
pair.pA$len_log10 = log10(pair.pA$len)

ggplot(data = pair.pA, aes(x = pair.pA$len_log10)) + geom_histogram(fill="blue", colour="blue", binwidth = 0.05)

# Calculate relative aUTR length
pair.pA = merge(pair.pA, unique(pA.df[, c("gene_symbol", "CDS_end", "CDS_size")]), all.x=T, sort=F)
pair.pA = subset(pair.pA, CDS_size > 0)
pair.pA = pair.pA[!is.na(pair.pA$CDS_end),]
pair.pA$cUTRlen = abs(pair.pA$CDS_end - pair.pA$proximal_pA_pos)
pair.pA$relative_len = log2(abs(pair.pA$CDS_end - pair.pA$distal_pA_pos)) - log2(pair.pA$cUTRlen)
pair.pA$len_bin = cut(pair.pA$len, breaks = quantile(pair.pA$len,probs = seq(0, 1, 0.2), na.rm = T), include.lowest = T)
pair.pA$relative_len_bin = cut(pair.pA$relative_len, breaks = quantile(pair.pA$relative_len,probs = seq(0, 1, 0.2), na.rm = T), include.lowest = T)

Num_pA = as.data.frame(table(pA3utr$gene_symbol))
names(Num_pA) = c("gene_symbol", "Num_pA")
pair.pA = merge(pair.pA, Num_pA, by = "gene_symbol", all.x = T, sort=F)
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
            col_name_1 =  paste(treatment, fractions[i], batch, "d2p_ratio", sep="_")
            col_name_2 =  paste(treatment, fractions[j], batch, "d2p_ratio", sep="_")
            if(all(c(col_name_1, col_name_2) %in% names(pair.pA))){
              pair.pA[, paste(treatment, fractions[j], "vs", fractions[i], batch, "RED", sep="_")] = pair.pA[, col_name_2]-pair.pA[,col_name_1] 
            }
          }
          if(length(batches) > 1){
            col_to_be_avaraged = grep(paste0(treatment, "_", fractions[j], "_vs_", fractions[i], ".+_RED"), names(pair.pA), value=T)
            col_name_avarage = paste0(treatment, "_", fractions[j], "_vs_", fractions[i], "_mean_RED")
            if(length(col_to_be_avaraged) > 1){
              pair.pA[, col_name_avarage] = rowMeans(pair.pA[, col_to_be_avaraged])
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
            col_name_1 =  paste(treatments[i], fraction, batch, "d2p_ratio", sep="_")
            col_name_2 =  paste(treatments[j], fraction, batch, "d2p_ratio", sep="_")
            if(all(c(col_name_1, col_name_2) %in% names(pair.pA))){
              pair.pA[, paste(treatments[j], "vs", treatments[i], fraction, batch, "RED", sep="_")] = pair.pA[, col_name_2]-pair.pA[,col_name_1]
            }
          }
          if(length(batches) > 1){
            col_to_be_avaraged = grep(paste0(treatments[j], "_vs_", treatments[i], "_", fraction,  ".+_RED"), names(pair.pA), value=T)
            col_name_avarage = paste0(treatments[j], "_vs_", treatments[i], "_", fraction, "_mean_RED")
            if(length(col_to_be_avaraged) > 1){
              pair.pA[, col_name_avarage] = rowMeans(pair.pA[, col_to_be_avaraged])
            }
          }
        }
      }
    }
  }
}


keeped_columns = c("gene_symbol", "uORF", "UTR5_size", "UTR5_GC", "CDS_size", "CDS_GC", "max_intron_size", "min_intron_size", 
                   "total_intron_size", "num_intron", "description")
pair.pA = merge(pair.pA, unique(pA3utr[, keeped_columns]), all.x = T, sort=F)
pair.pA = pair.pA[, -c(grep("^chr$|pA_pos$|strand|_log10$", names(pair.pA)))]

if(exists("additional_RED_pattern")){
  pair.pA = pair.pA[, c("gene_symbol","gene_id", "proximal_pA", "distal_pA", "UTR_length_proximal", "UTR_length_distal", "len", 
                        grep(paste0(c(sub("\\$$", "", additional_RED_pattern), "^exonic_", "_RED$", "_count_"), collapse = "|"), names(pair.pA), value=T), 
                        keeped_columns)]
}else{
  pair.pA = pair.pA[, c("gene_symbol","gene_id", "proximal_pA", "distal_pA", "UTR_length_proximal", "UTR_length_distal", "len", 
                        grep(paste0(c("^exonic_", "_RED$", "_count_"), collapse = "|"), names(pair.pA), value=T), 
                        keeped_columns)]
}

# Add 3'UTR GC
require(Biostrings)
cUTR_seq = DNAStringSet(pair.pA$exonic_3UTR_seq_proximal)
base_fraction = alphabetFrequency(cUTR_seq, baseOnly=T, as.prob=T)
pair.pA$cUTR_GC = base_fraction[, "C"] + base_fraction[, "G"]

fUTR_seq = DNAStringSet(pair.pA$exonic_3UTR_seq_distal)
base_fraction = alphabetFrequency(fUTR_seq, baseOnly=T, as.prob=T)
pair.pA$fUTR_GC = base_fraction[, "C"] + base_fraction[, "G"]

aUTR_seq = apply(pair.pA, 1, 
                 function(row) substr(row["exonic_3UTR_seq_distal"], 
                                      start = as.integer(row["UTR_length_proximal"])+1, 
                                      stop = row["UTR_length_distal"]))
aUTR_seq = DNAStringSet(aUTR_seq)
base_fraction = alphabetFrequency(aUTR_seq, baseOnly=T, as.prob=T)
pair.pA$aUTR_GC = base_fraction[, "C"] + base_fraction[, "G"]

# # Include CDS size
# gene_CDS_size = unique(pA.df[, c("gene_symbol", "CDS_size")]) %>%
#   filter(CDS_size > 3) %>%
#   group_by(gene_symbol) %>%
#   top_n(1, CDS_size) %>%
#   ungroup()
# pair.pA = left_join(pair.pA, gene_CDS_size)   


names(pair.pA) = sub("^len$", "aUTR_length", names(pair.pA))
write.csv(pair.pA, file.path(result_dir, "pair.pA.csv"), row.names = F)

rm(pA3utr)

