library(MASS) 
library(ggrepel) 
require(tidyr)
require(dplyr)
require(ggplot2)
theme_set(theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5))) 

# fisher's exact test: accept a row of 4 integers and return a p value. 
# Can be used like apply(x, 1, row.fisher), where x is a 4 column matrix
row.fisher = function(rowdat){
  m = matrix(rowdat, nrow = 2, byrow=T)
  fisher.test(m)$p.value
}

genewise.APA.fisher = function(gene.df){
  # gene.df is a dataframe with 4 columns containing integers
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
  # dat must contain "gene_symbol", two "_count" an two corresponding "_gene" columns
  # the "_gene" column contains total read counts for each gene
  # split the data 
  dat.lst = split(dat, dat$gene_symbol)
  # work on one gene each time
  p.vals = sapply(dat.lst, function(x) genewise.APA.fisher(x[, -1])) #the first column is not needed
  # recombind the data with results
  p.vals = unlist(p.vals)
}

row.chi = function(rowdat){
  #library(MASS)
  m = matrix(rowdat, nrow = 2, byrow=T)
  if(any(is.na(m)) | any(m < 0)){
    p_val = 1
  }else{
    chisq.test(m)$p.value
  }
}

genewise.chi = function(dat){
  # dat is a dataframe with 4 columns containing integers
  # all the cells in the 3rd/4th column contain the sum of the 1st/2nd column
  dat[, 3] = dat[, 3] - dat[, 1]
  dat[, 4] = dat[, 4] - dat[, 2]
  
  p_val = p.adjust(apply(dat, 1, row.chi), "BH")
}

## "Mixed Case" Capitalizing - toupper( every first letter of a word ) :
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}


########################################## align pAs of the same gene ##################
#### for a gene with n pA isoforms, return n-1 combination of neighboring pA sites
gene2neighorpAPairs = function(gene){
  if(nrow(gene) > 1){
    # set up the columns
  gene_id = rep(gene$gene_id[1],nrow(gene)-1) 
  gene_symbol = rep(gene$gene_symbol[1], nrow(gene)-1) 
  comparedPAi = rep(NA, length(gene_id)) 
  comparedPAj = rep(NA, length(gene_id)) 
  
  # Sort the pAs according to their position. 
  # This ensure that the pAi is always the proximal pA compared to pAj
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
  # This ensure that the pAi is always the proximal pA compared to pAj
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
  
  names(cmp) = sub("\\.x$", "_distal", names(cmp))
  names(cmp) = sub("\\.y$", "_proximal", names(cmp))
  names(cmp) = sub("comparedPAi", "proximal_pA", names(cmp))
  names(cmp) = sub("comparedPAj", "distal_pA", names(cmp))
  
  if(match_only == F){
    if(any(grep("^RAI", names(data)))){
      # calculate delta RAIs
      for(string in sub("_distal", "", grep("(RAI.*)_distal", names(cmp), value = T, perl = T))){
        distal.col.name = paste0(string, "_distal")
        proximal.col.name = paste0(string, "_proximal")
        cmp[, paste0("delt_", string)] = cmp[, distal.col.name] - cmp[, proximal.col.name]
      }
    }
    
    if(any(grep("^RPM", names(data)))){
      # calculate ratio of distal to proximal pA RPMs
      for(string in sub("_distal", "", grep("(RPM.*)_distal", names(cmp), value = T, perl = T))){
        
        distal.col.name = paste0(string, "_distal")
        proximal.col.name = paste0(string, "_proximal")
        cmp[, paste0("ratio_", sub("^RPM_", "", string))] = log2(cmp[, distal.col.name]) - log2(cmp[, proximal.col.name])
      }
    }
    
    if(any(grep("^num", names(data)))){
      # calculate ratio of distal to proximal pA read numbers
      for(string in sub("_distal", "", grep("(num.*)_distal", names(cmp), value = T, perl = T))){
        
        distal.col.name = paste0(string, "_distal")
        proximal.col.name = paste0(string, "_proximal")
        cmp[, paste0("ratio_", sub("^num_", "", string))] = log2(cmp[, distal.col.name]) - log2(cmp[, proximal.col.name])
      }
      #cmp = cmp[, -grep("^num", names(cmp))]
      #names(cmp) = sub("ratio_", "", names(cmp))
    }
    
    if(any(grep("count$", names(data)))){
      # calculate ratio of distal to proximal pA read numbers
      for(string in sub("_distal", "", grep("(count.*)_distal", names(cmp), value = T, perl = T))){
        
        distal.col.name = paste0(string, "_distal")
        proximal.col.name = paste0(string, "_proximal")
        cmp[, sub("counts?", "d2p_ratio", string)] = log2(cmp[, distal.col.name]) - log2(cmp[, proximal.col.name])
      }
    }
    
    if(any(grep("^HL\\.", names(data)))){
      # calculate ratio of distal to proximal pA halflifes
      for(string in sub("_distal", "", grep("(^HL\\..*)_distal", names(cmp), value = T, perl = T))){
        
        distal.col.name = paste0(string, "_distal")
        proximal.col.name = paste0(string, "_proximal")
        cmp[, paste0(string, "_ratio")] = log2(cmp[, distal.col.name]) - log2(cmp[, proximal.col.name])
        cmp[, paste0(string, "_difference")] = cmp[, distal.col.name] - cmp[, proximal.col.name]
      }
    }
  }
  cmp
}

get.pA.pairs = function(data, cols = "UTR_length", neighbor = F, toptwo = F, match_only=F){ 
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
  names(cmp) = sub("\\.x$", "_distal", names(cmp))
  names(cmp) = sub("\\.y$", "_proximal", names(cmp))
  names(cmp) = sub("comparedPAi", "proximal_pA", names(cmp))
  names(cmp) = sub("comparedPAj", "distal_pA", names(cmp))
  
  # calculate distal to proximal ratios within each sample
  if(match_only == F){
    if(any(grep("^RAI", names(data)))){
      # calculate delta RAIs
      for(string in sub("_distal", "", grep("(RAI.*)_distal", names(cmp), value = T, perl = T))){
        distal.col.name = paste0(string, "_distal")
        proximal.col.name = paste0(string, "_proximal")
        cmp[, paste0("delt_", string)] = cmp[, distal.col.name] - cmp[, proximal.col.name]
      }
    }
    
    if(any(grep("^RPM", names(data)))){
      # calculate ratio of distal to proximal pA RPMs
      for(string in sub("_distal", "", grep("(RPM.*)_distal", names(cmp), value = T, perl = T))){
        
        distal.col.name = paste0(string, "_distal")
        proximal.col.name = paste0(string, "_proximal")
        cmp[, paste0("ratio_", sub("^RPM_", "", string))] = log2(cmp[, distal.col.name]) - log2(cmp[, proximal.col.name])
      }
    }
    
    if(any(grep("^num", names(data)))){
      # calculate ratio of distal to proximal pA read numbers
      for(string in sub("_distal", "", grep("(num.*)_distal", names(cmp), value = T, perl = T))){
        
        distal.col.name = paste0(string, "_distal")
        proximal.col.name = paste0(string, "_proximal")
        cmp[, paste0("ratio_", sub("^num_", "", string))] = log2(cmp[, distal.col.name]) - log2(cmp[, proximal.col.name])
      }
      #cmp = cmp[, -grep("^num", names(cmp))]
      #names(cmp) = sub("ratio_", "", names(cmp))
    }
    
    if(any(grep("count$", names(data)))){
      # calculate ratio of distal to proximal pA read numbers
      for(string in sub("_distal", "", grep("(count.*)_distal", names(cmp), value = T, perl = T))){
        
        distal.col.name = paste0(string, "_distal")
        proximal.col.name = paste0(string, "_proximal")
        cmp[, sub("counts?", "d2p_ratio", string)] = log2(cmp[, distal.col.name]) - log2(cmp[, proximal.col.name])
      }
    }
    
    if(any(grep("^HL\\.", names(data)))){
      # calculate ratio of distal to proximal pA halflifes
      for(string in sub("_distal", "", grep("(^HL\\..*)_distal", names(cmp), value = T, perl = T))){
        
        distal.col.name = paste0(string, "_distal")
        proximal.col.name = paste0(string, "_proximal")
        cmp[, paste0(string, "_ratio")] = log2(cmp[, distal.col.name]) - log2(cmp[, proximal.col.name])
        cmp[, paste0(string, "_difference")] = cmp[, distal.col.name] - cmp[, proximal.col.name]
      }
    }
  }
  
  # return
  cmp = merge(cmp, unique(data[, c("gene_symbol", "description")]), all.x = T, sort=F)
}

get_red = function(data, cols = "UTR_length", neighbor = F, toptwo = F, match_only=F){ 
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
  names(cmp) = sub("\\.x$", "_distal", names(cmp))
  names(cmp) = sub("\\.y$", "_proximal", names(cmp))
  names(cmp) = sub("comparedPAi", "proximal_pA", names(cmp))
  names(cmp) = sub("comparedPAj", "distal_pA", names(cmp))
  
  # calculate distal to proximal ratios within each sample
  if(match_only == F){
    if(any(grep("^RPM", names(data)))){
      # calculate ratio of distal to proximal pA RPMs
      for(string in sub("_distal", "", grep("(RPM.*)_distal", names(cmp), value = T, perl = T))){
        distal.col.name = paste0(string, "_distal")
        proximal.col.name = paste0(string, "_proximal")
        cmp[, paste0("ratio_", sub("^RPM_", "", string))] = log2(cmp[, distal.col.name]) - log2(cmp[, proximal.col.name])
      }
    }
    if(any(grep("count$", names(data)))){
      # calculate ratio of distal to proximal pA read numbers
      for(string in sub("_distal", "", grep("(count.*)_distal", names(cmp), value = T, perl = T))){
        distal.col.name = paste0(string, "_distal")
        proximal.col.name = paste0(string, "_proximal")
        cmp[, sub("counts?", "d2p_ratio", string)] = log2(cmp[, distal.col.name]) - log2(cmp[, proximal.col.name])
      }
    }
  }
  
  # return
  cmp = merge(cmp, unique(data[, c("gene_symbol", "description")]), all.x = T, sort=F)
}


#### Convert pAid into a dataframe containing chr, pA_pos, and strand info
pAid2pos = function(pairwise.pA.df){
  pAid.i = pairwise.pA.df$proximal_pA
  pAid.j = pairwise.pA.df$distal_pA
  pA.pos.df = as.data.frame(do.call(rbind, strsplit(pAid.i, "[+-]")), 
                            stringsAsFactors = F)
  #pos.i.string = "proximal_pA_pos"
  names(pA.pos.df) = c("chr", "proximal_pA_pos")
  pA.pos.df$chr = as.character(pA.pos.df$chr)
  pA.pos.df$proximal_pA_pos = as.integer(pA.pos.df$proximal_pA_pos)
  pA.pos.df$strand = as.vector(sub(".*?([+-])\\d*", "\\1", pAid.i, perl = T))
  
  pA.pos.df$distal_pA_pos = as.data.frame(do.call(rbind, strsplit(pAid.j, "[+-]")), 
                                          stringsAsFactors = F)[,2]
  pA.pos.df$distal_pA_pos = as.integer(pA.pos.df$distal_pA_pos)
  
  cbind(pA.pos.df, pairwise.pA.df)
}

#### If the input is just pA.df (with only one pAid column)
split.pAid.in.pA.df = function(pA.df){
  pA.pos.df = as.data.frame(do.call(rbind, strsplit(pA.df$pAid, "[+-]")), 
                            stringsAsFactors = F)
  names(pA.pos.df) = c("chr", "pA_pos")
  pA.pos.df$chr = as.character(pA.pos.df$chr)
  pA.pos.df$pA_pos = as.integer(pA.pos.df$pA_pos)
  pA.pos.df$strand = as.vector(sub(".*?([+-])\\d*", "\\1", pA.df$pAid, perl = T))
  cbind(pA.pos.df, pA.df)
}

#### Match aUTRs with slightly different pA sites in two data frames allowing distance of "radius" nt
match.pA.pair = function(df1, df2, radius = 12, keep = "intersect", rm_dup_col = T){
  # "keep" can be one of "union", "intersect", "df1", "df2")
  # "rm_dup_col": should duplicated columns be removed?
  
  if(!keep %in% c("union", "intersect", "df1", "df2")){
    stop('"keep" can only be one of "union", "intersect", "df1", and "df2"')
  }
  
  # calculate chr, strand, and position
  if("proximal_pA" %in% names(df1) & "distal_pA" %in% names(df1)){
    df1 = pAid2pos(df1)
  }
  if("proximal_pA" %in% names(df2) & "distal_pA" %in% names(df2)){
    df2 = pAid2pos(df2)
  }
  
  # calculate pA position in case it is not provided (for backward compatibility, when pAid is not supplied)
  if(!"proximal_pA_pos" %in% names(df1)){
    df1$proximal_pA_pos = as.numeric(sub("chr.+[-+]", "", df1$proximal_pA))
  }
  if(!"distal_pA_pos" %in% names(df1)){
    df1$distal_pA_pos = as.numeric(sub("chr.+[-+]", "", df1$distal_pA))
  }
  if(!"proximal_pA_pos" %in% names(df2)){
    df2$proximal_pA_pos = as.numeric(sub("chr.+[-+]", "", df2$proximal_pA))
  }
  if(!"distal_pA_pos" %in% names(df2)){
    df2$distal_pA_pos = as.numeric(sub("chr.+[-+]", "", df2$distal_pA))
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
        if(abs(df1_intersect.lst[[i]][j,"proximal_pA_pos"] - df2_intersect.lst[[i]][k, "proximal_pA_pos"]) <= radius*2 &
           abs(df1_intersect.lst[[i]][j,"distal_pA_pos"] - df2_intersect.lst[[i]][k,"distal_pA_pos"]) <= radius*2){
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
  df1 = split.pAid.in.pA.df(as.data.frame(df1))
  df2$chr = NULL
  df2$strand = NULL
  df2$pA_pos = NULL
  df2 = split.pAid.in.pA.df(as.data.frame(df2))
  
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


#### An updated wrapper to do GO analysis using hypergeometric test
library(GOstats)
hg.go = function(entrezUniverse, 
                 selectedEntrezIds, 
                 annotation = "org.Mm.eg.db", 
                 adjust.p = F,
                 hgCutoff = 0.01, 
                 conditional = F, 
                 categorySize = 5,
                 simplify = T){
  params.bp <- new("GOHyperGParams",
                   geneIds = unique(as.character(selectedEntrezIds)),
                   universeGeneIds = unique(as.character(entrezUniverse)),
                   annotation = annotation, 
                   ontology="BP",
                   pvalueCutoff=hgCutoff,
                   conditional=conditional,        
                   testDirection="over")
  params.cc <- params.mf <- params.bp
  ontology(params.cc) = "CC"
  # ontology(params.mf) = "MF"
  
  ## BP
  hgOver.bp <- hyperGTest(params.bp)
  # How p values are calculated:
  # length(geneIds(hgOver.bp)) # 67, number of genes in the cluster
  # length(geneIdUniverse(hgOver.bp)[["GO:0042981"]]) # 562, number of background genes associated with the GO 
  # length(unique(unlist(geneIdUniverse(hgOver.bp)))) #6532, number of background genes  
  # p value for seeing 12 or more genes associated with GO:0042981 in the cluster of 67 genes: 1-phyper(11, 562, 6532-562, 67)
  df1 = summary(hgOver.bp, categorySize=categorySize) # categorySize is the minimum category size
  if(nrow(df1) == 0){
    df1 = summary(hgOver.bp, categorySize=2)
    if(nrow(df1) == 0){
      pvalueCutoff(params.bp) = 1
      hgOver.bp <- hyperGTest(params.bp)
      df1 = summary(hgOver.bp, categorySize=5)
    }
  }
  names(df1)[1] = "GOID"
  df1$Ontology = "BP"
  # Add gene symbols
  df1$GenesFound = ""
  geneIdList = geneIdsByCategory(hgOver.bp)[sigCategories(hgOver.bp)]
  if(annotation == "org.Mm.eg.db"){
    library(org.Mm.eg.db)
    for(i in 1:nrow(df1)){
      df1$GenesFound[i] = paste(AnnotationDbi::select(org.Mm.eg.db, keys = as.character(geneIdList[[df1$GOID[i]]]), columns = "SYMBOL")$SYMBOL, collapse = ", ")
    }
  }else if(annotation == "org.Hs.eg.db"){
    library(org.Hs.eg.db)
    for(i in 1:nrow(df1)){
      df1$GenesFound[i] = paste(AnnotationDbi::select(org.Hs.eg.db, keys = as.character(geneIdList[[df1$GOID[i]]]), columns = "SYMBOL")$SYMBOL, collapse = ", ")
    }
  }
  
  # Merge GO terms with significant overlaps
  if(nrow(df1) > 1 & simplify){
    # Step 1. Remove GO terms associated with > 1000 genes in the universe
    df1 = subset(df1, Size <= 1000 & Count > 1)
    
    # Step 2. Sort GO terms by p values (descending) 
    df1$filter = F
    df1 = df1[order(df1$Pvalue, decreasing = T),]
    
    # Step 3. From top to bottom, compare each GO (i) with each (j) of the GOs below it. If >= 75% genes in i overlap with genes in j,
    #         remove i. Otherwise, do nothing
    for(i in 1:(nrow(df1)-1)){
      for(j in (i+1):nrow(df1)){
        if(df1$filter[i]){next}
        overlap_i = length(intersect(geneIdList[[df1$GOID[i]]], geneIdList[[df1$GOID[j]]]))/length(geneIdList[[df1$GOID[i]]])
        if(overlap_i > 0.75){
          df1$filter[i] = T
        }
      }
    }
    df1 = df1[!df1$filter, ]
    df1$filter = NULL
  }
  df1 = df1[order(df1$Pvalue),]
  
  ## CC
  hgOver.cc <- hyperGTest(params.cc)
  df2 = summary(hgOver.cc, categorySize=categorySize)
  if(nrow(df2) == 0){
    df2 = summary(hgOver.cc, categorySize=2)
    if(nrow(df2) == 0){
      pvalueCutoff(params.cc) = 1
      hgOver.cc <- hyperGTest(params.cc)
      df2 = summary(hgOver.cc, categorySize=5)
    }
  }
  names(df2)[1] = "GOID"
  df2$Ontology = "CC"
  # Add gene symbols
  df2$GenesFound = ""
  geneIdList = geneIdsByCategory(hgOver.cc)[sigCategories(hgOver.cc)]
  if(annotation == "org.Mm.eg.db"){
    library(org.Mm.eg.db)
    for(i in 1:nrow(df2)){
      df2$GenesFound[i] = paste(AnnotationDbi::select(org.Mm.eg.db, keys = as.character(geneIdList[[df2$GOID[i]]]), columns = "SYMBOL")$SYMBOL, collapse = ", ")
    }
  }else if(annotation == "org.Hs.eg.db"){
    library(org.Hs.eg.db)
    for(i in 1:nrow(df2)){
      df2$GenesFound[i] = paste(AnnotationDbi::select(org.Hs.eg.db, keys = as.character(geneIdList[[df2$GOID[i]]]), columns = "SYMBOL")$SYMBOL, collapse = ", ")
    }
  }
  
  # Merge GO terms with significant overlaps
  if(nrow(df2) > 1 & simplify){
    # Step 1. Remove GO terms associated with > 1000 genes in the universe
    df2 = subset(df2, Size <= 1000 & Count > 1)
    # Step 2. Sort GO terms by p values (descending) 
    df2$filter = F
    df2 = df2[order(df2$Pvalue, decreasing = T),]
    
    # Step 3. From top to bottom, compare each GO (i) with each (j) of the GOs below it. If >= 75% genes in i overlap with genes in j,
    #         remove i. Otherwise, do nothing
    for(i in 1:(nrow(df2)-1)){
      for(j in (i+1):nrow(df2)){
        if(df2$filter[i]){next}
        overlap_i = length(intersect(geneIdList[[df2$GOID[i]]], geneIdList[[df2$GOID[j]]]))/length(geneIdList[[df2$GOID[i]]])
        if(overlap_i > 0.75){
          df2$filter[i] = T
        }
      }
    }
    df2 = df2[!df2$filter, ]
    df2$filter = NULL
  }
  df2 = df2[order(df2$Pvalue),]
  
  # Adjust p-values for multiple-testing?
  if(adjust.p){
    df1$Pvalue = p.adjust(df1$Pvalue)
    df2$Pvalue = p.adjust(df2$Pvalue)
  }
  
  df = rbind(df1, df2)
  
  # Reorder columns and return
  df[, c("Ontology", "GOID", "Pvalue", "OddsRatio",  "ExpCount", "Count", "Size", "GenesFound", "Term")]
  
}
# # How to use:
# go_result = hg.go.2(selectedEntrezIds = df$selected_gene_id),
#                    entrezUniverse = df$all_gene_id,
#                    annotation ="org.Hs.eg.db",
#                    adjust.p = F)


### a helper function to get genomic sequences upstream of the cleavage sites
get.pA.upstream = function(pA.df = pA.df, from = -149, to = -50){
  #   required column names: "chromosome", "strand", "pA_pos", "gene_symbol" 
  
  # the GRanges for the cleavage sites, used to determine if they are within introns or 3'UTR extentions
  require(GenomicRanges)
  pA = GRanges(seqnames = pA.df$chromosome, 
               ranges = IRanges(start = pA.df$pA_pos, 
                                end = pA.df$pA_pos, 
                                names = pA.df$gene_symbol),
               strand = pA.df$strand)
  # since the pA site might be a few nt outside of exons defined by refseq, expand pA by +/-12 nt
  pA = pA + 12
  
  require(GenomicFeatures)
  # get intron and exon definitions
  txdb = loadDb(file.path("~/projects/fud","mm9.refGene.txdb.sqlite"))
  introns = intronsByTranscript(txdb, use.names=F) #transcriptsBy() returns preRNAs containing introns.
  exons = exonsBy(txdb, by = "tx", use.names=F)
  #   #The last exon is still the last one in the GRanges when strand is "-"
  #   #last.exons = lapply(exons, function(x) x[length(x)]) # Too slow to run.
  #   num.exons = elementLengths(exons)
  #   #last.exons = mapply(function(x,i) x[i], exons, num.exons) # Too slow to run.
  #   last.exons = unlist(exons)
  #   last.exons = last.exons[cumsum(num.exons)] #fast!
  
  # Some pAs only overlap with exons 
  pA[pA %over% exons & !(pA %over% introns)] # case 1
  # Some pAs overlap with both exons and introns. They will be considered in exon. 
  pA[pA %over% exons & pA %over% introns] # case 2
  
  # Some pAs only overlap with introns. # case 3
  pA[!(pA %over% exons) & pA %over% introns] # no need to adjust ranges for the 3'ends of 3'UTRs, assuming the -150 to 0 region does not overlap with a different intron 
  # Some pAs overlap with neither exons nor introns (3' extention) # case 4
  pA[!(pA %over% exons) & !(pA %over% introns)] # no need to adjust ranges for the 3'ends of 3'UTRs, assuming the splicing pattern in the 3'extention region is unknown
  
  ########## In case 1 and 2, the 3'end of the 3'UTR needs to be adjusted
  exon.pA = pA[pA %over% exons] # pAs of both case 1 and 2
  # transcripts containing exons overlapping with pAs
  tx.with.pA.in.exons = subsetByOverlaps(exons, exon.pA) 
  
  # get the exonic genomic sequence upstream of the pA in case 1 and 2
  tx.seq = vector("character", length(exon.pA))
  for (i in 1:length(exon.pA)){
    # the transcript overlaping with the pA
    ol.tx = subsetByOverlaps(tx.with.pA.in.exons, exon.pA[i])[[1]] # >1 transcripts may overlap with one pA
    # the exon overlaping with the pA
    ol.exon = ol.tx[ol.tx %over% exon.pA[i]]
    # discard the exons downstream of the pA
    ol.tx = ol.tx[ol.tx$exon_rank <= ol.exon$exon_rank]
    # only keep the part of overlapping exon upstream of the pA site
    if (attributes(strand(exon.pA[i]))$values == "+"){
      end(ol.tx[ol.tx %over% exon.pA[i]]) = (end(exon.pA[i]) + start(exon.pA[i]))/2
    }else{
      start(ol.tx[ol.tx %over% exon.pA[i]]) = (end(exon.pA[i]) + start(exon.pA[i]))/2
    }
    # get genomic sequence for the exons 5' of the pA 
    require(BSgenome.Mmusculus.UCSC.mm9) 
    exon.seq = getSeq(Mmusculus, ol.tx)
    # combine the exons into one transcript
    exon.seq = unlist(exon.seq)
    # subset the sequence
    tx.seq[i] = substring(as.character(exon.seq), length(exon.seq) + from, length(exon.seq) + to)
  }
  exon.pA$seq = tx.seq # some returned sequences will be very short, because the pA is too close to the TSS
  
  ######### In case 3 and 4, the 3'end of the 3'UTR does not needs to be adjusted
  other.pA = pA[!(pA %over% exons)] # pA overlapping with introns or 3' extensions 
  # the GRanges for the 3'end of 3'UTRs
  gr = GRanges(seqnames = pA.df$chromosome, 
               # the start and end positions for the searched region
               ranges = IRanges(start = ifelse(pA.df$strand == "+", pA.df$pA_pos + from, pA.df$pA_pos - to), 
                                end = ifelse(pA.df$strand == "+", pA.df$pA_pos + to, pA.df$pA_pos - from), 
                                names = pA.df$gene_symbol),
               strand = pA.df$strand)
  # only care about the 3'UTRs in case 3 and 4  
  gr = gr[!(pA %over% exons)]
  # get genomic sequence 
  require(BSgenome.Mmusculus.UCSC.mm9) 
  gr.seq = getSeq(Mmusculus, gr)
  other.pA$seq = as.vector(gr.seq)
  
  # shrink the pA position to its original value
  exon.pA = exon.pA - 12
  other.pA = other.pA - 12
  
  # combine the two GRanges
  pA = c(exon.pA, other.pA)
  
  # convert to dataframe
  pA = data.frame(pA_pos = start(pA), seq = pA$seq)
  # merge
  pA.df = merge(pA.df, pA, by = "pA_pos", all = T, sort = F)
  
  # return the dataframe
  pA.df
}

### a helper function to test if the sequence is A-rich
is.Arich = function(window_seq, max = 6, window.size = 10){
  window_seq = toupper(window_seq)
  # slide a window of width window.size
  for (i in 1:(nchar(window_seq) - window.size + 1)){
    subseq = substr(window_seq, i, i+window.size-1)
    # count the number of "A"s in the window
    num_A = 0
    for (j in 1:window.size){
      if (substr(subseq, j, j) == "A"){
        num_A = num_A +1
      }
    }
    # if bigger than threshold, return True
    if (num_A >= max){
      return(T)
    }
  }
  # if reached here, the sequence is not A-rich
  return(F)
}

### a helper function to calculate GC%
GC.percent = function(window_seq){
  # count numbers of G, C ## could have tried oligonucleotideFrequency()
  window_seq = toupper(window_seq)
  # count numbers of G, C
  num_GC = 0
  for (i in 1:nchar(window_seq)){
    if (substr(window_seq, i, i) %in% c("C", "G")){
      num_GC = num_GC +1
    }
  }
  # calculate GC%
  num_GC/nchar(window_seq)
}

### a helper function to check if there are 1~3 G/C in the last 5 nucleotide
GC.clamp = function(window_seq){
  # get the last 5 nt
  window_seq = substr(window_seq, nchar(window_seq)-4, nchar(window_seq))
  window_seq = toupper(window_seq)
  # count numbers of G, C
  num_GC = 0
  for (i in 1:nchar(window_seq)){
    if (substr(window_seq, i, i) %in% c("C", "G")){
      num_GC = num_GC +1
    }
  }
  # check num_GC is in the right range
  if (num_GC >= 1 & num_GC <= 3){
    return(T)
  }else{
    return(F)
  }
}

# in silico PCR tp check if the primers are specific
iPCR = function(fprimers=pA.df$forward_primer,rprimers=pA.df$reverse_primer, genome = "mm9",max.size=4000){
  # trim the primers to allow for mismatches at the 5' end
  fprimers = substr(fprimers, start=6, stop=nchar(fprimers))
  rprimers = substr(rprimers, start=6, stop=nchar(rprimers))
  # vector to save result
  uniqueProduct = rep(1, length(fprimers))
  
  if (genome == "mm9"){
    require(BSgenome.Mmusculus.UCSC.mm9)
    require(Biostrings)
    require(GenomicRanges)
    
    for(i in 1:length(fprimers)){
      f.gr = vmatchPattern(fprimers[i], Mmusculus)
      r.gr = vmatchPattern(rprimers[i], Mmusculus)
      strand(r.gr) = ifelse(strand(r.gr) == "+","-","+") #flip strand
      for(j in 1:length(r.gr)){
        if(as.vector(strand(r.gr)[j]) == "+"){ #if there is any PCR product, the product is on the + strand
          start(r.gr)[j] = start(r.gr)[j] - max.size #expand r.gr to upstream region
        }
        else{ #if there is any PCR product, the product is on the - strand
          end(r.gr)[j] = end(r.gr)[j] + max.size #expand r.gr to downstream region
        }
      }
      
      if(sum(countOverlaps(f.gr, r.gr)) > 1){
        uniqueProduct[i] = 0
      }
    }
  }
  uniqueProduct
}


#### Improved version of the ggcorr function in GGally
ggcorr2 <- function(data, method = "pairwise", palette = "Reds", label.size = 2) {
  require(ggplot2)
  require(reshape2)
  M <- cor(data[1:ncol(data)], use = method)
  M <- M * lower.tri(M)
  M <- as.data.frame(M)
  M <- data.frame(row = names(data), M)
  M <- melt(M)
  M$value[M$value == 0] <- NA
  s <- seq(-1, 1, by = 0.1)
  M$value <- cut(M$value, breaks = s, include.lowest = TRUE, label = cut(s, breaks = s)[-1])
  M$row <- factor(M$row, levels = (unique(as.character(M$variable))))
  diag <- subset(M, row == variable)
  
  po.nopanel <- list(theme(panel.background = element_blank(), panel.grid.minor = element_blank(), 
                           panel.grid.major = element_blank(), axis.text.x = element_text(angle = -90)))
  
  ggplot(M, aes(row, variable)) + scale_fill_brewer(palette = palette, name = "Correlation\ncoefficient") + 
    geom_tile(aes(fill = value), colour = "white") + geom_text(data = diag, 
                                                               aes(label = variable, hjust = 0.5), size = label.size) + scale_x_discrete(breaks = NULL) + 
    scale_y_discrete(breaks = NULL) + labs(x = "", y = "") + po.nopanel
}

#### Count number of motifs upstream of a pA 
countMotif_mm9 = function(pA.df, search_from = 0, search_len = 50, motif = motif){
  # pA.df must contain "chr", "pA_pos", and "strand" columns
  require(BSgenome.Mmusculus.UCSC.mm9)
  require(GenomicRanges)
  require(Biostrings)
  
  upseq.gr = GRanges(seqnames = pA.df$chr, 
                     ranges = IRanges(start = pA.df$pA_pos, end = pA.df$pA_pos),
                     strand = pA.df$strand)
  start(upseq.gr[strand(upseq.gr) == "+", ]) = start(upseq.gr[strand(upseq.gr) == "+", ]) - search_from - search_len + 1
  end(upseq.gr[strand(upseq.gr) == "+", ]) = end(upseq.gr[strand(upseq.gr) == "+", ]) - search_from + 1
  end(upseq.gr[strand(upseq.gr) == "-", ]) = end(upseq.gr[strand(upseq.gr) == "-", ]) + search_from + search_len - 1
  start(upseq.gr[strand(upseq.gr) == "-", ]) = start(upseq.gr[strand(upseq.gr) == "-", ]) + search_from - 1
  
  upseq = getSeq(Mmusculus, upseq.gr)
  
  motif_counts = oligonucleotideFrequency(upseq, width = nchar(motif), step=1,
                                          as.prob=FALSE, as.array=F,
                                          fast.moving.side="right", with.labels=TRUE,
                                          simplify.as="matrix")
  
  motif_counts[, motif]
}


countAllMotif = function(pA.df, geno = "mm9", search_from = 0, search_len = 50, motif_width = 4){
  # pA.df must contain "chr", "pA_pos", and "strand" columns
  if(grepl("^mm", geno)){
    require(BSgenome.Mmusculus.UCSC.mm9)
  }else if(grepl("^hg", geno)){
    require(BSgenome.Hsapiens.UCSC.hg19)
  }
  
  require(GenomicRanges)
  require(Biostrings)
  
  upseq.gr = GRanges(seqnames = pA.df$chr, 
                     ranges = IRanges(start = pA.df$pA_pos, end = pA.df$pA_pos),
                     strand = pA.df$strand)
  start(upseq.gr[strand(upseq.gr) == "+", ]) = start(upseq.gr[strand(upseq.gr) == "+", ]) - search_from - search_len + 1
  end(upseq.gr[strand(upseq.gr) == "+", ]) = end(upseq.gr[strand(upseq.gr) == "+", ]) - search_from + 1
  end(upseq.gr[strand(upseq.gr) == "-", ]) = end(upseq.gr[strand(upseq.gr) == "-", ]) + search_from + search_len - 1
  start(upseq.gr[strand(upseq.gr) == "-", ]) = start(upseq.gr[strand(upseq.gr) == "-", ]) + search_from - 1
  
  
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


## Get and write aUTR sequence to a file
aUTR2fasta = function(tmp=pair.pA, geno="mm9", key="AS_CN", direction = "<", threshold = -0.5, output_dir = result_dir){
  if(geno=="mm9"){
    require(BSgenome.Mmusculus.UCSC.mm9) 
    require(GenomicRanges)
    require(Biostrings)
    direction = paste0(" ", direction, " ") # otherwise "AS_CN<-0.5" will look like assigning
    tmp = subset(tmp, eval(parse(text=paste0(key, direction, threshold))))
    tmp = tmp[order(tmp[, key]),]
      
    gr = GRanges(seqnames = tmp$chr, 
                 ranges = IRanges(start = ifelse(tmp$strand == "+", tmp$proximal_pA_pos, tmp$distal_pA_pos), 
                                  end = ifelse(tmp$strand == "+", tmp$distal_pA_pos, tmp$proximal_pA_pos), 
                                  names = paste0(tmp$gene_symbol, tmp$proximal_pA, "|", tmp$distal_pA_pos)),
                 strand = tmp$strand)
    tmp = getSeq(Mmusculus, gr)
    # write FASTA files
    comp = ifelse(grepl(">", direction), "greater_than", "less_than")
    output_file = paste0(key, "_", comp, "_", threshold, ".fasta")
    writeXStringSet(tmp, file=file.path(output_dir, output_file), format="fasta", width=80)
  }
}

# Get mature transcript sequence before the pA and calculate number of each nucleotide
txpA2fasta = function(pA.df, geno = "mm9", save_fasta = F){
  # pA.df should contain the "gene_id", "strand", "chr", and "pA_pos" columns
  # geno should be either "mm9" or "hg19" for now 
  
  # load pakages for the right genome
  if(geno == "mm9"){
    require(BSgenome.Mmusculus.UCSC.mm9) 
    library(GenomicFeatures)
    txdb = loadDb("../../../fud/mm9.201509.refGene.txdb.sqlite")
  }else if(geno == "hg19"){
    require(BSgenome.Hsapiens.UCSC.hg19) 
    library(GenomicFeatures)
    txdb = loadDb("../../../fud/hg19.refGene.txdb.sqlite")
  }else{
    stop("geno should be either mm9 or hg19 for now.")
  }
  exons = exonsBy(txdb, by="gene")
  exons = reduce(exons)
  
  # Make sure that gene_id is character (instead of integer, which will break exons[[pA.df$gene_id[i]]])
  if(!is.character(pA.df$gene_id)){
    pA.df$gene_id = as.character(pA.df$gene_id)
  }
  
  # set up
  pA.df$A_counts = NA
  pA.df$C_counts = NA
  pA.df$G_counts = NA
  pA.df$T_counts = NA
  file_name = paste0("transcript_seq_", gsub(" |:|-", "_", Sys.time()), ".fasta")
  
  # go through each row of pA.df to get mature transcript sequence
  i = 1
  for(i in 1:nrow(pA.df)){
    if(pA.df$gene_id[i] %in% names(exons)){ #names(exons) are Entrez Gene IDs
      gene = exons[[pA.df$gene_id[i]]] 
      # make sure that the chromosome and strand are correct
      if(as.vector(seqnames(gene))[1] == as.vector(pA.df[i, ]$chr) & 
         as.vector(strand(gene))[1] == as.vector(pA.df[i, ]$strand)){
        
        # for genes on the "+" strand
        if(as.vector(strand(gene))[1] == "+"){
          # relative position of pA to the starts and ends of exons
          right_of_starts = sign(pA.df[i, ]$pA_pos - start(gene))
          left_of_ends = sign(end(gene) - pA.df[i, ]$pA_pos)
          # loop through the exons
          for(j in 1:length(gene)){
            if(right_of_starts[j] %in% c(-1, 0) & left_of_ends[j] == 1){# pA is in intron before exon j
              end(gene[j-1]) = pA.df[i, ]$pA_pos
              gene = gene[1:(j-1)]
              break
            }
            if(right_of_starts[j] %in% c(1, 0) & left_of_ends[j] %in% c(1, 0)){# pA is in exon j
              end(gene[j]) = pA.df[i, ]$pA_pos
              gene = gene[1:j]
              break
            }
            if(j == length(gene) & right_of_starts[j] == 1 & left_of_ends[j] %in% c(1, 0)){ # pA is on the right of the last exon
              end(gene[j]) = pA.df[i, ]$pA_pos
            }
          }
        }
        # for genes on the "-" strand
        if(as.vector(strand(gene))[1] == "-"){
          # relative position of pA to the starts and ends of exons
          right_of_starts = sign(pA.df[i, ]$pA_pos - start(gene))
          left_of_ends = sign(end(gene) - pA.df[i, ]$pA_pos)
          # loop through the exons
          for(j in 1:length(gene)){
            if(j == 1 & right_of_starts[j] %in% c(-1,0) & left_of_ends[j] == 1){ # pA is on the left of the last exon
              start(gene[j]) = pA.df[i, ]$pA_pos
            }
            if(right_of_starts[j] %in% c(1, 0) & left_of_ends[j] %in% c(1, 0)){# pA is in exon j
              start(gene[j]) = pA.df[i, ]$pA_pos
              gene = gene[j:length(gene)]
              break
            }
            if(right_of_starts[j] %in% c(-1, 0) & left_of_ends[j] == 1){# pA is in intron before exon j
              start(gene[j]) = pA.df[i, ]$pA_pos
              gene = gene[j:length(gene)]
              break
            }
          }
        }
      }
      
      # get exonic sequence
      tx_seq = unlist(getSeq(Mmusculus, gene))
      # count nucleotide
      counts = oligonucleotideFrequency(tx_seq, width=1)
      pA.df[i,]$T_counts = counts["T"]
      pA.df[i,]$A_counts = counts["A"]
      pA.df[i,]$C_counts = counts["C"]
      pA.df[i,]$G_counts = counts["G"]
      
      
      # wite to fasta files
      if(save_fasta == T){
        line = paste0(">", pA.df[i,]$pAid, "\n", toString(tx_seq))
        write.table(line, file = file.path(result_dir, file_name), append = T, 
                    col.names=F, row.names=F, quote=F)
      }
      
    }
    # print out the progress
    if(!i%%100){cat(round(i/nrow(pA.df), 2), "of transcripts analyzed\n")}
    i = i + 1
  }
  # save the pA.df for reuse, to avoid redo the calculation 
  write.csv(pA.df, file.path(result_dir, "pA.df.csv"), row.names = F)
  rm(txdb)
  rm(exons)
  return(pA.df)
}

#### get mature transcript sequence before the pA and calculate number of each nucleotide to be used in foreach()
txpA2fasta.par = function(df, db_file = "../../../fud/mm9.201509.refGene.txdb.sqlite"){
  # df is a sub data frame of pA.df
  # df should contain the "gene_id", "strand", "chr", and "pA_pos" columns
  # geno should be either "mm9" or "hg19" for now 
  
  ### Example of using this function:
  # library(foreach)
  # library(doParallel)
  # cl=makeCluster(8)
  # registerDoParallel(cl)
  # pA.lst = split(pA.df, cut(1:nrow(pA.df), 8))
  # db_file = "../../../fud/mm9.201509.refGene.txdb.sqlite"
  # pA.df2 = foreach(k=1:length(pA.lst), .packages=c("BSgenome.Mmusculus.UCSC.mm9", "GenomicFeatures"), .combine = rbind) %dopar% {
  #   txpA2fasta.par(df=pA.lst[[k]], db_file)
  # }
  
  txdb = loadDb(db_file)
  exons = exonsBy(txdb, by="gene")
  exons = reduce(exons)
  
  # Make sure that gene_id is character (instead of integer, which will break exons[[df$gene_id[i]]])
  if(!is.character(df$gene_id)){
    df$gene_id = as.character(df$gene_id)
  }
  
  # set up
  df$A_counts = NA
  df$C_counts = NA
  df$G_counts = NA
  df$T_counts = NA
  
  # go through each row of df to get mature transcript sequence
  i = 1
  for(i in 1:nrow(df)){
    if(df$gene_id[i] %in% names(exons)){ #names(exons) are Entrez Gene IDs
      gene = exons[[df$gene_id[i]]]
      # make sure that the chromosome and strand are correct
      if(as.vector(seqnames(gene))[1] == as.vector(df[i, ]$chr) & 
         as.vector(strand(gene))[1] == as.vector(df[i, ]$strand)){
        
        # for genes on the "+" strand
        if(as.vector(strand(gene))[1] == "+"){
          # relative position of pA to the starts and ends of exons
          right_of_starts = sign(df[i, ]$pA_pos - start(gene))
          left_of_ends = sign(end(gene) - df[i, ]$pA_pos)
          # loop through the exons
          for(j in 1:length(gene)){
            if(right_of_starts[j] %in% c(-1, 0) & left_of_ends[j] == 1){# pA is in intron before exon j
              end(gene[j-1]) = df[i, ]$pA_pos
              gene = gene[1:(j-1)]
              break
            }
            if(right_of_starts[j] %in% c(1, 0) & left_of_ends[j] %in% c(1, 0)){# pA is in exon j
              end(gene[j]) = df[i, ]$pA_pos
              gene = gene[1:j]
              break
            }
            if(j == length(gene) & right_of_starts[j] == 1 & left_of_ends[j] %in% c(1, 0)){ # pA is on the right of the last exon
              end(gene[j]) = df[i, ]$pA_pos
            }
          }
        }
        # for genes on the "-" strand
        if(as.vector(strand(gene))[1] == "-"){
          # relative position of pA to the starts and ends of exons
          right_of_starts = sign(df[i, ]$pA_pos - start(gene))
          left_of_ends = sign(end(gene) - df[i, ]$pA_pos)
          # loop through the exons
          for(j in 1:length(gene)){
            if(j == 1 & right_of_starts[j] %in% c(-1,0) & left_of_ends[j] == 1){ # pA is on the left of the last exon
              start(gene[j]) = df[i, ]$pA_pos
            }
            if(right_of_starts[j] %in% c(1, 0) & left_of_ends[j] %in% c(1, 0)){# pA is in exon j
              start(gene[j]) = df[i, ]$pA_pos
              gene = gene[j:length(gene)]
              break
            }
            if(right_of_starts[j] %in% c(-1, 0) & left_of_ends[j] == 1){# pA is in intron before exon j
              start(gene[j]) = df[i, ]$pA_pos
              gene = gene[j:length(gene)]
              break
            }
          }
        }
      }
      
      # get exonic sequence
      tx_seq = unlist(getSeq(Mmusculus, gene))
      # count nucleotide
      counts = oligonucleotideFrequency(tx_seq, width=1)
      df[i,]$T_counts = counts["T"]
      df[i,]$A_counts = counts["A"]
      df[i,]$C_counts = counts["C"]
      df[i,]$G_counts = counts["G"]
    }
  }
  df
}

#### Calculate exonic UTR sequence and length by parallel computing
get_exonic_3UTR = function(df = pA.df, threeUTRs = threeUTRs, geno = "mm9"){
  # df is a sub data frame of pA.df. It should contain the "strand", "chr", "pA_pos", "cds_end", and "pAid" columns
  # threeUTRs is a GRange object containing 3'UTR definition
  # geno should be either "mm9" or "hg19" for now 
  
  ### Example of using this function:
  # df = get_exonic_3UTR(df, threeUTRs, "mm9")
  # Two columns (exonic_3UTR_seq, UTR_length) will be appended to df. The columns may contain NA.
  
  library(foreach)
  library(doParallel)
  cl=makeCluster(8)
  registerDoParallel(cl)
  
  if(geno=="mm9"){
    require(BSgenome.Mmusculus.UCSC.mm9) 
    geno = Mmusculus
  }else if(geno=="hg19"){
    require("BSgenome.Hsapiens.UCSC.hg19")
    geno = Hsapiens
  }
  
  require(GenomicRanges)
  require(Biostrings)
  
  
  pA_3UTR_df = df[df$region == "3UTR", ]
  pA_3UTR_df$exonic_3UTR_seq = NA
  pA_3UTR_df$UTR_length = NA
  
  pA_3UTR = GRanges(seqnames = pA_3UTR_df$chr, 
                    ranges = IRanges(start = pA_3UTR_df$pA_pos, 
                                     end = pA_3UTR_df$pA_pos, 
                                     names = pA_3UTR_df$pAid),
                    strand = pA_3UTR_df$strand,
                    cds_end = pA_3UTR_df$cds_end)
  
  
  olp = findOverlaps(pA_3UTR+3, threeUTRs)
  # # debugging
  # for (i in 1:length(olp)){
  #   if(length(threeUTRs[[subjectHits(olp[i,])]]) > 1) {
  #     print(i)
  #     break
  #   }
  # }
  
  # Parallel computing
  pA_3UTR_df[queryHits(olp), ]$exonic_3UTR_seq = foreach(i=1:length(olp), 
                                                         .packages = "BSgenome",
                                                         .combine = "c") %dopar%{
                                                           # i = 8
                                                           this_pA = pA_3UTR[queryHits(olp[i])]
                                                           this_3UTR = threeUTRs[[subjectHits(olp[i])]]
                                                           
                                                           #if(length(strand(this_3UTR)@values) > 1 || is.na(this_pA$cds_end)){ # skip incorrect 3'UTR definitions
                                                           if(length(strand(this_3UTR)@values) > 1){ # TODO: delete this line
                                                            NA
                                                           }else{
                                                             if(strand(this_pA)@values == "+"){
                                                               # Sometimes the there are multiple CDS ends in the last exon
                                                               # for example, check chr1:193919462-193925703 on mm9 (UCSC genome browser)
                                                               #this_3UTR = this_3UTR[start(this_3UTR) >= this_pA$cds_end] 
                                                               if(length(this_3UTR) < 1){
                                                                 NA
                                                               }else{
                                                                 # check if the overlap is true without moving the pA +/- 24 nt
                                                                 pA_in_3UTR = this_pA %over% this_3UTR 

                                                                 if(end(this_pA) >= min(start(this_3UTR))){
                                                                   start(this_pA) = min(start(this_3UTR))
                                                                 }else{ # to avoid negative width in some cases
                                                                   end(this_pA) = min(start(this_3UTR))
                                                                   start(this_pA) = min(start(this_3UTR))
                                                                 }
                                                                 
                                                                 # if the pA is not exactly in the annotated 3'UTR, modify the 3'UTR annotation
                                                                 if(!pA_in_3UTR){
                                                                   # calculate the index of the 3'UTR exon that needs to be extended
                                                                   exon_index = max(which(end(this_3UTR) - end(this_pA) < 0))
                                                                   # extend the exon to the pA position
                                                                   end(this_3UTR)[exon_index] = end(this_pA)
                                                                 }
                                                                 # update threeUTR 
                                                                 # threeUTRs[[subjectHits(olp[i,])]] = this_3UTR 
                                                                 this_3UTR = GenomicRanges::intersect(this_pA, this_3UTR)
                                                                 as.character(unlist(getSeq(geno, this_3UTR)))
                                                               }
                                                             }else if(strand(this_pA)@values == "-"){
                                                               # Sometimes the there are multiple CDS ends in the last exon: 
                                                               # for example, check chr1:193919462-193925703 on mm9 (UCSC genome browser)
                                                               #this_3UTR = this_3UTR[end(this_3UTR) <= this_pA$cds_end]
                                                               if(length(this_3UTR) < 1){ # TODO: delete this line
                                                                 NA
                                                               }else{
                                                                 # check if the overlap is true without moving the pA +/- 24 nt
                                                                 pA_in_3UTR = this_pA %over% this_3UTR 
                                                                 
                                                                 #end(this_pA) = max(end(this_3UTR))
                                                                 if(start(this_pA) <= max(end(this_3UTR))){
                                                                   end(this_pA) = max(end(this_3UTR))
                                                                 }else{ # to avoid negative width in some cases
                                                                   start(this_pA) = max(end(this_3UTR))
                                                                   end(this_pA) = max(end(this_3UTR))
                                                                 }
                                                                 # if the pA is not exactly in the annotated 3'UTR, modify the 3'UTR annotation
                                                                 if(!pA_in_3UTR){
                                                                   # calculate the index of the 3'UTR exon that needs to be extended
                                                                   exon_index = min(which(start(this_3UTR) - start(this_pA) > 0))
                                                                   # extend the exon to the pA position
                                                                   start(this_3UTR)[exon_index] = start(this_pA)
                                                                 }
                                                                 # update threeUTR 
                                                                 # threeUTRs[[subjectHits(olp[i,])]] = this_3UTR 
                                                                 this_3UTR = GenomicRanges::intersect(this_pA, this_3UTR)
                                                                 # when the 3'UTR has >1 exons, their sequences need to be re-ordered:
                                                                 as.character(unlist(getSeq(geno, this_3UTR)[length(this_3UTR):1, ]))
                                                               }
                                                             }
                                                           }
                                                         }
  
  pA_3UTR_df[queryHits(olp), ]$UTR_length = nchar(pA_3UTR_df[queryHits(olp), ]$exonic_3UTR_seq)
  
  merge(df, pA_3UTR_df[, c("pAid", "exonic_3UTR_seq", "UTR_length")], all.x = T, sort = F)
}
## use the function like so:
#df = get_exonic_3UTR(df = pA.df, threeUTRs, "mm9")


#### A function to make scatter plots between two RPM ratios from two data frames, 
######### after filtering out low reads # pAs
rpm_ratio_scatter_plot_from_different_pA_df = function(df1_file="/home/dinghai/sirius/projects/percy/result/batch05_06_09_13/PM_seperated/pA.df.csv", 
                                                       df2_file="/home/dinghai/sirius/projects/rip/result/batch12/NT_all/pA.df.csv", 
                                                       pattern1="_P_M.*rpm_ratio$", 
                                                       pattern2="IgG.*rpm_ratio$", 
                                                       lowest_counts = c(5, 10, 20, 50, 100),
                                                       UTR_only = T,
                                                       output_prefix = "PM_RIP(IgG)"){
  ### inputs: 
  # df1, df2: data frames that have pA site identifier, read counts, rpm ratio
  # pattern1, pattern2: strings for finding the column names of rpm ratios
  # lowest read counts used for filtering
  
  require(xlsx)
  
  # read pA.df files
  df1 = read.csv(df1_file)
  df2 = read.csv(df2_file)
  
  # confirm that the columns can be found
  stopifnot(length(grep(pattern1, names(df1))) >= 1, length(grep(pattern2, names(df2))) >= 1,
            length(grep("_count$", names(df1))) >= 2, length(grep("_count$", names(df2))) >= 2)
  
  # merge data frame based on pAs
  df1 = df1[, c("chr", "strand", "pA_pos", "UTR_length", grep("_count$", names(df1), value=T),
                grep(pattern1, names(df1), value=T))]
  df2 = df2[, c("chr", "strand", "pA_pos", "UTR_length", grep("_count$", names(df2), value=T),
                grep(pattern2, names(df2), value=T))]
  df = match.pAs(df1, df2)
  # for some plots, only mRNAs are meaningful
  if(UTR_only){ df = subset(df, UTR_length > 0)}
  
  rm(df1)
  rm(df2)
  
  # filter pAs based on read #
  for(n in lowest_counts){
    p.list = list() # for saving plots
    k = 1
    correlations = vector("numeric", length=length(grep(pattern1, names(df)))*length(grep(pattern2, names(df))))
    for(ratio1 in sort(grep(pattern1, names(df), value=T))){
      # for each rpm ratio, calculate the column names that contain read counts
      count_1_1 = paste0(paste(strsplit(ratio1, "_")[[1]][c(1,2,4)], collapse="_"), "_count")
      count_1_2 = paste0(paste(strsplit(ratio1, "_")[[1]][c(1,3,4)], collapse="_"), "_count")
      for(ratio2 in sort(grep(pattern2, names(df), value=T))){
        # for each rpm ratio, calculate the column names that contain read counts
        count_2_1 = paste0(paste(strsplit(ratio2, "_")[[1]][c(1,2,4)], collapse="_"), "_count")
        count_2_2 = paste0(paste(strsplit(ratio2, "_")[[1]][c(1,3,4)], collapse="_"), "_count")
        # filter
        clean.df = df[rowSums(df[, c(count_1_1, count_1_2, count_2_1, count_2_2)] >= n, na.rm = T) == 4, ]
        clean.df = clean.df[abs(clean.df[,ratio1]) < Inf & abs(clean.df[,ratio2]) < Inf,]
        
        # calculate correlation
        c1 = round(cor(clean.df[,ratio1], clean.df[,ratio2]), 2)
        c2 = round(cor(clean.df[,ratio1], clean.df[,ratio2], method = "spearman"), 2)   
        correlations[k] = c2
        #     c1.group = round(sapply(split(clean.df, f=clean.df$len_range), 
        #                             function(x) cor(x[,x_string], x[,y_string])), 2)
        #     c2.group = round(sapply(split(clean.df, f=clean.df$len_range), 
        #                             function(x) cor(x[,x_string], x[,y_string],method = "spearman")), 2)
        
        #c = cor(clean.df$UTR_length, clean.df[, y_string])
        p_val = cor.test(clean.df[,ratio1], clean.df[,ratio2])$p.value
        p.list[[k]] = ggplot(data = clean.df, aes_string(x= ratio2, y = ratio1)) + 
          geom_point(shape = 1, alpha = 0.3, size = 1) + 
          geom_smooth(method="rlm") +
          ggtitle(paste0("# of UTRs: ", nrow(clean.df), " | Read# >= ", n, "\n",
                         "Pearson: ", c1, " | ", "Spearman: ", c2, "\n",
                         "p = ", sprintf("%.2E", p_val)))
        
        k = k + 1
      }
    }
    p.list[["nrow"]] = length(grep(pattern1, names(df)))
    png(file.path(result_dir, paste0(output_prefix, "_", n, ".png")), width=300*length(grep(pattern2, names(df))), height=350*length(grep(pattern1, names(df))))
    do.call(grid.arrange, p.list)
    dev.off()
    
    correlations = matrix(correlations, nrow = length(grep(pattern1, names(df))), byrow = T,
                          dimnames = list(sort(grep(pattern1, names(df), value=T)), sort(grep(pattern2, names(df), value=T))))
    write.xlsx(correlations, file.path(result_dir, paste0(output_prefix, "_correlations.xlsx")), 
               sheetName = paste0("Read# cutoff ", n),
               row.names = T,
               append = T)
  }
}

#http://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter7.html
findPotentialStartsAndStops <- function(sequence){
  library("Biostrings")
  # Define a vector with the sequences of potential start and stop codons
  codons <- c("ATG", "TAA", "TAG", "TGA")
  sequence = toupper(sequence)
  # Find the number of occurrences of each type of potential start or stop codon
  for(i in 1:4){
    codon <- codons[i]
    # Find all occurrences of codon "codon" in sequence "sequence"
    occurrences <- matchPattern(codon, sequence)
    # Find the start positions of all occurrences of "codon" in sequence "sequence"
    codonpositions <- start(occurrences)
    # Find the total number of potential start and stop codons in sequence "sequence"
    numoccurrences <- length(codonpositions)
    if (i == 1){
      # Make a copy of vector "codonpositions" called "positions"
      positions <- codonpositions
      # Make a vector "types" containing "numoccurrences" copies of "codon"
      types <- rep(codon, numoccurrences)
    }else{
      # Add the vector "codonpositions" to the end of vector "positions":
      positions   <- append(positions, codonpositions, after=length(positions))
      # Add the vector "rep(codon, numoccurrences)" to the end of vector "types":
      types       <- append(types, rep(codon, numoccurrences), after=length(types))
    }
  }
  # Sort the vectors "positions" and "types" in order of position along the input sequence:
  indices <- order(positions)
  positions <- positions[indices]
  types <- types[indices]
  # Return a list variable including vectors "positions" and "types":
  mylist <- list(positions,types)
  return(mylist)
}
# s1 <- "aaaatgcagtaacccatgccc"
# findPotentialStartsAndStops(s1)

findORFsinSeq <- function(sequence){
  require(Biostrings)
  # Make vectors "positions" and "types" containing information on the positions of ATGs in the sequence:
  mylist <- findPotentialStartsAndStops(sequence)
  positions <- mylist[[1]]
  types <- mylist[[2]]
  # Make vectors "orfstarts" and "orfstops" to store the predicted start and stop codons of ORFs
  orfstarts <- numeric()
  orfstops <- numeric()
  # Make a vector "orflengths" to store the lengths of the ORFs
  orflengths <- numeric()
  # Print out the positions of ORFs in the sequence:
  # Find the length of vector "positions"
  numpositions <- length(positions)
  # There must be at least one start codon and one stop codon to have an ORF.
  if (numpositions >= 2){
    for (i in 1:(numpositions-1)){
      posi <- positions[i]
      typei <- types[i]
      found <- 0
      while (found == 0){
        for (j in (i+1):numpositions){
          posj  <- positions[j]
          typej <- types[j]
          posdiff <- posj - posi
          posdiffmod3 <- posdiff %% 3
          # Add in the length of the stop codon
          orflength <- posj - posi + 3
          if (typei == "ATG" && (typej == "TAA" || typej == "TAG" || typej == "TGA") && posdiffmod3 == 0){
            # Check if we have already used the stop codon at posj+2 in an ORF
            numorfs <- length(orfstops)
            usedstop <- -1
            if (numorfs > 0){
              for (k in 1:numorfs){
                orfstopk <- orfstops[k]
                if (orfstopk == (posj + 2)) { usedstop <- 1 }
              }
            }
            if (usedstop == -1){
              orfstarts <- append(orfstarts, posi, after=length(orfstarts))
              orfstops <- append(orfstops, posj+2, after=length(orfstops)) # Including the stop codon.
              orflengths <- append(orflengths, orflength, after=length(orflengths))
            }
            found <- 1
            break
          }
          if (j == numpositions) { found <- 1 }
        }
      }
    }
  }
  # Sort the final ORFs by start position:
  indices <- order(orfstarts)
  orfstarts <- orfstarts[indices]
  orfstops <- orfstops[indices]
  # Find the lengths of the ORFs that we have
  orflengths <- numeric()
  numorfs <- length(orfstarts)
  for (i in 1:numorfs)
  {
    orfstart <- orfstarts[i]
    orfstop <- orfstops[i]
    orflength <- orfstop - orfstart + 1
    orflengths <- append(orflengths,orflength,after=length(orflengths))
  }
  mylist <- list(orfstarts, orfstops, orflengths)
  names(mylist) = c("starts", "ends", "lengths")
  return(mylist)
}
s1 <- "aaaatgcagtaacccatgccc"
findORFsinSeq(s1)

motifFisher = function(df = pA.df, keyPattern = "tail_length", minReadCount = 10, 
                       sequenceRegion = "3UTR", widths = c(4,6), geno = "mm9"){
  # df: either pA.df or pair.pA
  # keyPattern: the pattern used to select columns for sorting sequences
  # sequenceRegion (from, to): either 3UTR(cds_end, pA_pos), aUTR(proximal_pA_pos, distal_pA_pos)
  
  
  ### define helper functions
  fisher_p = function(thisrow, sums){
    m = rbind(thisrow, sums - thisrow)
    sign = ifelse(m[1,1]/m[2,1] > m[1,2]/m[2,2], -1, 1) # postive: more enriched in upper
    fisher.test(m)$p.value * sign
  }
  oligofisher = function(df=df, key = "NT_IgG_C_1", minReadCount = 10, from = from, to = to, 
                         width=4, geno = "mm9", frac = 0.2, max.padj = 0.1, file_name = file_name){
    require(GenomicRanges)
    require(Biostrings)
    require(xlsx)
    # remove NA and +/-Inf
    df = as.data.frame(df[abs(df[, key]) != Inf & !is.na(df[, key]), ])
    ## filter low values with read counts
    colPatterns = unique(unlist(strsplit(sub(keyPattern, "", key), "_")))
    countColumns = grep(paste0(colPatterns, ".*_count", collapse = "|"), names(df))
    df = df[rowSums(df[, countColumns] >= minReadCount, na.rm = T) == length(countColumns), ]
    
    gr = GRanges(seqnames = df$chr, 
                 ranges = IRanges(start = pmin(df[,from], df[,to]), 
                                  end = pmax(df[,from], df[,to])),
                 strand = df$strand)
    
    if(geno == "mm9"){
      require(BSgenome.Mmusculus.UCSC.mm9) 
      gr = getSeq(Mmusculus, gr)
    }else if(geno == "hg19"){
      require("BSgenome.Hsapiens.UCSC.hg19") 
      gr = getSeq(Hsapiens, gr)
    }
    
    counts = oligonucleotideFrequency(gr, width = width, step = 1)
    #df_with_oligo_counts = cbind(df,counts)
    
    lower_index = df[, key] < quantile(df[, key], probs = seq(0,1, frac))[2]
    upper_index = df[, key] > quantile(df[, key], probs = seq(0,1, frac))[5]
    
    counts = t(rbind(colSums(counts[lower_index,]), colSums(counts[upper_index,])))
    colnames(counts) = c("lower_counts", "upper_counts")
    
    sums = colSums(counts)
    
    fisher = apply(counts, 1, fisher_p, sums) # postive: more enriched in upper
    
    counts = as.data.frame(cbind(counts, fisher))
    counts$padj = p.adjust(abs(counts$fisher), "BH")
    #counts = subset(counts, padj < max.padj)
    
    #hist(counts$fisher)
    # postive: more enriched in upper
    counts$enriched = ifelse(counts$padj > max.padj, "Middle", ifelse(counts$fisher > 0, "Top", "Bottom")) 
    
    counts = counts[order(counts$enriched, counts$padj),]
    
    rownames(counts) = gsub("T", "U", rownames(counts))
    
    counts$sequence = rownames(counts)
    
    rownames(counts) = NULL
    
    #filename = paste0("motif_fisher_", width, "_", key, "_", frac, ".csv")
    #write.csv(subset(counts, padj < 0.05), file.path(result_dir, filename), row.names = F)
    counts2 = subset(counts, padj < max.padj)
    counts2 = counts2[, c("sequence", "fisher")]
    # creat dfs of the same number of rows
    if(sum(counts2$fisher > 0) < sum(counts2$fisher < 0)){
      left = rbind(subset(counts2, fisher > 0), cbind(sequence=rep("", sum(counts2$fisher < 0) - sum(counts2$fisher > 0)), 
                                                      fisher=rep("", sum(counts2$fisher < 0) - sum(counts2$fisher > 0))))
      right = subset(counts2, fisher < 0)
    }else{
      left = subset(counts2, fisher > 0)
      right = rbind(subset(counts2, fisher < 0), cbind(sequence=rep("", sum(counts2$fisher > 0) - sum(counts2$fisher < 0)), 
                                                       fisher=rep("", sum(counts2$fisher > 0) - sum(counts2$fisher < 0))))
    }
    counts2 = cbind(left, right)
    names(counts2) = c("Sequence_Top", "Fisher_Top", "Sequence_Bottom", "Fisher_Bottom")
    ####### as.numeric() is necessary because counts2$Fisher_Top is character
    counts2$Fisher_Top = format(as.numeric(counts2$Fisher_Top), digits= 2, scientific=T)
    counts2$Fisher_Bottom = format(as.numeric(counts2$Fisher_Bottom), digits= 2, scientific=T)
    ########### remove NA
    counts2 = apply(counts2, 2, function(col) sub("NA", "", col)) 
    write.xlsx2(counts2, file.path(result_dir, file_name), 
                sheetName = paste0(key, "_", width, "mer"),
                row.names = F,
                append = T)
    
    counts 
  }
  
  ## calculate from, to
  if(sequenceRegion == "3UTR"){
    from = "cds_end"
    to = "pA_pos"
    df = subset(df, !is.na(cds_end)) # remove ncRNA
  }else if(sequenceRegion == "aUTR"){
    from = "proximal_pA_pos"
    to = "distal_pA_pos"
  }
  
  ## calculate keys
  keys = grep(keyPattern, names(df), value=T)
   
  ## do the calculation
  fisher.ls = list()
  # Sys.time() was used to avoid appending the files with the same name
  file_name = paste0(sequenceRegion, "_", keyPattern, "_minReadCount", minReadCount,  "_Motif_Fisher_", gsub(" |-|:", "_", Sys.time()), ".xlsx")
  for(width in widths){
    for(key in keys){
      fisher.ls[[paste0(key, "_", width,"mer")]] = oligofisher(df=df, key = key, minReadCount = 10, width=width, from=from, to=to,
                                                               geno = "mm9", frac = 0.2, max.padj = 1,
                                                               file_name = file_name)
    }
  }
  
  ## ploting
  if(length(keys) > 1){
    require(dplyr)
    require(tidyr)
    require(ggplot2)
    
    p.list = list()
    k = 1
    for(width in widths){
      for(i in 1:length(keys[-1])){
        for(j in (i+1):length(keys)){
          fisher.ls.name_1 = paste0(keys[i], "_", width, "mer")
          fisher.ls.name_2 = paste0(keys[j], "_", width, "mer")
          #cmp.ls.name = paste0("cmp.", width, "mer")
          df = inner_join(fisher.ls[[fisher.ls.name_1]][, c("sequence", "fisher", "padj")], 
                          fisher.ls[[fisher.ls.name_2]][, c("sequence", "fisher", "padj")], 
                          by = "sequence")
          fisher.x = df$fisher.x
          fisher.y = df$fisher.y
          
          df$fisher.x = -1*sign(fisher.x)*log10(abs(fisher.x))
          df$fisher.y = -1*sign(fisher.y)*log10(abs(fisher.y))
          
          df$significance = c(0.001, 1, 1)[(df$padj.x < 0.05) + (df$padj.y < 0.05) + 1]
          df$color = as.factor((df$padj.x < 0.05) + (df$padj.y < 0.05) + 1)
          
          p.list[[k]] = ggplot(df, aes(x=fisher.x, y=fisher.y))+
            geom_text(aes(label = sequence, x=fisher.x, y=fisher.y, alpha = significance, colour=color), size=3) +
            scale_colour_manual(values = c("black", "orange", "red")) +
            xlab(paste0("SS for ", fisher.ls.name_1))+
            ylab(paste0("SS for ", fisher.ls.name_2))+ #guides(color = FALSE)+
            geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme(legend.position="none")
          k = k + 1
        }
      }
    }
    p.list[["nrow"]] = length(widths)
    bmp(file.path(result_dir, paste0(sequenceRegion, "_", keyPattern, "_minReadCount", minReadCount, "_motif_SS_fisher.bmp")), height=length(widths)*500, width=500*(k-1)/length(widths))
    do.call(grid.arrange, p.list)
    dev.off()
  }
}

inf.omit = function(df){
  index = apply(df, 1, function(row) sum(row == Inf | row == -Inf) == 0)
  df[index,]
}

UTR3_with_oligo_counts = function(df=utr3, from = "cds_end", to = "pA_pos", width=6, 
                                  geno = "mm9", keep_sequence = F){
  require(GenomicRanges)
  require(Biostrings)
  gr = GRanges(seqnames = df$chr, 
               ranges = IRanges(start = pmin(df[,from], df[,to]), 
                                end = pmax(df[,from], df[,to])),
               strand = df$strand)
  
  if(geno == "mm9"){
    require(BSgenome.Mmusculus.UCSC.mm9) 
    gr = getSeq(Mmusculus, gr)
  }else if(geno == "hg19"){
    require("BSgenome.Hsapiens.UCSC.hg19") 
    gr = getSeq(Hsapiens, gr)
  }
  
  counts = oligonucleotideFrequency(gr, width = width, step = 1)
  #df_with_oligo_counts = cbind(df, apply(counts, 2, log10))
  if(keep_sequence){
    sequence = as.character(gr)
    df_with_oligo_counts = cbind(df, counts, sequence)
  }else{
    df_with_oligo_counts = cbind(df, counts)
  }
  df_with_oligo_counts
}


####
fishermain = function(df = tmp_n, keys, sequenceRegions = "3UTR", 
                      widths=c(4,6), frac=0.2, geno="mm9", key_pairs = "all.combinations", filename=filename){
  # tmp_n: df like pA.df (with column pA_pos) or pair.pA 
  # sequenceRegion: "-100_-41", "-40_-1", "1_100", "aUTR", or "3UTR"
  # widths = c(4,6)
  # geno: currently either hg19 or mm9
  # key_pairs: "all.combinations" or vectors like c("1,2", "1,3") to plot the SS of motifs on key1 vs key2 and key1 vs key3 plots
  
  ############ define helper functions
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  fisher_row_sum_p = function(thisrow, sums){
    m = rbind(thisrow, sums - thisrow)
    sign = ifelse(m[1,1]/m[2,1] > m[1,2]/m[2,2], -1, 1) # postive: more enriched in upper
    fisher.test(m)$p.value * sign
  }
  
  oligofisher = function(df=df, key = "tail",  sequenceRegion = "aUTR", width=4, frac = 0.2, geno = "mm9", 
                         filename = filename){
    # sequenceRegion: "-100_-41", "-40_-1", "1_100", "aUTR", or "3UTR"
    # widths = c(4,6)
    # cat("key: ", key, "\n")
    # cat("colnames: ", names(df))
    
    require(GenomicRanges)
    require(Biostrings)
    require(xlsx)
    # remove NA and +/-Inf
    if(key %in% names(df)){
      df = as.data.frame(df[abs(df[, key]) != Inf & !is.na(df[, key]), ])
    }
    # calculate gr
    if(sequenceRegion == "aUTR"){
      from = "proximal_pA_pos" 
      to = "distal_pA_pos" 
      df = pAid2pos(df)
      gr = GRanges(seqnames = df$chr, 
                   ranges = IRanges(start = pmin(df[,from], df[,to]), 
                                    end = pmax(df[,from], df[,to])),
                   strand = df$strand)
    }else if(sequenceRegion == "3UTR"){
      df = subset(df, region=="3UTR")
      from = "cds_end" 
      to = "pA_pos" 
      gr = GRanges(seqnames = df$chr, 
                   ranges = IRanges(start = pmin(df[,from], df[,to]), 
                                    end = pmax(df[,from], df[,to])),
                   strand = df$strand)
    }else if(grepl("-?\\d+_-?\\d+", sequenceRegion)){ # sequenceRegion = "-100_-41"
      from = as.integer(strsplit(sequenceRegion, "_")[[1]][1])
      to = as.integer(strsplit(sequenceRegion, "_")[[1]][2])
      signs = ifelse(df$strand == "+", 1, -1)
      gr = GRanges(seqnames = df$chr, 
                   ranges = IRanges(start = pmin(df$pA_pos + from*signs, df$pA_pos + to*signs), 
                                    end = pmax(df$pA_pos + from*signs, df$pA_pos + to*signs)),
                   strand = df$strand)
    }
    # get sequence
    if(geno == "mm9"){
      require(BSgenome.Mmusculus.UCSC.mm9) 
      gr = getSeq(Mmusculus, gr)
    }else if(geno == "hg19"){
      require("BSgenome.Hsapiens.UCSC.hg19") 
      gr = getSeq(Hsapiens, gr)
    }
    # get counts
    counts = oligonucleotideFrequency(gr, width = width, step = 1)
    
    # limit the impact of a few sequences with large numbers of cis-elements
    ceilings = ceiling(apply(counts, 2, function(c) quantile(c, 0.95)))
    for(j in 1:length(ceilings)){
      counts[counts[, j] > ceilings[j], j] = ceilings[j]
    }
    
    #if(length(lower_index) == 1 && is.na(lower_index) && length(upper_index) == 1 && is.na(upper_index)){
    lower_index = df[, key] < quantile(df[, key], frac) 
    upper_index = df[, key] > quantile(df[, key], 1-frac) 
    #}
    
    counts = t(rbind(colSums(counts[lower_index,]), colSums(counts[upper_index,])))
    colnames(counts) = c("lower_counts", "upper_counts")
    
    sums = colSums(counts)
    
    fisher = apply(counts, 1, fisher_row_sum_p, sums)
    
    counts = as.data.frame(cbind(counts, fisher))
    counts$padj = p.adjust(abs(counts$fisher), "BH")
    
    rownames(counts) = gsub("T", "U", rownames(counts))
    
    counts$sequence = rownames(counts)
    
    rownames(counts) = NULL
    
    counts = counts[, c("sequence", names(counts)[1:4])]
    counts_original = counts
    # creat dfs of the same number of rows
    row_number_difference = sum(counts$fisher > 0) - sum(counts$fisher < 0)
    if(row_number_difference != 0){
      filling = rep("", abs(row_number_difference)) # fill the left or right part of the table, so that the two parts have the same number of rows
      filling = data.frame(sequence = filling,
                           lower_counts = filling,
                           upper_counts = filling,
                           fisher = 1,
                           padj = 1)
      if(row_number_difference < 0){
        left = rbind(subset(counts, fisher > 0), filling)
        right = subset(counts, fisher < 0)
      }else{
        right = rbind(subset(counts, fisher < 0), filling)
        left = subset(counts, fisher > 0)
      }
    }else{
      left = subset(counts, fisher > 0)
      right = subset(counts, fisher < 0)
    }
    
    left = left[order(left$fisher),]
    right = right[order(abs(right$fisher)),]
    counts = cbind(left, right)
    names(counts) = c("Sequence1", "Sequence1_in_Bottom",
                      "Sequence1_in_Top", "Sequence1_p", "Sequence1_padj",
                      "Sequence2", "Sequence2_in_Bottom",
                      "Sequence2_in_Top", "Sequence2_p", "Sequence2_padj")
    
    
    ####### as.numeric() is necessary because counts2$Fisher_Top is character
    counts$Sequence1_p = format(as.numeric(counts$Sequence1_p), digits= 2, scientific=T)
    counts$Sequence1_padj = format(as.numeric(counts$Sequence1_padj), digits= 2, scientific=T)
    counts$Sequence2_p = format(as.numeric(counts$Sequence2_p), digits= 2, scientific=T)
    counts$Sequence2_padj = format(as.numeric(counts$Sequence2_padj), digits= 2, scientific=T)
    ########### remove NA
    counts = apply(counts, 2, function(col) sub("NA", "", col)) 
    sheetName = paste0(sequenceRegion, "_", key, "_", width, "mer")
    write.xlsx2(counts, 
                filename, 
                sheetName = sheetName,
                row.names = F,
                append = T)
    
    return(counts_original) 
  }
  
  ## use the helper functions
  fisher.ls = list()
  for(sequenceRegion in sequenceRegions){#sequenceRegions = c("-100_-41", "-40_0", "1_100")
    for(key in keys){
      for(width in widths){
        fisher.ls[[paste0(sequenceRegion, "_", key, "_", width,"mer")]] = oligofisher(df=df, 
                                                                                      key = key, 
                                                                                      width=width, 
                                                                                      sequenceRegion=sequenceRegion,
                                                                                      geno = geno, #lower_index, upper_index,
                                                                                      frac = 0.2, #max.padj = max.padj,
                                                                                      filename = filename)
      }
    }
    print(names(fisher.ls))
    # ploting
    require(dplyr)
    require(tidyr)
    require(ggplot2)
    require(ggrepel)
    max.padj = 0.05
    
    p.list = list()
    k = 1
    for(width in widths){
      for(i in 1:length(keys[-1])){
        for(j in (i+1):length(keys)){
          if(key_pairs != "all.combinations"){
            if(!paste(i,j,sep=",") %in% key_pairs){
              next
            }
          }
          fisher.ls.name_1 = paste0(sequenceRegion, "_", keys[i], "_", width, "mer")
          fisher.ls.name_2 = paste0(sequenceRegion, "_", keys[j], "_", width, "mer")
          #cmp.ls.name = paste0("cmp.", width, "mer")
          df2 = inner_join(fisher.ls[[fisher.ls.name_1]][, c("sequence", "fisher", "padj")], 
                           fisher.ls[[fisher.ls.name_2]][, c("sequence", "fisher", "padj")], 
                           by = "sequence")
          fisher.x = df2$fisher.x
          fisher.y = df2$fisher.y
          
          df2$fisher.x = -1*sign(fisher.x)*log10(abs(fisher.x))
          df2$fisher.y = -1*sign(fisher.y)*log10(abs(fisher.y))
          
          df2$significance = c(0.001, 0.5, 1)[(df2$padj.x < max.padj) + (df2$padj.y < max.padj) + 1]
          df2$color = factor((df2$padj.x < max.padj) + (df2$padj.y < max.padj) + 1)
          
          p.list[[k]] = ggplot(df2, aes(x=fisher.x, y=fisher.y)) + geom_point(aes(alpha = significance, color=color)) +
            scale_colour_manual(values = c("black", "orange", "red")) +
            #geom_text(aes(label = sequence, x=fisher.x, y=fisher.y, alpha = significance), size = 4, nudge_y = 0.5) +
            geom_text(data=filter(df2, significance < 1), aes(label = sequence, x=fisher.x, y=fisher.y, alpha = significance), size = 4, nudge_y = 0.5) +
            geom_text_repel(data=filter(df2, significance ==1), aes(label = sequence, alpha = significance), size = 4.5) +
            xlab(paste0("SS for ", fisher.ls.name_1))+
            ylab(paste0("SS for ", fisher.ls.name_2))+ #guides(color = FALSE)+
            geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) + theme(legend.position="none",panel.background=element_rect(fill = "white")) 
          k = k + 1
        }
      }
    }
    p.list = c(p.list, nrow=length(widths))
    newname = sub("\\.xlsx", paste0("_", sequenceRegion, ".png"), filename)
    png(newname, height=length(widths)*800, width=800*(k-1)/length(widths))
    do.call(grid.arrange, p.list)
    dev.off()
  }
}

####
fishermain.2 = function(df = tmp_n, keys, sequenceRegions = "3UTR", 
                        widths=c(4,6), frac=0.2, geno="mm9", key_pairs = "all.combinations", filename=filename){
  # improvement compared to fishermain: use exonic sequences instead of genomic sequences
  # tmp_n: df like pA.df (with column pA_pos) or pair.pA, with exonic 3'UTR sequences 
  # sequenceRegion: "-100_-41", "-40_-1", "1_100", "aUTR", or "3UTR"
  # widths = c(4,6)
  # geno: currently either hg19 or mm9
  # key_pairs: "all.combinations" or vectors like c("1,2", "1,3") to plot the SS of motifs on key1 vs key2 and key1 vs key3 plots
  
  # Define helper functions
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  fisher_row_sum_p = function(thisrow, sums){
    m = rbind(thisrow, sums - thisrow)
    sign = ifelse(m[1,1]/m[2,1] > m[1,2]/m[2,2], -1, 1) # postive: more enriched in upper
    fisher.test(m)$p.value * sign
  }
  
  oligofisher = function(df=df, key = "tail",  sequenceRegion = "aUTR", width=4, frac = 0.2, geno = "mm9", 
                         filename = filename){
    # sequenceRegion: "-100_-41", "-40_-1", "1_100", "aUTR", or "3UTR"
    # widths = c(4,6)
    # cat("key: ", key, "\n")
    # cat("colnames: ", names(df))
    
    require(GenomicRanges)
    require(Biostrings)
    require(xlsx)
    # remove NA and +/-Inf
    if(key %in% names(df)){
      df = as.data.frame(df[abs(df[, key]) != Inf & !is.na(df[, key]), ])
    }
    
    # calculate gr
    if(sequenceRegion == "aUTR"){
      df = pAid2pos(df)
      seq_df = df[, c("exonic_3UTR_seq_proximal", "exonic_3UTR_seq_distal")]
      seq_df$start = nchar(seq_df$exonic_3UTR_seq_proximal) + 1
      seq_df$stop = nchar(seq_df$exonic_3UTR_seq_distal)
      seqs = apply(seq_df, 1, function(row) substr(row[2], start = row[3], stop = row[4]))
      seqs = DNAStringSet(seqs)
    }else if(sequenceRegion == "3UTR"){
      seqs = df$exonic_3UTR_seq
      seqs = DNAStringSet(seqs)
    }else if(grepl("-?\\d+_-?\\d+", sequenceRegion)){ # sequenceRegion = "-100_-41"
      from = as.integer(strsplit(sequenceRegion, "_")[[1]][1])
      to = as.integer(strsplit(sequenceRegion, "_")[[1]][2])
      signs = ifelse(df$strand == "+", 1, -1)
      gr = GRanges(seqnames = df$chr, 
                   ranges = IRanges(start = pmin(df$pA_pos + from*signs, df$pA_pos + to*signs), 
                                    end = pmax(df$pA_pos + from*signs, df$pA_pos + to*signs)),
                   strand = df$strand)
      
      # get genomic sequence
      if(geno == "mm9"){
        require(BSgenome.Mmusculus.UCSC.mm9) 
        seqs = getSeq(Mmusculus, gr)
      }else if(geno == "hg19"){
        require("BSgenome.Hsapiens.UCSC.hg19") 
        seqs = getSeq(Hsapiens, gr)
      }
    }
    
    # get counts
    counts = oligonucleotideFrequency(seqs, width = width, step = 1)
    
    # limit the impact of a few sequences with large numbers of cis-elements
    ceilings = ceiling(apply(counts, 2, function(c) quantile(c, 0.95)))
    for(j in 1:length(ceilings)){
      counts[counts[, j] > ceilings[j], j] = ceilings[j]
    }
    
    #if(length(lower_index) == 1 && is.na(lower_index) && length(upper_index) == 1 && is.na(upper_index)){
    lower_index = df[, key] < quantile(df[, key], frac) 
    upper_index = df[, key] > quantile(df[, key], 1-frac) 
    #}
    
    counts = t(rbind(colSums(counts[lower_index,]), colSums(counts[upper_index,])))
    colnames(counts) = c("lower_counts", "upper_counts")
    
    sums = colSums(counts)
    
    fisher = apply(counts, 1, fisher_row_sum_p, sums)
    
    counts = as.data.frame(cbind(counts, fisher))
    counts$padj = p.adjust(abs(counts$fisher), "BH")
    
    rownames(counts) = gsub("T", "U", rownames(counts))
    
    counts$sequence = rownames(counts)
    
    rownames(counts) = NULL
    
    counts = counts[, c("sequence", names(counts)[1:4])]
    counts_original = counts
    # creat dfs of the same number of rows
    row_number_difference = sum(counts$fisher > 0) - sum(counts$fisher < 0)
    if(row_number_difference != 0){
      filling = rep("", abs(row_number_difference)) # fill the left or right part of the table, so that the two parts have the same number of rows
      filling = data.frame(sequence = filling,
                           lower_counts = filling,
                           upper_counts = filling,
                           fisher = 1,
                           padj = 1)
      if(row_number_difference < 0){
        left = rbind(subset(counts, fisher > 0), filling)
        right = subset(counts, fisher < 0)
      }else{
        right = rbind(subset(counts, fisher < 0), filling)
        left = subset(counts, fisher > 0)
      }
    }else{
      left = subset(counts, fisher > 0)
      right = subset(counts, fisher < 0)
    }
    
    left = left[order(left$fisher),]
    right = right[order(abs(right$fisher)),]
    counts = cbind(left, right)
    names(counts) = c("Sequence1", "Sequence1_in_Bottom",
                      "Sequence1_in_Top", "Sequence1_p", "Sequence1_padj",
                      "Sequence2", "Sequence2_in_Bottom",
                      "Sequence2_in_Top", "Sequence2_p", "Sequence2_padj")
    
    
    ####### as.numeric() is necessary because counts2$Fisher_Top is character
    counts$Sequence1_p = format(as.numeric(counts$Sequence1_p), digits= 2, scientific=T)
    counts$Sequence1_padj = format(as.numeric(counts$Sequence1_padj), digits= 2, scientific=T)
    counts$Sequence2_p = format(as.numeric(counts$Sequence2_p), digits= 2, scientific=T)
    counts$Sequence2_padj = format(as.numeric(counts$Sequence2_padj), digits= 2, scientific=T)
    ########### remove NA
    counts = apply(counts, 2, function(col) sub("NA", "", col)) 
    sheetName = paste0(sequenceRegion, "_", key, "_", width, "mer")
    write.xlsx2(counts, 
                filename, 
                sheetName = sheetName,
                row.names = F,
                append = T)
    
    return(counts_original) 
  }
  
  ## use the helper functions
  fisher.ls = list()
  cat("Number of rows: ", nrow(df), "\n")
  for(sequenceRegion in sequenceRegions){#sequenceRegions = c("-100_-41", "-40_0", "1_100")
    for(key in keys){
      for(width in widths){
        fisher.ls[[paste0(sequenceRegion, "_", key, "_", width,"mer")]] = oligofisher(df=df, 
                                                                                      key = key, 
                                                                                      width=width, 
                                                                                      sequenceRegion=sequenceRegion,
                                                                                      geno = geno, #lower_index, upper_index,
                                                                                      frac = 0.2, #max.padj = max.padj,
                                                                                      filename = filename)
      }
    }
    
    if(length(keys) > 1){
      # Scatter plot
      require(dplyr)
      require(tidyr)
      require(ggplot2)
      require(ggrepel)
      max.padj = 0.05
      
      p.list = list()
      k = 1
      for(width in widths){
        for(i in 1:length(keys[-1])){
          for(j in (i+1):length(keys)){
            if(key_pairs != "all.combinations"){
              if(!paste(i,j,sep=",") %in% key_pairs){
                next
              }
            }
            fisher.ls.name_1 = paste0(sequenceRegion, "_", keys[i], "_", width, "mer")
            fisher.ls.name_2 = paste0(sequenceRegion, "_", keys[j], "_", width, "mer")
            #cmp.ls.name = paste0("cmp.", width, "mer")
            df2 = inner_join(fisher.ls[[fisher.ls.name_1]][, c("sequence", "fisher", "padj")], 
                             fisher.ls[[fisher.ls.name_2]][, c("sequence", "fisher", "padj")], 
                             by = "sequence")
            fisher.x = df2$fisher.x
            fisher.y = df2$fisher.y
            
            df2$fisher.x = -1*sign(fisher.x)*log10(abs(fisher.x))
            df2$fisher.y = -1*sign(fisher.y)*log10(abs(fisher.y))
            
            df2$significance = c(0.001, 0.5, 1)[(df2$padj.x < max.padj) + (df2$padj.y < max.padj) + 1]
            df2$color = factor((df2$padj.x < max.padj) + (df2$padj.y < max.padj) + 1)
            
            p.list[[k]] = ggplot(df2, aes(x=fisher.x, y=fisher.y)) + geom_point(aes(alpha = significance, color=color)) +
              scale_colour_manual(values = c("black", "orange", "red")) +
              #geom_text(aes(label = sequence, x=fisher.x, y=fisher.y, alpha = significance), size = 4, nudge_y = 0.5) +
              geom_text(data=filter(df2, significance < 1), aes(label = sequence, x=fisher.x, y=fisher.y, alpha = significance), size = 4, nudge_y = 0.5) +
              geom_text_repel(data=filter(df2, significance ==1), aes(label = sequence, alpha = significance), size = 4.5) +
              xlab(paste0("SS for ", fisher.ls.name_1))+
              ylab(paste0("SS for ", fisher.ls.name_2))+ #guides(color = FALSE)+
              geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) + theme(legend.position="none",panel.background=element_rect(fill = "white")) 
            k = k + 1
          }
        }
      }
      p.list = c(p.list, nrow=length(widths))
      newname = sub("\\.xlsx", paste0("_", sequenceRegion, ".png"), filename)
      png(newname, height=length(widths)*800, width=800*(k-1)/length(widths))
      do.call(grid.arrange, p.list)
      dev.off()
    }
    
  }
}
fishermain = fishermain.2
#### how to use fishermain() and fishermain.2() 
# ## cis-element analysis for aUTRs 
# keys = grep("EP_CP_._red", names(pair.pA), value=T)[1:3]
# filename = file.path(result_dir, paste0("delta_EP_CP_aUTR_Motif_Fisher_", gsub(" |-|:", "_", Sys.time()), ".xlsx"))
# 
# fishermain.2(df = pair.pA, keys=keys, sequenceRegions = "aUTR", 
#            widths=c(4,6), frac=0.2, geno="mm9", key_pairs = "all.combinations", filename=filename)
# 
# pA.df = read.csv("pA.df.csv", as.is = T)
# 
# ## motif around pAs:
# for(n in c(10, 20)){
#   # n = 20
#   names(pA.df)
#   tmp_n = subset(pA.df, rowSums(pA.df[, grep("_count$", names(pA.df))] > n) == length(grep("_count$", names(pA.df))))
#   keys = grep("_Cs_Cm_[12]_rpm_ratio", names(tmp_n), value=T)
#   sequenceRegions = c("-100_-41", "-40_-1", "1_100")
#   
#   filename = file.path(result_dir = ".", paste0("Around_pA_Motif_Fisher_n=", n,"_", gsub(" |-|:", "_", Sys.time()), ".xlsx"))
#   fishermain.2(df = tmp_n, keys=keys, sequenceRegions = sequenceRegions, widths=c(4,6), frac=0.2, geno="mm9", 
#              key_pairs=c("1,3", "2,4"), filename = filename)
# }
# 
# ## motif in 3'UTRs:
# for(n in c(10, 20)){
#   # n = 10
#   names(pA.df)
#   tmp_n = subset(pA.df, rowSums(pA.df[, grep("_count$", names(pA.df))] > n) == length(grep("_count$", names(pA.df))))
#   keys = grep("_Cs_Cm_[12]_rpm_ratio", names(tmp_n), value=T)
#   sequenceRegions = "3UTR"
#   
#   filename = file.path(result_dir=".", paste0("3UTR_Motif_Fisher_n=", n,"_", gsub(" |-|:", "_", Sys.time()), ".xlsx"))
#   fishermain.2(df = tmp_n, keys=keys, sequenceRegions = sequenceRegions, widths=c(4,6), frac=0.2, geno="mm9", 
#              key_pairs=c("1,3", "2,4"), filename = filename)
# }
####################### end of fishermain ##################################

####
fishermain.color = function(df = pair.pA, sequenceRegions = c("cUTR", "aUTR"), 
                            widths=c(2,4,6), geno="mm9", filename=filename){
  # pair.pA: must contain a color column 
  # sequenceRegion:  "aUTR" or "cUTR" 
  # widths = c(4,6)
  # geno: currently either hg19 or mm9
  if(!"color" %in% names(df)){
    stop("The input data frame must contain a column named 'color'")
  }
  
  ############ define helper functions
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  fisher_row_sum_p = function(thisrow, sums){
    m = rbind(thisrow, sums - thisrow)
    sign = ifelse(m[1,1]/m[2,1] > m[1,2]/m[2,2], 1, -1) # negative: more enriched in background (column 2)
    fisher.test(m)$p.value * sign
  }
  
  oligofisher = function(df=df, color = "blue", sequenceRegion = "cUTR", width=2, geno = "mm9", filename = filename){
    require(GenomicRanges)
    require(Biostrings)
    require(xlsx)
    
    # calculate gr
    if(sequenceRegion == "aUTR"){
      df = pAid2pos(df)
      seq_df = df[, c("exonic_3UTR_seq_proximal", "exonic_3UTR_seq_distal")]
      seq_df$start = nchar(seq_df$exonic_3UTR_seq_proximal) + 1
      seq_df$stop = nchar(seq_df$exonic_3UTR_seq_distal)
      seqs = apply(seq_df, 1, function(row) substr(row[2], start = row[3], stop = row[4]))
      seqs = DNAStringSet(seqs)
    }else if(sequenceRegion == "cUTR"){
      seqs = DNAStringSet(df$exonic_3UTR_seq_proximal)
    } 
    
    # get counts
    counts = oligonucleotideFrequency(seqs, width = width, step = 1)
    
    # limit the impact of a few sequences with large numbers of cis-elements
    ceilings = ceiling(apply(counts, 2, function(c) quantile(c, 0.95)))
    for(j in 1:length(ceilings)){
      counts[counts[, j] > ceilings[j], j] = ceilings[j]
    }
    
    fg_index = df$color == color 
    bg_index = df$color != color 
    
    counts = t(rbind(colSums(counts[fg_index,]), colSums(counts[bg_index,])))
    colnames(counts) = c("fg_counts", "bg_counts")
    
    sums = colSums(counts)
    
    fisher = apply(counts, 1, fisher_row_sum_p, sums)
    
    counts = as.data.frame(cbind(counts, fisher))
    #counts$padj = p.adjust(abs(counts$fisher), "BH")
    
    rownames(counts) = gsub("T", "U", rownames(counts))
    
    counts$sequence = rownames(counts)
    
    rownames(counts) = NULL
    
    counts = counts[, c("sequence", names(counts)[1:3])]
    counts_original = counts
    # creat dfs of the same number of rows
    row_number_difference = sum(counts$fisher > 0) - sum(counts$fisher < 0)
    if(row_number_difference != 0){
      filling = rep("", abs(row_number_difference)) # fill the left or right part of the table, so that the two parts have the same number of rows
      filling = data.frame(sequence = filling,
                           fg_counts = filling,
                           bg_counts = filling,
                           fisher = 1)
      if(row_number_difference < 0){
        left = rbind(subset(counts, fisher > 0), filling)
        right = subset(counts, fisher < 0)
      }else{
        right = rbind(subset(counts, fisher < 0), filling)
        left = subset(counts, fisher > 0)
      }
    }else{
      left = subset(counts, fisher > 0)
      right = subset(counts, fisher < 0)
    }
    
    left = left[order(left$fisher), ]
    right = right[order(abs(right$fisher)), ]
    counts = cbind(left, right)
    
    sheetName = paste0(color, "_", sequenceRegion, "_", width, "mer")
    write.xlsx2(counts, 
                filename, 
                sheetName = sheetName,
                row.names = F,
                append = T)
    
    return(counts_original) 
  }
  
  ## use the helper functions
  fisher.ls = list()
  for(width in widths){
    for(color in c("blue", "red")){
      for(sequenceRegion in sequenceRegions){ 
        fisher.ls[[paste0(color, "_", sequenceRegion, "_", width, "mer")]] = oligofisher(df=df, 
                                                                                         color = color, 
                                                                                         sequenceRegion=sequenceRegion,
                                                                                         width=width,
                                                                                         geno = geno,  
                                                                                         filename = filename)
      }
    }
    
    
  }
}

fishermain.fg.vs.bg = function(df = pair.pA, sequenceRegions = c("cUTR", "aUTR"), widths=c(2,4,6), 
                               fg_color = "red", bg_color = "blue", geno="mm9", filename=filename){
  # pair.pA: must contain a color column 
  # sequenceRegion:  "aUTR" or "cUTR" 
  # widths = c(4,6)
  # geno: currently either hg19 or mm9
  if(!"color" %in% names(df)){
    stop("The input data frame must contain a column named 'color'")
  }
  
  # define helper functions
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  fisher_row_sum_p = function(thisrow, sums){
    m = rbind(thisrow, sums - thisrow)
    sign = ifelse(m[1,1]/m[2,1] > m[1,2]/m[2,2], 1, -1) # negative: more enriched in background (column 2)
    fisher.test(m)$p.value * sign
  }
  
  oligofisher = function(df=df, sequenceRegion = "cUTR", width=2, geno = "mm9", filename = filename){
    require(GenomicRanges)
    require(Biostrings)
    require(xlsx)
    
    # calculate gr
    if(sequenceRegion == "aUTR"){
      df = pAid2pos(df)
      seq_df = df[, c("exonic_3UTR_seq_proximal", "exonic_3UTR_seq_distal")]
      seq_df$start = nchar(seq_df$exonic_3UTR_seq_proximal) + 1
      seq_df$stop = nchar(seq_df$exonic_3UTR_seq_distal)
      seqs = apply(seq_df, 1, function(row) substr(row[2], start = row[3], stop = row[4]))
      seqs = DNAStringSet(seqs)
    }else if(sequenceRegion == "cUTR"){
      seqs = DNAStringSet(df$exonic_3UTR_seq_proximal)
    } 
    
    # get counts
    counts = oligonucleotideFrequency(seqs, width = width, step = 1)
    
    # limit the impact of a few sequences with large numbers of cis-elements
    ceilings = ceiling(apply(counts, 2, function(c) quantile(c, 0.95)))
    for(j in 1:length(ceilings)){
      counts[counts[, j] > ceilings[j], j] = ceilings[j]
    }
    
    fg_index = df$color == fg_color 
    bg_index = df$color == bg_color 
    
    counts = t(rbind(colSums(counts[fg_index,]), colSums(counts[bg_index,])))
    colnames(counts) = c("fg_counts", "bg_counts")
    
    sums = colSums(counts)
    
    fisher = apply(counts, 1, fisher_row_sum_p, sums)
    
    counts = as.data.frame(cbind(counts, fisher))
    #counts$padj = p.adjust(abs(counts$fisher), "BH")
    
    rownames(counts) = gsub("T", "U", rownames(counts))
    
    counts$sequence = rownames(counts)
    
    rownames(counts) = NULL
    
    counts = counts[, c("sequence", names(counts)[1:3])]
    counts_original = counts
    # creat dfs of the same number of rows
    row_number_difference = sum(counts$fisher > 0) - sum(counts$fisher < 0)
    if(row_number_difference != 0){
      filling = rep("", abs(row_number_difference)) # fill the left or right part of the table, so that the two parts have the same number of rows
      filling = data.frame(sequence = filling,
                           fg_counts = filling,
                           bg_counts = filling,
                           fisher = 1)
      if(row_number_difference < 0){
        left = rbind(subset(counts, fisher > 0), filling)
        right = subset(counts, fisher < 0)
      }else{
        right = rbind(subset(counts, fisher < 0), filling)
        left = subset(counts, fisher > 0)
      }
    }else{
      left = subset(counts, fisher > 0)
      right = subset(counts, fisher < 0)
    }
    
    left = left[order(left$fisher), ]
    right = right[order(abs(right$fisher)), ]
    counts = cbind(left, right)
    
    sheetName = paste0(sequenceRegion, "_", width, "mer")
    write.xlsx2(counts, 
                filename, 
                sheetName = sheetName,
                row.names = F,
                append = T)
    
    return(counts_original) 
  }
  
  ## use the helper functions
  fisher.ls = list()
  for(width in widths){
    for(sequenceRegion in sequenceRegions){ 
      fisher.ls[[paste0(sequenceRegion, "_", width, "mer")]] = oligofisher(df=df, 
                                                                           sequenceRegion=sequenceRegion,
                                                                           width=width,
                                                                           geno = geno,  
                                                                           filename = filename)
    }
  }
}
fishermain.color = fishermain.fg.vs.bg


cUTR.aUTR.Zsw = function(df = tmp_n, keys, 
                            widths=c(4,6), group = "group", frac=0.2, geno="mm9", 
                            filename=filename){
  # improvement compared to fishermain: use exonic sequences instead of genomic sequences
  # tmp_n: df like pA.df (with column pA_pos) or pair.pA, with exonic 3'UTR sequences 
  # sequenceRegion: "-100_-41", "-40_-1", "1_100", "aUTR", or "3UTR"
  # widths = c(4,6)
  # geno: currently either hg19 or mm9
  
  Zsw = function(rowdat){
    # rowdat: "lower_cUTR_counts", "upper_cUTR_counts", "lower_aUTR_counts", "upper_aUTR_counts"
    # ns: # of k-mer in "strong" set
    # nw: # of k-mer in "weak" set
    # Ns: total number of k-mers in "strong" set
    # Nw: total number of k-mers in "weak" set
    ns = rowdat[1]
    nw = rowdat[2]
    Ns = rowdat[3]
    Nw = rowdat[4]
    p = (ns + nw)/(Ns + Nw)
    (ns/Ns - nw/Nw)/sqrt((1/Ns + 1/Nw)*p*(1-p))
  }
  
  oligofisher = function(df=df, key = "tail", 
                         width = 4, group = "group", frac = 0.2, geno = "mm9", 
                         filename = filename){
    require(GenomicRanges)
    require(Biostrings)
    require(xlsx)
    
    if(key %in% names(df)){
      df = as.data.frame(df[abs(df[, key]) != Inf & !is.na(df[, key]), ])
    }
    
    df = pAid2pos(df)
    seq_df = df[, c("exonic_3UTR_seq_proximal", "exonic_3UTR_seq_distal")]
    cUTR_seqs = DNAStringSet(df$exonic_3UTR_seq_proximal)
    cUTR_counts = oligonucleotideFrequency(cUTR_seqs, width = width, step = 1)
    # limit the impact of a few sequences with large numbers of cis-elements
    ceilings = ceiling(apply(cUTR_counts, 2, function(c) quantile(c, 0.95)))
    for(j in 1:length(ceilings)){
      cUTR_counts[cUTR_counts[, j] > ceilings[j], j] = ceilings[j]
    }
    
    
    seq_df$start = nchar(seq_df$exonic_3UTR_seq_proximal) + 1
    seq_df$stop = nchar(seq_df$exonic_3UTR_seq_distal)
    aUTR_seqs = apply(seq_df, 1, function(row) substr(row[2], start = row[3], stop = row[4]))
    aUTR_seqs = DNAStringSet(aUTR_seqs)
    aUTR_counts = oligonucleotideFrequency(aUTR_seqs, width = width, step = 1)
    # limit the impact of a few sequences with large numbers of cis-elements
    ceilings = ceiling(apply(aUTR_counts, 2, function(c) quantile(c, 0.95)))
    for(j in 1:length(ceilings)){
      aUTR_counts[aUTR_counts[, j] > ceilings[j], j] = ceilings[j]
    }
    
    if(group %in% names(df) && all(c("Top", "Bottom") %in% unique(df$group))){
      lower_index = df[, group] == "Bottom"
      upper_index = df[, group] == "Top"
    }else{
      lower_index = df[, key] < quantile(df[, key], frac) 
      upper_index = df[, key] > quantile(df[, key], 1-frac) 
    }
    
    cUTR_counts = t(rbind(colSums(cUTR_counts[lower_index,]), colSums(cUTR_counts[upper_index,])))
    colnames(cUTR_counts) = c("lower_cUTR_counts", "upper_cUTR_counts")
    
    aUTR_counts = t(rbind(colSums(aUTR_counts[lower_index,]), colSums(aUTR_counts[upper_index,])))
    colnames(aUTR_counts) = c("lower_aUTR_counts", "upper_aUTR_counts")
    
    counts = cbind(cUTR_counts, aUTR_counts)
    counts = as.data.frame(counts)
    
    counts$cUTR_Ns = sum(counts$upper_cUTR_counts)
    counts$cUTR_Nw = sum(counts$lower_cUTR_counts)
    counts$aUTR_Ns = sum(counts$upper_aUTR_counts)
    counts$aUTR_Nw = sum(counts$lower_aUTR_counts)
    
    counts$cUTR_Zsw = apply(counts[, c("upper_cUTR_counts", "lower_cUTR_counts", "cUTR_Ns", "cUTR_Nw")], 1, Zsw)
    counts$aUTR_Zsw = apply(counts[, c("upper_aUTR_counts", "lower_aUTR_counts", "aUTR_Ns", "aUTR_Nw")], 1, Zsw)
    
    counts$Top_Zac = apply(counts[, c("upper_aUTR_counts", "upper_cUTR_counts", "aUTR_Ns", "cUTR_Ns")], 1, Zsw)
    counts$Bottom_Zac = apply(counts[, c("lower_aUTR_counts", "lower_cUTR_counts", "aUTR_Nw", "cUTR_Nw")], 1, Zsw)
    
    rownames(counts) = gsub("T", "U", rownames(counts))
    counts$sequence = rownames(counts)
    counts = counts[, c("sequence", "lower_cUTR_counts", "upper_cUTR_counts", "lower_aUTR_counts", "upper_aUTR_counts", 
                        "cUTR_Zsw", "aUTR_Zsw", "Top_Zac", "Bottom_Zac")]
    
    counts = apply(counts, 2, function(col) sub("NA", "", col)) 
    sheetName = paste0(key, "_", width, "mer")
    write.xlsx2(counts, 
                filename, 
                sheetName = sheetName,
                row.names = F,
                append = T)
    
  }
  
  ## use the helper function
  cat("Number of rows: ", nrow(df), "\n")
  for(key in keys){
    # key = keys[1]
    for(width in widths){
      # width = widths[1]
      oligofisher(df=df, 
                  key = key, 
                  width=width, 
                  geno = geno, 
                  group = "group",
                  frac = 0.2, 
                  filename = filename)
      
    }
  }
}

#### A helper function to plot contour and smooth line in ggpairs
contour.smooth <- function(data, mapping, method="loess", ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point(shape=1, alpha=0.5) + 
    stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') +
    scale_fill_continuous(low="green",high="red") +
    geom_smooth(method=method, ...)
  p
}
## example:
# n = 10
# p = ggpairs(df, 
#             upper = list(continuous = wrap(contour.smooth, method="lm")),
#             lower = list(continuous = "cor"),
#             title = paste0("Read # >= ", n, " | Transcript # = ", nrow(df)))
# print(p)


#####
align_polysome_profile_with_reference = function(pp, ref, profile_name = "NT_Cytosol"){
  ## Arguments:
  # pp: a data frame with two columns ("Second", "A254"). This is the polysome profile to be aligned with reference (ref) profile
  # ref: The reference polysome profile, which will not be truncated
  # profile_names: names for the polysome profile pp
  ## Output: 
  # pp aligned to ref. Part of the beginning of pp may be deleted or padded with 0. "ref" profile will not change
  
  # slide pp up with respect to ref (delete the first best_off_set_up rows of pp) 
  max_cor_up = 0
  best_off_set_up = 0
  for(off_set in 0:round(0.2*nrow(pp))){ # shift pp no more than 1/5 of pp length
    tmp_cor = cor(pp[(off_set+1):nrow(pp), ]$A254, ref[1:(nrow(pp) - off_set), ]$A254)
    if(tmp_cor > max_cor_up){
      max_cor_up = tmp_cor
      best_off_set_up = off_set
    }
  }
  # slide pp down with respect to ref (padd the first best_off_set_down rows of pp with 0s) 
  max_cor_down = 0
  best_off_set_down = 0
  for(off_set in 0:round(0.2*nrow(pp))){ # shift pp no more than 1/5 of pp length
    tmp_cor = cor(pp[1:min(nrow(pp), nrow(ref) - off_set), ]$A254, ref[(off_set + 1):min(nrow(pp) + off_set, nrow(ref)), ]$A254)
    if(tmp_cor > max_cor_down){
      max_cor_down = tmp_cor
      best_off_set_down = off_set
    }
  }
  # Use the offset with maximum correlation to align data
  if(max_cor_up >= max_cor_down){
    pp_new = pp[(best_off_set_up + 1):nrow(pp), ] # only allowing discarding data points at beginning
    #ref_new = ref[1:(nrow(pp) - best_off_set_up), ] # discarding useful data, bad!
  }else{
    padding = matrix(rep(NA, best_off_set_down*2), ncol = 2)
    colnames(padding) = colnames(pp)
    pp_new = as.data.frame(rbind(padding, pp)) # pad pp, so that we don't need to discard data points in ref
  }
  # add sample name
  pp_new$Sample = profile_name
  # Recalculate time
  ref$Second = ref$Second - ref$Second[1]
  pp_new$Second = ref$Second[1:nrow(pp_new)]
  
  # remove padded rows
  pp_new = subset(pp_new, !is.na(A254))
  
  # smooth data
  tmp = c(pp_new$A254, rep(0, 10))
  for(i in 1:10){
    new_col = c(rep(0, i), pp_new$A254, rep(0, 10-i))
    tmp = cbind(tmp, new_col)
  }
  pp_new$A254 = rowMeans(tmp)[6:(nrow(tmp)-5)]
  
  # return
  pp_new 
}

###### Sometimes ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl") will fail
robust_useMart = function(biomart = "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl"){
  tryCatch(useMart(biomart = biomart, dataset=dataset),
           #warning = function(w) print(w),
           error = function(e) {
             #print("Retrying ....")
             robust_useMart(biomart = biomart, dataset=dataset)}
  )
}
# used like: ensembl = robust_useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

#### Add mouse protein localization annotation without combining different locations
add_mouse_protein_localization = function(df = pA.df, 
                                    annotation_file = "/home/dinghai/projects/fud/mouse_protein_locations.csv"){
  # Data: "MetazSecKB.all_metazoa_protein_subloc.txt" downloaded from http://proteomics.ysu.edu/secretomes/animal/index.php
  require(dplyr)
  
  protein_locations = read.csv(annotation_file, as.is = T)
  protein_locations$localizations = NULL  
  
  names(protein_locations) = paste0("YSU.", names(protein_locations)) 

  names(protein_locations) = sub("YSU.mgi_symbol", "gene_symbol", names(protein_locations)) 
  
  # The genes not found in protein_locations will be NA
  df = left_join(df, protein_locations) 
  
  # When df is pA.df, pAs in introns, UA, etc should not be annotated with protein localization
  if("region" %in% names(df)){ 
    for(loc in names(protein_locations)[-1]){
      df[df$region != "3UTR", loc] = NA
    }
  }
  
  df
}

#### Add protein localization annotation 
add_protein_localization = function(df = pA.df, 
                                    annotation_file = "/home/dinghai/projects/fud/mouse_protein_locations.csv"){
  require(dplyr)
  
  protein_locations = read.csv(annotation_file, as.is = T) 
  protein_locations = protein_locations %>%  
    mutate(locate2membrane = apply(protein_locations[, c("ER", "Golgi", "Membrane", "Lysosome", "Peroxisome")], 1, any)) %>%
    mutate(secreted_or_highlylikely_secreted = apply(protein_locations[, c("Secreted", "Secreted_highlylikely")], 1, any)) 
  
  names(protein_locations) = sub("mgi_symbol", "gene_symbol", names(protein_locations)) 
  
  # The MoS of genes not found in protein_locations will be NA
  df = left_join(df, protein_locations[, c("gene_symbol", "locate2membrane", "secreted_or_highlylikely_secreted")]) 
  
  # When df is pA.df, pAs in introns, UA, etc should not be annotated with protein localization
  if("region" %in% names(df)){ 
    df = mutate(df, locate2membrane = ifelse(region != "3UTR", NA, locate2membrane)) %>%
      mutate(secreted_or_highlylikely_secreted = ifelse(region != "3UTR", NA, secreted_or_highlylikely_secreted))
  }
  
  df
}


##### Add human protein localization annotation 
add_human_protein_localization = function(df = pA.df, 
                                          annotation_file = "/home/dinghai/projects/fud/cell_atlas/subcellular_location.tsv"){
  # Using data downloaded from https://www.proteinatlas.org/download/subcellular_location.tsv.zip
  # See the paper titled "A subcellular map of the human proteome"
  #    In the Cell Atlas, we provide a reliability score for every annotated location and
  # protein on a four-tiered scale: “validated (enhanced?),” “supported,” “approved,” and “uncertain.” Locations
  # obtained the score “validated” if the antibody was validated according to one of the validation “pillars”
  # proposed by an international working group (36) as suitable for IF: (i) genetic methods using
  # short interfering RNA (siRNA) silencing (37) or CRISPR-Cas9 knockout, (ii) expression of a fluorescent protein–tagged 
  # protein at endogenous levels (38), or (iii) independent antibodies targeting different epitopes (see fig. S4 for examples).
  # The second tier, “supported” locations, is defined by agreement with external experimental data
  # from the UniProt database. An “approved” location score indicates a lack of external experimental information about 
  # the protein location. Last, an “uncertain” location is contradictory to complementary information, such as literature
  # or transcriptomics data. “Uncertain” locations are only shown when it cannot be ruled out that the data are correct.
  
  require(dplyr)
  
  protein_locations = read.table(annotation_file, header=T, sep="\t", as.is = T) # The data are not tid
  
  locations_1 = unique(strsplit(paste(unique(protein_locations$Enhanced), collapse = ";"), ";")[[1]])
  locations_2 = unique(strsplit(paste(unique(protein_locations$Supported), collapse = ";"), ";")[[1]])
  locations_3 = unique(strsplit(paste(unique(protein_locations$Approved), collapse = ";"), ";")[[1]])
  locations_4 = unique(strsplit(paste(unique(protein_locations$Uncertain), collapse = ";"), ";")[[1]])
  
  locations = unique(c(locations_1,locations_2,locations_3,locations_4)) 
  locations = locations[locations != ""]
  
  # Loop through each location and determine if the value should be "Enhanced", "Supported", etc for each gene
  for(loc in locations){
    # loc = "Vesicles"
    protein_locations[, loc] = NA
    for(reliability in c("Enhanced", "Supported", "Approved", "Uncertain")){
      # reliability = "Enhanced"
      protein_locations[, loc] = ifelse(grepl(loc, protein_locations[, reliability]), reliability, protein_locations[, loc])
    }
  }
  
  names(protein_locations) = gsub(" ", ".", names(protein_locations))
  names(protein_locations) = gsub("&", "and", names(protein_locations))
  names(protein_locations) = paste0("Atlas.", names(protein_locations))
  locations = gsub(" ", ".", locations)
  locations = gsub("&", "and", locations)
  locations = paste0("Atlas.", locations)
  
  protein_locations = protein_locations %>%  
    dplyr::rename(gene_symbol = Atlas.Gene.name) %>%
    mutate(gene_symbol = Hmisc::capitalize(tolower(gene_symbol)))
  
  # Genes not found in protein_locations will be NA
  df = left_join(df, protein_locations[, c("gene_symbol", locations)])
  
  # When df is pA.df, pAs in introns, UA, etc should not be annotated with protein localization
  if("region" %in% names(df)){ 
    for(loc in locations){
      df[, loc] = ifelse(df$region != "3UTR", NA, df[, loc])
    }
  }
  
  df
}

#### Filter RED values based on read numbers
filterRED = function(input_df = pair.pA, red_type = "C2C12D_EP_vs_EM_2_RED", read_num_cutoff = 5){
  # Returns a subset of input_df containing genes whose read numbers involved in calculating the red_type are all >= the read_num_cutoff
  input_df = unique(input_df)
  pieces  = setdiff(strsplit(red_type, "_")[[1]], "vs")[1:3]
  
  col_candidates = grep("_count_(distal|proximal)", names(input_df), value=T)
  col_selected = vector()
  for(candidate in col_candidates){
    candidate_pieces = strsplit(candidate, "_")[[1]][1:2]
    if(all(candidate_pieces %in% pieces)){
      col_selected = c(col_selected, candidate)
    }
  }
  
  print("Columns used for filtering data: ")
  print(col_selected)
  
  expressed = rowSums(input_df[, col_selected] >= read_num_cutoff) == length(col_selected)
  input_df[expressed, ]
}


#### Filter RE values based on read numbers
filterRE = function(input_df = pA.df, re_type = "stau1_egfp_IP_mean_rpm_ratio", read_num_cutoff = 5){
  # Returns a subset of input_df containing genes whose read numbers involved in calculating the re_type are all >= the read_num_cutoff
  input_df = unique(input_df)
  pieces  = setdiff(strsplit(re_type, "_")[[1]], "vs")[1:3]
  
  col_candidates = grep("_count$", names(input_df), value=T)
  col_selected = vector()
  for(candidate in col_candidates){
    candidate_pieces = strsplit(candidate, "_")[[1]][1:2]
    if(all(candidate_pieces %in% pieces)){
      col_selected = c(col_selected, candidate)
    }
  }
  
  print("Columns used for filtering data: ")
  print(col_selected)
  
  expressed = rowSums(input_df[, col_selected] >= read_num_cutoff) == length(col_selected)
  input_df[expressed, ]
}

###########################
filter_RED = function(input_df = pair.pA, red = "C2C12D_EP_vs_EM_2_RED", batch = 2, read_num_cutoff = 5){
  # Returns a subset of input_df containing genes whose read numbers involved in calculating the red_type are all >= the read_num_cutoff
  input_df = unique(input_df)
  pieces  = setdiff(strsplit(red, "_")[[1]], "vs")
  
  col_candidates = grep("_count_(distal|proximal)", names(input_df), value=T)
  col_selected = vector()
  for(candidate in col_candidates){
    candidate_pieces = strsplit(candidate, "_")[[1]][1:2]
    
    if(batch != "all"){
      if(!batch %in% pieces) stop("'batch' should be consistent with 'red'!")
      candidate_pieces = c(candidate_pieces, batch)
    }
    
    if(all(candidate_pieces %in% pieces)){
      col_selected = c(col_selected, candidate)
    }
  }
  
  print("Columns used for filtering data: ")
  print(col_selected)
  
  expressed = rowSums(input_df[, col_selected] >= read_num_cutoff) == length(col_selected)
  input_df[expressed, ]
}

###########################
filter_RE = function(input_df = pA.df, re = "stau1_egfp_IP_mean_rpm_ratio", batch = "all", read_num_cutoff = 5){
  # Returns a subset of input_df containing genes whose read numbers involved in calculating the re_type are all >= the read_num_cutoff
  input_df = unique(input_df)
  pieces  = setdiff(strsplit(re, "_")[[1]], "vs")
  
  col_candidates = grep("_count$", names(input_df), value=T)
  col_selected = vector()
  for(candidate in col_candidates){
    candidate_pieces = strsplit(candidate, "_")[[1]][1:2]
    
    if(batch != "all"){
      if(!batch %in% pieces) stop("'batch' should be consistent with 're'!")
      candidate_pieces = c(candidate_pieces, batch)
    }
    
    if(all(candidate_pieces %in% pieces)){
      col_selected = c(col_selected, candidate)
    }
  }
  
  print("Columns used for filtering data: ")
  print(col_selected)
  
  expressed = rowSums(input_df[, col_selected] >= read_num_cutoff) == length(col_selected)
  input_df[expressed, ]
}


####### Heatmap with bins
one_heatmap = function(df1=pair.pA, bin_x = "CDS_size", bin_y = "aUTR_length", bin_num = 5, 
                       fill = "C2C12D_EP_EM_1_RED", read_num_cutoff = 5, require_CDS = T, one_tx_per_gene=T){
  # A function to generate one heatmap at a time
  require(dplyr)
  require(ggplot2)
  require(Hmisc)
  
  # Step 1: filter records with low read number 
  if(require_CDS){
    df1 = df1[df1$CDS_size > 3, ]
  }
  if(grepl("_RED", fill)){
    fill_pieces = strsplit(fill, "_")[[1]][1:3]
    col_candidates = grep("_count_(distal|proximal)", names(df1), value=T)
    col_selected = vector()
    for(candidate in col_candidates){
      candidate_pieces = strsplit(candidate, "_")[[1]][1:2]
      if(all(candidate_pieces %in% fill_pieces)){
        col_selected = c(col_selected, candidate)
      }
    }
    expressed = rowSums(df1[, col_selected] >= read_num_cutoff) == length(col_selected)
    df1 = df1[expressed,]
    df1 = df1[, c(bin_x, bin_y, fill)]
  }
  if(grepl("_RE$|_rpm_ratio", fill)){
    fill_pieces = strsplit(fill, "_")[[1]][1:3]
    col_candidates = grep("_count$", names(df1), value=T)
    col_selected = vector()
    for(candidate in col_candidates){
      candidate_pieces = strsplit(candidate, "_")[[1]][1:2]
      if(all(candidate_pieces %in% fill_pieces)){
        col_selected = c(col_selected, candidate)
      }
    }
    expressed = rowSums(df1[, col_selected] >= read_num_cutoff) == length(col_selected)
    df1 = df1[expressed, ]
    if(one_tx_per_gene){
      df1 = df1 %>%
        mutate(total_count = rowSums(df1[, col_selected])) %>%
        group_by(gene_symbol) %>% 
        top_n(1, total_count) %>%
        ungroup()
    }
    df1 = df1[, c(bin_x, bin_y, fill)]
  }
  
  df1 = df1 %>%
    na.omit() 

  # Step 2: Calculations
  df1 = as.data.frame(df1)
  df1[, bin_x] = cut2(df1[, bin_x], g = bin_num)
  df1[, bin_y] = cut2(df1[, bin_y], g = bin_num)
  df1 = df1 %>%
    group_by_at(vars(bin_x, bin_y)) %>%
    mutate(count = n()) %>%
    summarize_at(vars(fill, "count"), funs(mean, median)) %>%
    ungroup() %>%
    mutate_at(vars(paste0(fill, "_median")), funs(scale)) %>%
    mutate_at(vars(paste0(fill, "_median")), funs(ifelse(. > 2, 2, .))) %>% # maximum: 2
    mutate_at(vars(paste0(fill, "_median")), funs(ifelse(. < -2, -2, .)))  # minimum: -2
  
  total_count = sum(df1$count_median)
  
  main_title = sub("EP_EM", "PM(ER)", fill)
  main_title = sub("CP_CM", "PM(Cyto)", main_title)
  main_title = sub("EM_CM", "EC(M)", main_title)
  main_title = sub("EP_CP", "EC(P)", main_title)

  if(one_tx_per_gene & grepl("_RE$|_rpm_ratio", fill)){
    main_title = paste0("Normalized median of ", main_title, "\nTotal count: ", total_count, " | One transcript/gene")
  }else{
    main_title = paste0("Normalized median of ", main_title, "\nTotal count: ", total_count)
  }
  
  # Step 3: Plotting
  ggplot(df1, aes_string(x=bin_x, y=bin_y, fill = paste0(fill, "_median"))) + 
    geom_tile() +
    #scale_fill_distiller(palette = "RdYlBu", guide = guide_legend(title = "median"), limits = c(-0.8, 0.8)) +
    scale_fill_gradientn(colours = rainbow(100, start=0, end = 0.7)[100:1], limits = c(-2, 2)) +
    ggtitle(main_title) +
    geom_text(aes(label=count_median)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title=element_blank(), plot.title=element_text(hjust = 0.5))
}


#### A function that generate multiple heatmaps and save them
multi_heatmap = function(df = pair.pA, bin_x = "CDS_size", bin_y = "aUTR_length", bin_nums = c(5, 10),
                         fill_pattern = "_RED$", read_num_cutoffs = c(2, 5), require_CDS = T, one_tx_per_gene=T, 
                         width_per_plot = 650, height_per_plot = 600, nrow = 2){
  names(df) = sub("maxCDS", "CDS_size", names(df))
  for(bin_num in bin_nums){
    for(read_num_cutoff in read_num_cutoffs){
      plts = list()
      k = 1
      for(fill in grep(fill_pattern, names(df), value=T)){
        plts[[k]] = one_heatmap(df1 = df, bin_x = bin_x, bin_y = bin_y, bin_num = bin_num, 
                                fill = fill, read_num_cutoff = read_num_cutoff, require_CDS = require_CDS, 
                                one_tx_per_gene=one_tx_per_gene)
        k = k + 1
      }
      plts[["nrow"]] = nrow
      png(file.path(result_dir, paste0(gsub("_|\\$", "", fill_pattern), "_", bin_x, "_", bin_y, "_heatmap_read_num_cutoff=", 
                                       read_num_cutoff, "_bin_num=", bin_num, ".png")), 
          width=width_per_plot*(k-1)/nrow, height=height_per_plot*nrow)
      do.call(grid.arrange, plts)
      dev.off()
    }
  }
}
# Usage:
# multi_heatmap(df = pA.df, bin_x = "CDS_size", bin_y = "UTR_length", bin_nums = c(5, 10),
#               fill_pattern = "_rpm_ratio$", read_num_cutoffs = c(5, 10), require_CDS = T, one_tx_per_gene=T, nrow = 2)
# multi_heatmap(df = pair.pA, bin_x = "CDS_size", bin_y = "aUTR_length", bin_nums = c(5, 10),
#               fill_pattern = "_RED$", read_num_cutoffs = c(2, 5), require_CDS = T, nrow = 2)


#### A function to plot correlation between REDs and output gene lists
## To be polished
RED_cor = function(input_df=pair.pA.combined, read_num_cutoff=5, 
                   compared_idxs=list(c(5, 1), c(5, 3), c(6,2), c(6,4)), # compared_idxs: indexes of REDs to be compared.
                   treatments = c("egfp", "stau1", "koStau1", "WT"),
                   fractions = c("input", "IP", "T", "N", "C", "P"),
                   nrow = 2, file_prefix = "",
                   x_threshold_1 = -0.5, y_threshold_1 = -0.5, # z-scores thresholds used for selecting genes
                   x_threshold_2 = 0.5, y_threshold_2 = 0.5,
                   highlight = F){
  
  input_df = unique(input_df)
  
  plts = list()
  k = 1
  select_genes = data.frame()
  for(idx in compared_idxs){
    print(idx)
    # idx = compared_idxs[[2]]
    red1 = names(input_df)[idx[1]]
    red2 = names(input_df)[idx[2]]
    
    name_pieces1 = strsplit(red1, "_")[[1]]
    name_pieces1 = intersect(name_pieces1, c(treatments, fractions))
    name_pieces2 = strsplit(red2, "_")[[1]]
    name_pieces2 = intersect(name_pieces2, c(treatments, fractions))
    col_candidates = grep("_count_", names(input_df), value=T)
    col_selected = vector()
    for(candidate in col_candidates){
      # candidate = col_candidates[1]
      candidate_pieces = strsplit(candidate, "_")[[1]][1:2]
      if(all(candidate_pieces %in% name_pieces1) | all(candidate_pieces %in% name_pieces2)){
        col_selected = c(col_selected, candidate)
      }
    }
    expressed = rowSums(input_df[, col_selected] >= read_num_cutoff) == length(col_selected)
    df = input_df[expressed, ]
    
    x_z = scale(df[, red1])
    y_z = scale(df[, red2])
    
    # df[, red1] = x_z
    # df[, red2] = y_z
    
    df$highlighted = F
    x_min = min(x_threshold_1, x_threshold_2)
    x_max = max(x_threshold_1, x_threshold_2)
    y_min = min(y_threshold_1, y_threshold_2)
    y_max = max(y_threshold_1, y_threshold_2)
    if(x_threshold_1*y_threshold_1 > 0 & x_threshold_2*y_threshold_2 > 0){
      df[x_z <= x_min & y_z <= y_min, "highlighted"] = T
      df[x_z >= x_max & y_z >= y_max, "highlighted"] = T
    }else if(x_threshold_1*y_threshold_1 < 0 & x_threshold_2*y_threshold_2 < 0){
      df[x_z <= x_min & y_z >= y_max, "highlighted"] = T
      df[x_z >= x_max & y_z <= y_min, "highlighted"] = T
    }
    
    x = red1
    y = red2
    PCC = round(cor(df[, x], df[, y]), 2)
    
    main_title0 = paste0(y, " vs ", x)
    main_title0 = sub("EP_EM", "PM(ER)", main_title0)
    main_title0 = sub("CP_CM", "PM(Cyto)", main_title0)
    main_title0 = sub("EM_CM", "EC(M)", main_title0)
    main_title0 = sub("EP_CP", "EC(P)", main_title0)
    main_title0 = sub("E_C", "EC(M+P)", main_title0)
    
    main_title = paste0(main_title0, "\nRead#>=", read_num_cutoff, " | # of points: ", nrow(df), " | PCC: ", PCC)
    
    df$comparison = main_title0
    
    if(highlight){
      p = ggplot(df, aes_string(x=x, y=y)) + 
        geom_point(aes(color=highlighted, shape=highlighted, alpha=highlighted)) +
        scale_color_manual(values = c("black", "red")) +
        scale_alpha_manual(values = c(0.5, 1)) +
        scale_shape_manual(values = c(1, 19)) +
        geom_smooth(method = "rlm") +
        ggtitle(main_title) +
        theme(plot.title=element_text(hjust = 0.5))
      
      if(sum(df$highlighted) > 200){
        plts[[k]] = p
      }else{
        plts[[k]] = p + geom_text_repel(data=subset(df, highlighted), aes(label = gene_symbol), size = 4, color = "blue")
      }
    }else{
      plts[[k]] = ggplot(df, aes_string(x=x, y=y)) + 
        geom_point(shape=1) +
        geom_smooth(method = "rlm") +
        ggtitle(main_title) +
        theme(plot.title=element_text(hjust = 0.5))
    }
    
    
    k = k + 1
    
    
    df = df[, !grepl("_count_|_3UTR_seq_|^chr$|^strand$|_pA_pos$|gene_id", names(df))]
    df = data.frame(df[, !grepl("description", names(df))], description = df$description)
    if("CDS" %in% names(df)){
      new_cols = c(1:3, match("CDS", names(df)), 4:(match("CDS", names(df))-1), ((match("CDS", names(df))+1):length(names(df))))
      df = df[, new_cols]
    }
    select_genes = rbind(select_genes, subset(df, highlighted))
  }
  
  plts[["nrow"]] = nrow
  png(file.path(result_dir, paste0(file_prefix, "_RED_cor_n=", read_num_cutoff, ".png")), width=650*(k-1)/nrow, height=600*nrow)
  do.call(grid.arrange, plts)
  dev.off()
  
  write.csv(select_genes, 
            file.path(result_dir, paste0(file_prefix, "_selected_genes_RED_cor_n=", read_num_cutoff, ".csv")),
            row.names = F)
}


RED_cor_IRAlu = function(input_df=pair.pA_combined, read_num_cutoff=5, 
                   compared_idxs=list(c(5, 1), c(5, 3), c(6,2), c(6,4)), # compared_idxs: indexes of REDs to be compared.
                   treatments = c("egfp", "stau1", "koStau1", "WT"),
                   fractions = c("input", "IP", "T", "N", "C", "P"),
                   nrow = 2, file_prefix = ""){
  
  input_df = unique(input_df)
  
  plts = list()
  k = 1
  select_genes = data.frame()
  for(idx in compared_idxs){
    print(idx)
    # idx = compared_idxs[[2]]
    red1 = names(input_df)[idx[1]]
    red2 = names(input_df)[idx[2]]
    
    name_pieces1 = strsplit(red1, "_")[[1]]
    name_pieces1 = intersect(name_pieces1, c(treatments, fractions))
    name_pieces2 = strsplit(red2, "_")[[1]]
    name_pieces2 = intersect(name_pieces2, c(treatments, fractions))
    col_candidates = grep("_count_", names(input_df), value=T)
    col_selected = vector()
    for(candidate in col_candidates){
      # candidate = col_candidates[1]
      candidate_pieces = strsplit(candidate, "_")[[1]][1:2]
      if(all(candidate_pieces %in% name_pieces1) | all(candidate_pieces %in% name_pieces2)){
        col_selected = c(col_selected, candidate)
      }
    }
    expressed = rowSums(input_df[, col_selected] >= read_num_cutoff) == length(col_selected)
    df = input_df[expressed, ]
    
    x_z = scale(df[, red1])
    y_z = scale(df[, red2])
    
    # df[, red1] = x_z
    # df[, red2] = y_z
    
    
    
    x = red1
    y = red2
    PCC_dIRAlu_0 = round(cor(df[!df$delta_num_IRAlu_positive, x], df[!df$delta_num_IRAlu_positive, y]), 2)
    PCC_dIRAlu_positive = round(cor(df[df$delta_num_IRAlu_positive, x], df[df$delta_num_IRAlu_positive, y]), 2)
    
    main_title0 = paste0(y, " vs ", x)
    main_title0 = sub("EP_EM", "PM(ER)", main_title0)
    main_title0 = sub("CP_CM", "PM(Cyto)", main_title0)
    main_title0 = sub("EM_CM", "EC(M)", main_title0)
    main_title0 = sub("EP_CP", "EC(P)", main_title0)
    main_title0 = sub("E_C", "EC(M+P)", main_title0)
    
    main_title = paste0(main_title0, "\nRead#>=", read_num_cutoff, " | # of points: ", nrow(df), " | highlighted: ", sum(df$delta_num_IRAlu_positive)," | PCC: ", PCC_dIRAlu_0, ",", PCC_dIRAlu_positive)
    
    df$comparison = main_title0
    
    plts[[k]] = ggplot(df, aes_string(x=x, y=y)) + 
      geom_point(aes(color=delta_num_IRAlu_positive, size=delta_num_IRAlu_positive, alpha=delta_num_IRAlu_positive)) +
      scale_color_manual(values = c("black", "red")) +
      scale_alpha_manual(values = c(0.25, 1)) +
      scale_size_manual(values = c(0.5, 2)) +
      geom_smooth(aes(color = delta_num_IRAlu_positive), method = "rlm") +
      ggtitle(main_title) +
      theme(plot.title=element_text(hjust = 0.5)) 
    
    k = k + 1
    
    # 
    # df = df[, !grepl("_count_|_3UTR_seq_|^chr$|^strand$|_pA_pos$|gene_id", names(df))]
    # df = data.frame(df[, !grepl("description", names(df))], description = df$description)
    # if("CDS" %in% names(df)){
    #   new_cols = c(1:3, match("CDS", names(df)), 4:(match("CDS", names(df))-1), ((match("CDS", names(df))+1):length(names(df))))
    #   df = df[, new_cols]
    # }
    # select_genes = rbind(select_genes, subset(df, highlighted))
  }
  
  plts[["nrow"]] = nrow
  png(file.path(result_dir, paste0(file_prefix, "_RED_cor_n=", read_num_cutoff, ".png")), width=750*(k-1)/nrow, height=600*nrow)
  do.call(grid.arrange, plts)
  dev.off()
  
  # write.csv(select_genes, 
  #           file.path(result_dir, paste0(file_prefix, "_selected_genes_RED_cor_n=", read_num_cutoff, ".csv")),
  #           row.names = F)
}


RED_cor_Alu_IRAlu = function(input_df=pair.pA_combined, read_num_cutoff=5, 
                         compared_idxs=list(c(5, 1), c(5, 3), c(6,2), c(6,4)), # compared_idxs: indexes of REDs to be compared.
                         treatments = c("egfp", "stau1", "koStau1", "WT"),
                         fractions = c("input", "IP", "T", "N", "C", "P"),
                         nrow = 2, file_prefix = ""){
  
  input_df = unique(input_df)
  
  plts = list()
  k = 1
  select_genes = data.frame()
  for(idx in compared_idxs){
    print(idx)
    # idx = compared_idxs[[2]]
    red1 = names(input_df)[idx[1]]
    red2 = names(input_df)[idx[2]]
    
    name_pieces1 = strsplit(red1, "_")[[1]]
    name_pieces1 = intersect(name_pieces1, c(treatments, fractions))
    name_pieces2 = strsplit(red2, "_")[[1]]
    name_pieces2 = intersect(name_pieces2, c(treatments, fractions))
    col_candidates = grep("_count_", names(input_df), value=T)
    col_selected = vector()
    for(candidate in col_candidates){
      # candidate = col_candidates[1]
      candidate_pieces = strsplit(candidate, "_")[[1]][1:2]
      if(all(candidate_pieces %in% name_pieces1) | all(candidate_pieces %in% name_pieces2)){
        col_selected = c(col_selected, candidate)
      }
    }
    expressed = rowSums(input_df[, col_selected] >= read_num_cutoff) == length(col_selected)
    df = input_df[expressed, ]
    
    x_z = scale(df[, red1])
    y_z = scale(df[, red2])
    
    # df[, red1] = x_z
    # df[, red2] = y_z
    
    
    
    x = red1
    y = red2
    # 
    PCC_dAlu_0 = cor.test(df[df$delta_Alu_IRAlu == "dAlu_0", x], df[df$delta_Alu_IRAlu == "dAlu_0", y], method = "pearson")
    PCC_dAlu_positive = cor.test(df[df$delta_Alu_IRAlu == "dAlu_positive", x], df[df$delta_Alu_IRAlu == "dAlu_positive", y], method = "pearson")
    PCC_dIRAlu_positive = cor.test(df[df$delta_Alu_IRAlu == "dIRAlu_positive", x], df[df$delta_Alu_IRAlu == "dIRAlu_positive", y], method = "pearson")

    main_title0 = paste0(y, " vs ", x)
    main_title0 = sub("EP_EM", "PM(ER)", main_title0)
    main_title0 = sub("CP_CM", "PM(Cyto)", main_title0)
    main_title0 = sub("EM_CM", "EC(M)", main_title0)
    main_title0 = sub("EP_CP", "EC(P)", main_title0)
    main_title0 = sub("E_C", "EC(M+P)", main_title0)
    
    main_title = paste0(main_title0, " | Read#>=", read_num_cutoff, "\n# of points: ", 
                        paste(table(df$delta_Alu_IRAlu), collapse = ", "), 
                        " | r: ", paste(round(c(PCC_dAlu_0$estimate, PCC_dAlu_positive$estimate, PCC_dIRAlu_positive$estimate), 2), collapse = ", "),
                        " | p: ", paste(sprintf("%.2E", c(PCC_dAlu_0$p.value, PCC_dAlu_positive$p.value, PCC_dIRAlu_positive$p.value)), collapse = ", "))
    
    df$comparison = main_title0
    
    plts[[k]] = ggplot(df, aes_string(x=x, y=y)) + 
      geom_point(aes(color=delta_Alu_IRAlu, alpha=delta_Alu_IRAlu)) +
      geom_smooth(aes(color = delta_Alu_IRAlu), method = "rlm") +
      scale_color_manual(values = c("black", "limegreen", "tomato"), name = "Alu Type") +
      scale_alpha_manual(values = c(0.1, 1, 1)) +
      facet_grid(~delta_Alu_IRAlu) +
      ggtitle(main_title) +
      theme(plot.title=element_text(hjust = 0.5)) 
    
    k = k + 1
    
  }
  
  plts[["nrow"]] = nrow
  png(file.path(result_dir, paste0(file_prefix, "_RED_cor_n=", read_num_cutoff, ".png")), width=1500*(k-1)/nrow, height=500*nrow)
  do.call(grid.arrange, plts)
  dev.off()
  
}

#### Add description to gene symbols
add_gene_description = function(df, geno){
  # Input requirements:
  # 1. df is a dataframe that must have a gene_symbol column
  # 2. geno must begin with "mm" or "hg"
  # Output: modified df with a gene "description" column
  
  library("biomaRt")
  if(grepl("^mm", geno)){
    ensembl = robust_useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
    bm <- biomaRt::select(ensembl, keys = df$gene_symbol,
                          keytype = "mgi_symbol",
                          columns = c(
                            "mgi_symbol",
                            "description"))
    bm$description = sub(" \\[Source.*\\]", "", bm$description) 
    names(bm) = sub("mgi_", "gene_", names(bm))
    
    # get additional descriptions from an earlier database
    missing_genes = setdiff(df$gene_symbol, bm$gene_symbol)
    if(length(missing_genes) > 0){
      ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host = "jul2015.archive.ensembl.org")
      bm2 <- biomaRt::select(ensembl, keys = missing_genes,
                             keytype = "mgi_symbol",
                             columns = c("mgi_symbol", "description"))
      bm2$description = sub(" \\[Source.*\\]", "", bm2$description) 
      names(bm2) = sub("mgi_", "gene_", names(bm2))
      bm = rbind(bm, bm2)
    }
    
    
  }else if(grepl("^hg", geno)){
    ensembl = robust_useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
    bm <- biomaRt::select(ensembl, keys = df$gene_symbol,
                          keytype = "hgnc_symbol",
                          columns = c("hgnc_symbol", "description"))
    bm$description = sub(" \\[Source.*\\]", "", bm$description) 
    names(bm) = sub("hgnc_", "gene_", names(bm))
    
    # get additional descriptions from an earlier database
    missing_genes = setdiff(df$gene_symbol, bm$gene_symbol)
    if(length(missing_genes) > 0){
      ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl", host = "jul2015.archive.ensembl.org")
      bm2 <- biomaRt::select(ensembl, keys = missing_genes,
                             keytype = "hgnc_symbol",
                             columns = c("hgnc_symbol", "description"))
      bm2$description = sub(" \\[Source.*\\]", "", bm2$description) 
      names(bm2) = sub("hgnc_", "gene_", names(bm2))
      bm = rbind(bm, bm2)
    }
  }
  # Deal with gene_symbols with >1 descriptions
  deduplicated_bm = data.frame()
  for(this_symbol in bm$gene_symbol[duplicated(bm$gene_symbol)]){
    duplicated_bm = bm[bm$gene_symbol == this_symbol, ]
    # guess gene symbol based on the description
    duplicated_bm$gussed_symbol = gsub("(.)[^\\w]+?(\\s|\\b)", "\\1", duplicated_bm$description)
    duplicated_bm$gussed_symbol = gsub("\\-|\\/|\\(|\\)|\\:|\\|\\+,", "",  duplicated_bm$gussed_symbol)
    # calculate string distance between gene symbol and guessed gene symbol
    string_distance = vector("numeric")
    for(i in 1:nrow(duplicated_bm)){
      string_distance[i] = adist(duplicated_bm$gene_symbol[i], duplicated_bm$gussed_symbol[i], ignore.case = T)
    }
    # pick the description with the closest distance
    deduplicated_bm = rbind(deduplicated_bm, duplicated_bm[which.min(string_distance),])
  }
  deduplicated_bm$gussed_symbol = NULL
  deduplicated_bm = deduplicated_bm[!duplicated(deduplicated_bm$gene_symbol), ] # in case the above method faild to pick the best description
  bm = rbind(bm[!bm$gene_symbol %in% bm$gene_symbol[duplicated(bm$gene_symbol)], ], deduplicated_bm)
  
  merge(df, bm, by = "gene_symbol", all.x = T, sort = F)
}

#### A function for extracting REDs from multiple folders under root_dir
extract_RED = function(root_dir = "~/projects", project = "*", exclude = "", RED_type, gene_symbols, geno = "mm9"){
  # project = "stress"
  # RED_type = "H_vs_M"
  # gene_symbols = "GPX3, DIO2, SELENOK, SELENOS, GPX2, GPX6, TXNRD1, SELENOH, SELENOM, MSRB1"
  # gene_symbols = "GPX3 DIO2   SELENOK SELENOS GPX2 GPX6 TXNRD1 SELENOH SELENOM         MSRB1"
  # gene_symbols = "GPX3, DIO2.   SELENOK; SELENOS GPX2 GPX6 TXNRD1 SELENOH SELENOM         MSRB1"
  # gene_symbols = c("GPX3", "DIO2", "SELENOK")
  if(length(gene_symbols) == 1){
    gene_symbols = strsplit(gene_symbols, ";|\\.|,|\\s+")[[1]]
    gene_symbols = gene_symbols[gene_symbols != ""]
  }
  gene_symbols = Hmisc::capitalize(tolower(gene_symbols))
  
  output = data.frame()
  fnames = Sys.glob(file.path(root_dir, project, "result", "*", "*", paste0("blue_red_*", RED_type, "*.csv")))
  for(fname in fnames){
    if(exclude == "" | (exclude != "" && !grepl(exclude, fname))){
      # fname = fnames[1]
      df = read.csv(fname, as.is = T)
      cat("Reading ", fname, "\n")
      df = subset(df, gene_symbol %in% gene_symbols, c("gene_symbol", "proximal_pA", "distal_pA", 
                                                       "UTR_length_proximal", "UTR_length_distal",
                                                       grep(paste0(RED_type,".+_(Prox|Dis)$"), names(df), value=T),
                                                       "fisher", "colors"))
      if(nrow(df)){
        #print(dim(df))
        RED_name = grep(paste0(RED_type,".+_Prox$"), names(df), value=T)
        RED_name = sub("_Prox", "", RED_name)
        df$RED_name = RED_name
        
        df$RED_value = df[, paste0(RED_name,"_Dis")] - df[, paste0(RED_name,"_Prox")]
        
        #df$fname = fname
        
        output = rbind(output, df[, -grep(".+_(Prox|Dis)$", names(df))])
      }
    }
    
  }
  output$aUTR_size = output$UTR_length_distal - output$UTR_length_proximal
  output = output[, c("gene_symbol", "proximal_pA", "distal_pA", "UTR_length_proximal", "UTR_length_distal", "aUTR_size",
                "RED_name", "RED_value", "fisher", "colors")]
  output$fisher = sprintf("%.2E", output$fisher)
  names(output) = sub("fisher", "p-value", names(output))
  names(output) = sub("colors", "color_in_scatter_plot", names(output))
  
  output = output[order(output$gene_symbol, output$RED_name),]
  rownames(output) = NULL
  output$RED_value = round(output$RED_value, 2)
  output = unique(output)
  output = add_gene_description(output, geno)
  
  output
}
# Test: 
# df1 = extract_RED(project = "stress", RED_type = "P_vs_M", gene_symbols = "GPX3, DIO2, SELENOK, SELENOS, GPX2, GPX6, TXNRD1, SELENOH, SELENOM, MSRB1")
# df2 = extract_RED(project = "stress", RED_type = "H_vs_M", gene_symbols = "GPX3, DIO2, SELENOK, SELENOS, GPX2, GPX6, TXNRD1, SELENOH, SELENOM, MSRB1")
# df = rbind(df1, df2)


#### A function for extracting REDs from multiple folders under root_dir
extract_median_RED = function(root_dir = "~/projects", project = "*", exclude = "", RED_type){
  # project = "stress"
  # exclude = "190130_hg19"
  # RED_type = "C2C12D_C2C12P"
  output = vector()
  fnames = Sys.glob(file.path(root_dir, project, "result", "*", "*", "pair.pA.csv"))
  for(fname in fnames){
    if(exclude == "" | (exclude != "" && !grepl(exclude, fname))){
      # fname = fnames[1]
      cat("Reading ", fname, "\n")
      df = read.csv(fname, as.is = T)
      df = df[, grep("_red$|_RED$", names(df), value=T)]
      df = df[, grep(RED_type, names(df), value=T)]
      if(!is.null(ncol(df)) && ncol(df) > 0 && !is.null(nrow(df)) && nrow(df) > 50){
        output = c(output, apply(df, 2, sd))
      }
    }
  }
  as.data.frame(output[grep(RED_type, names(output))])
}
# Test: 
# df1 = extract_RED(project = "stress", RED_type = "P_vs_M", gene_symbols = "GPX3, DIO2, SELENOK, SELENOS, GPX2, GPX6, TXNRD1, SELENOH, SELENOM, MSRB1")
# df2 = extract_RED(project = "stress", RED_type = "H_vs_M", gene_symbols = "GPX3, DIO2, SELENOK, SELENOS, GPX2, GPX6, TXNRD1, SELENOH, SELENOM, MSRB1")
# df = rbind(df1, df2

#### A function for extracting gene expression fold changes from multiple folders under root_dir
extract_gene_FC = function(root_dir = "~/projects", project = "*", exclude = "", 
                           FC_type, gene_symbols, geno = "mm9",
                           output_fname = paste0(FC_type, "_gene_FC.xlsx")){
  require(openxlsx)
  # project = "stress"
  # FC_type = "AS_vs_NT_T"
  # gene_symbols = "GPX3, DIO2, SELENOK, SELENOS, GPX2, GPX6, TXNRD1, SELENOH, SELENOM, MSRB1"
  # gene_symbols = "GPX3 DIO2   SELENOK SELENOS GPX2 GPX6 TXNRD1 SELENOH SELENOM         MSRB1"
  # gene_symbols = "GPX3, DIO2.   SELENOK; SELENOS GPX2 GPX6 TXNRD1 SELENOH SELENOM         MSRB1"
  # gene_symbols = c("GPX3", "DIO2", "SELENOK")
  if(length(gene_symbols) == 1){
    gene_symbols = strsplit(gene_symbols, ";|\\.|,|\\s+")[[1]]
    gene_symbols = gene_symbols[gene_symbols != ""]
  }
  gene_symbols = Hmisc::capitalize(tolower(gene_symbols))
  
  output = createWorkbook()
  s = 1
  #root_dir = sub("~", "/home/dinghai", root_dir)
  require(openxlsx)
  fnames = Sys.glob(file.path(root_dir, project, "result", "*", "*", paste0("Gene-list.xlsx")))
  for(fname in fnames){
    if(exclude == "" | (exclude != "" && !grepl(exclude, fname))){
      # fname = fnames[1]
      cat("Reading ", fname, "\n")
      df = read.xlsx(fname)
      names(df) = sub("Gene_sym", "gene_symbol", names(df))
      df = subset(df, gene_symbol %in% gene_symbols, c("gene_symbol", grep(paste0("Gene_(FC|SS)_.*", FC_type), names(df), value=T)))
      if(nrow(df) && ncol(df) >= 3){
        df$project = sub(".+projects/(.+?)/.+", "\\1", fname)
        df$experiment = sub(".+result/(.+?)/.+", "\\1", fname)
        df = add_gene_description(df, geno)
        addWorksheet(output, sheetName = s)
        writeData(output, s, df)
        s = s + 1
      }
    }
  }
  saveWorkbook(output, output_fname, overwrite = T)
}
# # Testing:
# extract_gene_FC(root_dir = "~/projects", project = "stress", exclude = "", FC_type = "AS_vs_NT_T", geno = "mm9",
#                 gene_symbols = "GPX3, DIO2, SELENOK, SELENOS, GPX2, GPX6, TXNRD1, SELENOH, SELENOM, MSRB1",
#                 output_fname = "AS_vs_NT_T_gene_FC.xlsx")

####
create_3UTR_from_pAs = function(pA.df){
  required_columns = c("region", "gene_symbol", "pAid", "cds_start", "cds_end")
  if(!all(required_columns %in% names(pA.df))) {
    stop(cat("The input dataframe must contain the following columns: ", paste(required_columns, collapse=", ")))
  }
  df = subset(pA.df, region == "3UTR", c("gene_symbol", "pAid", "cds_start", "cds_end"))
  df$chr = sub("(chr.+)[+-]\\d+$", "\\1", df$pAid)
  df$strand = sub("chr.+([+-])\\d+$", "\\1", df$pAid)
  df$pos = as.numeric(sub("chr.+[+-](\\d+)$", "\\1", df$pAid))
  df$start = as.numeric(ifelse(df$strand == "+", df$cds_end + 1, df$pos))
  df$end = as.numeric(ifelse(df$strand == "+", df$pos, df$cds_end - 1))
  df = df[df$end >= df$start, ]
  df[, c("pos", "cds_end", "cds_start")] = NULL
  df = df[complete.cases(df),]
  
  makeGRangesFromDataFrame(df, keep.extra.columns=T)
}

####
create_cUTR_aUTR_from_pAs = function(pair.pA, pA.df){
  require(GenomicRanges)
  
  pair.pA = merge(pair.pA, unique(na.omit(pA.df[, c("gene_symbol", "cds_end")])), sort = F)
  df = pair.pA[, c("gene_symbol", "distal_pA", "cds_end")]
  df$chr = sub("(chr.+)[+-]\\d+$", "\\1", df$distal_pA)
  df$strand = sub("chr.+([+-])\\d+$", "\\1", df$distal_pA)
  df$pos = sub("chr.+[+-](\\d+)$", "\\1", df$distal_pA)
  df$start = ifelse(df$strand == "+", df$cds_end + 1, df$pos)
  df$end = ifelse(df$strand == "+", df$pos, df$cds_end - 1)
  df = df[df$end >= df$start, ]
  df$pos = NULL
  df$pAid = df$distal_pA
  df$distal_pA = NULL
  # dUTR: UTR using the "Distal" pA
  dUTR = makeGRangesFromDataFrame(df, keep.extra.columns=T)
  
  df = pair.pA[, c("gene_symbol", "proximal_pA", "cds_end")]
  df$chr = sub("(chr.+)[+-]\\d+$", "\\1", df$proximal_pA)
  df$strand = sub("chr.+([+-])\\d+$", "\\1", df$proximal_pA)
  df$pos = sub("chr.+[+-](\\d+)$", "\\1", df$proximal_pA)
  df$start = ifelse(df$strand == "+", df$cds_end + 1, df$pos)
  df$end = ifelse(df$strand == "+", df$pos, df$cds_end - 1)
  df = df[df$end >= df$start, ]
  df$pos = NULL
  df$pAid = df$proximal_pA
  df$proximal_pA = NULL
  cUTR = makeGRangesFromDataFrame(df, keep.extra.columns=T)
  
  df = pair.pA[, c("gene_symbol", "cds_end", "proximal_pA", "distal_pA")]
  df$chr = sub("(chr.+)[+-]\\d+$", "\\1", df$proximal_pA)
  df$strand = sub("chr.+([+-])\\d+$", "\\1", df$proximal_pA)
  df$pos1 = sub("chr.+[+-](\\d+)$", "\\1", df$proximal_pA)
  df$pos2 = sub("chr.+[+-](\\d+)$", "\\1", df$distal_pA)
  df$start = ifelse(df$strand == "+", df$pos1, df$pos2)
  df$end = ifelse(df$strand == "+", df$pos2, df$pos1)
  df = df[df$end >= df$start, ]
  df$pos1 = NULL
  df$pos2 = NULL
  # df$distal_pA = NULL
  # df$proximal_pA = NULL
  aUTR = makeGRangesFromDataFrame(df, keep.extra.columns=T)
  aUTR = aUTR[aUTR$gene_symbol %in% cUTR$gene_symbol]
  
  cUTR = cUTR[cUTR$pAid %in% aUTR$proximal_pA]
  dUTR = dUTR[dUTR$pAid %in% aUTR$distal_pA]
  
  result = list(cUTR, aUTR, dUTR)
  names(result) = c("cUTR", "aUTR", "dUTR")
  
  return(result)
}

####
create_conserved_cUTR_aUTR_from_polyA_DB = function(geno = "mm9"){
  # geno can only be "mm9" or "hg19"
  require(dplyr)
  require(tidyr)
  require(GenomicRanges)
  
  # Read data
  if(geno == "mm9"){
    # cds_ends = sapply(cds, function(x) ifelse(strand(x)[1] == "-", start(x)[1], 
    #                                           ifelse(strand(x)[1] == "+", end(x)[length(x)], "NA")))
    # write.csv(cds_ends, "../../../fud/mm9.cds_ends.csv")
    cds_ends_df = read.csv("../../../fud/mm9.cds_ends.csv", header=T)
    padb_raw = read.csv("~/projects/fud/lab/mouse.cental.decompact.txt", header=T, as.is=T, sep="\t")
  }else if(geno == "hg19"){
    # cds_ends = sapply(cds, function(x) ifelse(strand(x)[1] == "-", start(x)[1], 
    #                                           ifelse(strand(x)[1] == "+", end(x)[length(x)], "NA")))
    # write.csv(cds_ends, "../../../fud/hg19.cds_ends.csv")
    cds_ends_df = read.csv("../../../fud/hg19.cds_ends.csv", header=T)
    padb_raw = read.csv("~/projects/fud/lab/human.cental.decompact.txt", header=T, as.is=T, sep="\t") 
  }
  
  # Add cds_end column 
  cds_ends = cds_ends_df[,2]
  names(cds_ends) = cds_ends_df[,1]
  rm(cds_ends_df)
  padb_raw$cds_end = cds_ends[padb_raw$gene_symbol]
  # Only keep genes with cds, otherwise cUTR is undefined. 
  padb_raw = subset(padb_raw, !is.na(cds_end))
  
  # Create aUTR and cUTR GRanges
  padb = padb_raw %>%
    # Only keep pAs in 3'UTRs, otherwise cUTR is undefined.
    filter(coserveLV != "na" & gene_symbol != "na" & region == "3UTR") %>%
    group_by(gene_symbol) %>%
    mutate(num_chrom = length(unique(Chrom)), num_strand = length(unique(Strand))) %>%
    filter(num_chrom == 1 & num_strand ==1) %>%
    mutate(signed_position = ifelse(Strand[1] == "+", Pos, -1*Pos)) %>%
    mutate(signed_cds_end = ifelse(Strand[1] == "+", cds_end, -1*cds_end)) %>%
    # Remove potential 3'UTRs located upstream of the CDS end
    filter(signed_position >= signed_cds_end) %>% 
    mutate(num_pA = n()) %>%
    filter(num_pA > 1) %>%
    arrange(gene_symbol, Pos) %>% 
    filter(Pos %in% c(min(Pos), max(Pos))) %>%
    mutate(Pos_type = c("start", "end")) %>%
    ungroup() %>%
    dplyr::select(gene_symbol, Chrom, Strand, cds_end, Pos, Pos_type) %>%
    spread(key = Pos_type, value = Pos)
  
  names(padb) = c("gene_symbol", "chr", "strand", "cds_end", "end", "start") 
  
  # Some cds_end values are abnormal
  padb_pos = subset(padb, strand == "+" & start >= cds_end)
  # aUTRs start after the proximal pA on "+" strand
  padb_pos$start = padb_pos$start + 1
  padb_neg = subset(padb, strand == "-" & end < cds_end)
  # aUTRs end before the proximal pA on "-" strand
  padb_neg$end = padb_neg$end - 1
  padb = rbind(padb_pos, padb_neg)
  # Create the aUTR GRanges
  aUTR = makeGRangesFromDataFrame(padb, keep.extra.columns =T)
  # Re-calculate start and end for cUTRs
  fist_pA = ifelse(padb$strand == "+", padb$start-1, padb$end+1)
  padb$start = pmin(padb$cds_end+1, fist_pA)
  padb$end = pmax(padb$cds_end-1, fist_pA)
  # Create the cUTR GRanges
  cUTR = makeGRangesFromDataFrame(padb, keep.extra.columns =T)
  
  # Create 3'UTR GRanges
  padb = padb_raw %>%
    # Only keep pAs in 3'UTRs, otherwise cUTR is undefined.
    filter(coserveLV %in% c("H", "H.C", "H.R", "H.R.C", "R", "R.C") & gene_symbol != "na" & region == "3UTR") %>%
    group_by(gene_symbol) %>%
    mutate(num_chrom = length(unique(Chrom)), num_strand = length(unique(Strand))) %>%
    filter(num_chrom == 1 & num_strand ==1) %>%
    mutate(signed_position = ifelse(Strand[1] == "+", Pos, -1*Pos)) %>%
    mutate(signed_cds_end = ifelse(Strand[1] == "+", cds_end, -1*cds_end)) %>%
    # Remove potential 3'UTRs located upstream of the CDS end
    filter(signed_position >= signed_cds_end) %>% 
    # Only keep the last pA
    filter(Pos == max(signed_position)) %>%
    ungroup() %>%
    dplyr::select(gene_symbol, Chrom, Strand, cds_end, Pos) 
  
  names(padb) = c("gene_symbol", "chr", "strand", "cds_end", "last_pA") 
  
  padb$start = pmin(padb$cds_end+1, padb$last_pA)
  padb$end = pmax(padb$cds_end-1, padb$last_pA)
  padb$last_pA = NULL
  # dUTR: UTR using the "Distal" pA
  dUTR = makeGRangesFromDataFrame(padb, keep.extra.columns =T)
  
  result = list(cUTR, aUTR, dUTR)
  names(result) = c("cUTR", "aUTR", "dUTR")
  
  return(result)
}

#### Add gene id using gene_symbol
gene_id_from_symbol = function(geno = "mm9", input_df = pA.df){
  input_df$gene_symbol = as.character(input_df$gene_symbol)
  if(grepl("^mm", geno)){
    require(org.Mm.eg.db)
    sym2id = AnnotationDbi::select(org.Mm.eg.db, keys = input_df$gene_symbol, keytype = "SYMBOL", 
                                   columns = "ENTREZID")
  }else if(grepl("^hg", geno)){
    require(org.Hs.eg.db)
    sym2id = AnnotationDbi::select(org.Hs.eg.db, keys = toupper(input_df$gene_symbol), keytype = "SYMBOL", 
                                   columns = "ENTREZID")
  }
  #sum(!is.na(sym2id$ENTREZID)) 
  names(sym2id) = c("gene_symbol", "gene_id")
  sym2id = unique(sym2id) ################## Never blindly trust downloaded data !!!!!!!!!!
  sym2id = na.omit(sym2id) ################## Never blindly trust downloaded data !!!!!!!!!!
  # some gene symbols have multiple gene id:
  # sum(duplicated(sym2id$gene_symbol))
  sym2id = sym2id[!duplicated(sym2id$gene_symbol), ]
  
  merge(input_df, sym2id, by = "gene_symbol", all.x = T, sort = F) ######### The gene_id is stored as character, not integer!
}

#### Add Alu annotation (hg19)
add_hg19_Alu = function(pA.df, ignore.intron.Alu = T, keep.alu.direction = F){
  # Download the repeat makser file from https://genome.ucsc.edu/cgi-bin/hgTables:
  # Group: Repeats; track: RepeatMasker; table: rmsk; region: genome; output format: all fields from selected table; output file: hg19_rmsk
  
  # Schema for RepeatMasker:
  
  # field   	example	      SQL type	            description
  # bin	      585	          smallint(5) unsigned	Indexing field to speed chromosome range queries.
  # swScore	  1504	        int(10) unsigned	    Smith Waterman alignment score
  # milliDiv	13	          int(10) unsigned	    Base mismatches in parts per thousand
  # milliDel	4	            int(10) unsigned	    Bases deleted in parts per thousand
  # milliIns	13	          int(10) unsigned	    Bases inserted in parts per thousand
  # genoName	chr1	        varchar(255)	        Genomic sequence name
  # genoStart	10000	        int(10) unsigned	    Start in genomic sequence
  # genoEnd	  10468	        int(10) unsigned	    End in genomic sequence
  # genoLeft	-249240153	  int(11)	              -#bases after match in genomic sequence
  # strand	  +	            char(1)	              Relative orientation + or -
  # repName	  (CCCTAA)n	    varchar(255)	        Name of repeat
  # repClass	Simple_repeat	varchar(255)	        Class of repeat
  # repFamily	Simple_repeat	varchar(255)	        Family of repeat
  # repStart	1	            int(11)	              Start (if strand is +) or -#bases after match (if strand is -) in repeat sequence
  # repEnd	  463	          int(11)	              End in repeat sequence
  # repLeft	  0	            int(11)	              -#bases after match (if strand is +) or start (if strand is -) in repeat sequence
  # id	      1	            char(1)	              First digit of id field in RepeatMasker .out file. Best ignored.
  
  
  input_file = "/home/dinghai/projects/fud/hg19_rmsk"
  rmsk = read.table(input_file, header = T, sep = "\t", as.is = T, comment.char = "")
  alu = subset(rmsk, repFamily == "Alu", select = names(rmsk)[c(2:11, 14:16)])
  # "repStart" and "repLeft" have opposite signs. Iif "strand" is "+", "repStart" is positive.
  # Iif "strand" is "+": the sequence between "genoStart" and "genoEnd" on the positive strand is similar to the "repName"
  # Iif "strand" is "-": the sequence between "genoStart" and "genoEnd" on the minus strand is similar to the "repName"
  
  # If a 3'UTR is on the "+" strand and it contains an Alu element defined by ("genoName", "genoStart", "genoEnd", "strand"):
  #    if the Alu "strand" is "+": the Alu in 3'UTR is "S" (sense)
  #    if the Alu "strand" is "-": the Alu in 3'UTR is "AS" (anti-sense)
  # If a 3'UTR is on the "-" strand and it contains an Alu element defined by ("genoName", "genoStart", "genoEnd", "strand"):
  #    if the Alu "strand" is "+": the Alu in 3'UTR is "AS" (anti-sense)
  #    if the Alu "strand" is "-": the Alu in 3'UTR is "S" (sense)
  
  alu$alu.id = 1:nrow(alu)
  
  ############# Map Alu to 3'UTRs
  require(GenomicRanges)
  pA.gr = create_3UTR_from_pAs(pA.df)
  # utr3 = subset(pA.df, region == "3UTR" & cds_end > 0, c("gene_symbol", "pAid", "chr", "strand", "cds_start", "cds_end", "pA_pos"))
  # idx = ifelse(utr3$strand == "+", utr3$cds_end < utr3$pA_pos, utr3$cds_end > utr3$pA_pos)
  # utr3 = utr3[idx,]
  # idx = ifelse(utr3$strand == "+", utr3$cds_end > utr3$cds_start, utr3$cds_end < utr3$cds_start)
  # utr3 = utr3[idx,]
  # 
  # pA.pos = makeGRangesFromDataFrame(subset(utr3, strand == "+"), keep.extra.columns = T, start.field="cds_end", end.field="pA_pos")
  # # Wrong:
  # # pA.neg = makeGRangesFromDataFrame(subset(utr3, strand == "-"), keep.extra.columns = T, end.field="cds_end", start.field="pA_pos")
  # pA.neg = makeGRangesFromDataFrame(subset(utr3, strand == "-"), keep.extra.columns = T, end.field="cds_start", start.field="pA_pos")
  # pA.gr = c(pA.pos, pA.neg)
  
  alu.s = alu[, c("genoName", "genoStart", "genoEnd", "strand", "repName", "repStart", "repEnd", "alu.id")]
  names(alu.s) = c("chr", "start", "end", "strand", "repName", "repStart", "repEnd", "alu.id")
  alu.as = alu.s
  alu.as$strand = ifelse(alu.as$strand == "+", "-", "+")
  
  alu.s.gr = makeGRangesFromDataFrame(alu.s, keep.extra.columns = T)
  alu.as.gr = makeGRangesFromDataFrame(alu.as, keep.extra.columns = T)
  
  # Remove alu in introns first
  if(ignore.intron.Alu){
    txdb = loadDb("../../../fud/hg19.refGene.txdb.sqlite")
    introns = intronsByTranscript(txdb)
    alu.s.gr = subsetByOverlaps(alu.s.gr, introns, type = "within", invert = T)
    alu.as.gr = subsetByOverlaps(alu.as.gr, introns, type = "within", invert = T)
  }
  
  pA.alu.s = mergeByOverlaps(alu.s.gr, pA.gr, type="within")
  pA.alu.s$direction = "s"
  names(pA.alu.s)[1] = "alu.gr"
  
  pA.alu.as = mergeByOverlaps(alu.as.gr, pA.gr, type="within")
  pA.alu.as$direction = "as"
  names(pA.alu.as)[1] = "alu.gr"
  
  pA.alu = rbind(pA.alu.s, pA.alu.as)
  names(pA.alu) = sub("^pA.gr$", "utr3.gr", names(pA.alu))
  pA.alu$alu.gr = as.character(pA.alu$alu.gr)     
  pA.alu$utr3.gr = as.character(pA.alu$utr3.gr)  
  
  require(dplyr)
  require(tidyr)
  df = as.data.frame(pA.alu) %>% 
    mutate(signed.start = as.numeric(sub(".+?:(\\d+?)-\\d+:([+-])", "\\2\\1", alu.gr))) %>% 
    group_by(gene_symbol, pAid, utr3.gr) %>%
    arrange(signed.start) %>%
    summarize(num_Alu = n(), 
              alu.directions = paste(direction, collapse=";"), 
              repNames = paste(repName, collapse=";"), 
              alu.ids = paste(alu.id, collapse=";"),
              alu.grs = paste(alu.gr, collapse=";")) %>%
    dplyr::filter(num_Alu < 100)
  
  
  df$alu.s = !is.na(df$alu.directions) & !grepl("as", df$alu.directions)
  df$alu.as = !is.na(df$alu.directions) & grepl("as", df$alu.directions) & !grepl("^s", df$alu.directions) & !grepl("[^a]s", df$alu.directions)
  df$alu.as.s = !is.na(df$alu.directions) & grepl("as;s", df$alu.directions)
  df$alu.s.as = !is.na(df$alu.directions) & grepl("[^a]s;as", df$alu.directions) 
  
  require(stringr)
  # df$num_IRAlu = str_count(df$alu.directions, pattern = "as;s") + str_count(df$alu.directions, pattern = "[^a]s;as") - str_count(df$alu.directions, pattern = "as;s;as") -   str_count(df$alu.directions, pattern = "[^a]s;as;s")
  df$num_IRAlu = pmin(str_count(df$alu.directions, pattern = "[^a]s|^s"), str_count(df$alu.directions, pattern = "as"))
  IRAlu = df$alu.as.s | df$alu.s.as
  
  df$Alu_type = ifelse(IRAlu, "IRAlu", df$num_Alu)
  
  df$Alu_type[!df$Alu_type %in% c("1", "2", "IRAlu")] = ">=3"
  
  if(keep.alu.direction){
    df = df[, c("pAid", "Alu_type", "num_Alu", "num_IRAlu", "alu.s", "alu.as", "alu.as.s", "alu.s.as", "alu.directions")]
  }else{
    df = df[, c("pAid", "Alu_type", "num_Alu", "num_IRAlu", "alu.s", "alu.as", "alu.as.s", "alu.s.as")]
  }
  pA.df = merge(pA.df, df, by = "pAid", all.x=T)
  pA.df$Alu_type[is.na(pA.df$Alu_type)] = "0"
  pA.df$num_Alu[is.na(pA.df$num_Alu)] = 0
  pA.df$num_IRAlu[is.na(pA.df$num_IRAlu)] = 0
  
  pA.df
}

add_hg19_RIPiT_Stau1 = function(pA.df){
  library(rtracklayer)
  library(AnnotationHub)
  require(dplyr)
  require(tidyr)
  
  # list.files("~/projects/fud/RIPiT/")
  stau1_wt = import("~/projects/fud/RIPiT/GSE52447_Stau1-WT_Sonication_Formaldehyde.bw", format = "BigWig") # hg18
  # summary(stau1_wt$score)
  # hist(log10(stau1_wt$score), breaks=50)
  
  stau1_mt = import("~/projects/fud/RIPiT/GSE52447_Stau1-MUT_Sonication_Formaldehyde.bw", format = "BigWig") # hg18
  # summary(stau1_mt$score)
  # hist(log10(stau1_mt$score), breaks=50)
  
  ahub <- AnnotationHub()
  # table(ahub$rdataclass)
  ahub.chain <- subset(ahub, rdataclass == "ChainFile" & species == "Homo sapiens")
  # query(ahub.chain, c("hg18", "hg19"))
  
  chain <- ahub.chain[ahub.chain$title == "hg18ToHg19.over.chain.gz"]
  chain <- chain[[1]]
  stau1_wt.hg19 <- liftOver(stau1_wt, chain)
  stau1_wt.hg19 = unlist(stau1_wt.hg19)
  stau1_mt.hg19 <- liftOver(stau1_mt, chain)
  stau1_mt.hg19 = unlist(stau1_mt.hg19)
  
  UTR3 = create_3UTR_from_pAs(pA.df)
  
  UTR.wt = as.data.frame(mergeByOverlaps(UTR3, stau1_wt.hg19)[, c("pAid", "score")]) %>%
    group_by(pAid) %>%
    summarise(stau1WT_RIPiT_count = sum(score))
  UTR.mt = as.data.frame(mergeByOverlaps(UTR3, stau1_mt.hg19)[, c("pAid", "score")]) %>%
    group_by(pAid) %>%
    summarise(stau1MT_RIPiT_count = sum(score))
  
  pA.df = merge(pA.df, UTR.wt, all.x = T, sort = F)  
  pA.df = merge(pA.df, UTR.mt, all.x = T, sort = F)   
  
  pA.df[is.na(pA.df$stau1MT_RIPiT_count), "stau1MT_RIPiT_count"] = 0
  pA.df[is.na(pA.df$stau1WT_RIPiT_count), "stau1WT_RIPiT_count"] = 0
  
  pA.df %<>% replace_na(list(stau1MT_RIPiT_count = 0, stau1WT_RIPiT_count = 0)) %>%
    mutate(RIPiT_WT_vs_MT_rpm_ratio = 
             log2((stau1WT_RIPiT_count/sum(stau1WT_RIPiT_count))/(stau1MT_RIPiT_count/sum(stau1MT_RIPiT_count)))
    )
  
  pA.df
}


#### Length-controlled CDF of RED grouped by "group_by"
aUTR_size_controlled_cdf = function(pair.pA, group_by, bg_level = FALSE, red_pattern, file_name, result_dir = result_dir){
  # group_by = "aUTR_cUTR_duplex"
  # red_pattern = "stau1_vs_egfp_IP_mean_RED"
  # n_row = 2
  # file_name = "aUTR_size_controlled_dSBI_grouped_by_aUTR_cUTR_duplex.png"
  
  plts = list()
  k = 1
  n_row = length(grep(red_pattern, names(pair.pA)))
  for(red in grep(red_pattern, names(pair.pA), value=T)){
    # red = "stau1_vs_egfp_IP_mean_RED"
    df = filterRED(input_df = pair.pA, red_type = red, read_num_cutoff=5)
    df[, group_by] = factor(df[, group_by])
    df[, group_by] = relevel(df[, group_by], ref = paste0(bg_level))
    
    df0 = df[df[, group_by] == bg_level, ]
    
    for(other_group in setdiff(levels(df[, group_by]), bg_level)){
      # other_group = "TRUE"
      sub.df = df[df[, group_by] == other_group, ]
      UTR_breaks = quantile(sub.df$aUTR_length, probs = seq(0, 1, 0.1))
      fdr = vector()
      for(t in 1:100){
        idx = vector()
        for(i in 2:(length(UTR_breaks)-1)){
          # i = 2
          # print(UTR_breaks[i-1])
          # print(UTR_breaks[i])
          candidate_idx = which(df0$aUTR_length >= UTR_breaks[i-1] & df0$aUTR_length < UTR_breaks[i])
          # summary(df0[candidate_idx, ]$UTR_length)
          idx = c(idx, sample(candidate_idx, sum(sub.df$aUTR_length >= UTR_breaks[i-1] & sub.df$aUTR_length < UTR_breaks[i]), replace=T))
        }
        i = length(UTR_breaks) # When i is reaches length(UTR_breaks), the ranges are slightly different
        candidate_idx = which(df0$aUTR_length >= UTR_breaks[i-1] & df0$aUTR_length <= UTR_breaks[i])
        idx = c(idx, sample(candidate_idx, sum(sub.df$aUTR_length >= UTR_breaks[i-1] & sub.df$aUTR_length <= UTR_breaks[i]), replace=T))
        # mixed.df = rbind(df0[idx, ], sub.df)[, c("other_group", "WT_C_N_2_rpm_ratio")]
        #sampled.df = cbind(sampled.df, sort(df0[idx, "WT_C_N_2_rpm_ratio"]))
        if(t == 1){
          # p = ggplot(mixed.df, aes(x=WT_C_N_2_rpm_ratio, color = other_group)) + stat_ecdf()
          sampled.m = sort(df0[idx, red])
        }else{
          # p = p + stat_ecdf(aes(x=WT_C_N_2_rpm_ratio, color = other_group), data = mixed.df)
          sampled.m = cbind(sampled.m, sort(df0[idx, red]))
        }
        fdr = c(fdr, 1-mean(sort(sub.df[,red]) >= sort(df0[idx, red])))
      }
      # p
      sampled.m = t(apply(sampled.m, 1, quantile, probs=c(0.05, 0.5, 0.95)))
      sampled.df = as.data.frame(sampled.m) %>%
        gather(Value_type, Value) %>%
        mutate(Value_type = paste0("Sampled ", Value_type)) %>%
        mutate(group = bg_level) 
      mixed.df = data.frame(Value_type = "Observed", Value = sub.df[, red], group = other_group)
      mixed.df$group = as.character(mixed.df$group)
      mixed.df = rbind(mixed.df, sampled.df)
      # cols = c(paste0(bg_level) = "black", paste0(other_group) = "red")
      w.p = round(wilcox.test(sub.df[, red], sampled.m[, 2])$p.value, 3)
      names(mixed.df) = sub("group", group_by, names(mixed.df)) 
      plts[[k]] = ggplot(mixed.df, aes_string(x="Value", linetype="Value_type", size="Value_type", color = group_by)) + 
        stat_ecdf() +
        scale_color_manual(values = c("black", "red")) +
        scale_linetype_manual(values = c(1, 2, 1, 2)) +
        scale_size_manual(values = c(2, 1, 1, 1)) +
        xlab(paste0("Size-controlled RED:", red)) +
        ylab("Cumulative Fraction") +
        ggtitle(paste0("# of genes: ", sum(mixed.df[, group_by] == other_group), "; Mean FDR: ", round(mean(fdr),2), "\nWilcox.p-val: ", w.p))
      
      k = k + 1
    }
  }
  n_col = k - 1
  plts[["ncol"]] = n_col
  png(file.path(result_dir, file_name), 750*n_col, 600)
  do.call(grid.arrange, plts)
  dev.off()
  
}


#### Length-controlled CDF of RED or RE, grouped by "group_by"
size_controlled_cdf = function(df = pair.pA, group_by, bg_level, x_pattern = "stau1_vs_egfp_IP_mean_RED", 
                               size_by = "aUTR_length", file_name, result_dir = result_dir){
  # group_by = "aUTR_cUTR_duplex"
  # x_pattern = "stau1_vs_egfp_IP_mean_RED"
  # n_row = 2
  # file_name = "aUTR_size_controlled_dSBI_grouped_by_aUTR_cUTR_duplex.png"
  
  plts = list()
  k = 1
  n_row = length(grep(x_pattern, names(df)))
  for(x_string in grep(x_pattern, names(df), value=T)){
    # x_string = "CP_vs_CM_mean_dRED"
    if(grepl("_RED$", x_string)){
      df = filterRED(input_df = df, red_type = x_string, read_num_cutoff=5)
    }
    if(grepl("(_RE$)|(_rpm_ratio$)", x_string)){
      df = filterRE(input_df = df, re_type = x_string, read_num_cutoff=5)
    }
    df = as.data.frame(df)
    df[, group_by] = factor(df[, group_by])
    df[, group_by] = relevel(df[, group_by], ref = paste0(bg_level))
    
    df0 = df[df[, group_by] == bg_level, ]
    
    for(other_group in setdiff(levels(df[, group_by]), bg_level)){
      # other_group = "[0.4186,8.1234]"
      sub.df = df[df[, group_by] == other_group, ]
      UTR_breaks = quantile(sub.df[, size_by], probs = seq(0, 1, 0.1))
      fdr = vector()
      for(t in 1:1000){
        idx = vector()
        for(i in 2:(length(UTR_breaks)-1)){
          # i = 2
          # print(UTR_breaks[i-1])
          # print(UTR_breaks[i])
          candidate_idx = which(df0[, size_by] >= UTR_breaks[i-1] & df0[, size_by] < UTR_breaks[i])
          # summary(df0[candidate_idx, ]$UTR_length)
          if(length(candidate_idx)){
            idx = c(idx, sample(candidate_idx, sum(sub.df[, size_by] >= UTR_breaks[i-1] & sub.df[, size_by] < UTR_breaks[i]), replace=T))
          }else{
            sub.df = sub.df[(sub.df[, size_by] < UTR_breaks[i-1]) | (sub.df[, size_by] >= UTR_breaks[i]), ] 
          }
        }
        i = length(UTR_breaks) # When i is reaches length(UTR_breaks), the ranges are slightly different
        candidate_idx = which(df0[, size_by] >= UTR_breaks[i-1] & df0[, size_by] <= UTR_breaks[i])
        if(length(candidate_idx)){
          idx = c(idx, sample(candidate_idx, sum(sub.df[, size_by] >= UTR_breaks[i-1] & sub.df[, size_by] <= UTR_breaks[i]), replace=T))
        }else{
          sub.df = sub.df[sub.df[, size_by] < UTR_breaks[i-1], ] 
        }
        
        # mixed.df = rbind(df0[idx, ], sub.df)[, c("other_group", "WT_C_N_2_rpm_ratio")]
        # sampled.df = cbind(sampled.df, sort(df0[idx, "WT_C_N_2_rpm_ratio"]))
        if(t == 1){
          # p = ggplot(mixed.df, aes(x=WT_C_N_2_rpm_ratio, color = other_group)) + stat_ecdf()
          sampled.m = sort(df0[idx, x_string])
        }else{
          # p = p + stat_ecdf(aes(x=WT_C_N_2_rpm_ratio, color = other_group), data = mixed.df)
          sampled.m = cbind(sampled.m, sort(df0[idx, x_string]))
        }
        fdr = c(fdr, 1-mean(sort(sub.df[,x_string]) >= sort(df0[idx, x_string])))
      }
      # p
      sampled.m = t(apply(sampled.m, 1, quantile, probs=c(0.05, 0.5, 0.95)))
      sampled.df = as.data.frame(sampled.m) %>%
        gather(Value_type, Value) %>%
        mutate(Value_type = paste0("Sampled ", Value_type)) %>%
        mutate(group = bg_level) 
      mixed.df = data.frame(Value_type = "Observed", Value = sub.df[, x_string], group = other_group)
      mixed.df$group = as.character(mixed.df$group)
      mixed.df = rbind(mixed.df, sampled.df)
      # cols = c(paste0(bg_level) = "black", paste0(other_group) = "x_string")
      w.p = sprintf("%.2E", wilcox.test(sub.df[, x_string], sampled.m[, 2])$p.value)
      names(mixed.df) = sub("group", group_by, names(mixed.df)) 
      mixed.df[, group_by] = factor(mixed.df[, group_by], levels = c(bg_level, other_group))
      plts[[k]] = ggplot(mixed.df, aes_string(x="Value", linetype="Value_type", size="Value_type", color = group_by)) + 
        stat_ecdf() +
        scale_color_manual(values = c("black", "red")) +
        scale_linetype_manual(values = c(1, 2, 1, 2)) +
        scale_size_manual(values = c(2, 1, 1, 1)) +
        xlab(paste0("Size-controlled: ", x_string)) +
        ylab("Cumulative Fraction") +
        ggtitle(paste0("# of genes: ", sum(mixed.df[, group_by] == other_group), "; Mean FDR: ", round(mean(fdr),2), "\nWilcox.p-val: ", w.p))
      
      k = k + 1
    }
  }
  n_col = k - 1
  plts[["ncol"]] = n_col
  png(file.path(result_dir, file_name), 750*n_col, 600)
  do.call(grid.arrange, plts)
  dev.off()
}

#### Length-controlled CDF of RED or RE, grouped by "group_by"
size_controlled_cdf_2 = function(input_df = df2, # Input dataframe
                                 group_by = "C2C12DR_HuR_RE_bin", # A column for grouping data. Each group will have its own CDF curve.
                                 control_group = "[0.0000,0.0372)", # One of the values in "group_by". Data in this group will be plotted as black CDF curves.
                                 value_pattern = "EP_vs_EM_mean_dRED", # A regex pattern for selecting columns (>= 1) used as the x-axis values for CDF curves.  
                                 size_column = "aUTR_length", # Either "aUTR_length" or "UTR_length"
                                 read_num_cutoff = 5, # Minimum read number for filtering data
                                 batch = "all", # Batch of data used for filtering and plotting. If mean RED or RE is plotted, use "all". Otherwise use a batch number.
                                 plot_width = 800, # Width for each small plot
                                 plot_height = 600, # Height for each small plot
                                 file_name = "", 
                                 result_dir = "."){
  
  if(!group_by %in% names(input_df)){
    stop("The 'group_by' column is not in the input data frame!")
  }else(
    if(!control_group %in% input_df[, group_by]){
      stop("The 'control_group' is not a value of the 'group_by' column!")
    }
  )
  
  if(!size_column %in% names(input_df)){
    stop("The 'size_column' is not in the input data frame!")
  }
  
  if(!length(grep(value_pattern, names(input_df)))){
    stop("Failed to select columns using 'value_pattern'!")
  }
  
  # A list container for saving the plots
  plts = list()
  k = 1
  
  # Loop through columns matching value_pattern to make plots
  for(value_column in grep(value_pattern, names(input_df), value=T)){
    # Filter data using read number cutoff
    if(grepl("_RED$", value_column)){
      df = filter_RED(input_df = input_df, red  = value_column, batch = batch, read_num_cutoff = read_num_cutoff)
    }else if(grepl("(_RE$)|(_rpm_ratio$)", value_column)){
      df = filter_RE(input_df = input_df, re = value_column, batch = batch, read_num_cutoff = read_num_cutoff)
    }else{
      df = input_df
      warning(paste0("After filtering ", value_column, ", row number of the input data frame remains the same. \n  You may want to filter the data manually before calling size_controlled_cdf()."))
    }

    # Prepare the group_by column
    df = unique(as.data.frame(df))
    df[, group_by] = factor(df[, group_by])
    df[, group_by] = relevel(df[, group_by], ref = control_group)
    
    # Prepare the control data frame
    df_control = df[df[, group_by] == control_group, ]
    
    # Loop through the test groups and compare each of them to the control group
    for(test_group in setdiff(levels(df[, group_by]), control_group)){
      # Calculate the quantiles of aUTR_length or UTR_length 
      df_test = df[df[, group_by] == test_group, ]
      size_quantiles = quantile(df_test[, size_column], probs = seq(0, 1, 0.05))
      
      # Bootstrap 1000 times
      for(t in 1:1000){
        # Set up a container for sampled indexes
        idx = vector()
        
        # Make a copy of df_test, so that it can be modified in the for-loop below
        df_test_copy = df_test
        
        # Loop through each size range and randomly sample (with replacement) the control data for the size range
        for(i in 2:length(size_quantiles)){
          if(i < length(size_quantiles)){
            candidate_idx = which((df_control[, size_column] >= size_quantiles[i-1]) & (df_control[, size_column] < size_quantiles[i]))
            data_points_in_df_test = sum((df_test[, size_column] >= size_quantiles[i-1]) & (df_test[, size_column] < size_quantiles[i]))
          }else{ # When i is equal to length(size_quantiles), the ranges are slightly different
            candidate_idx = which((df_control[, size_column] >= size_quantiles[i-1]) & (df_control[, size_column] <= size_quantiles[i]))
            data_points_in_df_test = sum((df_test[, size_column] >= size_quantiles[i-1]) & (df_test[, size_column] <= size_quantiles[i]))
          }
          
          # If no data points are found in the size range of the control group, 
          # remove the data points in the size range of the test group (df_test)
          if(!length(candidate_idx)){ # No data points are found in the size range of the control group
            # print(c(t, i))
            if(i < length(size_quantiles)){ # When the size range is not the last size range
              # Remove the data points in the size range of the test group (df_test)
              df_test = df_test[(df_test[, size_column] < size_quantiles[i-1]) | (df_test[, size_column] >= size_quantiles[i]), ] 
            }else{ # When the size range is the last size range
              # Remove the data points in the size range of the test group (df_test)
              df_test = df_test[df_test[, size_column] < size_quantiles[i-1], ]
            }
          }else{ 
            # Simply save the sampled index
            idx = c(idx, sample(candidate_idx, data_points_in_df_test, replace=T))
          }
        }
        
        # Retrieve the sampled data using sampled indexes, sort each batch of sampled data, and put them in a column of a matrix
        if(t == 1){
          sampled_value_matrix = sort(df_control[idx, value_column]) # sorted
        }else{
          sampled_value_matrix = cbind(sampled_value_matrix, sort(df_control[idx, value_column])) # sorted
        }
        
        # Retrieve the aUTR or 3'UTR sizes of sampled data
        if(t == 1){
          sampled_size_matrix = sort(df_control[idx, size_column]) # sorted
        }else{
          sampled_size_matrix = cbind(sampled_size_matrix, sort(df_control[idx, size_column])) # sorted
        }
      }
      
      # Calculate row-wise quantiles of the sampled data
      sampled_value_summary = t(apply(sampled_value_matrix, 1, quantile, probs=c(0.05, 0.5, 0.95)))
      sampled_size_summary = t(apply(sampled_size_matrix, 1, quantile, probs=c(0.05, 0.5, 0.95)))
      
      # Compare the observed data with the row-wise median of sampled data. Calculate p value using Wilcoxon rank sum test
      v.p = sprintf("%.2E", wilcox.test(df_test[, value_column], sampled_value_summary[, 2])$p.value)
      s.p = sprintf("%.2E", wilcox.test(df_test[, size_column], sampled_size_summary[, 2])$p.value)
      
      # Compute FDR
      v.fdr = mean(sampled_value_summary[, 2] >= median(df_test[, value_column]))
      v.fdr = min(v.fdr, 1 - v.fdr)
      
      s.fdr = mean(sampled_size_summary[, 2] >= median(df_test[, size_column]))
      s.fdr = min(s.fdr, 1 - s.fdr)
      
      require(tidyr)
      require(dplyr)
      # Reshape the sampled data for plotting
      value_sampled = as.data.frame(sampled_value_summary) %>%
        gather(Value_type, Value) %>%
        mutate(Value_type = paste0("Sampled ", Value_type)) %>%
        mutate(group = control_group) 
      size_sampled = as.data.frame(sampled_size_summary) %>%
        gather(Value_type, Value) %>%
        mutate(Value_type = paste0("Sampled ", Value_type)) %>%
        mutate(group = control_group) 
      
      # Reshape the observed test data for plotting
      value_test = data.frame(Value_type = "Observed", Value = df_test[, value_column], group = test_group)
      value_test$group = as.character(value_test$group)
      size_test = data.frame(Value_type = "Observed", Value = df_test[, size_column], group = test_group)
      size_test$group = as.character(size_test$group)
      
      # Combine the sampled data with observed data in the test group
      value_mixed = rbind(value_test, value_sampled)
      names(value_mixed) = sub("group", group_by, names(value_mixed)) 
      value_mixed[, group_by] = factor(value_mixed[, group_by], levels = c(control_group, test_group))
      
      size_mixed = rbind(size_test, size_sampled)
      names(size_mixed) = sub("group", group_by, names(size_mixed)) 
      size_mixed[, group_by] = factor(size_mixed[, group_by], levels = c(control_group, test_group))
      
      # Plotting CDF
      require(ggplot2)
      require(gridExtra)
      require(ggpubr)
      size_cdf = ggplot(size_mixed, aes_string(x="Value", linetype="Value_type", size="Value_type", color = group_by)) + 
        stat_ecdf() +
        scale_color_manual(values = c("black", "red")) +
        scale_linetype_manual(values = c(1, 2, 1, 2)) +
        scale_size_manual(values = c(2, 1, 1, 1)) +
        scale_x_log10() +
        xlab(paste0(size_column)) +
        ylab("Cumulative Fraction") +
        ggtitle(paste0("# of genes: ", sum(size_mixed[, group_by] == test_group), " | FDR = ", s.fdr,  " | " , "Wilcox.p-val: ", s.p))
      
      value_cdf = ggplot(value_mixed, aes_string(x="Value", linetype="Value_type", size="Value_type", color = group_by)) + 
        stat_ecdf() +
        scale_color_manual(values = c("black", "red")) +
        scale_linetype_manual(values = c(1, 2, 1, 2)) +
        scale_size_manual(values = c(2, 1, 1, 1)) +
        xlab(paste0("Size-controlled ", value_column)) +
        ylab("Cumulative Fraction") +
        ggtitle(paste0("# of genes: ", sum(value_mixed[, group_by] == test_group), " | FDR = ", v.fdr,  " | " , "Wilcox.p-val: ", v.p))
      
      plts[[k]] = ggarrange(size_cdf, value_cdf)
      
      k = k + 1
    }
  }
  
  # Calculate output file name if not provided
  if(file_name == "") {
    file_name = paste(size_column, "controlled", value_pattern, "grouped_by", group_by, 
                      "using", control_group, "as_control.png", sep = "_")
  }
  
  # Calcualte output image size and create the image
  plts[["nrow"]] = length(grep(value_pattern, names(input_df)))
  png(file.path(result_dir, file_name), plot_width*2*(k-1)/plts[["nrow"]], plot_height*plts[["nrow"]])
  do.call(grid.arrange, plts)
  dev.off()
}

#### Length-controlled CDF of RED or RE, grouped by "group_by"
size_controlled_cdf_4 = function(input_df = df2, # Input dataframe
                                 group_by = "C2C12DR_HuR_RE_bin", # A column for grouping data. Each group will have its own CDF curve.
                                 control_group = "[0.0000,0.0372)", # One of the values in "group_by". Data in this group will be plotted as black CDF curves.
                                 value_pattern = "EP_vs_EM_mean_dRED", # A regex pattern for selecting columns (>= 1) used as the x-axis values for CDF curves.  
                                 size_column = "aUTR_length", # Either "aUTR_length" or "UTR_length"
                                 read_num_cutoff = 5, # Minimum read number for filtering data
                                 batch = "all", # Batch of data used for filtering and plotting. If mean RED or RE is plotted, use "all". Otherwise use a batch number.
                                 plot_width = 800, # Width for each small plot
                                 plot_height = 600, # Height for each small plot
                                 file_name = "", 
                                 result_dir = "."){
  
  if(!group_by %in% names(input_df)){
    stop("The 'group_by' column is not in the input data frame!")
  }else(
    if(!control_group %in% input_df[, group_by]){
      stop("The 'control_group' is not a value of the 'group_by' column!")
    }
  )
  
  if(!size_column %in% names(input_df)){
    stop("The 'size_column' is not in the input data frame!")
  }
  
  if(!length(grep(value_pattern, names(input_df)))){
    stop("Failed to select columns using 'value_pattern'!")
  }
  
  # A list container for saving the plots
  plts = list()
  k = 1
  
  # Loop through columns matching value_pattern to make plots
  for(value_column in grep(value_pattern, names(input_df), value=T)){
    # value_column = grep(value_pattern, names(input_df), value=T)[1]
    # Filter data using read number cutoff
    if(grepl("_RED$", value_column)){
      df = filter_RED(input_df = input_df, red  = value_column, batch = batch, read_num_cutoff = read_num_cutoff)
    }else if(grepl("(_RE$)|(_rpm_ratio$)", value_column)){
      df = filter_RE(input_df = input_df, re = value_column, batch = batch, read_num_cutoff = read_num_cutoff)
    }else{
      df = input_df
      warning(paste0("After filtering ", value_column, ", row number of the input data frame remains the same. \n  You may want to filter the data manually before calling size_controlled_cdf()."))
    }
    
    # Prepare the group_by column
    df = unique(as.data.frame(df))
    df[, group_by] = factor(df[, group_by])
    df[, group_by] = relevel(df[, group_by], ref = control_group)
    
    # Prepare the candidate control data, which will be filtered to remove data points using size_column
    df_control_candidate = df[df[, group_by] == control_group, ]
    control_size_range = range(df_control_candidate[, size_column])
    
    # Loop through the test groups and compare each of them to the control group
    for(test_group in setdiff(levels(df[, group_by]), control_group)){
      # test_group = setdiff(levels(df[, group_by]), control_group)[2]
      # Prepare the candidate test data, which will be filtered to remove data points using size_column
      df_test_candidate = df[df[, group_by] == test_group, ]
      test_size_range = range(df_test_candidate[, size_column])
      dim(df_test_candidate)
      
      # Remove data points in the test group beyond the range of control size_column 
      df_test = df_test_candidate[df_test_candidate[, size_column] >= control_size_range[1], ]
      df_test = df_test[df_test[, size_column] <= control_size_range[2], ]
      # dim(df_test)
      # range(df_test[, size_column])
      
      # Remove data points in the control group beyond the range of test size_column
      df_control = df_control_candidate[df_control_candidate[, size_column] >= test_size_range[1], ]
      df_control = df_control[df_control[, size_column] <= test_size_range[2], ]
      # dim(df_control)
      
      # Calculate the quantiles of aUTR_length or UTR_length 
      size_quantiles = quantile(df_test[, size_column], probs = seq(0, 1, 0.1))
      
      df_test = df_test[df_test[, size_column] <= size_quantiles[length(size_quantiles) - 1], ]
      
      # Bootstrap 1000 times
      for(t in 1:1000){
        # Set up a container for sampled indexes
        idx = vector()
        
        # Loop through each size range and randomly sample (with replacement) the control data for the size range
        for(i in 1:(length(size_quantiles) - 2)){
          candidate_idx = which((df_control[, size_column] >= size_quantiles[i+1]) & (df_control[, size_column] < size_quantiles[i+2]))
          data_points_in_df_test = sum((df_test[, size_column] >= size_quantiles[i]) & (df_test[, size_column] < size_quantiles[i+1]))
          
          # If no data points are found in the size range of the control group, 
          # remove the data points in the size range of the test group (df_test)
          if(!length(candidate_idx)){ # No data points are found in the size range of the control group
            df_test = df_test[(df_test[, size_column] < size_quantiles[i]) | (df_test[, size_column] >= size_quantiles[i+1]), ] 
          }else{ 
            # Simply save the sampled index
            idx = c(idx, sample(candidate_idx, data_points_in_df_test, replace=T))
          }
        }
        
        # Retrieve the sampled data using sampled indexes, sort each batch of sampled data, and put them in a column of a matrix
        if(t == 1){
          sampled_value_matrix = sort(df_control[idx, value_column]) # sorted
        }else{
          sampled_value_matrix = cbind(sampled_value_matrix, sort(df_control[idx, value_column])) # sorted
        }
        
        # Retrieve the aUTR or 3'UTR sizes of sampled data
        if(t == 1){
          sampled_size_matrix = sort(df_control[idx, size_column]) # sorted
        }else{
          sampled_size_matrix = cbind(sampled_size_matrix, sort(df_control[idx, size_column])) # sorted
        }
      }
      
      # Calculate row-wise quantiles of the sampled data
      sampled_value_summary = t(apply(sampled_value_matrix, 1, quantile, probs=c(0.05, 0.5, 0.95)))
      sampled_size_summary = t(apply(sampled_size_matrix, 1, quantile, probs=c(0.05, 0.5, 0.95)))
      
      # Compare the observed data with the row-wise median of sampled data. Calculate p value using Wilcoxon rank sum test
      v.p = sprintf("%.2E", wilcox.test(df_test[, value_column], sampled_value_summary[, 2])$p.value)
      s.p = sprintf("%.2E", wilcox.test(df_test[, size_column], sampled_size_summary[, 2])$p.value)
      
      require(tidyr)
      require(dplyr)
      # Reshape the sampled data for plotting
      value_sampled = as.data.frame(sampled_value_summary) %>%
        gather(Value_type, Value) %>%
        mutate(Value_type = paste0("Sampled ", Value_type)) %>%
        mutate(group = control_group) 
      size_sampled = as.data.frame(sampled_size_summary) %>%
        gather(Value_type, Value) %>%
        mutate(Value_type = paste0("Sampled ", Value_type)) %>%
        mutate(group = control_group) 
      
      # Reshape the observed test data for plotting
      value_test = data.frame(Value_type = "Observed", Value = df_test[, value_column], group = test_group)
      value_test$group = as.character(value_test$group)
      size_test = data.frame(Value_type = "Observed", Value = df_test[, size_column], group = test_group)
      size_test$group = as.character(size_test$group)
      
      # Combine the sampled data with observed data in the test group
      value_mixed = rbind(value_test, value_sampled)
      names(value_mixed) = sub("group", group_by, names(value_mixed)) 
      value_mixed[, group_by] = factor(value_mixed[, group_by], levels = c(control_group, test_group))
      
      size_mixed = rbind(size_test, size_sampled)
      names(size_mixed) = sub("group", group_by, names(size_mixed)) 
      size_mixed[, group_by] = factor(size_mixed[, group_by], levels = c(control_group, test_group))
      
      # Plotting CDF
      require(ggplot2)
      require(gridExtra)
      require(ggpubr)
      size_cdf = ggplot(size_mixed, aes_string(x="Value", linetype="Value_type", size="Value_type", color = group_by)) + 
        stat_ecdf() +
        scale_color_manual(values = c("black", "red")) +
        scale_linetype_manual(values = c(1, 2, 1, 2)) +
        scale_size_manual(values = c(2, 1, 1, 1)) +
        scale_x_log10() +
        xlab(paste0(size_column)) +
        ylab("Cumulative Fraction") +
        ggtitle(paste0("# of genes: ", sum(size_mixed[, group_by] == test_group), " | " , "Wilcox.p-val: ", s.p))
      
      value_cdf = ggplot(value_mixed, aes_string(x="Value", linetype="Value_type", size="Value_type", color = group_by)) + 
        stat_ecdf() +
        scale_color_manual(values = c("black", "red")) +
        scale_linetype_manual(values = c(1, 2, 1, 2)) +
        scale_size_manual(values = c(2, 1, 1, 1)) +
        xlab(paste0("Size-controlled ", value_column)) +
        ylab("Cumulative Fraction") +
        ggtitle(paste0("# of genes: ", sum(value_mixed[, group_by] == test_group), " | " , "Wilcox.p-val: ", v.p))
      
      plts[[k]] = ggarrange(size_cdf, value_cdf)
      
      k = k + 1
    }
  }
  
  # Calculate output file name if not provided
  if(file_name == ""){
    file_name = paste(size_column, "controlled", value_pattern, "grouped_by", group_by, 
                      "using", control_group, "as_control.png", sep = "_")
  }
  
  # Calcualte output image size and create the image
  plts[["nrow"]] = length(grep(value_pattern, names(input_df)))
  png(file.path(result_dir, file_name), plot_width*2*(k-1)/plts[["nrow"]], plot_height*plts[["nrow"]])
  do.call(grid.arrange, plts)
  dev.off()
}

#### Length-controlled CDF of RED or RE, grouped by "group_by"
size_controlled_cdf_5 = function(input_df = df2, # Input dataframe
                                 group_by = "C2C12DR_HuR_RE_bin", # A column for grouping data. Each group will have its own CDF curve.
                                 control_group = "[0.0000,0.0372)", # One of the values in "group_by". Data in this group will be plotted as black CDF curves.
                                 value_pattern = "EP_vs_EM_mean_dRED", # A regex pattern for selecting columns (>= 1) used as the x-axis values for CDF curves.  
                                 size_column = "aUTR_length", # Either "aUTR_length" or "UTR_length"
                                 read_num_cutoff = 5, # Minimum read number for filtering data
                                 batch = "all", # Batch of data used for filtering and plotting. If mean RED or RE is plotted, use "all". Otherwise use a batch number.
                                 plot_width = 800, # Width for each small plot
                                 plot_height = 600, # Height for each small plot
                                 file_name = "", 
                                 result_dir = "."){
  
  if(!group_by %in% names(input_df)){
    stop("The 'group_by' column is not in the input data frame!")
  }else(
    if(!control_group %in% input_df[, group_by]){
      stop("The 'control_group' is not a value of the 'group_by' column!")
    }
  )
  
  if(!size_column %in% names(input_df)){
    stop("The 'size_column' is not in the input data frame!")
  }
  
  if(!length(grep(value_pattern, names(input_df)))){
    stop("Failed to select columns using 'value_pattern'!")
  }
  
  # A list container for saving the plots
  plts = list()
  k = 1
  
  # Loop through columns matching value_pattern to make plots
  for(value_column in grep(value_pattern, names(input_df), value=T)){
    # value_column = grep(value_pattern, names(input_df), value=T)[1]
    # Filter data using read number cutoff
    if(grepl("_RED$", value_column)){
      df = filter_RED(input_df = input_df, red  = value_column, batch = batch, read_num_cutoff = read_num_cutoff)
    }else if(grepl("(_RE$)|(_rpm_ratio$)", value_column)){
      df = filter_RE(input_df = input_df, re = value_column, batch = batch, read_num_cutoff = read_num_cutoff)
    }else{
      df = input_df
      warning(paste0("After filtering ", value_column, ", row number of the input data frame remains the same. \n  You may want to filter the data manually before calling size_controlled_cdf()."))
    }
    
    # Prepare the group_by column
    df = unique(as.data.frame(df))
    df[, group_by] = factor(df[, group_by])
    df[, group_by] = relevel(df[, group_by], ref = control_group)
    
    # Prepare the candidate control data, which will be filtered to remove data points using size_column
    df_control_candidate = df[df[, group_by] == control_group, ]
    control_size_range = range(df_control_candidate[, size_column])
    
    # Loop through the test groups and compare each of them to the control group
    for(test_group in setdiff(levels(df[, group_by]), control_group)){
      # test_group = setdiff(levels(df[, group_by]), control_group)[2]
      # Prepare the candidate test data, which will be filtered to remove data points using size_column
      df_test_candidate = df[df[, group_by] == test_group, ]
      test_size_range = range(df_test_candidate[, size_column])
      dim(df_test_candidate)
      
      # Remove data points in the test group beyond the range of control size_column 
      df_test = df_test_candidate[df_test_candidate[, size_column] >= control_size_range[1], ]
      df_test = df_test[df_test[, size_column] <= control_size_range[2], ]
      # dim(df_test)
      # range(df_test[, size_column])
      
      # Remove data points in the control group beyond the range of test size_column
      df_control = df_control_candidate[df_control_candidate[, size_column] >= test_size_range[1], ]
      df_control = df_control[df_control[, size_column] <= test_size_range[2], ]
      # dim(df_control)
      
      # Gradually increase quantile step 
      for(quantile_step in seq(0.1, 0.01, -0.01)){
        # Calculate the quantiles of aUTR_length or UTR_length 
        size_quantiles = quantile(df_test[, size_column], probs = seq(0, 1, quantile_step))
        
        df_test = df_test[df_test[, size_column] <= size_quantiles[length(size_quantiles) - 1], ]
        
        # Bootstrap 1000 times
        for(t in 1:1000){
          # Set up a container for sampled indexes
          idx = vector()
          
          # Loop through each size range and randomly sample (with replacement) the control data for the size range
          for(i in 1:(length(size_quantiles) - 2)){
            candidate_idx = which((df_control[, size_column] >= size_quantiles[i+1]) & (df_control[, size_column] < size_quantiles[i+2]))
            data_points_in_df_test = sum((df_test[, size_column] >= size_quantiles[i]) & (df_test[, size_column] < size_quantiles[i+1]))
            
            # If no data points are found in the size range of the control group, 
            # remove the data points in the size range of the test group (df_test)
            if(!length(candidate_idx)){ # No data points are found in the size range of the control group
              df_test = df_test[(df_test[, size_column] < size_quantiles[i]) | (df_test[, size_column] >= size_quantiles[i+1]), ] 
            }else{ 
              # Simply save the sampled index
              idx = c(idx, sample(candidate_idx, data_points_in_df_test, replace=T))
            }
          }
          
          # Retrieve the sampled data using sampled indexes, sort each batch of sampled data, and put them in a column of a matrix
          if(t == 1){
            sampled_value_matrix = sort(df_control[idx, value_column]) # sorted
          }else{
            sampled_value_matrix = cbind(sampled_value_matrix, sort(df_control[idx, value_column])) # sorted
          }
          
          # Retrieve the aUTR or 3'UTR sizes of sampled data
          if(t == 1){
            sampled_size_matrix = sort(df_control[idx, size_column]) # sorted
          }else{
            sampled_size_matrix = cbind(sampled_size_matrix, sort(df_control[idx, size_column])) # sorted
          }
        }
        
        # Calculate row-wise quantiles of the sampled data
        sampled_value_summary = t(apply(sampled_value_matrix, 1, quantile, probs=c(0.05, 0.5, 0.95)))
        sampled_size_summary = t(apply(sampled_size_matrix, 1, quantile, probs=c(0.05, 0.5, 0.95)))
        
        # Compare the observed data with the row-wise median of sampled data. Calculate p value using Wilcoxon rank sum test
        v.p = wilcox.test(df_test[, value_column], sampled_value_summary[, 2])$p.value
        s.p = wilcox.test(df_test[, size_column], sampled_size_summary[, 2])$p.value
        
        # Compute FDRs
        v.fdr = mean(sampled_value_summary[, 2] >= median(df_test[, value_column]))
        v.fdr = min(v.fdr, 1 - v.fdr)
        
        s.fdr = mean(sampled_size_summary[, 2] >= median(df_test[, size_column]))
        s.fdr = min(s.fdr, 1 - s.fdr)
        
        if(s.p <= 0.05){
          cat("Number of bins:", length(size_quantiles) - 1, "\n")
          cat("p-value for comparing", size_column, ":",  s.p, "\n")
          cat("Skipping ...")
          next
        }
        
        # if(s.fdr <= 0.45){
        #   cat("Number of bins:", length(size_quantiles) - 1, "\n")
        #   cat("FDR of", size_column, ":",  s.fdr, "\n")
        #   cat("Skipping ...\n")
        #   next
        # }
        
        s.p = sprintf("%.2E", s.p)
        v.p = sprintf("%.2E", v.p)
        
        
        
        require(tidyr)
        require(dplyr)
        # Reshape the sampled data for plotting
        value_sampled = as.data.frame(sampled_value_summary) %>%
          gather(Value_type, Value) %>%
          mutate(Value_type = paste0("Sampled ", Value_type)) %>%
          mutate(group = control_group) 
        size_sampled = as.data.frame(sampled_size_summary) %>%
          gather(Value_type, Value) %>%
          mutate(Value_type = paste0("Sampled ", Value_type)) %>%
          mutate(group = control_group) 
        
        # Reshape the observed test data for plotting
        value_test = data.frame(Value_type = "Observed", Value = df_test[, value_column], group = test_group)
        value_test$group = as.character(value_test$group)
        size_test = data.frame(Value_type = "Observed", Value = df_test[, size_column], group = test_group)
        size_test$group = as.character(size_test$group)
        
        # Combine the sampled data with observed data in the test group
        value_mixed = rbind(value_test, value_sampled)
        names(value_mixed) = sub("group", group_by, names(value_mixed)) 
        value_mixed[, group_by] = factor(value_mixed[, group_by], levels = c(control_group, test_group))
        
        size_mixed = rbind(size_test, size_sampled)
        names(size_mixed) = sub("group", group_by, names(size_mixed)) 
        size_mixed[, group_by] = factor(size_mixed[, group_by], levels = c(control_group, test_group))
        
        # Plotting CDF
        require(ggplot2)
        require(gridExtra)
        require(ggpubr)
        size_cdf = ggplot(size_mixed, aes_string(x="Value", linetype="Value_type", size="Value_type", color = group_by)) + 
          stat_ecdf() +
          scale_color_manual(values = c("black", "red")) +
          scale_linetype_manual(values = c(1, 2, 1, 2)) +
          scale_size_manual(values = c(2, 1, 1, 1)) +
          scale_x_log10() +
          xlab(paste0(size_column)) +
          ylab("Cumulative Fraction") +
          ggtitle(paste0("# of genes: ", sum(size_mixed[, group_by] == test_group), " | ", 
                         "# of bins: ", length(size_quantiles) - 1, 
                         " | Wilcox.p = ", s.p, " | FDR = ", round(s.fdr, 3)))
        
        value_cdf = ggplot(value_mixed, aes_string(x="Value", linetype="Value_type", size="Value_type", color = group_by)) + 
          stat_ecdf() +
          scale_color_manual(values = c("black", "red")) +
          scale_linetype_manual(values = c(1, 2, 1, 2)) +
          scale_size_manual(values = c(2, 1, 1, 1)) +
          xlab(paste0("Size-controlled ", value_column)) +
          ylab("Cumulative Fraction") +
          ggtitle(paste0("# of genes: ", sum(value_mixed[, group_by] == test_group), " | " ,
                         "# of bins: ", length(size_quantiles) - 1, 
                         " | Wilcox.p = ", v.p, " | FDR = ", round(v.fdr, 3)))
        
        break
        
      }
      
      plts[[k]] = ggarrange(size_cdf, value_cdf)
      
      k = k + 1
    }
  }
  
  # Calculate output file name if not provided
  if(file_name == ""){
    file_name = paste(size_column, "controlled", value_pattern, "grouped_by", group_by, 
                      "using", control_group, "as_control.png", sep = "_")
  }
  
  # Calcualte output image size and create the image
  plts[["nrow"]] = length(grep(value_pattern, names(input_df)))
  png(file.path(result_dir, file_name), plot_width*2*(k-1)/plts[["nrow"]], plot_height*plts[["nrow"]])
  do.call(grid.arrange, plts)
  dev.off()
}

####
design.APA.primers = function(aUTR_df,
                              species = "Mm", geno = "mm9",
                              txdb_path = "/home/dinghai/projects/fud/mm9.201509.refGene.txdb.sqlite",
                              primer3_path = "/home/dinghai/tools/primer3/src/primer3_core",
                              PRIMER_MISPRIMING_LIBRARY = "/home/dinghai/projects/fud/mm9.201509.refGene.fasta",
                              result_dir = "/home/dinghai/projects/ER/result/190411_mm9/PM",
                              primer_file_name = "primers_for_selected_aUTRs_designed_by_primer3.csv"
){
  #### aUTR_df must have the following columns: gene_symbol, proximal_pA, distal_pA, aUTR_length ################
  
  # Step 1. Create mispriming library unless it has been created before
  get.mispriming.library = function(species = "Mm", geno = "mm9",
                                    txdb_path = "/home/dinghai/projects/fud/mm9.201509.refGene.txdb.sqlite",
                                    mispriming_lib_path = "/home/dinghai/projects/fud/mm9.201509.refGene.fasta"){ 
    require(GenomicFeatures)
    txdb = loadDb(txdb_path)
    exons = exonsBy(txdb, by = "tx", use.names=T)
    
    # Add gene symbol to accession number
    db = paste0("org.", species, ".eg.db")
    require(db, character.only = T) #### Cool!
    acc2sym = AnnotationDbi::select(eval(parse(text = db)), keys = names(exons), keytype =  "ACCNUM", columns =  "SYMBOL") #### Cool
    acc2sym$SYMBOL_ACCNUM = paste(acc2sym$SYMBOL, acc2sym$ACCNUM, sep = "_")
    names(exons) = acc2sym$SYMBOL_ACCNUM
    
    # Get exon sequences
    if(species == "Mm"){
      BSgeno = paste0("BSgenome.Mmusculus.UCSC.", geno)
      BSspec = "Mmusculus"
    }else if(species == "Hs"){
      BSgeno = paste0("BSgenome.Hsapiens.UCSC.", geno)
      BSspec = "Hsapiens"
    }
    require(BSgeno, character.only = T) 
    exon.seq = getSeq(eval(parse(text = BSspec)), exons)
    
    # Combine exon sequences 
    tx.seq = lapply(exon.seq, function(x) unlist(x))
    tx.seq = DNAStringSet(tx.seq)
    
    # Write fasta file
    writeXStringSet(tx.seq, mispriming_lib_path, format = 'fasta')
  }
  
  if(!file.exists(PRIMER_MISPRIMING_LIBRARY)){
    get.mispriming.library(species = species, 
                           geno = geno, 
                           refseq_path = refseq_path,
                           mispriming_lib_path = PRIMER_MISPRIMING_LIBRARY)
  }
  
  # Step 2. Get template sequences for different target regions of each gene
  get.template = function(species = "Mm", 
                          geno = "mm9",
                          aUTR_df,
                          txdb_path = "/home/dinghai/projects/fud/mm9.201509.refGene.txdb.sqlite"){
    
    if(species == "Mm"){
      BSgeno = paste0("BSgenome.Mmusculus.UCSC.", geno)
      BSspec = "Mmusculus"
    }else if(species == "Hs"){
      BSgeno = paste0("BSgenome.Hsapiens.UCSC.", geno)
      BSspec = "Hsapiens"
    }
    require(BSgeno, character.only = T) 
    require(dplyr)
    require(tidyr)
    
    df = aUTR_df[, c("gene_symbol", "proximal_pA", "distal_pA", "aUTR_length")]
    df$chromosome = sub("(chr.+)([-+])(\\d+)$", "\\1", df$proximal_pA)
    df$strand = sub("(chr.+)([-+])(\\d+)$", "\\2", df$proximal_pA)
    df$p = as.integer(sub("(chr.+)([-+])(\\d+)$", "\\3", df$proximal_pA))
    df$d = as.integer(sub("(chr.+)([-+])(\\d+)$", "\\3", df$distal_pA))
    df$proximal_pA = NULL
    df$distal_pA = NULL
    df = gather(df, key = "p_or_d", value = "pA_pos", -gene_symbol, -chromosome, -strand, -aUTR_length)
    
    # Create the GRanges for the cleavage sites. 
    require(GenomicRanges)
    pA = GRanges(seqnames = df$chromosome, 
                 strand = df$strand,
                 ranges = IRanges(start = df$pA_pos, 
                                  end = df$pA_pos, 
                                  names = df$gene_symbol))
    # Since the pA site might be a few nt outside of exons defined by refseq (microheterogeniety), expand pA by +/-12 nt
    pA = pA + 12
    # Get intron and exon definitions
    require(GenomicFeatures)
    txdb = loadDb(txdb_path)
    introns = intronsByTranscript(txdb, use.names=F) #transcriptsBy() returns preRNAs containing introns.
    exons = exonsBy(txdb, by = "tx", use.names=F)
    
    ### Get UTR sequences upstream of the pA using the refseq database
    from = -325
    to = -25
    ## In case 1 (pAs only overlap with exons) and 2 (pAs overlap with both exons and introns)
    exon.pA = pA[pA %over% exons] 
    tx.seq = vector("character", length(exon.pA))
    for(i in 1:length(exon.pA)){
      # Transcript overlaping with the pA
      ol.tx = subsetByOverlaps(exons, exon.pA[i])[[1]] # >1 transcripts may overlap with one pA
      # Exons overlaping with the pA
      ol.exon = ol.tx[ol.tx %over% exon.pA[i]]
      # Discard the exons downstream of the pA
      ol.tx = ol.tx[ol.tx$exon_rank <= ol.exon$exon_rank]
      # Only keep the part of overlapping exon upstream of the pA site
      if (attributes(strand(exon.pA[i]))$values == "+"){
        end(ol.tx[ol.tx %over% exon.pA[i]]) = (end(exon.pA[i]) + start(exon.pA[i]))/2
      }else{
        start(ol.tx[ol.tx %over% exon.pA[i]]) = (end(exon.pA[i]) + start(exon.pA[i]))/2
      }
      # Get exonic sequence 5' of the pA 
      exon.seq = getSeq(eval(parse(text = BSspec)), ol.tx)
      # Combine the exons into one transcript
      exon.seq = unlist(exon.seq)
      # Subset the sequence
      tx.seq[i] = substring(as.character(exon.seq), length(exon.seq) + from, length(exon.seq) + to)
    }
    exon.pA$UTRseq = tx.seq # some returned sequences will be very short, because the pA is too close to the TSS
    ## In case 3 (pAs only overlap with introns) and 4 (pAs overlap with neither exons nor introns (3' extention)), 
    ## the genomic sequence upstream of the pAs will be returned
    other.pA = pA[!(pA %over% exons)] 
    # GRanges for the 3'end of 3'UTRs
    gr = GRanges(seqnames = df$chromosome, 
                 # the start and end positions for the searched region
                 ranges = IRanges(start = ifelse(df$strand == "+", df$pA_pos + from, df$pA_pos - to), 
                                  end = ifelse(df$strand == "+", df$pA_pos + to, df$pA_pos - from), 
                                  names = df$gene_symbol),
                 strand = df$strand)
    # Only keep the 3'UTRs in case 3 and 4  
    gr = gr[!(pA %over% exons)]
    # Get genomic sequence 
    gr.seq = getSeq(eval(parse(text = BSspec)), gr)
    other.pA$UTRseq = as.vector(gr.seq)
    
    # Shrink the pA position to its original value
    exon.pA = exon.pA - 12
    other.pA = other.pA - 12
    # Combine the two GRanges
    pA = c(exon.pA, other.pA)
    
    # Convert to dataframe and merge
    UTRseq = data.frame(pA_pos = start(pA), UTRseq = pA$UTRseq)
    df = merge(df, UTRseq, by = "pA_pos", all = T, sort = F)
    
    ### Get genomic sequences surrounding the pAs (The primers may be in introns)
    pA = pA + 150
    pA$FlankSeq = getSeq(eval(parse(text = BSspec)), pA)
    pA = pA - 150
    FlankSeq = data.frame(pA_pos = start(pA), FlankSeq = pA$FlankSeq)
    df = merge(df, FlankSeq, by = "pA_pos", all = T, sort = F)
    
    # Return the dataframe
    gather(df, key = "seq_type", value = "seq", -pA_pos, -gene_symbol, -chromosome, -strand, -p_or_d, -aUTR_length)
  }
  
  df = get.template(species = species, geno = geno, aUTR_df = aUTR_df, txdb_path = txdb_path)
  
  # Step 3. Create the input file for Primer3
  SEQUENCE_ID = with(df, paste0("SEQUENCE_ID=", gene_symbol, "_", p_or_d, "_", seq_type, "\n"))
  SEQUENCE_TEMPLATE = paste0("SEQUENCE_TEMPLATE=", df$seq,"\n")
  PRIMER_MISPRIMING_LIBRARY = paste0("PRIMER_MISPRIMING_LIBRARY=", PRIMER_MISPRIMING_LIBRARY, "\n")
  # Make sure that the reverse primer for the distal pA is located in the aUTR region
  SEQUENCE_PRIMER_PAIR_OK_REGION_LIST = ifelse(df$p_or_d == "d" & df$seq_type == "UTRseq", 
                                               paste0("SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=", ",,", pmax(300-df$aUTR_length, 1), ",", pmin(df$aUTR_length-25, 300),"\n"), 
                                               ifelse(df$seq_type == "FlankSeq", "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,150,150,150\n", "") 
  )
  
  lines = paste0(SEQUENCE_ID, SEQUENCE_TEMPLATE, SEQUENCE_PRIMER_PAIR_OK_REGION_LIST, PRIMER_MISPRIMING_LIBRARY,
                 "PRIMER_TASK=pick_detection_primers\n", 
                 "PRIMER_MIN_SIZE=15\n",
                 "PRIMER_MAX_SIZE=30\n", 
                 "PRIMER_PRODUCT_SIZE_RANGE=75-300\n", 
                 "PRIMER_EXPLAIN_FLAG=1\n",
                 "PRIMER_PAIR_NUM_RETURNED=1\n", 
                 "PRIMER_MAX_LIBRARY_MISPRIMING=999\n",
                 "PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=999\n",
                 "=\n")
  
  # Write "lines" into a file as input for primer3
  write.table(lines, file = file.path(result_dir, "primer3_input.txt"), col.names=F, row.names=F, quote=F, eol = "")
  
  # Step 4. Call Primer3
  cmd = paste(primer3_path, "<", file.path(result_dir, "primer3_input.txt"),
              ">", file.path(result_dir, "primer3_output.txt"), sep = " ")
  print("Calling primer3 using the following command:")
  print(cmd)
  system(cmd)
  
  # Step 5. Parse Primer3 output
  primers = data.frame()
  for(line in readLines(file.path(result_dir, "primer3_output.txt"))){
    if(grepl("SEQUENCE_ID=", line)){
      target_name = sub("SEQUENCE_ID=", "", line)
      gene_symbol = sub("(.+?)_.+", "\\1", target_name)
      forward_primers = vector()
      reverse_primers = vector()
      misprimings = vector()
      amplicon_sizes = vector()
    }else if(line == "="){
      target_summary = data.frame(target_name = target_name, forward_primer = forward_primers, 
                                  reverse_primer = reverse_primers, mispriming = misprimings, 
                                  amplicon_size = amplicon_sizes)
      target_summary$amplicon_size = as.integer(as.character(target_summary$amplicon_size))
      target_summary = target_summary[order(target_summary$amplicon_size), ]
      if(sum(!misprimings)){
        target_summary = subset(target_summary, !mispriming)
      }
      primers = rbind(primers, target_summary)
    }else{
      if(grepl("PRIMER_LEFT_[1234]_SEQUENCE=", line)){
        forward_primers = c(forward_primers, sub("PRIMER_LEFT_[1234]_SEQUENCE=", "",  line))
      }else if(grepl("PRIMER_RIGHT_[1234]_SEQUENCE=", line)){
        reverse_primers = c(reverse_primers, sub("PRIMER_RIGHT_[1234]_SEQUENCE=", "",  line))
      }else if(grepl("PRIMER_PAIR_[1234]_LIBRARY_MISPRIMING=", line)){
        off_target = sub("PRIMER_PAIR_[1234]_LIBRARY_MISPRIMING=.+, ", "",  line)
        off_target = sub("(.+?)_.+", "\\1", off_target)
        if(grepl("reverse", off_target) | toupper(off_target) == toupper(gene_symbol)){
          misprimings = c(misprimings, F)
        }else{
          misprimings = c(misprimings, T)
        }
      }else if(grepl("PRIMER_PAIR_[1234]_PRODUCT_SIZE=", line)){
        amplicon_sizes = c(amplicon_sizes, sub("PRIMER_PAIR_[1234]_PRODUCT_SIZE=", "", line))
      }
    } 
  }
  
  df$target_name = paste(df$gene_symbol, df$p_or_d, df$seq_type, sep = "_")
  df = merge(df, primers, all.x = T)
  df$seq = NULL
  
  write.csv(df, file.path(result_dir, primer_file_name), row.names = F)
}

add_num_upstream_3UTR_pA = function(pA.df){
  require(dplyr)
  require(tidyr)
  df = pA.df %>%
    dplyr::filter(region == "3UTR") %>%
    arrange(UTR_length) %>%
    group_by(gene_symbol) %>%
    mutate(num_upstream_pA = min_rank(UTR_length) - 1) 
  
  merge(pA.df, df[, c("pAid", "num_upstream_pA")], all.x=T)
}

add_num_conserved_upstream_3UTR_pA = function(pA.df){
  require(dplyr)
  require(tidyr)
  df = pA.df %>%
    dplyr::filter(region == "3UTR") %>%
    arrange(UTR_length) %>%
    group_by(gene_symbol) %>%
    mutate(num_conserved_upstream_pA = cumsum(!is.na(conserveLV)) - (!is.na(conserveLV))) 
  
  merge(pA.df, df[, c("pAid", "num_conserved_upstream_pA")], all.x=T)
}


get_markov_samples = function(input_sequence, num_output_sequence){
  # A function that use Markov chain to sample new sequences
  require(Biostrings)
  seq2tm = function(input_sequence = "AATTGGCCC"){
    # From sequence to transition matrix 
    transition_counts = oligonucleotideFrequency(DNAString(input_sequence), 2)
    
    count_matrix = matrix(0, nrow = 4, ncol = 4,
                          dimnames=list(c("A", "C", "G", "T"),
                                        c("A", "C", "G", "T")))
    
    for(i in c("A", "C", "G", "T")){
      for(j in c("A", "C", "G", "T")){
        count_matrix[i, j] = transition_counts[paste0(i,j)]
      }
    }
    
    # If all 0, set equal probability to each state
    count_matrix[rowSums(count_matrix) == 0, ] = c(1, 1, 1, 1)
    
    count_matrix/rowSums(count_matrix)
 
  }
  
  tm2seq = function(m, first, seqlen){
    # Sample from transition matrix 
    seq = character(seqlen)
    
    seq[1] = ifelse(first %in% rownames(m), first, 
                    sample(rownames(m)[rowSums(is.na(m)) == 0], 1))
    
    for(k in 1:(seqlen - 1)){
      seq[k+1] = sample(rownames(m), replace=T, 1, m[seq[k],])}
    
    paste0(seq, collapse="")
  }
  
  
  m = seq2tm(input_sequence)
  
  output_sequences = lapply(1:num_output_sequence, 
                            function(i) tm2seq(m, substr(input_sequence[1], 1, 1), nchar(input_sequence)))
  
  do.call(cbind, output_sequences)
}

########### A function to map unstranded 3'end-seq data to polyA_Db ###################
map_3_prime_seq_to_polyA_Db = function(reads_file, polyA_Db_file, min_mapq = 10, pA_window = 24, max_read_length = 250){
  library(GenomicAlignments)
  library(dplyr)
  # Read polyA_DB 
  pADb = read.table(polyA_Db_file, sep = "\t", header = F, stringsAsFactors = F,
                    col.names = c("chr", "start", "end", "pAid", "gene_symbol", "gene_id", "pA_type", "strand", "LOCATION", "genicLOC", "region")) %>% 
    dplyr::select(pAid, chr, strand, start, end) %>% 
    mutate(pos = as.integer(sub("chr.+:(\\d+):[+-]", "\\1", pAid))) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = T)
  
  # Define 3'UTR GRanges by extending the pA position upstream max_read_length nt and downstream pA_window nt
  pADb_p = pADb[strand(pADb) == "+"]
  start(pADb_p) = pADb_p$pos - max_read_length
  end(pADb_p) = pADb_p$pos + pA_window
  
  pADb_m = pADb[strand(pADb) == "-"]
  start(pADb_m) = pADb_m$pos - pA_window
  end(pADb_m) = pADb_m$pos + max_read_length
  
  # Read the reads file
  reads = read.table(reads_file, sep = "\t", header = F, stringsAsFactors = F, 
                     col.names = c("chr", "start", "end", "name", "score", "strand")) 
  raw_read_num = nrow(reads)
  reads = dplyr::filter(reads, strand %in% c("+", "-") & score >= min_mapq) %>% mutate(strand = "*") %>% select(-name)
  mapq_read_num = nrow(reads)
  reads = makeGRangesFromDataFrame(reads, keep.extra.columns = T)
  
  # Find overlaps, requiring that the reads are within the 3'UTR 
  o_p = findOverlaps(reads, pADb_p, type="within")
  o_m = findOverlaps(reads, pADb_m, type="within")
  
  # For reads mapped to the + strand, require end(read) - end(3'UTR) is between [-48, 0];
  dis = end(reads[queryHits(o_p)]) - end(pADb_p[subjectHits(o_p)])
  o_p = o_p[dis <= 0 & dis >= -pA_window*2]
  
  # For reads mapped to the - strand, require start(read) - start(3'UTR) is between [0, 48]
  dis = start(reads[queryHits(o_m)]) - start(pADb_m[subjectHits(o_m)])
  o_m = o_m[dis >= 0 & dis <= pA_window*2]
  
  total_PASS_num = length(o_p) + length(o_m) 
  
  # Remove ambiguious reads that are mapped to >1 pAs on the same strand or different strands
  ambiguious_PASS = unique(c(queryHits(o_p), queryHits(o_m))[duplicated(c(queryHits(o_p), queryHits(o_m)))])
  o_p = o_p[!queryHits(o_p) %in% ambiguious_PASS]
  o_m = o_m[!queryHits(o_m) %in% ambiguious_PASS]
  
  ambiguious_PASS_num = total_PASS_num - length(o_p) - length(o_m) 
  
  # Combine reads and pADb
  pA_p = cbind(as.data.frame(pADb_p[subjectHits(o_p)]), mcols(reads[queryHits(o_p)]))
  pA_m = cbind(as.data.frame(pADb_m[subjectHits(o_m)]), mcols(reads[queryHits(o_m)]))
  
  # Combine + and - strand
  pA = as.data.frame(rbind(pA_p, pA_m)) %>% group_by(pAid) %>% summarise(count = n())
  names(pA) = sub("^count$", sub(".+\\/(.+)\\.bed", "\\1", reads_file), names(pA)) 
  # Return pA, and input/output read num
  return(list(PASS_count = as.data.frame(pA), 
              stats = data.frame(raw_read = raw_read_num, mapq_read = mapq_read_num, 
                                 total_PASS = total_PASS_num, ambiguious_PASS = ambiguious_PASS_num, 
                                 PASS = total_PASS_num - ambiguious_PASS_num)))
}
# Example:
# reads_file = "../data/PC1.bed"
# polyA_Db_file = "../data/human.cental.polyAdb.V3.3.bed"
# outputs = map_3_prime_seq_to_polyA_Db(reads_file, polyA_Db_file, pA_window = 24, max_read_length = 250)

