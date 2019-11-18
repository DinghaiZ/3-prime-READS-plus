#### Plot rpm ratio between fractions for combined upstream pA vs downstream pA 
p.list = list()
d.list = list()
k = 1
if(length(treatments) > 1){
  for(i in 1:(length(treatments)-1)){
    for (j in (i+1):length(treatments)){
      # only keep desirable comparisons
      if(!exists("comparisons") | (exists("comparisons") && paste(treatments[j], treatments[i], sep="_") %in% comparisons)){
        for(fraction in fractions){
          for(batch in batches){
            # focus on these two columns
            col1_count = paste(treatments[i], fraction, batch, "count", sep="_")
            col2_count = paste(treatments[j], fraction, batch, "count", sep="_")
            col1_rpm = paste(treatments[i], fraction, batch, "rpm", sep="_")
            col2_rpm = paste(treatments[j], fraction, batch, "rpm", sep="_")
            
            if(all(c(col1_count, col2_count, col1_rpm, col2_rpm) %in% names(pA.df))){
              tmp = pA.df[, c(names(pA.df)[1:6], c(col1_count, col2_count, col1_rpm, col2_rpm, "pA_type"))]
              
              # ds: downstream; us: upstream
              tmp$region_group = ifelse(tmp$pA_type %in% c("S", "F", "M", "L"), "ds", 
                                        ifelse(tmp$pA_type %in% c("UA"), "ua", "us")
                                        )
              # only keep genes with both upstream and downstream pAs
              tmp = subset(tmp, gene_symbol %in% tmp[tmp$region_group == "us",]$gene_symbol)
              tmp = subset(tmp, gene_symbol %in% tmp[tmp$region_group == "ds",]$gene_symbol)
              tmp = tmp[tmp$region_group != "ua",]
              
              # combine read counts in upstream and downstream regions
              tmp.list  = split(tmp, f = tmp$gene_symbol)
              
              for(k2 in 1:length(tmp.list)){
                subtotal = aggregate(tmp.list[[k2]][,c(col1_count, col2_count, col1_rpm, col2_rpm)], list(region_group = tmp.list[[k2]]$region_group), sum) 
                subtotal = cbind(tmp.list[[k2]][1:2, c("gene_symbol", "chr", "strand")], subtotal)
                # fisher's exact test
                p_val = fisher.test(subtotal[, c(col1_count, col2_count)])$p.value
                # calculate values for the scatter plot
                x_val = log2(subtotal[subtotal$region_group == "us", col2_rpm]/subtotal[subtotal$region_group == "us", col1_rpm])
                y_val = log2(subtotal[subtotal$region_group == "ds", col2_rpm]/subtotal[subtotal$region_group == "ds", col1_rpm])
                
                if(k2 == 1){
                  tmp2 = data.frame(gene_symbol = tmp.list[[k2]]$gene_symbol[1],
                                    x_val = x_val,
                                    y_val = y_val,
                                    p_val = p_val, stringsAsFactors = F)
                }else{
                  tmp2 = rbind(tmp2, c(tmp.list[[k2]]$gene_symbol[1], x_val, y_val, p_val))
                }
              }
              for(col_num in 2:4){
                tmp2[, col_num] = as.numeric(tmp2[, col_num])
              }
              tmp2 = tmp2[abs(rowSums(tmp2[, c("x_val", "y_val")])) < Inf, ]
              tmp2 = na.omit(tmp2)
              rm(tmp.list)
              
              tmp2$colors = ifelse(tmp2$p_val >= 0.05, "grey", 
                                   ifelse(tmp2$y_val > tmp2$x_val, "red", "blue"))
              
              blue_over_red = ifelse(sum(tmp2$colors == "red") == 0, 
                                     log2(sum(tmp2$colors == "blue")+1) - log2(sum(tmp2$colors == "red")+1),
                                     log2(sum(tmp2$colors == "blue")) - log2(sum(tmp2$colors == "red")))
              
              log2_ratio = paste("log2(Blue/Red) = ", round(blue_over_red, 2), sep = "")
              
              d.list[[paste(fraction, treatments[j], "vs", treatments[i], batch, sep = "_")]] = tmp2
              
              # scatter plot ###plot colored points later
              layer1 = subset(tmp2, colors == "grey")
              layer2 = subset(tmp2, colors != "grey")
              print(paste0("Plotting ", paste(fraction, paste0(treatments[j], "/", treatments[i], ", Upstream vs Downstream"), batch, sep = "_")))
  
              p.list[[k]] = ggplot() +
                geom_point(data=layer1, aes_string(x="x_val", y="y_val"), colour = "grey", shape = 1) + 
                geom_point(data=layer2, aes_string(x="x_val", y="y_val"), colour = layer2$colors, shape = 1) +
                ggtitle(paste0("Grey: ", sum(tmp2$colors == "grey"), " | ",
                               "Red: ", sum(tmp2$colors == "red"), " | ",
                               "Blue: ", sum(tmp2$colors == "blue"), "\n",
                               log2_ratio))  +
                xlim(range(tmp2$x_val)) + ylim(range(tmp2$y_val)) +
                xlab(paste0("Upstream Region\n", paste(treatments[j], "vs", treatments[i], fraction, batch, sep = "_"))) +
                ylab(paste0(paste(treatments[j], "vs", treatments[i], fraction, batch, sep = "_"), "\n3'UTR"))
              
              k = k + 1
            }
          }
        }
      }
    }
  }
  if(k > 1){
    p.list = c(p.list, ncol = length(fractions)) 
    png(file.path(result_dir, "Upstream_downstream_treatment_ratios_p=0.05.png"), height=425*(k-1)/length(fractions), width=400*length(fractions))
    do.call(grid.arrange, p.list)
    dev.off()
  }
}


#### Plot rpm ratio between fractions for combined upstream pA vs downstream pA 
p.list = list()
d.list = list()
k = 1
if(length(fractions) > 1){
  for(i in 1:(length(fractions)-1)){
    for (j in (i+1):length(fractions)){
      # only keep desirable comparisons
      if((exists("comparisons") && paste(fractions[j], fractions[i], sep="_") %in% comparisons) | 
         !exists("comparisons")){
        for (treatment in treatments){
          for (batch in batches){
            # focus on these two columns
            col1_count = paste(treatment, fractions[i], batch, "count", sep="_")
            col2_count = paste(treatment, fractions[j], batch, "count", sep="_")
            col1_rpm = paste(treatment, fractions[i], batch, "rpm", sep="_")
            col2_rpm = paste(treatment, fractions[j], batch, "rpm", sep="_")
            
            if(all(c(col1_count, col2_count, col1_rpm, col2_rpm) %in% names(pA.df))){
              tmp = pA.df[, c(names(pA.df)[1:6], c(col1_count, col2_count, col1_rpm, col2_rpm, "pA_type"))]
              
              # ds: downstream; us: upstream
              tmp$region_group = ifelse(tmp$pA_type %in% c("S", "F", "M", "L"), "ds", 
                                        ifelse(tmp$pA_type %in% c("UA"), "ua", "us")
                                        ) 
              # only keep genes with both upstream and downstream pAs
              tmp = subset(tmp, gene_symbol %in% tmp[tmp$region_group == "us",]$gene_symbol)
              tmp = subset(tmp, gene_symbol %in% tmp[tmp$region_group == "ds",]$gene_symbol)
              tmp = tmp[tmp$region_group != "ua",]
              
              # combine read counts in upstream and downstream regions
              tmp.list  = split(tmp, f = tmp$gene_symbol)
              
              for(k2 in 1:length(tmp.list)){
                subtotal = aggregate(tmp.list[[k2]][,c(col1_count, col2_count, col1_rpm, col2_rpm)], list(region_group = tmp.list[[k2]]$region_group), sum) 
                subtotal = cbind(tmp.list[[k2]][1:2, c("gene_symbol", "chr", "strand")], subtotal)
                # fisher's exact test
                p_val = fisher.test(subtotal[, c(col1_count, col2_count)])$p.value
                # calculate values for the scatter plot
                x_val = log2(subtotal[subtotal$region_group == "us", col2_rpm]/subtotal[subtotal$region_group == "us", col1_rpm])
                y_val = log2(subtotal[subtotal$region_group == "ds", col2_rpm]/subtotal[subtotal$region_group == "ds", col1_rpm])
                
                if(k2 == 1){
                  tmp2 = data.frame(gene_symbol = tmp.list[[k2]]$gene_symbol[1],
                                    x_val = x_val,
                                    y_val = y_val,
                                    p_val = p_val, stringsAsFactors = F)
                }else{
                  tmp2 = rbind(tmp2, c(tmp.list[[k2]]$gene_symbol[1], x_val, y_val, p_val))
                }
              }
              for(col_num in 2:4){
                tmp2[, col_num] = as.numeric(tmp2[, col_num])
              }
              tmp2 = tmp2[abs(rowSums(tmp2[, c("x_val", "y_val")])) < Inf, ]
              tmp2 = na.omit(tmp2)
              rm(tmp.list)
              
              tmp2$colors = ifelse(tmp2$p_val >= 0.05, "grey", 
                                   ifelse(tmp2$y_val > tmp2$x_val, "red", "blue"))
              
              blue_over_red = ifelse(sum(tmp2$colors == "red") == 0, 
                                     log2(sum(tmp2$colors == "blue")+1) - log2(sum(tmp2$colors == "red")+1),
                                     log2(sum(tmp2$colors == "blue")) - log2(sum(tmp2$colors == "red")))
              
              log2_ratio = paste("log2(Blue/Red) = ", round(blue_over_red, 2), sep = "")
              
              d.list[[paste(treatment, fractions[j], "vs", fractions[i], batch, sep = "_")]] = tmp2
              
              # scatter plot ###plot colored points later
              layer1 = subset(tmp2, colors == "grey")
              layer2 = subset(tmp2, colors != "grey")
              print(paste0("Plotting ", col2_rpm, "/", col1_rpm, ", Upstream vs Downstream"))
              p.list[[k]] = ggplot() +
                geom_point(data=layer1, aes_string(x="x_val", y="y_val"), colour = "grey", shape = 1) + 
                geom_point(data=layer2, aes_string(x="x_val", y="y_val"), colour = layer2$colors, shape = 1) +
                ggtitle(paste0("Grey: ", sum(tmp2$colors == "grey"), " | ",
                               "Red: ", sum(tmp2$colors == "red"), " | ",
                               "Blue: ", sum(tmp2$colors == "blue"), "\n",
                               log2_ratio))  +
                xlim(range(tmp2$x_val)) + ylim(range(tmp2$y_val)) +
                xlab(paste0("Upstream Region\n", paste(treatment, fractions[j], "vs", fractions[i], batch, sep = "_"))) +
                ylab(paste0(paste(treatment, fractions[j], "vs", fractions[i], batch, sep = "_"), "\n3'UTR"))

              k = k + 1
            }
          }
        }
      }
    }
  }
  if(k > 1){
    p.list[["nrow"]] = length(treatments)
    png(file.path(result_dir, "Upstream_downstream_fraction_ratios_p=0.05.png"), height=425*length(treatments), width=400*(k-1)/length(treatments))
    do.call(grid.arrange, p.list)
    dev.off()
  }
}


#### Scatter plot, distal pA vs proximal pA
require(grid)
require(ggplot2)
p.list = list() # list for plots
d.list = list() # list for data
p2.list = list() # centered scatter plot
k = 1
names(pA.df) = sub("Isoform_abn_(.+)", "\\1_usage", names(pA.df))
### calculate and plot ratio between fraction1/fraction2 counts for each pA
if(length(fractions) > 1){
  for(treatment in treatments){
    #for(i in c(1)){
    for(i in 1:(length(fractions)-1)){
      for (j in (i+1):length(fractions)){
        # only keep desirable comparisons
        if((exists("comparisons") && paste(fractions[j], fractions[i], sep="_") %in% comparisons) | 
           !exists("comparisons")){
          for(batch in batches){
            # focus on these two columns
            col1_count = paste(treatment, fractions[i], batch, "count", sep="_")
            col2_count = paste(treatment, fractions[j], batch, "count", sep="_")
            col1_usage = paste(treatment, fractions[i], batch, "usage", sep="_")
            col2_usage = paste(treatment, fractions[j], batch, "usage", sep="_")
            col1_rpm = paste(treatment, fractions[i], batch, "rpm", sep="_")
            col2_rpm = paste(treatment, fractions[j], batch, "rpm", sep="_")
            
            if(all(c(col1_count, col2_count) %in% names(pA.df))){
              tmp = pA.df[, c(c("gene_symbol", "chr", "strand", "tss_position", "pA_pos", "pAid", "region", "CDS_end", "UTR_length"), c(col1_count, col2_count, col1_usage, 
                                                   col2_usage, col1_rpm, col2_rpm, grep("(HuR|PTB)_(sites$|counts$)", names(pA.df), value=T), "pA_type", "gene_id", "description", "Num_pA"))]
              #colnames(tmp) = sub("position", "pA_pos", colnames(tmp))
              tmp = tmp[pA.df$pA_type %in% c("F", "M", "L"),]
              # remove isoforms that have low read counts in the two samples
              tmp = subset(tmp, (tmp[,col1_count] + tmp[,col2_count] >= 20) & tmp[,col1_count] > 0 & tmp[,col2_count] > 0)
              # pseudocounts
              #tmp[,c(col1_count, col2_count)][tmp[,c(col1_count, col2_count)] ==0] = 1
              #tmp_unormalized = tmp
              # normalize (RPM)
              #tmp[, col1_count] = round(tmp[, col1_rpm])
              #tmp[, col2_count] = round(tmp[, col2_rpm])
              # the rpm columns are needed for choosing the top two isoforms by get.pA.pairs()
              #             tmp[, sub("count", "rpm", col1_count)] = tmp[, col1_count] 
              #             tmp[, sub("count", "rpm", col2_count)] = tmp[, col2_count] 
              #             tmp_unormalized[, col1_rpm] = tmp[, col1_count] 
              #             tmp_unormalized[, col2_rpm] = tmp[, col2_count] 
              
              # reorganize data to get pairwise ratios
              #pair.pA = pairwise.pAs(data = tmp, neighbor = F)
              pair.pA = get.pA.pairs(data = tmp, cols = c(col1_usage, col2_usage, col1_rpm, col2_rpm, 
                                                          grep("(HuR|PTB)_(sites$|counts$)", names(pA.df), value=T), "UTR_length", "description", "Num_pA"), neighbor = neighbor, toptwo = toptwo, match_only=T)
              #pair.pA_unormalized = pairwise.pAs(data = tmp_unormalized, neighbor = F)
              #pair.pA_unormalized = get.pA.pairs(data = tmp_unormalized, 
              #                                   neighbor = neighbor, toptwo = toptwo, match_only=T)
              # discard the calculated Distal/Proximal ratio
              #pair.pA = pair.pA[, -grep("_ratio", names(pair.pA))]
              #pair.pA_unormalized = pair.pA_unormalized[, -grep("_ratio", names(pair.pA_unormalized))]
              # calculate the ratio accross fractions
              x=paste(treatment, fractions[j], "vs", fractions[i], batch, "Prox", sep = "_")
              y=paste(treatment, fractions[j], "vs", fractions[i], batch, "Dis", sep = "_")
              pair.pA[, x] = log2(pair.pA[,paste(treatment, fractions[j], batch, "rpm_proximal", sep = "_")]) - log2(pair.pA[,paste(treatment, fractions[i], batch, "rpm_proximal", sep = "_")])
              pair.pA[, y] = log2(pair.pA[,paste(treatment, fractions[j], batch, "rpm_distal", sep = "_")]) - log2(pair.pA[,paste(treatment, fractions[i], batch, "rpm_distal", sep = "_")])
              
              # calculate p values using fisher's exact test
              #pair.pA$fisher = p.adjust(apply(pair.pA[, grep("_count_", names(pair.pA))], 1, row.fisher))
              pair.pA$fisher = apply(pair.pA[, grep("_count_", names(pair.pA))], 1, row.fisher)
              # calculate colors
              pair.pA$colors = "grey"
              pair.pA$colors[(pair.pA$fisher < p_threshold) & (pair.pA[, y] > pair.pA[, x]) & 
                               (pair.pA[, paste0(col2_usage, "_distal")] - pair.pA[, paste0(col1_usage, "_distal")] >= 5)] = "red"
              pair.pA$colors[(pair.pA$fisher < p_threshold) & (pair.pA[, y] < pair.pA[, x]) &
                               ((pair.pA[, paste0(col2_usage, "_proximal")] - pair.pA[, paste0(col1_usage, "_proximal")] >= 5))] = "blue"
              #             pair.pA$colors = ifelse(pair.pA$fisher >= p_threshold, "grey", 
              #                                     ifelse(pair.pA[, y] > pair.pA[, x], "red", "blue"))
              
              # direction of regulation, not considering p values
              pair.pA$Num_pA = pair.pA$Num_pA_distal
              pair.pA$Num_pA_distal = NULL
              pair.pA$Num_pA_proximal = NULL
              pair.pA$direction = ifelse(pair.pA[, y] > pair.pA[, x], "up", "down")
              pair.pA$description = pair.pA$description_distal
              pair.pA$description_distal = NULL
              pair.pA$description_proximal = NULL
              
              d.list[[paste(treatment, fractions[j], "vs", fractions[i], batch, sep = "_")]] = pair.pA
              # calculate annotation text
              # red_over_blue = ifelse(sum(pair.pA$colors == "blue") == 0, 
              #                        log2(sum(pair.pA$colors == "red")+1) - log2(sum(pair.pA$colors == "blue")+1),
              #                        log2(sum(pair.pA$colors == "red")) - log2(sum(pair.pA$colors == "blue")))
              blue_over_red = ifelse(sum(pair.pA$colors == "red") == 0, 
                                     log2(sum(pair.pA$colors == "blue")+1) - log2(sum(pair.pA$colors == "red")+1),
                                     log2(sum(pair.pA$colors == "blue")) - log2(sum(pair.pA$colors == "red")))

              log2_ratio = paste("log2(Blue/Red) = ", round(blue_over_red, 2), sep = "")
              
              # scatter plot ###plot colored points later
              layer1 = subset(pair.pA, colors == "grey")
              layer2 = subset(pair.pA, colors != "grey")
              print(paste0("Plotting ", y, " ~ ", x))
              p.list[[k]] = ggplot() +
                geom_point(data=layer1, aes_string(x=x, y=y), colour = "grey", shape = 1) + 
                geom_point(data=layer2, aes_string(x=x, y=y), colour = layer2$colors, shape = 1) +
                ggtitle(paste0("Grey:", sum(pair.pA$colors == "grey"), "|",
                               "Red:", sum(pair.pA$colors == "red"), "|",
                               "Blue:", sum(pair.pA$colors == "blue"), "\n",
                               log2_ratio)) #+
              #xlim(-3,3) + ylim(-3,3)
              #theme(panel.grid = element_blank())
              
              # centered scatter plot
              p2.list[[k]] = p.list[[k]] + xlim(-3,3) + ylim(-3,3) #+
              #theme(panel.grid = element_blank())
              
              k = k + 1
            }
          }
        }
        
      }
    }
  }
  if(k > 1){ 
    #   data_range = c(min(sapply(d.list, function(x) min(x[,grep("vs", names(x))]))), 
    #                  max(sapply(d.list, function(x) max(x[,grep("vs", names(x))]))))
    #   p.list = lapply(p.list, function(p) p+xlim(data_range)+ylim(data_range))
    p.list[["nrow"]] = length(treatments) 
    png(file.path(result_dir, paste0("Dis_Prox_fraction_ratios_p=", p_threshold, ".png")), width=450*(k-1)/length(treatments), height=500*length(treatments))
    do.call(grid.arrange, p.list)
    dev.off()
    
    # centered scatter plot
    p2.list = c(p2.list, nrow = length(treatments)) 
    png(file.path(result_dir, paste0("Dis_Prox_fraction_ratios_p=", p_threshold, "_centered.png")), width=450*(k-1)/length(treatments), height=500*length(treatments))
    do.call(grid.arrange, p2.list)
    dev.off()
    
    # save data
    for(list_name in names(d.list)){
      write.csv(d.list[[list_name]], file.path(result_dir, paste0("blue_red_", list_name, ".csv")), row.names = F)
    }
  }
}


###### scatter plot, distal pA vs proximal pA, based on mean of replicates
##### using d.list calculated earlier
if(length(batches) > 1){
  require(grid)
  require(ggplot2)
  p.list = list() # list for plots
  p2.list = list() # centered scatter plot
  v.list = list()
  k = 1
  # merge data in the d.list calculated earlier
  #comparisons2 = unique(sub("\\d+[a-z]?$", "", names(d.list)))
  comparisons2 = sub("\\d+[a-z]?$", "", names(d.list))
  comparisons2 = unique(comparisons2[duplicated(comparisons2)])
  if(length(comparisons2) >= 1){
    for(comparison in comparisons2){
      #df = d.list[[paste0(comparison, batches[1])]]
      list.names = grep(comparison, names(d.list), value=T)
      df = d.list[[list.names[1]]]
      
      for(list.name in list.names[-1]){
        tmp = d.list[[list.name]]
        df = match.pA.pair(df, tmp[,-grep("gene_id", names(tmp))], rm_dup_col = F)
      }
      
      # combine colors in different samples
      # if relaxed criteria are used to call blue or red aUTRs:
      if(exists("relaxed_blue_red_criteria")){ 
        if(relaxed_blue_red_criteria){# p < 0.05 in >= half of the replicates and with the same regulation directions in all replicates
          df$colors = ifelse(rowSums(df[, grep("direction", names(df))] == "up",na.rm=T) == length(batches) & rowSums(df[, grep("colors", names(df))] == "red",na.rm=T) >= length(batches)/2, "red",
                             ifelse(rowSums(df[, grep("direction", names(df))] == "down",na.rm=T) == length(batches) & rowSums(df[, grep("colors", names(df))] == "blue", na.rm=T) >= length(batches)/2, "blue", 
                                    "grey"))
        }else{
          # only keep consistant colors
          df$colors = ifelse(rowSums(df[, grep("colors", names(df))] == "red",na.rm=T) == length(list.names), "red",
                             ifelse(rowSums(df[, grep("colors", names(df))] == "blue",na.rm=T) == length(list.names), "blue", "grey"))
        }
      }else{
        # only keep consistant colors
        df$colors = ifelse(rowSums(df[, grep("colors", names(df))] == "red",na.rm=T) == length(list.names), "red",
                           ifelse(rowSums(df[, grep("colors", names(df))] == "blue",na.rm=T) == length(list.names), "blue", "grey"))
      }
      
      # average the ratios
      x = sub("_$", "_Prox", comparison)
      y = sub("_$", "_Dis", comparison)
      df[, x] = rowMeans(df[, grep(paste0(comparison, ".+_Prox"), names(df))], na.rm=T)
      df[, y] = rowMeans(df[, grep(paste0(comparison, ".+_Dis"), names(df))], na.rm=T)
      
      # save data for venn diagram
      tmp = df[df$colors=="red",] 
      write.csv(tmp, file.path(result_dir, paste0(comparison, "red.csv")), row.names = F)
      v.list[[paste0(comparison, "Red")]] = paste0(tmp$proximal_pA, tmp$distal_pA)
      tmp = df[df$colors=="blue",] 
      write.csv(tmp, file.path(result_dir, paste0(comparison, "blue.csv")), row.names = F)
      v.list[[paste0(comparison, "Blue")]] = paste0(tmp$proximal_pA, tmp$distal_pA)
      rm(tmp)
      
      # scatter plot ###plot colored points later
      layer1 = subset(df, colors == "grey")
      layer2 = subset(df, colors != "grey")
      
      # calculate annotation text
      # red_over_blue = ifelse(sum(df$colors == "blue") == 0, 
      #                        log2(sum(df$colors == "red")+1) - log2(sum(df$colors == "blue")+1),
      #                        log2(sum(df$colors == "red")) - log2(sum(df$colors == "blue")))
      blue_over_red = ifelse(sum(df$colors == "red") == 0, 
                             log2(sum(df$colors == "blue")+1) - log2(sum(df$colors == "red")+1),
                             log2(sum(df$colors == "blue")) - log2(sum(df$colors == "red")))
      
      # log2_ratio = ifelse(red_over_blue >= 0,
      #                     paste("Red/Blue=", round(red_over_blue,2), sep = ""), 
      #                     paste("Blue/Red=", round(blue_over_red,2), sep = ""))
      log2_ratio = paste("log2(Blue/Red) = ", round(blue_over_red, 2), sep = "")
      
      print(paste0("Plotting ", y, " ~ ", x))
      p.list[[k]] = ggplot() +
        geom_point(data=layer1, aes_string(x=x, y=y), colour = "grey", shape = 1) + 
        geom_point(data=layer2, aes_string(x=x, y=y), colour = layer2$colors, shape = 1) +
        ggtitle(paste0("Grey:", sum(df$colors == "grey"), "|",
                       "Red:", sum(df$colors == "red"), "|",
                       "Blue:", sum(df$colors == "blue"), "\n",
                       log2_ratio)) #+ xlim(-2,2) + ylim(-2,2) #+
      
      # centered scatter plot
      p2.list[[k]] = p.list[[k]] + xlim(-3,3) + ylim(-3,3) #+
      #theme(panel.grid = element_blank())
      
      k = k + 1
    }
    
    # scatter plot
    p.list = c(p.list, nrow = length(comparisons2)) 
    png(file.path(result_dir, paste0("Dis_Prox_mean_fraction_ratios_p=", p_threshold, ".png")), height=550*length(comparisons2), width=500*(k-1)/length(comparisons2))
    do.call(grid.arrange, p.list)
    dev.off()
    
    # centered scatter plot
    p2.list = c(p2.list, nrow = length(comparisons2)) 
    png(file.path(result_dir, paste0("Dis_Prox_mean_fraction_ratios_p=", p_threshold, "_centered.png")), height=550*length(comparisons2), width=500*(k-1)/length(comparisons2))
    do.call(grid.arrange, p2.list)
    dev.off()
    # venn diagram
    if(length(v.list) < 6 & length(v.list) > 2){
      require(VennDiagram)
      venn.diagram(
        x = v.list,
        filename = file.path(result_dir, paste0("Dis_Prox_fraction_ratios_venn_p=", p_threshold, ".png")),
        height=2000, width=2000,
        imagetype="png",
        cat.cex = 0.5,
        margin = 0.1,
        col = "transparent",
        fill = c("cornflowerblue", "green", "yellow", "darkorchid1")
      )
    }
  }
}


#### Calculate and plot ratio between treatment1/treatment2 counts for each pA
p.list = list()
d.list = list()
p2.list = list() 
k = 1
if(length(treatments) > 1){
  for(i in 1:(length(treatments)-1)){
    for (j in (i+1):length(treatments)){
      # only keep desirable comparisons
      if((exists("comparisons") && paste(treatments[j], treatments[i], sep="_") %in% comparisons) | 
         !exists("comparisons")){
        for(fraction in fractions){
          for (batch in batches){
            # focus on these two columns
            col1_count = paste(treatments[i], fraction, batch, "count", sep="_")
            col2_count = paste(treatments[j], fraction, batch, "count", sep="_")
            col1_usage = paste(treatments[i], fraction, batch, "usage", sep="_")
            col2_usage = paste(treatments[j], fraction, batch, "usage", sep="_")
            col1_rpm = paste(treatments[i], fraction, batch, "rpm", sep="_")
            col2_rpm = paste(treatments[j], fraction, batch, "rpm", sep="_")
            
            if(all(c(col1_count, col2_count) %in% names(pA.df))){
              tmp = pA.df[, c(c("gene_symbol", "chr", "strand", "tss_position", "pA_pos", "pAid", "region", "CDS_end", "UTR_length"), c(col1_count, col2_count, col1_usage, col2_usage, col1_rpm, col2_rpm,
                                                   grep("(HuR|PTB)_(sites$|counts$)", names(pA.df), value=T), "pA_type", "gene_id", "description", "Num_pA"))]
              tmp = tmp[pA.df$pA_type %in% c("F", "M", "L"),]
              # remove isoforms that have low read counts in the two samples
              tmp = subset(tmp, (tmp[,col1_count] + tmp[,col2_count] >= 20) & tmp[,col1_count] > 0 & tmp[,col2_count] > 0)
              
              pair.pA = get.pA.pairs(data = tmp, cols = c(col1_usage, col2_usage, col1_rpm, col2_rpm, grep("(HuR|PTB)_(sites$|counts$)", names(pA.df), value=T), "UTR_length", "description", "Num_pA"), 
                                     neighbor = neighbor, toptwo = toptwo, match_only=T)
              
              # calculate the ratio accross treatments
              x=paste(treatments[j], "vs", treatments[i], fraction, batch, "Prox", sep = "_")
              y=paste(treatments[j], "vs", treatments[i], fraction, batch, "Dis", sep = "_")
              pair.pA[, x] = log2(pair.pA[,paste(treatments[j], fraction, batch, "rpm_proximal", sep = "_")]) - log2(pair.pA[,paste(treatments[i], fraction, batch, "rpm_proximal", sep = "_")])
              pair.pA[, y] = log2(pair.pA[,paste(treatments[j], fraction, batch, "rpm_distal", sep = "_")]) - log2(pair.pA[,paste(treatments[i], fraction, batch, "rpm_distal", sep = "_")])
              
              # calculate p values using fisher's exact test
              pair.pA$fisher = apply(pair.pA[, grep("_count_", names(pair.pA))], 1, row.fisher)
              # calculate colors
              pair.pA$colors = "grey"
              pair.pA$colors[(pair.pA$fisher < p_threshold) & (pair.pA[, y] > pair.pA[, x]) & 
                               (pair.pA[, paste0(col2_usage, "_distal")] - pair.pA[, paste0(col1_usage, "_distal")] >= 5)] = "red"
              pair.pA$colors[(pair.pA$fisher < p_threshold) & (pair.pA[, y] < pair.pA[, x]) &
                               ((pair.pA[, paste0(col2_usage, "_proximal")] - pair.pA[, paste0(col1_usage, "_proximal")] >= 5))] = "blue"
              
              pair.pA$Num_pA = pair.pA$Num_pA_distal
              pair.pA$Num_pA_distal = NULL
              pair.pA$Num_pA_proximal = NULL
              # direction of regulation, not considering p values
              pair.pA$direction = ifelse(pair.pA[, y] > pair.pA[, x], "up", "down")
              pair.pA$description = pair.pA$description_distal
              pair.pA$description_distal = NULL
              pair.pA$description_proximal = NULL
              
              
              # save data for later steps
              d.list[[paste(fraction, treatments[j], "vs", treatments[i], batch, sep = "_")]] = pair.pA
              
              # calculate annotation text
              blue_over_red = ifelse(sum(pair.pA$colors == "red") == 0, 
                                     log2(sum(pair.pA$colors == "blue")+1) - log2(sum(pair.pA$colors == "red")+1),
                                     log2(sum(pair.pA$colors == "blue")) - log2(sum(pair.pA$colors == "red")))
              
              log2_ratio = paste("log2(Blue/Red) = ", round(blue_over_red, 2), sep = "")
              
              # scatter plot ###plot colored points later
              layer1 = subset(pair.pA, colors == "grey")
              layer2 = subset(pair.pA, colors != "grey")
              print(paste0("Plotting ", y, " ~ ", x))
              p.list[[k]] = ggplot() +
                geom_point(data=layer1, aes_string(x=x, y=y), colour = "grey", shape = 1) + 
                geom_point(data=layer2, aes_string(x=x, y=y), colour = layer2$colors, shape = 1) +
                ggtitle(paste0("Grey:", sum(pair.pA$colors == "grey"), "|",
                               "Red:", sum(pair.pA$colors == "red"), "|",
                               "Blue:", sum(pair.pA$colors == "blue"), "\n",
                               log2_ratio)) 
              
              # centered scatter plot
              p2.list[[k]] = p.list[[k]] + xlim(-3,3) + ylim(-3,3) 
              
              k = k + 1
            }
          }
        }
      }
    }
  }
  if(k > 1){
    p.list[["ncol"]] = length(fractions)  
    png(file.path(result_dir, paste0("Dis_Prox_treatment_ratios_p=", p_threshold, ".png")), width = 450*length(fractions), height= 500*(k-1)/length(fractions))
    do.call(grid.arrange, p.list)
    dev.off()
    
    # centered scatter plot
    p2.list[["ncol"]] = length(fractions) 
    png(file.path(result_dir, paste0("Dis_Prox_treatment_ratios_p=", p_threshold, "_centered.png")), width = 450*length(fractions), height= 500*(k-1)/length(fractions))
    do.call(grid.arrange, p2.list)
    dev.off()
    
    # save data
    for(list_name in names(d.list)){
      write.csv(d.list[[list_name]], file.path(result_dir, paste0("blue_red_", list_name, ".csv")), row.names = F)
    }
  }
}

#### scatter plot, distal pA vs proximal pA, based on the mean of replicates
# using d.list calculated earlier
if(length(batches) > 1){
  require(grid)
  require(ggplot2)
  p.list = list() # list for plots
  p2.list = list() # list for centered plots
  v.list = list()
  k = 1
  # merge data in the d.list calculated earlier
  comparisons2 = sub("\\d+[a-z]?$", "", names(d.list))
  comparisons2 = unique(comparisons2[duplicated(comparisons2)])
  if(length(comparisons2) >= 1){
    for(comparison in comparisons2){
      list.names = grep(comparison, names(d.list), value=T)
      df = d.list[[list.names[1]]]
      
      for(list.name in list.names[-1]){
        tmp = d.list[[list.name]]
        df = match.pA.pair(df, tmp[,-grep("gene_id", names(tmp))], rm_dup_col = F)
      }
      
      # combine colors in different samples
      # if relaxed criteria are used to call blue or red aUTRs:
      if(exists("relaxed_blue_red_criteria")){ 
        if(relaxed_blue_red_criteria){# p < 0.05 in >= half of the replicates and in the same regulation direction in all replicates
          df$colors = ifelse(rowSums(df[, grep("direction", names(df))] == "up",na.rm=T) == length(batches) & rowSums(df[, grep("colors", names(df))] == "red",na.rm=T) >= length(batches)/2, "red",
                             ifelse(rowSums(df[, grep("direction", names(df))] == "down",na.rm=T) == length(batches) & rowSums(df[, grep("colors", names(df))] == "blue", na.rm=T) >= length(batches)/2, "blue", 
                                    "grey"))
        }else{
          # only keep consistant colors
          df$colors = ifelse(rowSums(df[, grep("colors", names(df))] == "red",na.rm=T) == length(list.names), "red",
                             ifelse(rowSums(df[, grep("colors", names(df))] == "blue",na.rm=T) == length(list.names), "blue", "grey"))
        }
      }else{
        # only keep consistant colors
        df$colors = ifelse(rowSums(df[, grep("colors", names(df))] == "red",na.rm=T) == length(list.names), "red",
                           ifelse(rowSums(df[, grep("colors", names(df))] == "blue",na.rm=T) == length(list.names), "blue", "grey"))
      }
      
      # save data for venn diagram
      tmp = df[df$colors=="red",] 
      write.csv(tmp, file.path(result_dir, paste0(comparison, "red.csv")), row.names = F)
      v.list[[paste0(comparison, "Red")]] = paste0(tmp$proximal_pA, tmp$distal_pA)
      tmp = df[df$colors=="blue",] 
      write.csv(tmp, file.path(result_dir, paste0(comparison, "blue.csv")), row.names = F)
      v.list[[paste0(comparison, "Blue")]] = paste0(tmp$proximal_pA, tmp$distal_pA)
      rm(tmp)
      
      # average the ratios
      x = sub("_$", "_Prox", comparison)
      y = sub("_$", "_Dis", comparison)
      df[, x] = rowMeans(df[, grep("\\d_Prox", names(df))], na.rm=T)
      df[, y] = rowMeans(df[, grep("\\d_Dis", names(df))], na.rm=T)
      
      # numbers at beginning of columns have to be removed
      names(df) = sub("^\\d+", "", names(df))
      x = sub("^\\d+", "", x)
      y = sub("^\\d+", "", y)
      # scatter plot ###plot colored points later
      layer1 = subset(df, colors == "grey")
      layer2 = subset(df, colors != "grey")
      
      # calculate annotation text
      blue_over_red = ifelse(sum(df$colors == "red") == 0, 
                             log2(sum(df$colors == "blue")+1) - log2(sum(df$colors == "red")+1),
                             log2(sum(df$colors == "blue")) - log2(sum(df$colors == "red")))
      
      log2_ratio = paste("log2(Blue/Red) = ", round(blue_over_red, 2), sep = "")
      
      print(paste0("Plotting ", y, " ~ ", x))
      p.list[[k]] = ggplot() +
        geom_point(data=layer1, aes_string(x=x, y=y), colour = "grey", shape = 1) + 
        geom_point(data=layer2, aes_string(x=x, y=y), colour = layer2$colors, shape = 1) +
        ggtitle(paste0("Grey:", sum(df$colors == "grey"), "|",
                       "Red:", sum(df$colors == "red"), "|",
                       "Blue:", sum(df$colors == "blue"), "\n",
                       log2_ratio)) 
      
      # centered scatter plot
      p2.list[[k]] = p.list[[k]] + xlim(-3,3) + ylim(-3,3) 
      
      k = k + 1
    }
    # scatter plot
    p.list = c(p.list, nrow = length(comparisons2)) 
    png(file.path(result_dir, paste0("Dis_Prox_mean_treatment_ratios_p=", p_threshold, ".png")), height=550*length(comparisons2), width=500*(k-1)/length(comparisons2))
    do.call(grid.arrange, p.list)
    dev.off()
    # centered scatter plot
    p2.list = c(p2.list, nrow = length(comparisons2)) 
    png(file.path(result_dir, paste0("Dis_Prox_mean_treatment_ratios_p=", p_threshold, "_centered.png")), height=550*length(comparisons2), width=500*(k-1)/length(comparisons2))
    do.call(grid.arrange, p2.list)
    dev.off()
    
    # venn diagram
    if(length(v.list) < 6 & length(v.list) > 2){
      require(VennDiagram)
      venn.diagram(
        x = v.list,
        filename = file.path(result_dir, paste0("Dis_Prox_treatment_ratios_venn_p=", p_threshold, ".png")),
        height=2000, width=2000,
        imagetype="png",
        cat.cex = 0.5,
        margin = 0.1,
        col = "transparent",
        fill = c("cornflowerblue", "green", "yellow", "darkorchid1")
      )
    }
  }
}


#### Treatment RED in each fraction vs aUTR len 
# Straight line and loess
if(length(treatments) > 1){
  p1.list = list()
  p2.list = list()
  p3.list = list()
  
  k = 1
  for(batch in batches){
    for(fraction in fractions){
      for(i in 1:(length(treatments)-1)){
        for(j in (i+1):length(treatments)){
          # only keep desirable comparisons
          if((exists("comparisons") && paste(treatments[j], treatments[i], sep="_") %in% comparisons) | 
             !exists("comparisons")){
            y_string = paste(treatments[j], "vs", treatments[i], fraction, batch, "RED", sep="_")
            if(y_string %in% names(pair.pA)){
              
              # calculate column names
              counts_1 = paste0(paste(strsplit(y_string, "_")[[1]][c(1,4,5)], collapse="_"), "_count")
              counts_2 = paste0(paste(strsplit(y_string, "_")[[1]][c(3,4,5)], collapse="_"), "_count")
              # filter pAs with low read numbers
              clean.df = pair.pA[rowSums(pair.pA[, grep(paste0(counts_1, "|", counts_2), names(pair.pA))] >= 5, na.rm = T) == 4, ]
              
              c = cor(clean.df[,"len_log10"], clean.df[, y_string])
              p_val = coefficients(summary(lm(clean.df[,"len"] ~ clean.df[, y_string])))[2,4]
              #paste0("r = ", round(c,2), ", p = ", sprintf("%.3G", p_val))
              
              c2 = cor(log10(clean.df[clean.df[,"CDS_size"] > 0,"CDS_size"]), clean.df[clean.df[,"CDS_size"] > 0, y_string])
              p_val2 = coefficients(summary(lm(log10(clean.df[, "CDS_size"]) ~ clean.df[, y_string])))[2,4]
              
              p1.list[[k]] = ggplot(data = clean.df, aes_string(x= "len", y = y_string)) + 
                geom_point(shape = 19, alpha = 0.3, size = 1, color = "blue") + 
                geom_smooth(method="lm", color = "red") + 
                scale_x_log10() + xlab("aUTR length") +
                ggtitle(paste0("r = ", round(c,2), ", p = ", sprintf("%.3G", p_val))) 
              
              p2.list[[k]] = ggplot(data = clean.df, aes_string(x= "len", y = y_string)) + 
                geom_point(shape = 19, alpha = 0.3, size = 1, color = "blue") + 
                geom_smooth(method="loess", color = "red") + 
                scale_x_log10() + xlab("aUTR length") +
                ggtitle(paste0("r = ", round(c,2), ", p = ", sprintf("%.3G", p_val))) 
              
              p3.list[[k]] = ggplot(data = clean.df, aes_string(x= "CDS_size", y = y_string)) + 
                geom_point(shape = 19, alpha = 0.3, size = 1, color = "blue") + 
                geom_smooth(method="loess", color = "red") + 
                scale_x_log10() + xlab("CDS Size") +
                ggtitle(paste0("r = ", round(c2, 2), ", p = ", sprintf("%.3G", p_val2))) 
              
              k = k + 1
            }
          }
        }
      }
    }
  }
  if(k > 1){
    nrow = ifelse(k <= 5, 1, (k-1) %/% 4 + 1)
    p1.list = c(p1.list, nrow = nrow) 
    png(file.path(result_dir, "APA_vs_len_line.png"), width=400*(k-1)/nrow, height=400*nrow)
    do.call(grid.arrange, p1.list)
    dev.off()
    p2.list = c(p2.list, nrow = nrow) 
    png(file.path(result_dir, "APA_vs_len_loess.png"), width=400*(k-1)/nrow, height=400*nrow)
    do.call(grid.arrange, p2.list)
    dev.off()
    
    p3.list = c(p3.list, nrow = nrow)
    png(file.path(result_dir, "APA_vs_CDS_size_loess.png"), width=400*(k-1)/nrow, height=400*nrow)
    do.call(grid.arrange, p3.list)
    dev.off()
  }
  
  # Treatment mean RED vs aUTR len bin
  require(tidyr)
  require(dplyr)
  require(magrittr)
  p.list = list()
  d.list = list()
  k = 1
  for(batch in batches){
    for(fraction in fractions){
      for(i in 1:(length(treatments)-1)){
        for (j in (i+1):length(treatments)){
          # only keep desirable comparisons
          if((exists("comparisons") && paste(treatments[j], treatments[i], sep="_") %in% comparisons) | 
             !exists("comparisons")){
            y_string =  paste(treatments[j], "vs", treatments[i], fraction, batch, "RED", sep="_")
            
            if(y_string %in% names(pair.pA)){
              
              # calculate column names
              counts_1 = paste0(paste(strsplit(y_string, "_")[[1]][c(1,4,5)], collapse="_"), "_count")
              counts_2 = paste0(paste(strsplit(y_string, "_")[[1]][c(3,4,5)], collapse="_"), "_count")
              # filter pAs with low read numbers
              df = pair.pA[rowSums(pair.pA[, grep(paste0(counts_1, "|", counts_2), names(pair.pA))] >= 5, na.rm = T) == 4, ]
              df$len_bin = cut(df$len, breaks = quantile(df$len,probs = seq(0, 1, 0.2), na.rm = T), include.lowest = T)
              
              ps = vector()
              for(bin in levels(df$len_bin)){
                p = wilcox.test(df[df$len_bin == bin, y_string], df[df$len_bin == levels(df$len_bin)[1], y_string])$p.value
                ps = c(ps, sprintf("%.2E", p))
              }
              
              df = df[, c(y_string, "len_bin")]
              names(df) = sub("len_bin", "aUTR_Size", names(df))
              
              df = gather(df, key=ScoreType, value=RED, -aUTR_Size) 
              df %<>% filter(!is.na(aUTR_Size), abs(RED) != Inf) %>%
                group_by(aUTR_Size, ScoreType) %>%
                dplyr::summarise(Mean=mean(RED), SD=sd(RED), SEM = sd(RED)/sqrt(n()), Num = n())
              
              # save min and max values for plotting
              if(k==1){
                min.val = min(df$Mean - df$SEM)
                max.val = max(df$Mean + df$SEM)
              }else{
                min.val = min(min.val, min(df$Mean - df$SEM))
                max.val = max(max.val, max(df$Mean + df$SEM))
              }
              
              # save data
              d.list[[y_string]] = df
              # plotting
              p.list[[k]] = ggplot(df, aes(x=aUTR_Size, y=Mean, group=ScoreType, color=ScoreType)) + 
                geom_line(size=2) +
                geom_point(size=8) + 
                theme(legend.position = "bottom") +
                xlab("aUTR size") + 
                geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=.4,
                              position=position_dodge(0.1)) +
                guides(colour=guide_legend(ncol=2)) +
                ggtitle(paste0("aUTR#: ", paste(df$Num, collapse="|"), "\np = ", paste0(ps, collapse = "|")))
              
              k = k + 1
              
            }
          }
        }
      }
    }
  }
  if(k > 1){
    p.list = lapply(p.list, function(p) p+ylim(min.val, max.val))
    nrow = ifelse(k <= 5, 1, length(batches))
    p.list = c(p.list, nrow = nrow)
    png(file.path(result_dir, "treatment_REDs_vs_len_spaghetti.png"), width=600*(k-1)/nrow, height=650*nrow)
    do.call(grid.arrange, p.list)
    dev.off()
    df = do.call(rbind, d.list)
    df = df[order(df$ScoreType),]
    write.csv(df, file.path(result_dir, "treatment_REDs_vs_len_spaghetti.csv"), row.names=F)
  }
  # Treatment mean mean RED vs aUTR bin spaghetti
  if(length(batches) > 1){
    p.list = list()
    d.list = list()
    k = 1
    for(fraction in fractions){
      for(i in 1:(length(treatments)-1)){
        for (j in (i+1):length(treatments)){
          # only keep desirable comparisons
          if((exists("comparisons") && paste(treatments[j], treatments[i], sep="_") %in% comparisons) | 
             !exists("comparisons")){
            y_strings =  paste(treatments[j], "vs", treatments[i], fraction, batches, "RED", sep="_")
            y_strings = y_strings[y_strings %in% names(pair.pA)]
            
            if(length(y_strings) > 1){
              # calculate column names
              count_columns_1 = sapply(strsplit(y_strings, "_"), function(x) paste0(paste(x[c(1,4,5)], collapse="_"), "_count"))
              count_columns_2 = sapply(strsplit(y_strings, "_"), function(x) paste0(paste(x[c(3,4,5)], collapse="_"), "_count"))
              pattern = c(count_columns_1, count_columns_2)
              count_columns = grep(paste0(pattern, collapse = "|"), names(pair.pA), value=T)
              # filter pAs with low read numbers
              df = pair.pA[rowSums(pair.pA[, count_columns] >= 5, na.rm = T) == length(count_columns), ]
              df$len_bin = cut(df$len, breaks = quantile(df$len,probs = seq(0, 1, 0.2), na.rm = T), include.lowest = T)
              
              y_string = sub(paste0("(", paste(batches, collapse="|"), ")_RED"), "mean_RED", y_strings[1])
              
              ps = vector()
              for(bin in levels(df$len_bin)){
                p = wilcox.test(df[df$len_bin == bin, y_string], df[df$len_bin == levels(df$len_bin)[1], y_string])$p.value
                ps = c(ps, sprintf("%.2E", p))
              }
              
              df = df[, c(y_string, "len_bin")]
              names(df) = sub("len_bin", "aUTR_Size", names(df))
              
              df = gather(df, key=ScoreType, value=RED, -aUTR_Size) 
              
              df %<>% group_by(ScoreType, aUTR_Size) %>%
                dplyr::summarise(Mean=mean(RED), SD=sd(RED), SEM = sd(RED)/sqrt(n()), Num = n())
              
              # save min and max values for plotting
              if(k==1){
                min.val = min(df$Mean - df$SEM)
                max.val = max(df$Mean + df$SEM)
              }else{
                min.val = min(min.val, min(df$Mean - df$SEM))
                max.val = max(max.val, max(df$Mean + df$SEM))
              }
              
              # save data
              d.list[[y_string]] = df
              # plotting
              p.list[[k]] = ggplot(df, aes(x=aUTR_Size, y=Mean, group=ScoreType, color=ScoreType)) + 
                geom_line(size=2) +
                geom_point(size=8) + 
                theme(legend.position = "bottom") +
                xlab("aUTR size") + 
                geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=.4,
                              position=position_dodge(0.1)) +
                guides(colour=guide_legend(ncol=2)) +
                ggtitle(paste0("aUTR#: ", paste(df$Num, collapse="|"), "\np = ", paste0(ps, collapse = "|")))
              
              k = k + 1
              
            }
          }
        }
      }
    }
    if(k > 1){
      p.list = lapply(p.list, function(p) p+ylim(min.val, max.val))
      nrow = ifelse(k <= 5, 1, length(batches))
      p.list = c(p.list, nrow = nrow)
      png(file.path(result_dir, "treatment_mean_REDs_vs_len_spaghetti.png"), width=600*(k-1)/nrow, height=675*nrow)
      do.call(grid.arrange, p.list)
      dev.off()
      df = do.call(rbind, d.list)
      df = df[order(df$ScoreType),]
      write.csv(df, file.path(result_dir, "treatment_mean_REDs_vs_len_spaghetti.csv"), row.names=F)
    }
  }
 
  # Treatment median RED vs aUTR len bin
  require(tidyr)
  require(dplyr)
  require(magrittr)
  p.list = list()
  d.list = list()
  k = 1
  for(batch in batches){
    for(fraction in fractions){
      for(i in 1:(length(treatments)-1)){
        for (j in (i+1):length(treatments)){
          # only keep desirable comparisons
          if((exists("comparisons") && paste(treatments[j], treatments[i], sep="_") %in% comparisons) | 
             !exists("comparisons")){
            y_string =  paste(treatments[j], "vs", treatments[i], fraction, batch, "RED", sep="_")
            
            if(y_string %in% names(pair.pA)){
              
              # calculate column names
              counts_1 = paste0(paste(strsplit(y_string, "_")[[1]][c(1,4,5)], collapse="_"), "_count")
              counts_2 = paste0(paste(strsplit(y_string, "_")[[1]][c(3,4,5)], collapse="_"), "_count")
              # filter pAs with low read numbers
              df = pair.pA[rowSums(pair.pA[, grep(paste0(counts_1, "|", counts_2), names(pair.pA))] >= 5, na.rm = T) == 4, ]
              df$len_bin = cut(df$len, breaks = quantile(df$len,probs = seq(0, 1, 0.2), na.rm = T), include.lowest = T)
              
              ps = vector()
              for(bin in levels(df$len_bin)){
                p = wilcox.test(df[df$len_bin == bin, y_string], df[df$len_bin == levels(df$len_bin)[1], y_string])$p.value
                ps = c(ps, sprintf("%.2E", p))
              }
              
              df = df[, c(y_string, "len_bin")]
              names(df) = sub("len_bin", "aUTR_Size", names(df))
              
              df = gather(df, key=ScoreType, value=RED, -aUTR_Size) 
              df %<>% filter(!is.na(aUTR_Size), abs(RED) != Inf) %>%
                group_by(aUTR_Size, ScoreType) %>%
                dplyr::summarise(Median=median(RED), MAD=mad(RED), Num = n())
              
              # save min and max values for plotting
              if(k==1){
                min.val = min(df$Median - df$MAD)
                max.val = max(df$Median + df$MAD)
              }else{
                min.val = min(min.val, min(df$Median - df$MAD))
                max.val = max(max.val, max(df$Median + df$MAD))
              }
              
              # save data
              d.list[[y_string]] = df
              # plotting
              p.list[[k]] = ggplot(df, aes(x=aUTR_Size, y=Median, group=ScoreType, color=ScoreType)) + 
                geom_line(size=2) +
                geom_point(size=8) + 
                theme(legend.position = "bottom") +
                xlab("aUTR size") + 
                geom_errorbar(aes(ymin=Median-MAD, ymax=Median+MAD), width=.4,
                              position=position_dodge(0.1)) +
                guides(colour=guide_legend(ncol=2)) +
                ggtitle(paste0("aUTR#: ", paste(df$Num, collapse="|"), "\np = ", paste0(ps, collapse = "|")))
              
              k = k + 1
              
            }
          }
        }
      }
    }
  }
  if(k > 1){
    p.list = lapply(p.list, function(p) p+ylim(min.val, max.val))
    nrow = ifelse(k <= 5, 1, length(batches))
    p.list = c(p.list, nrow = nrow)
    png(file.path(result_dir, "treatment_median_REDs_vs_len_spaghetti.png"), width=600*(k-1)/nrow, height=650*nrow)
    do.call(grid.arrange, p.list)
    dev.off()
    df = do.call(rbind, d.list)
    df = df[order(df$ScoreType),]
    write.csv(df, file.path(result_dir, "treatment_median_REDs_vs_len_spaghetti.csv"), row.names=F)
  }
  
  # Treatment RED vs relative aUTR relative_len bin
  require(tidyr)
  require(dplyr)
  require(magrittr)
  p.list = list()
  d.list = list()
  k = 1
  for(batch in batches){
    for(fraction in fractions){
      for(i in 1:(length(treatments)-1)){
        for (j in (i+1):length(treatments)){
          # only keep desirable comparisons
          if((exists("comparisons") && paste(treatments[j], treatments[i], sep="_") %in% comparisons) | 
             !exists("comparisons")){
            y_string =  paste(treatments[j], "vs", treatments[i], fraction, batch, "RED", sep="_")
            
            if(y_string %in% names(pair.pA)){
              # calculate column names
              counts_1 = paste0(paste(strsplit(y_string, "_")[[1]][c(1,4,5)], collapse="_"), "_count")
              counts_2 = paste0(paste(strsplit(y_string, "_")[[1]][c(3,4,5)], collapse="_"), "_count")
              # filter pAs with low read numbers
              df = pair.pA[rowSums(pair.pA[, grep(paste0(counts_1, "|", counts_2), names(pair.pA))] >= 5, na.rm = T) == 4, ]
              df$relative_len_bin = cut(df$relative_len, breaks = quantile(df$relative_len,probs = seq(0, 1, 0.2), na.rm = T), include.lowest = T)
              
              df = df[, c(y_string, "relative_len_bin")]
              names(df) = sub("relative_len_bin", "aUTR_Size", names(df))
              
              
              df = gather(df, key=ScoreType, value=RED, -aUTR_Size) 
              df %<>% filter(!is.na(aUTR_Size), abs(RED) != Inf) %>%
                group_by(aUTR_Size, ScoreType) %>%
                dplyr::summarise(Mean=mean(RED), SD=sd(RED), SEM = sd(RED)/sqrt(n()), Num = n())
              
              # save min and max values for plotting
              if(k==1){
                min.val = min(df$Mean - df$SEM)
                max.val = max(df$Mean + df$SEM)
              }else{
                min.val = min(min.val, min(df$Mean - df$SEM))
                max.val = max(max.val, max(df$Mean + df$SEM))
              }
              
              # save data
              d.list[[y_string]] = df
              # plotting
              p.list[[k]] = ggplot(df, aes(x=aUTR_Size, y=Mean, group=ScoreType, color=ScoreType)) + 
                geom_line(size=2) +
                geom_point(size=8) + 
                theme(legend.position = "bottom") +
                xlab("log2(aUTR size/cUTR size)") +
                geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=.4,
                              position=position_dodge(0.1)) +
                guides(colour=guide_legend(ncol=2)) +
                ggtitle(paste0("aUTR#: ", paste(df$Num, collapse="|")))
              
              k = k + 1
              
            }
          }
        }
      }
    }
  }
  if(k > 1){
    p.list = lapply(p.list, function(p) p+ylim(min.val, max.val))
    nrow = ifelse(k <= 5, 1, length(batches))
    p.list = c(p.list, nrow = nrow)
    png(file.path(result_dir, "treatment_REDs_vs_relative_len_spaghetti.png"), width=600*(k-1)/nrow, height=650*nrow)
    do.call(grid.arrange, p.list)
    dev.off()
    df = do.call(rbind, d.list)
    df = df[order(df$ScoreType),]
    write.csv(df, file.path(result_dir, "treatment_REDs_vs_relative_len_spaghetti.csv"), row.names=F)
  }
  # Treatment mean RED vs relative aUTR bin spaghetti
  if(length(batches) > 1){
    p.list = list()
    d.list = list()
    k = 1
    for(fraction in fractions){
      for(i in 1:(length(treatments)-1)){
        for (j in (i+1):length(treatments)){
          # only keep desirable comparisons
          if((exists("comparisons") && paste(treatments[j], treatments[i], sep="_") %in% comparisons) | 
             !exists("comparisons")){
            y_strings =  paste(treatments[j], "vs", treatments[i], fraction, batches, "RED", sep="_")
            y_strings = y_strings[y_strings %in% names(pair.pA)]
            
            if(length(y_strings) > 1){
              
              # calculate column names
              count_columns_1 = sapply(strsplit(y_strings, "_"), function(x) paste0(paste(x[c(1,4,5)], collapse="_"), "_count"))
              count_columns_2 = sapply(strsplit(y_strings, "_"), function(x) paste0(paste(x[c(3,4,5)], collapse="_"), "_count"))
              pattern = c(count_columns_1, count_columns_2)
              count_columns = grep(paste0(pattern, collapse = "|"), names(pair.pA), value=T)
              # filter pAs with low read numbers
              df = pair.pA[rowSums(pair.pA[, count_columns] >= 5, na.rm = T) == length(count_columns), ]
              df$relative_len_bin = cut(df$relative_len, breaks = quantile(df$relative_len,probs = seq(0, 1, 0.2), na.rm = T), include.lowest = T)
              
              y_string = sub(paste0("(", paste(batches, collapse="|"), ")_RED"), "mean_RED", y_strings[1])
              df = df[, c(y_string, "relative_len_bin")]
              names(df) = sub("relative_len_bin", "aUTR_Size", names(df))
              
              df = gather(df, key=ScoreType, value=RED, -aUTR_Size) 
              
              df %<>% group_by(ScoreType, aUTR_Size) %>%
                dplyr::summarise(Mean=mean(RED), SD=sd(RED), SEM = sd(RED)/sqrt(n()), Num = n())
              
              # save min and max values for plotting
              if(k==1){
                min.val = min(df$Mean - df$SEM)
                max.val = max(df$Mean + df$SEM)
              }else{
                min.val = min(min.val, min(df$Mean - df$SEM))
                max.val = max(max.val, max(df$Mean + df$SEM))
              }
              
              # save data
              d.list[[y_string]] = df
              # plotting
              p.list[[k]] = ggplot(df, aes(x=aUTR_Size, y=Mean, group=ScoreType, color=ScoreType)) + 
                geom_line(size=2) +
                geom_point(size=8) + 
                theme(legend.position = "bottom") +
                xlab("log2(aUTR size/cUTR size)") + 
                geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=.4,
                              position=position_dodge(0.1)) +
                guides(colour=guide_legend(ncol=2)) +
                ggtitle(paste0("aUTR#: ", paste(df$Num, collapse="|")))
              
              k = k + 1
              
            }
          }
        }
      }
    }
    if(k > 1){
      p.list = lapply(p.list, function(p) p+ylim(min.val, max.val))
      nrow = ifelse(k <= 5, 1, length(batches))
      p.list = c(p.list, nrow = nrow)
      png(file.path(result_dir, "treatment_mean_REDs_vs_relative_len_spaghetti.png"), width=600*(k-1)/nrow, height=650*nrow)
      do.call(grid.arrange, p.list)
      dev.off()
      df = do.call(rbind, d.list)
      df = df[order(df$ScoreType),]
      write.csv(df, file.path(result_dir, "treatment_mean_REDs_vs_relative_len_spaghetti.csv"), row.names=F)
    }
  }
}


#### Fraction RED across fractions in each treatment vs aUTR len 
# Straight line and loess
if(length(fractions) > 1){
  p1.list = list()
  p2.list = list()
  p3.list = list()
  k = 1
  for(batch in batches){
    for(treatment in treatments){
      for(i in 1:(length(fractions)-1)){
        for (j in (i+1):length(fractions)){
          # only keep desirable comparisons
          if((exists("comparisons") && paste(fractions[j], fractions[i], sep="_") %in% comparisons) | 
             !exists("comparisons")){
            y_string =  paste(treatment, fractions[j], "vs", fractions[i], batch, "RED", sep="_")
            
            if(y_string %in% names(pair.pA)){
              
              # calculate column names
              counts_1 = paste0(paste(strsplit(y_string, "_")[[1]][c(1,2,5)], collapse="_"), "_count")
              counts_2 = paste0(paste(strsplit(y_string, "_")[[1]][c(1,4,5)], collapse="_"), "_count")
              # filter pAs with low read numbers
              clean.df = pair.pA[rowSums(pair.pA[, grep(paste0(counts_1, "|", counts_2), names(pair.pA))] >= 5, na.rm = T) == 4, ]
              
              c = cor(clean.df[,"len_log10"], clean.df[, y_string], use="pairwise.complete.obs")
              p_val = coefficients(summary(lm(clean.df[,"len"] ~ clean.df[, y_string])))[2,4]
              
              c2 = cor(log10(clean.df[clean.df[,"CDS_size"] > 0,"CDS_size"]), clean.df[clean.df[,"CDS_size"] > 0, y_string])
              p_val2 = coefficients(summary(lm(log10(clean.df[, "CDS_size"]) ~ clean.df[, y_string])))[2,4]

              
              p1.list[[k]] = ggplot(data = clean.df, aes_string(x= "len", y = y_string)) + 
                geom_point(shape = 19, alpha = 0.3, size = 1, color = "blue") + #annotation_custom(my_grob) +
                geom_smooth(method="lm", color = "red") + 
                scale_x_log10() + xlab("aUTR length") +
                ggtitle(paste0("r = ", round(c,2), ", p = ", sprintf("%.3G", p_val))) 
              
              p2.list[[k]] = ggplot(data = clean.df, aes_string(x= "len", y = y_string)) + 
                geom_point(shape = 19, alpha = 0.3, size = 1, color = "blue") + 
                geom_smooth(method="loess", color = "red") + 
                scale_x_log10() + xlab("aUTR length") +
                ggtitle(paste0("r = ", round(c,2), ", p = ", sprintf("%.3G", p_val))) 
              
              p3.list[[k]] = ggplot(data = clean.df, aes_string(x= "CDS_size", y = y_string)) + 
                geom_point(shape = 19, alpha = 0.3, size = 1, color = "blue") + #annotation_custom(my_grob) +
                geom_smooth(method="loess", color = "red") + 
                scale_x_log10() + xlab("CDS Size") +
                ggtitle(paste0("r = ", round(c2, 2), ", p = ", sprintf("%.3G", p_val2))) 
              
              k = k + 1
            }
          }
        }
      }
    }
  }
  if(k > 1){
    nrow = ifelse(k <= 5, 1, (k-1) %/% 4 + 1)
    p1.list = c(p1.list, nrow = nrow) 
    png(file.path(result_dir, "RED_vs_len_line.png"), width=400*(k-1)/nrow, height=400*nrow)
    do.call(grid.arrange, p1.list)
    dev.off()
    #p2.list = lapply(p2.list, function(p) p+ylim(min.val, max.val))
    #p2.list = c(p2.list, nrow = length(treatments)) 
    p2.list = c(p2.list, nrow = nrow)
    png(file.path(result_dir, "RED_vs_len_loess.png"), width=400*(k-1)/nrow, height=400*nrow)
    do.call(grid.arrange, p2.list)
    dev.off()
    
    p3.list = c(p3.list, nrow = nrow)
    png(file.path(result_dir, "RED_vs_CDS_size_loess.png"), width=400*(k-1)/nrow, height=400*nrow)
    do.call(grid.arrange, p3.list)
    dev.off()
  }
  
  # Fraction RED vs aUTR len bin
  require(tidyr)
  require(dplyr)
  require(magrittr)
  p.list = list()
  d.list = list()
  k = 1
  for(batch in batches){
    for(treatment in treatments){
      for(i in 1:(length(fractions)-1)){
        for (j in (i+1):length(fractions)){
          # only keep desirable comparisons
          if((exists("comparisons") && paste(fractions[j], fractions[i], sep="_") %in% comparisons) | 
             !exists("comparisons")){
            y_string =  paste(treatment, fractions[j], "vs", fractions[i], batch, "RED", sep="_")
            
            if(y_string %in% names(pair.pA)){
              
              # calculate column names
              counts_1 = paste0(paste(strsplit(y_string, "_")[[1]][c(1,2,5)], collapse="_"), "_count")
              counts_2 = paste0(paste(strsplit(y_string, "_")[[1]][c(1,4,5)], collapse="_"), "_count")
              # filter pAs with low read numbers
              df = pair.pA[rowSums(pair.pA[, grep(paste0(counts_1, "|", counts_2), names(pair.pA))] >= 5, na.rm = T) == 4, ]
              df$len_bin = cut(df$len, breaks = quantile(df$len,probs = seq(0, 1, 0.2), na.rm = T), include.lowest = T)
              
              df = df[, c(y_string, "len_bin")]
              ps = vector()
              for(bin in levels(df$len_bin)){
                p = wilcox.test(df[df$len_bin == bin, y_string], df[df$len_bin == levels(df$len_bin)[1], y_string])$p.value
                ps = c(ps, sprintf("%.2E", p))
              }
              names(df) = sub("len_bin", "aUTR_Size", names(df))
              
              
              df = gather(df, key=ScoreType, value=RED, -aUTR_Size) 
              df %<>% filter(!is.na(aUTR_Size), abs(RED) != Inf) %>%
                group_by(aUTR_Size, ScoreType) %>%
                dplyr::summarise(Mean=mean(RED), SD=sd(RED), SEM = sd(RED)/sqrt(n()), Num = n())
              
              # save min and max values for plotting
              if(k==1){
                min.val = min(df$Mean - df$SEM)
                max.val = max(df$Mean + df$SEM)
              }else{
                min.val = min(min.val, min(df$Mean - df$SEM))
                max.val = max(max.val, max(df$Mean + df$SEM))
              }
              
              # save data
              d.list[[y_string]] = df
              # plotting
              p.list[[k]] = ggplot(df, aes(x=aUTR_Size, y=Mean, group=ScoreType, color=ScoreType)) + 
                geom_line(size=2) +
                geom_point(size=8) + 
                theme(legend.position = "bottom") +
                xlab("aUTR size") + 
                geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=.4,
                              position=position_dodge(0.1)) +
                guides(colour=guide_legend(ncol=2)) +
                ggtitle(paste0("aUTR#: ", paste(df$Num, collapse="|"), "\np = ", paste0(ps, collapse = "|")))
              
              k = k + 1
              
            }
          }
        }
      }
    }
  }
  if(k > 1){
    p.list = lapply(p.list, function(p) p+ylim(min.val, max.val))
    nrow = ifelse(k <= 5, 1, length(batches))
    p.list = c(p.list, nrow = nrow)
    png(file.path(result_dir, "fraction_REDs_vs_len_spaghetti.png"), width=600*(k-1)/nrow, height=650*nrow)
    do.call(grid.arrange, p.list)
    dev.off()
    df = do.call(rbind, d.list)
    df = df[order(df$ScoreType),]
    write.csv(df, file.path(result_dir, "fraction_REDs_vs_len_spaghetti.csv"), row.names=F)
    
    # replot, putting the same score type (fraction2_fraction1) in the same subfigure
    df = mutate(df, ScoreClass = sub(".+?_(.+)$", "\\1", ScoreType), Treatment = sub("(.+?)_.+$", "\\1", ScoreType))
    p = ggplot(df, aes(x=aUTR_Size, y=Mean, group=ScoreType, color=Treatment)) + 
      geom_line(size=2) +
      geom_point(size=8) + 
      theme(legend.position = "bottom") +
      xlab("aUTR size") + 
      geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=.4) +
      facet_grid(.~ScoreClass)
    png(file.path(result_dir, "fraction_REDs_vs_len_spaghetti_facet_by_ScoreClass.png"), 
        width=700*length(unique(df$ScoreClass)), height=750)
    print(p)
    dev.off()
  }
  # Fraction mean RED vs aUTR bin spaghetti
  if(length(batches) > 1){
    p.list = list()
    d.list = list()
    k = 1
    for(treatment in treatments){
      for(i in 1:(length(fractions)-1)){
        for (j in (i+1):length(fractions)){
          # only keep desirable comparisons
          if((exists("comparisons") && paste(fractions[j], fractions[i], sep="_") %in% comparisons) | 
             !exists("comparisons")){
            y_strings =  paste(treatment, fractions[j], "vs", fractions[i], batches, "RED", sep="_")
            y_strings = y_strings[y_strings %in% names(pair.pA)]
            
            if(length(y_strings) > 1){
              
              # calculate column names
              count_columns_1 = sapply(strsplit(y_strings, "_"), function(x) paste0(paste(x[c(1,2,5)], collapse="_"), "_count"))
              count_columns_2 = sapply(strsplit(y_strings, "_"), function(x) paste0(paste(x[c(1,4,5)], collapse="_"), "_count"))
              pattern = c(count_columns_1, count_columns_2)
              count_columns = grep(paste0(pattern, collapse = "|"), names(pair.pA), value=T)
              # filter pAs with low read numbers
              df = pair.pA[rowSums(pair.pA[, count_columns] >= 5, na.rm = T) == length(count_columns), ]
              df$len_bin = cut(df$len, breaks = quantile(df$len,probs = seq(0, 1, 0.2), na.rm = T), include.lowest = T)
              
              y_string = sub(paste0("(", paste(batches, collapse="|"), ")_RED"), "mean_RED", y_strings[1])
              ps = vector()
              for(bin in levels(df$len_bin)){
                p = wilcox.test(df[df$len_bin == bin, y_string], df[df$len_bin == levels(df$len_bin)[1], y_string])$p.value
                ps = c(ps, sprintf("%.2E", p))
              }
              df = df[, c(y_string, "len_bin")]
              names(df) = sub("len_bin", "aUTR_Size", names(df))
              
              df = gather(df, key=ScoreType, value=RED, -aUTR_Size) 
              
              df %<>% group_by(ScoreType, aUTR_Size) %>%
                dplyr::summarise(Mean=mean(RED), SD=sd(RED), SEM = sd(RED)/sqrt(n()), Num = n())
              
              # save min and max values for plotting
              if(k==1){
                min.val = min(df$Mean - df$SEM)
                max.val = max(df$Mean + df$SEM)
              }else{
                min.val = min(min.val, min(df$Mean - df$SEM))
                max.val = max(max.val, max(df$Mean + df$SEM))
              }
              
              # save data
              d.list[[y_string]] = df
              # plotting
              p.list[[k]] = ggplot(df, aes(x=aUTR_Size, y=Mean, group=ScoreType, color=ScoreType)) + 
                geom_line(size=2) +
                geom_point(size=8) + 
                theme(legend.position = "bottom") +
                xlab("aUTR size") + 
                geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=.4,
                              position=position_dodge(0.1)) +
                guides(colour=guide_legend(ncol=2)) +
                ggtitle(paste0("aUTR#: ", paste(df$Num, collapse="|"), "\np = ", paste0(ps, collapse = "|")))
              
              k = k + 1
              
            }
          }
        }
      }
    }
    if(k > 1){
      p.list = lapply(p.list, function(p) p+ylim(min.val, max.val))
      nrow = ifelse(k <= 5, 1, length(batches))
      p.list = c(p.list, nrow = nrow)
      png(file.path(result_dir, "fraction_mean_REDs_vs_len_spaghetti.png"), width=600*(k-1)/nrow, height=650*nrow)
      do.call(grid.arrange, p.list)
      dev.off()
      df = do.call(rbind, d.list)
      df = df[order(df$ScoreType),]
      write.csv(df, file.path(result_dir, "fraction_mean_REDs_vs_len_spaghetti.csv"), row.names=F)
    }
  }
  
  # Fraction median RED vs aUTR len bin
  require(tidyr)
  require(dplyr)
  require(magrittr)
  p.list = list()
  d.list = list()
  k = 1
  for(batch in batches){
    for(treatment in treatments){
      for(i in 1:(length(fractions)-1)){
        for (j in (i+1):length(fractions)){
          # only keep desirable comparisons
          if((exists("comparisons") && paste(fractions[j], fractions[i], sep="_") %in% comparisons) | 
             !exists("comparisons")){
            y_string =  paste(treatment, fractions[j], "vs", fractions[i], batch, "RED", sep="_")
            
            if(y_string %in% names(pair.pA)){
              
              # calculate column names
              counts_1 = paste0(paste(strsplit(y_string, "_")[[1]][c(1,2,5)], collapse="_"), "_count")
              counts_2 = paste0(paste(strsplit(y_string, "_")[[1]][c(1,4,5)], collapse="_"), "_count")
              # filter pAs with low read numbers
              df = pair.pA[rowSums(pair.pA[, grep(paste0(counts_1, "|", counts_2), names(pair.pA))] >= 5, na.rm = T) == 4, ]
              df$len_bin = cut(df$len, breaks = quantile(df$len,probs = seq(0, 1, 0.2), na.rm = T), include.lowest = T)
              
              df = df[, c(y_string, "len_bin")]
              ps = vector()
              for(bin in levels(df$len_bin)){
                p = wilcox.test(df[df$len_bin == bin, y_string], df[df$len_bin == levels(df$len_bin)[1], y_string])$p.value
                ps = c(ps, sprintf("%.2E", p))
              }
              names(df) = sub("len_bin", "aUTR_Size", names(df))
              
              
              df = gather(df, key=ScoreType, value=RED, -aUTR_Size) 
              df %<>% filter(!is.na(aUTR_Size), abs(RED) != Inf) %>%
                group_by(aUTR_Size, ScoreType) %>%
                dplyr::summarise(Median=median(RED), MAD=mad(RED), Num = n())
              
              # save min and max values for plotting
              if(k==1){
                min.val = min(df$Median - df$MAD)
                max.val = max(df$Median + df$MAD)
              }else{
                min.val = min(min.val, min(df$Median - df$MAD))
                max.val = max(max.val, max(df$Median + df$MAD))
              }
              
              # save data
              d.list[[y_string]] = df
              # plotting
              p.list[[k]] = ggplot(df, aes(x=aUTR_Size, y=Median, group=ScoreType, color=ScoreType)) + 
                geom_line(size=2) +
                geom_point(size=8) + 
                theme(legend.position = "bottom") +
                xlab("aUTR size") + 
                geom_errorbar(aes(ymin=Median-MAD, ymax=Median+MAD), width=.4,
                              position=position_dodge(0.1)) +
                guides(colour=guide_legend(ncol=2)) +
                ggtitle(paste0("aUTR#: ", paste(df$Num, collapse="|"), "\np = ", paste0(ps, collapse = "|")))
              
              k = k + 1
              
            }
          }
        }
      }
    }
  }
  if(k > 1){
    p.list = lapply(p.list, function(p) p+ylim(min.val, max.val))
    nrow = ifelse(k <= 5, 1, length(batches))
    p.list = c(p.list, nrow = nrow)
    png(file.path(result_dir, "fraction_median_REDs_vs_len_spaghetti.png"), width=600*(k-1)/nrow, height=650*nrow)
    do.call(grid.arrange, p.list)
    dev.off()
    df = do.call(rbind, d.list)
    df = df[order(df$ScoreType),]
    write.csv(df, file.path(result_dir, "fraction_median_REDs_vs_len_spaghetti.csv"), row.names=F)
    
    # replot, putting the same score type (fraction2_fraction1) in the same subfigure
    df = mutate(df, ScoreClass = sub(".+?_(.+)$", "\\1", ScoreType), Treatment = sub("(.+?)_.+$", "\\1", ScoreType))
    p = ggplot(df, aes(x=aUTR_Size, y=Median, group=ScoreType, color=Treatment)) + 
      geom_line(size=2) +
      geom_point(size=8) + 
      theme(legend.position = "bottom") +
      xlab("aUTR size") + 
      geom_errorbar(aes(ymin=Median-MAD, ymax=Median+MAD), width=.4) +
      facet_grid(.~ScoreClass)
    png(file.path(result_dir, "fraction_median_REDs_vs_len_spaghetti_facet_by_ScoreClass.png"), 
        width=700*length(unique(df$ScoreClass)), height=750)
    print(p)
    dev.off()
  }
  
  # Fraction mean RED vs relative aUTR bin spaghetti
  if(length(batches) > 1){
    p.list = list()
    d.list = list()
    k = 1
    for(treatment in treatments){
      for(i in 1:(length(fractions)-1)){
        for (j in (i+1):length(fractions)){
          # only keep desirable comparisons
          if((exists("comparisons") && paste(fractions[j], fractions[i], sep="_") %in% comparisons) | 
             !exists("comparisons")){
            y_strings =  paste(treatment, fractions[j], "vs", fractions[i], batches, "RED", sep="_")
            y_strings = y_strings[y_strings %in% names(pair.pA)]
            
            if(length(y_strings) > 1){
              # calculate column names
              count_columns_1 = sapply(strsplit(y_strings, "_"), function(x) paste0(paste(x[c(1,2,5)], collapse="_"), "_count"))
              count_columns_2 = sapply(strsplit(y_strings, "_"), function(x) paste0(paste(x[c(1,4,5)], collapse="_"), "_count"))
              pattern = c(count_columns_1, count_columns_2)
              count_columns = grep(paste0(pattern, collapse = "|"), names(pair.pA), value=T)
              # filter pAs with low read numbers
              df = pair.pA[rowSums(pair.pA[, count_columns] >= 5, na.rm = T) == length(count_columns), ]
              df$relative_len_bin = cut(df$relative_len, breaks = quantile(df$relative_len,probs = seq(0, 1, 0.2), na.rm = T), include.lowest = T)
              
              y_string = sub(paste0("(", paste(batches, collapse="|"), ")_RED"), "mean_RED", y_strings[1])
              df = df[, c(y_string, "relative_len_bin")]
              names(df) = sub("relative_len_bin", "aUTR_Size", names(df))
              
              df = gather(df, key=ScoreType, value=RED, -aUTR_Size) 
              
              df %<>% group_by(ScoreType, aUTR_Size) %>%
                dplyr::summarise(Mean=mean(RED), SD=sd(RED), SEM = sd(RED)/sqrt(n()), Num = n())
              
              # save min and max values for plotting
              if(k==1){
                min.val = min(df$Mean - df$SEM)
                max.val = max(df$Mean + df$SEM)
              }else{
                min.val = min(min.val, min(df$Mean - df$SEM))
                max.val = max(max.val, max(df$Mean + df$SEM))
              }
              
              # save data
              d.list[[y_string]] = df
              # plotting
              p.list[[k]] = ggplot(df, aes(x=aUTR_Size, y=Mean, group=ScoreType, color=ScoreType)) + 
                geom_line(size=2) +
                geom_point(size=8) + 
                theme(legend.position = "bottom") +
                xlab("log2(aUTR size/cUTR size)") +
                geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=.4,
                              position=position_dodge(0.1)) +
                guides(colour=guide_legend(ncol=2)) +
                ggtitle(paste0("aUTR#: ", paste(df$Num, collapse="|")))
              
              k = k + 1
              
            }
          }
        }
      }
    }
    if(k > 1){
      p.list = lapply(p.list, function(p) p+ylim(min.val, max.val))
      nrow = ifelse(k <= 5, 1, length(batches))
      p.list = c(p.list, nrow = nrow)
      png(file.path(result_dir, "fraction_mean_REDs_vs_relative_len_spaghetti.png"), width=600*(k-1)/nrow, height=650*nrow)
      do.call(grid.arrange, p.list)
      dev.off()
      df = do.call(rbind, d.list)
      df = df[order(df$ScoreType),]
      write.csv(df, file.path(result_dir, "fraction_mean_REDs_vs_relative_len_spaghetti.csv"), row.names=F)
    }
  }
}