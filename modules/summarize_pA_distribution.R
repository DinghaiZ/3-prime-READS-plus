# input: pA.df, read number cutoff
# outputs: 
# A table containing the number of pA sites for each sample, grouped by pA.df$region
# A table containing the number of pA sites for each sample, grouped by pA.df$pA_type
# A table containing the number of reads for each sample, grouped by pA.df$region
# A table containing the number of reads for each sample, grouped by pA.df$pA_type
# Stacked bar plots showing percentage for each of the above.

require(dplyr)
require(ggplot2)
require(scales)

min_count_ = 5


sample_counts = sort(grep("_count$", names(pA.df), value=T))
df = pA.df[, c("region", "pA_type", sample_counts)]

# Collapse pAs mapped to multiple features
df$region[df$region == "3UTR|CDS"] = "CDS"
df$region[df$region == "3UTR|intron"] = "intron"
df$region[df$region == "3UTR|intron|CDS"] = "CDS"
df$region[df$region == "intron|5UTR"] = "5UTR"
df$region[df$region == "intron|CDS"] = "CDS"
df$region[df$region == "CDS|5UTR"] = "CDS"
df$region[df$region == "intron|CDS|5UTR"] = "CDS"

df$region[df$region == "3UTR"] = "UTR3"
df$region[df$region == "5UTR"] = "UTR5"


# Initialize counting matrixes
region_summary = matrix(nrow = 0, ncol = length(names(table(df$region))))
colnames(region_summary) = names(table(df$region))

region_read_summary = matrix(ncol = 0, nrow = length(names(table(df$region))))
rownames(region_read_summary) = names(table(df$region))

pA_type_summary = matrix(nrow = 0, ncol = length(names(table(df$pA_type))))
colnames(pA_type_summary) = names(table(df$pA_type))

pA_type_read_summary = matrix(ncol = 0, nrow = length(names(table(df$pA_type))))
rownames(pA_type_read_summary) = names(table(df$pA_type))

# Go through each sample, filter out pAs with low read numbers, and summarize
for(sample_count in sample_counts){
  df_ = df[, c("region", "pA_type", sample_count)]
  df_ = df_[df_[, sample_count] >= min_count_, ]
  
  # summarize numbers of pA site
  region_summary = rbind(region_summary, table(df_$region))
  pA_type_summary = rbind(pA_type_summary, table(df_$pA_type))
  
  # summarize read numbers 
  col = group_by(df_, region) %>%
    summarize_at(vars(sample_count), sum)
  col = col[, 2]
  region_read_summary = cbind(region_read_summary, col)
  
  # summarize read numbers 
  col = group_by(df_, pA_type) %>%
    summarize_at(vars(sample_count), sum)
  col = col[, 2]
  pA_type_read_summary = cbind(pA_type_read_summary, col)
  
}

# Add row names
samples_ = sub("_count$", "", sample_counts)
rownames(region_summary) = samples_
rownames(pA_type_summary) = samples_

region_summary = t(region_summary)
pA_type_summary = t(pA_type_summary)

# Modify colum names
colnames(region_read_summary) = sub("_count", "_reads", colnames(region_read_summary))
colnames(pA_type_read_summary) = sub("_count", "_reads", colnames(pA_type_read_summary))

# Write to Excel
require(openxlsx)
wb <- createWorkbook(creator = "Dinghai")
addWorksheet(wb, "pA_region_type_summary")
writeData(wb, sheet = 1, region_summary, rowNames = T)
writeData(wb, sheet = 1, pA_type_summary, startRow = nrow(region_summary)+3, rowNames = T)
writeData(wb, sheet = 1, region_read_summary, startRow = nrow(pA_type_summary) + nrow(region_summary)+5, rowNames = T)
writeData(wb, sheet = 1, pA_type_read_summary, startRow = nrow(pA_type_summary) + nrow(region_summary) + nrow(region_read_summary) + 7, rowNames = T)
saveWorkbook(wb, file.path(result_dir, "pA_region_type_summary.xlsx"), overwrite = TRUE)


# Plotting
region_summary = as.data.frame(region_summary)
region_summary$region = rownames(region_summary)
region_summary = gather(region_summary, key=Sample, value=num_pA, -region)
region_summary$region = factor(region_summary$region, levels = c("intergenic", "ncRNA", "UA", "UTR5", "CDS", "intron", "UTR3"))
png(file.path(result_dir, "num_pA_region_summary.png"), 800, length(samples_)*50)
p = ggplot(data=region_summary, aes(x=Sample, y=num_pA, fill = region)) + 
  geom_col(position = "fill") + 
  scale_y_continuous(labels = percent_format()) +
  coord_flip()
print(p)
dev.off()

pA_type_summary = as.data.frame(pA_type_summary)
pA_type_summary$pA_type = rownames(pA_type_summary)
pA_type_summary = gather(pA_type_summary, key=Sample, value=num_pA, -pA_type)
pA_type_summary$pA_type = factor(pA_type_summary$pA_type, levels = c("Intergenic", "ncRNA", "UA", "E", "I", "S", "F", "M", "L"))
png(file.path(result_dir, "num_pA_type_summary.png"), 800, length(samples_)*50)
p = ggplot(data=pA_type_summary, aes(x=Sample, y=num_pA, fill = pA_type)) + 
  geom_col(position = "fill") + 
  scale_y_continuous(labels = percent_format()) + 
  coord_flip()
print(p)
dev.off()

region_read_summary = as.data.frame(region_read_summary)
region_read_summary$region = rownames(region_read_summary)
region_read_summary = gather(region_read_summary, key=Sample, value=num_read, -region)
region_read_summary$region = factor(region_read_summary$region, levels = c("intergenic", "ncRNA", "UA", "UTR5", "CDS", "intron", "UTR3"))
png(file.path(result_dir, "num_read_region_summary.png"), 800, length(samples_)*50)
p = ggplot(data=region_read_summary, aes(x=Sample, y=num_read, fill = region)) + 
  geom_col(position = "fill") + 
  scale_y_continuous(labels = percent_format()) +
  coord_flip()
print(p)
dev.off()


pA_type_read_summary = as.data.frame(pA_type_read_summary)
pA_type_read_summary$pA_type = rownames(pA_type_read_summary)
pA_type_read_summary = gather(pA_type_read_summary, key=Sample, value=num_read, -pA_type)
pA_type_read_summary$pA_type = factor(pA_type_read_summary$pA_type, levels = c("Intergenic", "ncRNA", "UA", "E", "I", "S", "F", "M", "L"))
png(file.path(result_dir, "num_read_pA_type_summary.png"), 800, length(samples_)*50)
p = ggplot(data=pA_type_read_summary, aes(x=Sample, y=num_read, fill = pA_type)) + 
  geom_col(position = "fill") + 
  scale_y_continuous(labels = percent_format()) +
  coord_flip()
print(p)
dev.off()

