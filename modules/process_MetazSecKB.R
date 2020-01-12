# This script will download the database from http://proteomics.ysu.edu/secretomes/animal/index.php
# See Meinken J, Walker G, Cooper CR, Min XJ (2015) 
# MetazSecKB: the human and animal secretome and subcellular proteome knowledgebase. 
# Database, bav077

database_file = file.path(SHARED_DATA_DIR, "all_metazoa_protein_subloc.txt")

if(!file.exists(database_file)){
  url = "http://proteomics.ysu.edu/publication/data/MetazSecKB/all_metazoa_protein_subloc.txt"
  download.file(url, database_file, quiet = F)
}

df = read.csv(database_file, sep = "\t", stringsAsFactors = F, header = F)
names(df) = c("uniprot_id", "localization")

library(biomaRt)

if(species == "mouse"){
  ensembl_dataset = useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
}else if(species == "human"){
  ensembl_dataset = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
}

# Query biomart in minibatches (nrow(df) is 4080818, too big for biomaRt in one batch) 
result = data.frame()
for(start in seq(1, nrow(df), 10000)){ 
  end = min(start + 10000 - 1, nrow(df))
  
  if(species == "mouse"){
    one_batch = biomaRt::select(ensembl_dataset, 
                         keytype = "uniprotswissprot",
                         keys = df$uniprot_id[start:end],
                         columns = c("uniprotswissprot", "mgi_symbol", "entrezgene"))
  }
  
  if(species == "human"){
    one_batch = biomaRt::select(ensembl_dataset, 
                         keytype = "uniprotswissprot",
                         keys = df$uniprot_id[start:end],
                         columns = c("uniprotswissprot", "hgnc_symbol", "entrezgene_id"))
  }
  
  result = rbind(result, one_batch)
}

names(result) = sub("mgi_symbol|hgnc_symbol", "gene_symbol", names(result))

# Process records
require(dplyr)
df = merge(df, result, by.x = "uniprot_id", by.y = "uniprotswissprot") %>%
  dplyr::select(gene_symbol, localization) %>%
  mutate(gene_symbol = Hmisc::capitalize(tolower(gene_symbol))) %>%
  distinct() %>%
  group_by(gene_symbol) %>%
  mutate(localizations = paste0(localization, collapse = ", ")) %>%
  dplyr::select(-localization) %>%
  ungroup() %>%
  mutate(ER = grepl("ER", localizations)) %>%
  mutate(Golgi = grepl("Golgi", localizations)) %>%
  mutate(Membrane = grepl("[Mm]embrane", localizations)) %>%
  mutate(Nucleus = grepl("[Nn]ucleus",localizations)
                 &!grepl("[Mm]embrane",localizations)
                 &!grepl("Secreted",localizations)) %>%
  mutate(Cytoplasm = grepl("[Cc]ytoplasm",localizations)
                   &!grepl("[Mm]embrane",localizations)
                   &!grepl("Secreted",localizations)) %>%
  mutate(Mitochondria = grepl("[Mm]itochondria", localizations)) %>%
  mutate(Lysosome = grepl("[Ll]ysosome", localizations)) %>%
  mutate(Peroxisome = grepl("[Pp]eroxisome", localizations)) %>%
  mutate(Cytoskeleton = grepl("[Cc]ytoskeleton",localizations)
                      &!grepl("[Mm]embrane",localizations)
                      &!grepl("Secreted",localizations)) %>%
  mutate(GPI_Anchored = grepl("GPI anchored", localizations)) %>%
  mutate(Secreted = grepl("Secreted \\(curated\\)", localizations)) %>%
  mutate(Secreted_likely = grepl("Secreted \\(likely\\)", localizations)) %>%
  mutate(Secreted_highlylikely = grepl("Secreted \\(highly likely\\)", localizations)) %>%
  mutate(Secreted_weaklylikely = grepl("Secreted \\(weakly likely\\)", localizations)) 

df[df == TRUE] = 1
df[df == FALSE] = 0

write.csv(df, file.path(SHARED_DATA_DIR, paste0(species, "_protein_localizations.csv")),
          row.names=F)
