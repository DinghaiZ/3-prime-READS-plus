# This script will download the database from http://proteomics.ysu.edu/secretomes/animal/index.php
# See Meinken J, Walker G, Cooper CR, Min XJ (2015) 
# MetazSecKB: the human and animal secretome and subcellular proteome knowledgebase. 
# Database, bav077

destfile = file.path(SHARED_DATA_DIR, "all_metazoa_protein_subloc.txt")

if(!file.exists(destfile)){
  url = "http://proteomics.ysu.edu/publication/data/MetazSecKB/all_metazoa_protein_subloc.txt"
  download.file(url, destfile, quiet = F)
}

df = read.csv(destfile, sep = "\t", stringsAsFactors = F, header = F)
names(df) = c("uniprot_id", "localization")

library(biomaRt)

if(GENOME == "mm9"){
  ensembl_dataset = useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  species = "mouse"
}

if(GENOME == "hg19"){
  ensembl_dataset = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  species = "human"
}

# # filters
# listFilters(ensembl_dataset)[grepl("symbol", listFilters(ensembl_dataset)[,1]),]
# listFilters(ensembl_dataset)[grepl("uniprot", listFilters(ensembl_dataset)[,1]),] 
# listFilters(ensembl_dataset)[grepl("wissprot", listFilters(ensembl_dataset)[,1]),] 
 
# # attributes
# listAttributes(ensembl_dataset)[,1][grepl("symbol", listAttributes(ensembl_dataset)[,1])] 
# listAttributes(ensembl_dataset)[,1][grepl("entrezgene", listAttributes(ensembl_dataset)[,1])] 
# listAttributes(ensembl_dataset)[grepl("description", listAttributes(ensembl_dataset)[,1]),] 

# Query biomart in minibatches (nrow(df) is 4080818, too big for biomaRt in one batch) 
result = data.frame()
for(start in seq(1, nrow(df), 10000)){ 
  print(start)
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

# Process records
require(dplyr)
df = merge(df, result, by.x = "uniprot_id", by.y = "uniprotswissprot") %>%
  dplyr::select(mgi_symbol, localization) %>%
  distinct() %>%
  group_by(mgi_symbol) %>%
  mutate(localizations = paste0(localization, collapse = ", ")) %>%
  dplyr::select(-localization) %>%
  ungroup() %>%
  mutate(ER = grepl("ER", localizations)) %>%
  mutate(Golgi = grepl("Golgi", localizations)) %>%
  mutate(Membrane = grepl("[Mm]embrane", localizations)) %>%
  mutate(Nucleus = grepl("Nucleus", localizations)) %>%
  mutate(Cytoplasm = grepl("Cytoplasm", localizations)) %>%
  mutate(Mitochondria = grepl("Mitochondria", localizations)) %>%
  mutate(Lysosome = grepl("Lysosome", localizations)) %>%
  mutate(Peroxisome = grepl("Peroxisome", localizations)) %>%
  mutate(Cytoskeleton = grepl("Cytoskeleton", localizations)) %>%
  mutate(Secreted = grepl("Secreted \\(curated\\)", localizations)) %>%
  mutate(Secreted_likely = grepl("Secreted \\(likely\\)", localizations)) %>%
  mutate(Secreted_highlylikely = grepl("Secreted \\(highly likely\\)", localizations)) %>%
  mutate(Secreted_weaklylikely = grepl("Secreted \\(weakly likely\\)", localizations)) %>%
  dplyr::select(-locations)

df[df == TRUE] = 1
df[df == FALSE] = 0

write.csv(df, file.path(SHARED_DATA_DIR, paste0(species, "_protein_locations.csv"), 
                        row.names=F) 
