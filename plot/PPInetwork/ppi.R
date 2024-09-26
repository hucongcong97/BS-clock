setwd('F:\\Methylation\\model\\train_data4\\plot_data\\PPInetwork\\')
# °²×° biomaRt °ü
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("biomaRt")
# BiocManager::install("STRINGdb")
# install.packages("visNetwork")

# load packages
library(biomaRt)
library(igraph)
library(STRINGdb)
library(visNetwork)
#-------------------------------------------------------------------------------
# Map gene name to Ensembl protein ID
#-------------------------------------------------------------------------------

# Connect to Ensembl database
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

tissue_list = c('Blood','Brain','Lung','Skin') # c('Blood','Brain','Lung','Skin')
for(tissue in tissue_list){
  print(tissue)
  data = read.csv(paste0('Manhattan_data_',tissue,'_new.csv'))
  gene_list <- unique(data$SYMBOL)
  
  # Get Ensembl protein ID
  df <- getBM(filters = "hgnc_symbol",
              attributes = c("hgnc_symbol", "ensembl_peptide_id"),
              values = gene_list,
              mart = mart)
  
  # Delete the row with empty string
  is_empty <- df == ""
  row_has_empty <- rowSums(is_empty) > 0
  df_clean <- df[!row_has_empty, ]
  
  # save result
  write.csv(df_clean, file = paste0("Manhattan_data_",tissue,"_covert.csv"), row.names = FALSE)
}


#-------------------------------------------------------------------------------
# Protein-protein interaction network
#-------------------------------------------------------------------------------
# read data
df1 = read.csv('Manhattan_data_Blood_new.csv')
df1 = df1[(abs(df1$Spearman.Correlation)>=0.45)&(df1$Spearman.P.value<1e-3),]
df2 = read.csv('Manhattan_data_Brain_new.csv')
df3 = read.csv('Manhattan_data_Lung_new.csv')
df3 = df1[df3$Spearman.P.value<1e-2,]
df4 = read.csv('Manhattan_data_Skin_new.csv')

# Merge gene names of all tissues
all_genes <- unique(c(df1$SYMBOL, df2$SYMBOL, df3$SYMBOL, df4$SYMBOL))

# Initialize STRINGdb object
string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=400, input_directory="")

# Map gene names to protein IDs of STRINGdb
mapped_genes <- string_db$map(data.frame(Gene=all_genes), "Gene")
mapped_genes = mapped_genes[complete.cases(mapped_genes),]

# Get protein interaction data
interactions <- string_db$get_interactions(mapped_genes$STRING_id)

# Remap back to gene names
mapped_interactions <- merge(interactions, mapped_genes[, c("STRING_id", "Gene")], by.x="from", by.y="STRING_id", all.x=TRUE)
mapped_interactions <- merge(mapped_interactions, mapped_genes[, c("STRING_id", "Gene")], by.x="to", by.y="STRING_id", all.x=TRUE)
colnames(mapped_interactions) <- c("to", "from","score","Gene_from","Gene_to")

# Add tissue source information to mapped_interactions data frame
mapped_interactions$source <- NA
mapped_interactions$source[mapped_interactions$Gene_from %in% df1$SYMBOL | mapped_interactions$Gene_to %in% df1$SYMBOL] <- "Blood"
mapped_interactions$source[mapped_interactions$Gene_from %in% df2$SYMBOL | mapped_interactions$Gene_to %in% df2$SYMBOL] <- "Brain"
mapped_interactions$source[mapped_interactions$Gene_from %in% df3$SYMBOL | mapped_interactions$Gene_to %in% df3$SYMBOL] <- "Lung"
mapped_interactions$source[mapped_interactions$Gene_from %in% df4$SYMBOL | mapped_interactions$Gene_to %in% df4$SYMBOL] <- "Skin"

# Only keep columns related to genes and tissue source information
final_interactions <- mapped_interactions[, c("Gene_from", "score","Gene_to", "source")]
colnames(final_interactions) <- c("protein1", "score","protein2", "source")

# Save the data frame as a CSV file
write.csv(final_interactions, file="interaction_all_temp.csv", row.names=FALSE)
