setwd('F:\\Methylation\\model\\train_data4\\plot_data\\PPInetwork\\')
# 安装 biomaRt 包
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("biomaRt")
# BiocManager::install("STRINGdb")
# install.packages("visNetwork")

# 加载R包
library(biomaRt)
library(igraph)
library(STRINGdb)
library(visNetwork)
#-------------------------------------------------------------------------------
# 将gene名映射 Ensembl 蛋白质 ID
#-------------------------------------------------------------------------------

# 连接到 Ensembl 数据库
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

tissue_list = c('Blood','Brain','Lung','Skin') # c('Blood','Brain','Lung','Skin')
for(tissue in tissue_list){
  print(tissue)
  data = read.csv(paste0('Manhattan_data_',tissue,'_new.csv'))
  gene_list <- unique(data$SYMBOL)
  
  # 获取 Ensembl 蛋白质 ID
  df <- getBM(filters = "hgnc_symbol",
              attributes = c("hgnc_symbol", "ensembl_peptide_id"),
              values = gene_list,
              mart = mart)
  
  # 删除空字符串所在的行
  is_empty <- df == ""
  row_has_empty <- rowSums(is_empty) > 0
  df_clean <- df[!row_has_empty, ]
  
  # 保存转换结果
  write.csv(df_clean, file = paste0("Manhattan_data_",tissue,"_covert.csv"), row.names = FALSE)
}


#-------------------------------------------------------------------------------
# 蛋白质-蛋白质互作网络
#-------------------------------------------------------------------------------
# 读取数据
df1 = read.csv('Manhattan_data_Blood_new.csv')
df1 = df1[(abs(df1$Spearman.Correlation)>=0.45)&(df1$Spearman.P.value<1e-3),]
df2 = read.csv('Manhattan_data_Brain_new.csv')
df3 = read.csv('Manhattan_data_Lung_new.csv')
df3 = df1[df3$Spearman.P.value<1e-2,]
df4 = read.csv('Manhattan_data_Skin_new.csv')

# 合并所有组织的基因名
all_genes <- unique(c(df1$SYMBOL, df2$SYMBOL, df3$SYMBOL, df4$SYMBOL))

# 初始化 STRINGdb 对象
string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=400, input_directory="")

# 映射基因名到 STRINGdb 的蛋白质 ID
mapped_genes <- string_db$map(data.frame(Gene=all_genes), "Gene")
mapped_genes = mapped_genes[complete.cases(mapped_genes),]

# 获取蛋白质互作数据
interactions <- string_db$get_interactions(mapped_genes$STRING_id)

# 重新映射回基因名
mapped_interactions <- merge(interactions, mapped_genes[, c("STRING_id", "Gene")], by.x="from", by.y="STRING_id", all.x=TRUE)
mapped_interactions <- merge(mapped_interactions, mapped_genes[, c("STRING_id", "Gene")], by.x="to", by.y="STRING_id", all.x=TRUE)
colnames(mapped_interactions) <- c("to", "from","score","Gene_from","Gene_to")

# 将组织来源信息添加到 mapped_interactions 数据框中
mapped_interactions$source <- NA
mapped_interactions$source[mapped_interactions$Gene_from %in% df1$SYMBOL | mapped_interactions$Gene_to %in% df1$SYMBOL] <- "Blood"
mapped_interactions$source[mapped_interactions$Gene_from %in% df2$SYMBOL | mapped_interactions$Gene_to %in% df2$SYMBOL] <- "Brain"
mapped_interactions$source[mapped_interactions$Gene_from %in% df3$SYMBOL | mapped_interactions$Gene_to %in% df3$SYMBOL] <- "Lung"
mapped_interactions$source[mapped_interactions$Gene_from %in% df4$SYMBOL | mapped_interactions$Gene_to %in% df4$SYMBOL] <- "Skin"

# 只保留与基因相关的列和组织来源信息
final_interactions <- mapped_interactions[, c("Gene_from", "score","Gene_to", "source")]
colnames(final_interactions) <- c("protein1", "score","protein2", "source")

# 将数据框保存为 CSV 文件
write.csv(final_interactions, file="interaction_all_temp.csv", row.names=FALSE)
