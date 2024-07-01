library(ComplexHeatmap)
library(circlize)
#500 CpG
# 读取normal_pos_data文件
normal_pos_data <- read.table("normal_positive_500.tsv", header = F, sep = "\t")
normal_pos_data <- normal_pos_data[,c(1,2,3,4,5,20)]
colnames(normal_pos_data) <- c("Ontology","ID","Term Name","Rank","P-Value","Total Genes")

# 读取disease_pos_data文件
disease_pos_data <- read.table("disease_positive_500.tsv", header = F, sep = "\t")
disease_pos_data <- disease_pos_data[,c(1,2,3,4,5,20)]
colnames(disease_pos_data) <- c("Ontology","ID","Term Name","Rank","P-Value","Total Genes")

# 读取normal_neg_data文件
normal_neg_data <- read.table("normal_negative_500.tsv", header = F, sep = "\t")
normal_neg_data <- normal_neg_data[,c(1,2,3,4,5,20)]
colnames(normal_neg_data) <- c("Ontology","ID","Term Name","Rank","P-Value","Total Genes")

# 读取disease_neg_data文件
disease_neg_data <- read.table("disease_negative_500.tsv", header = F, sep = "\t")
disease_neg_data <- disease_neg_data[,c(1,2,3,4,5,20)]
colnames(disease_neg_data) <- c("Ontology","ID","Term Name","Rank","P-Value","Total Genes")



########筛选GO Biological Process
normal_pos_GOBP_data <- normal_pos_data[normal_pos_data[,1] == "GO Biological Process", ]

disease_pos_GOBP_data <- disease_pos_data[disease_pos_data[,1] == "GO Biological Process", ]

normal_neg_GOBP_data <- normal_neg_data[normal_neg_data[,1] == "GO Biological Process",]

disease_neg_GOBP_data <- disease_neg_data[disease_neg_data[,1] == "GO Biological Process", ]


#四个数据框中第二列相同的数据
common_GOBP_data <- Reduce(intersect, list(normal_pos_GOBP_data[, 2], disease_pos_GOBP_data[, 2], normal_neg_GOBP_data[, 2]))

normal_pos_GOBP_data <- normal_pos_data[normal_pos_data[,2] %in% c("GO:0051572","GO:0000415"),]
normal_pos_GOBP_data <- normal_pos_GOBP_data[,c(1,3,5,6)]

disease_pos_GOBP_data <- disease_pos_data[disease_pos_data[,2] %in% c("GO:0051572","GO:0000415"), ]
#disease_pos_GOBP_data <- disease_pos_GOBP_data[c(2,1),]
disease_pos_GOBP_data <- disease_pos_GOBP_data[,c(5,6)]

normal_neg_GOBP_data <- normal_neg_data[normal_neg_data[,2] %in% c("GO:0051572","GO:0000415"),]
normal_neg_GOBP_data <- normal_neg_GOBP_data[,c(5,6)]

#disease_neg_GOBP_data <- disease_neg_data[disease_neg_data[,2] %in% c("GO:0032940","GO:0018394"), ]

#整理到一�?
GOBP_data <- cbind(normal_pos_GOBP_data,disease_pos_GOBP_data,normal_neg_GOBP_data)
GOBP_data <- GOBP_data[,c(1,2,3,4,5,7)]
colnames(GOBP_data) <- c("ontology","term name","normal_pos_p_value","total genes","disease_pos_p_value","normal_neg_p_value")
GOBP_data$disease_neg_p_value <- c(1,1)
GOBP_data <- GOBP_data[,c(1,2,3,5,6,7,4)]





########筛选GO Cellular Component
normal_pos_GOCC_data <- normal_pos_data[normal_pos_data[,1] == "GO Cellular Component", ]

disease_pos_GOCC_data <- disease_pos_data[disease_pos_data[,1] == "GO Cellular Component", ]

normal_neg_GOCC_data <- normal_neg_data[normal_neg_data[,1] == "GO Cellular Component",]

disease_neg_GOCC_data <- disease_neg_data[disease_neg_data[,1] == "GO Cellular Component", ]

#四个数据框中第二列相同的数据
common_GOCC_data <- Reduce(intersect, list(disease_pos_GOCC_data[, 2], normal_neg_GOCC_data[, 2],disease_neg_GOCC_data[, 2]))


#筛选交集数�?              data[data[, 1] %in% c(12, 13), ]

#normal_pos_GOCC_data <- normal_pos_data[normal_pos_data[,2] %in% c("GO:0044437","GO:0005765"),]
#normal_pos_GOCC_data <- normal_pos_GOCC_data[,c(1,3,5,6)]

disease_pos_GOCC_data <- disease_pos_data[disease_pos_data[,2] %in% common_GOCC_data, ]

disease_pos_GOCC_data <- disease_pos_data[disease_pos_data[,2] %in% c("GO:0030667","GO:0101003"), ]
#disease_pos_GOCC_data <- disease_pos_GOCC_data[c(2,1),]
disease_pos_GOCC_data <- disease_pos_GOCC_data[,c(1,3,5,6)]

normal_neg_GOCC_data <- normal_neg_data[normal_neg_data[,2] %in% c("GO:0030667","GO:0101003"),]
normal_neg_GOCC_data <- normal_neg_GOCC_data[c(2,1),]
normal_neg_GOCC_data <- normal_neg_GOCC_data[,c(5,6)]

disease_neg_GOCC_data <- disease_neg_data[disease_neg_data[,2] %in% c("GO:0030667","GO:0101003"), ]
disease_neg_GOCC_data <- disease_neg_GOCC_data[c(2,1),]
disease_neg_GOCC_data <- disease_neg_GOCC_data[,c(5,6)]

#整理到一�?
GOCC_data <- cbind(disease_pos_GOCC_data,normal_neg_GOCC_data,disease_neg_GOCC_data)
GOCC_data <- GOCC_data[,c(1,2,3,4,5,7)]
colnames(GOCC_data) <- c("ontology","term name","disease_pos_p_value","total genes","normal_neg_p_value","disease_neg_p_value")
GOCC_data$normal_pos_p_value <- c(1,1)

GOCC_data <- GOCC_data[,c(1,2,7,3,5,6,4)]




########筛选GO Molecular Function
normal_pos_GOMF_data <- normal_pos_data[normal_pos_data[,1] == "GO Molecular Function", ]

disease_pos_GOMF_data <- disease_pos_data[disease_pos_data[,1] == "GO Molecular Function", ]

normal_neg_GOMF_data <- normal_neg_data[normal_neg_data[,1] == "GO Molecular Function",]

disease_neg_GOMF_data <- disease_neg_data[disease_neg_data[,1] == "GO Molecular Function", ]

#四个数据框中第二列相同的数据
common_GOMF_data <- Reduce(intersect, list(normal_pos_GOMF_data[, 2], disease_pos_GOMF_data[, 2], normal_neg_GOMF_data[, 2],disease_neg_GOMF_data[, 2]))

#筛选交集数�?              data[data[, 1] %in% c(12, 13), ]
normal_pos_GOMF_data <- normal_pos_data[normal_pos_data[,2] %in% c("GO:0019900"),]

#normal_pos_GOMF_data <- normal_pos_data[normal_pos_data[,2] %in% c("GO:0003690","GO:1990837"),]
normal_pos_GOMF_data <- normal_pos_GOMF_data[,c(1,3,5,6)]

disease_pos_GOMF_data <- disease_pos_data[disease_pos_data[,2] %in% c("GO:0019900"), ]
disease_pos_GOMF_data <- disease_pos_GOMF_data[,c(5,6)]

normal_neg_GOMF_data <- normal_neg_data[normal_neg_data[,2] %in% c("GO:0019900"),]
#normal_neg_GOMF_data <- normal_neg_GOMF_data[c(2,1),]
normal_neg_GOMF_data <- normal_neg_GOMF_data[,c(5,6)]

disease_neg_GOMF_data <- disease_neg_data[disease_neg_data[,2] %in% c("GO:0019900"), ]
#disease_neg_GOMF_data <- disease_neg_GOMF_data[c(2,1),]
disease_neg_GOMF_data <- disease_neg_GOMF_data[,c(5,6)]

#整理到一�?
GOMF_data <- cbind(normal_pos_GOMF_data,disease_pos_GOMF_data,normal_neg_GOMF_data,disease_neg_GOMF_data)
GOMF_data <- GOMF_data[,c(1,2,3,4,5,7,9)]
colnames(GOMF_data) <- c("ontology","term name","normal_pos_p_value","total genes","disease_pos_p_value","normal_neg_p_value","disease_neg_p_value")
GOMF_data <- GOMF_data[,c(1,2,3,5,6,7,4)]



########筛选Human Phenotype
normal_pos_HP_data <- normal_pos_data[normal_pos_data[,1] == "Human Phenotype", ]

disease_pos_HP_data <- disease_pos_data[disease_pos_data[,1] == "Human Phenotype", ]

normal_neg_HP_data <- normal_neg_data[normal_neg_data[,1] == "Human Phenotype",]

disease_neg_HP_data <- disease_neg_data[disease_neg_data[,1] == "Human Phenotype", ]

#四个数据框中第二列相同的数据
common_HP_data <- Reduce(intersect, list(normal_pos_HP_data[, 2], normal_neg_HP_data[, 2], disease_neg_HP_data[, 2]))


#筛选交集数�?              data[data[, 1] %in% c(12, 13), ]
normal_pos_HP_data <- normal_pos_data[normal_pos_data[,2] %in% common_HP_data,]

normal_pos_HP_data <- normal_pos_data[normal_pos_data[,2] %in% c("HP:0005952","HP:0008404"),]
normal_pos_HP_data <- normal_pos_HP_data[,c(1,3,5,6)]

#disease_pos_HP_data <- disease_pos_data[disease_pos_data[,2] %in% c("HP:0005952","HP:0008404"), ]
#disease_pos_HP_data <- disease_pos_HP_data[c(2,1),]
#disease_pos_HP_data <- disease_pos_HP_data[,c(5,6)]

normal_neg_HP_data <- normal_neg_data[normal_neg_data[,2] %in% c("HP:0005952","HP:0008404"),]
#normal_neg_HP_data <- normal_neg_HP_data[c(2,1),]
normal_neg_HP_data <- normal_neg_HP_data[,c(5,6)]

disease_neg_HP_data <- disease_neg_data[disease_neg_data[,2] %in% c("HP:0005952","HP:0008404"), ]
#disease_neg_HP_data <- disease_neg_HP_data[c(2,1),]
disease_neg_HP_data <- disease_neg_HP_data[,c(5,6)]

#整理到一�?
HP_data <- cbind(normal_pos_HP_data,normal_neg_HP_data,disease_neg_HP_data)
HP_data <- HP_data[,c(1,2,3,4,5,7)]
colnames(HP_data) <- c("ontology","term name","normal_pos_p_value","total genes","normal_neg_p_value","disease_neg_p_value")
HP_data$disease_pos_p_value <- c(1,1)
HP_data <- HP_data[,c(1,2,3,7,5,6,4)]



#将上面的数据整理到一�?
heatmap_data <- rbind(GOBP_data,GOCC_data,GOMF_data,HP_data)
total_genes <- heatmap_data$`total genes`
heatmap_data <- heatmap_data[,-1]
rownames(heatmap_data) <- heatmap_data[,1]
heatmap_data <- heatmap_data[,-1]
heatmap_data <- heatmap_data[,-5]



heatmap_data_1 <- -log10(heatmap_data)
heatmap_data_1$type <- c(rep("GO Biological Process",2),rep("GO Cellular Component",2),rep("GO Molecular Function",1),rep("Human Phenotype",2))
heatmap_data_1$type <- factor(heatmap_data_1$type)
# 假设你的数据框名为data，需要将�?4列转换为数值型
heatmap_data_1[, 1:4] <- lapply(heatmap_data_1[, 1:4], as.numeric)
heatmap_data_2 <- as.matrix(heatmap_data_1[,1:4])


#首先，自定义颜色�?
col_fun = colorRamp2(c(0,10), c("white", "red"))

#继续尝试绘制热图，这里取消列聚类，添加cell边框�?
Heatmap(heatmap_data_2[,1:4],col = col_fun,
        name = "-log10(p_value)",
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        rect_gp = gpar(col= "white"),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        row_names_side = "left")


# �?5个柱子一个颜色，后面4个一个颜色，后面2个一个颜色，最�?3个一个颜�?
colors <- c(rep("#1f77b4", 2), rep("#ff7f0e", 2), rep("#9467bd", 1), rep("#2ca02c", 2))
#为type设置分组颜色                                                
colors_type <- c("#1f77b4", "#ff7f0e", "#9467bd", "#2ca02c")
names(colors_type) <- unique(heatmap_data_1$type)



bar2 = rowAnnotation(
  sum2 = anno_barplot(
    total_genes,
    bar_width = 0.9,
    #gp = gpar(col = "white", fill = "orange"),
    gp = gpar(col = "white", fill = colors),
    border = T,
    axis_param = list(at = c(0,324,647),
                      labels = c("0","324","647")),
    width = unit(1, "cm")), show_annotation_name = F)





###绘制热图
Heatmap(heatmap_data_2[,1:4],col = col_fun,
        width = unit(8, "cm"),
        height = unit(8, "cm"),
        name = "-log10(p_value)",
        cluster_columns = F,
        cluster_rows = F,
        show_row_dend = F,
        show_column_dend = F,
        rect_gp = gpar(col= "white",lwd = 1),
        #top_annotation = bar1,
        left_annotation =bar2,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        column_title = "Increased with age    Decreased with age",
        column_title_gp = gpar(fontsize = 12),
        cell_fun = function(j, i, x, y, width, height, fill, data=10^(-heatmap_data_2[,1:4])){
          value <- data[i, j]  # 使用行列名称来访问数�?
          formatted_value <- format(value, scientific = TRUE)
          # 将科学计数法字符串解析为数字
          scientific_notation_value <- as.numeric(formatted_value)
          # 四舍五入整数部分
          integer_part <- round(scientific_notation_value)
          # 提取小数部分
          decimal_part <- scientific_notation_value - integer_part
          # 将整数部分和小数部分重新格式化为字符�?
          formatted_value <- paste0(format(decimal_part, digits = 2))
          grid.text(formatted_value, x, y, gp = gpar(fontsize = 10))
        })+
  Heatmap(heatmap_data_1$type,name = "Ontology",width = unit(3, "mm"),col = colors_type)





#1000 CpG
# 读取normal_pos_data文件
normal_pos_data <- read.table("data/GREAT_heatmap/normal_positive_1000.tsv", header = F, sep = "\t")
normal_pos_data <- normal_pos_data[,c(1,2,3,4,5,20)]
colnames(normal_pos_data) <- c("Ontology","ID","Term Name","Rank","P-Value","Total Genes")

# 读取disease_pos_data文件
disease_pos_data <- read.table("data/GREAT_heatmap/disease_positive_1000.tsv", header = F, sep = "\t")
disease_pos_data <- disease_pos_data[,c(1,2,3,4,5,20)]
colnames(disease_pos_data) <- c("Ontology","ID","Term Name","Rank","P-Value","Total Genes")

# 读取normal_neg_data文件
normal_neg_data <- read.table("data/GREAT_heatmap/normal_negative_1000.tsv", header = F, sep = "\t")
normal_neg_data <- normal_neg_data[,c(1,2,3,4,5,20)]
colnames(normal_neg_data) <- c("Ontology","ID","Term Name","Rank","P-Value","Total Genes")

# 读取disease_neg_data文件
disease_neg_data <- read.table("data/GREAT_heatmap/disease_negative_1000.tsv", header = F, sep = "\t")
disease_neg_data <- disease_neg_data[,c(1,2,3,4,5,20)]
colnames(disease_neg_data) <- c("Ontology","ID","Term Name","Rank","P-Value","Total Genes")



########筛选GO Biological Process
normal_pos_GOBP_data <- normal_pos_data[normal_pos_data[,1] == "GO Biological Process", ]

disease_pos_GOBP_data <- disease_pos_data[disease_pos_data[,1] == "GO Biological Process", ]

normal_neg_GOBP_data <- normal_neg_data[normal_neg_data[,1] == "GO Biological Process",]

disease_neg_GOBP_data <- disease_neg_data[disease_neg_data[,1] == "GO Biological Process", ]


#四个数据框中第二列相同的数据
common_GOBP_data <- Reduce(intersect, list(normal_pos_GOBP_data[, 2], disease_pos_GOBP_data[, 2],normal_neg_GOBP_data[,2]))

#筛选交集数�? data[data[, 1] %in% c(12, 13), ]
normal_pos_GOBP_data <- normal_pos_data[normal_pos_data[,2] %in% common_GOBP_data,]

normal_pos_GOBP_data <- normal_pos_data[normal_pos_data[,2] %in% c("GO:0051572","GO:0000415"),]
normal_pos_GOBP_data <- normal_pos_GOBP_data[,c(1,3,5,6)]

disease_pos_GOBP_data <- disease_pos_data[disease_pos_data[,2] %in% c("GO:0051572","GO:0000415"), ]
#disease_pos_GOBP_data <- disease_pos_GOBP_data[c(2,1),]
disease_pos_GOBP_data <- disease_pos_GOBP_data[,c(5,6)]

normal_neg_GOBP_data <- normal_neg_data[normal_neg_data[,2] %in% c("GO:0051572","GO:0000415"),]
normal_neg_GOBP_data <- normal_neg_GOBP_data[,c(5,6)]

#disease_neg_GOBP_data <- disease_neg_data[disease_neg_data[,2] %in% c("GO:0032940","GO:0018394"), ]

#整理到一�?
GOBP_data <- cbind(normal_pos_GOBP_data,disease_pos_GOBP_data,normal_neg_GOBP_data)
GOBP_data <- GOBP_data[,c(1,2,3,4,5,7)]
colnames(GOBP_data) <- c("ontology","term name","normal_pos_p_value","total genes","disease_pos_p_value","normal_neg_p_value")
#GOBP_data$normal_neg_p_value <- c(1,1)
GOBP_data$disease_neg_p_value <- c(1,1)
GOBP_data <- GOBP_data[,c(1,2,3,5,6,7,4)]



########筛选GO Cellular Component
normal_pos_GOCC_data <- normal_pos_data[normal_pos_data[,1] == "GO Cellular Component", ]

disease_pos_GOCC_data <- disease_pos_data[disease_pos_data[,1] == "GO Cellular Component", ]

normal_neg_GOCC_data <- normal_neg_data[normal_neg_data[,1] == "GO Cellular Component",]

disease_neg_GOCC_data <- disease_neg_data[disease_neg_data[,1] == "GO Cellular Component", ]

#四个数据框中第二列相同的数据
common_GOCC_data <- Reduce(intersect, list(normal_pos_GOCC_data[, 2],normal_neg_GOCC_data[,2],disease_neg_GOCC_data[, 2]))


#筛选交集数�?              data[data[, 1] %in% c(12, 13), ]
normal_pos_GOCC_data <- normal_pos_data[normal_pos_data[,2] %in% common_GOCC_data,]

normal_pos_GOCC_data <- normal_pos_data[normal_pos_data[,2] %in% c("GO:0071818","GO:0031519"),]
normal_pos_GOCC_data <- normal_pos_GOCC_data[,c(1,3,5,6)]

#disease_pos_GOCC_data <- disease_pos_data[disease_pos_data[,2] %in% c("GO:0044437","GO:0005765"), ]
#disease_pos_GOCC_data <- disease_pos_GOCC_data[c(2,1),]
#disease_pos_GOCC_data <- disease_pos_GOCC_data[,c(5,6)]

normal_neg_GOCC_data <- normal_neg_data[normal_neg_data[,2] %in% c("GO:0071818","GO:0031519"),]
#normal_neg_GOCC_data <- normal_neg_GOCC_data[c(2,1),]
normal_neg_GOCC_data <- normal_neg_GOCC_data[,c(5,6)]

disease_neg_GOCC_data <- disease_neg_data[disease_neg_data[,2] %in% c("GO:0071818","GO:0031519"), ]
#disease_neg_GOCC_data <- disease_neg_GOCC_data[c(2,1),]
disease_neg_GOCC_data <- disease_neg_GOCC_data[,c(5,6)]

#整理到一�?
GOCC_data <- cbind(normal_pos_GOCC_data,normal_neg_GOCC_data,disease_neg_GOCC_data)
GOCC_data <- GOCC_data[,c(1,2,3,4,5,7)]
colnames(GOCC_data) <- c("ontology","term name","normal_pos_p_value","total genes","normal_neg_p_value","disease_neg_p_value")
GOCC_data$disease_pos_p_value <- c(1,1)
#GOCC_data$normal_neg_p_value <- c(1,1)
GOCC_data <- GOCC_data[,c(1,2,3,7,5,6,4)]








########筛选GO Molecular Function
normal_pos_GOMF_data <- normal_pos_data[normal_pos_data[,1] == "GO Molecular Function", ]

disease_pos_GOMF_data <- disease_pos_data[disease_pos_data[,1] == "GO Molecular Function", ]

normal_neg_GOMF_data <- normal_neg_data[normal_neg_data[,1] == "GO Molecular Function",]

disease_neg_GOMF_data <- disease_neg_data[disease_neg_data[,1] == "GO Molecular Function", ]

#四个数据框中第二列相同的数据
common_GOMF_data <- Reduce(intersect, list(normal_neg_GOMF_data[, 2], disease_neg_GOMF_data[, 2]))

#筛选交集数�?              data[data[, 1] %in% c(12, 13), ]
#normal_pos_GOMF_data <- normal_pos_data[normal_pos_data[,2] %in% common_GOMF_data,]

#normal_pos_GOMF_data <- normal_pos_data[normal_pos_data[,2] %in% c("GO:0008134","GO:0047192"),]
#normal_pos_GOMF_data <- normal_pos_GOMF_data[,c(1,3,5,6)]

#disease_pos_GOMF_data <- disease_pos_data[disease_pos_data[,2] %in% c("GO:0008134","GO:0047192"), ]
#disease_pos_GOMF_data <- disease_pos_GOMF_data[c(2,1),]
#disease_pos_GOMF_data <- disease_pos_GOMF_data[,c(5,6)]


normal_neg_GOMF_data <- normal_neg_data[normal_neg_data[,2] %in% common_GOMF_data,]

normal_neg_GOMF_data <- normal_neg_data[normal_neg_data[,2] %in% c("GO:0042813","GO:0005484"),]
#normal_neg_GOMF_data <- normal_neg_GOMF_data[c(2,1),]
normal_neg_GOMF_data <- normal_neg_GOMF_data[,c(1,3,5,6)]

disease_neg_GOMF_data <- disease_neg_data[disease_neg_data[,2] %in% c("GO:0042813","GO:0005484"), ]
#disease_neg_GOMF_data <- disease_neg_GOMF_data[c(2,1),]
disease_neg_GOMF_data <- disease_neg_GOMF_data[,c(5,6)]

#整理到一�?
GOMF_data <- cbind(normal_neg_GOMF_data,disease_neg_GOMF_data)
GOMF_data <- GOMF_data[,c(1,2,3,4,5)]
colnames(GOMF_data) <- c("ontology","term name","normal_neg_p_value","total genes","disease_neg_p_value")
GOMF_data$normal_pos_p_value <- c(1,1)
GOMF_data$disease_pos_p_value <- c(1,1)
GOMF_data <- GOMF_data[,c(1,2,6,7,3,5,4)]








########筛选Human Phenotype
normal_pos_HP_data <- normal_pos_data[normal_pos_data[,1] == "Human Phenotype", ]

disease_pos_HP_data <- disease_pos_data[disease_pos_data[,1] == "Human Phenotype", ]

normal_neg_HP_data <- normal_neg_data[normal_neg_data[,1] == "Human Phenotype",]

disease_neg_HP_data <- disease_neg_data[disease_neg_data[,1] == "Human Phenotype", ]

#四个数据框中第二列相同的数据
common_HP_data <- Reduce(intersect, list(disease_pos_HP_data[, 2], normal_neg_HP_data[, 2]))


#筛选交集数�?              data[data[, 1] %in% c(12, 13), ]
#normal_pos_HP_data <- normal_pos_data[normal_pos_data[,2] %in% c("HP:0002205","HP:0001633"),]
#normal_pos_HP_data <- normal_pos_HP_data[,c(1,3,5,6)]


disease_pos_HP_data <- disease_pos_data[disease_pos_data[,2] %in% common_HP_data, ]

disease_pos_HP_data <- disease_pos_data[disease_pos_data[,2] %in% c("HP:0100614","HP:0002715"), ]
#disease_pos_HP_data <- disease_pos_HP_data[c(2,1),]
disease_pos_HP_data <- disease_pos_HP_data[,c(1,3,5,6)]

normal_neg_HP_data <- normal_neg_data[normal_neg_data[,2] %in% c("HP:0100614","HP:0002715"),]
#normal_neg_HP_data <- normal_neg_HP_data[c(2,1),]
normal_neg_HP_data <- normal_neg_HP_data[,c(5,6)]



#整理到一�?
HP_data <- cbind(disease_pos_HP_data,normal_neg_HP_data)
HP_data <- HP_data[,c(1,2,3,4,5)]
colnames(HP_data) <- c("ontology","term name","disease_pos_p_value","total genes","normal_neg_p_value")
HP_data$normal_pos_p_value <- c(1,1)
HP_data$disease_neg_p_value <- c(1,1)
HP_data <- HP_data[,c(1,2,6,3,5,7,4)]


#将上面的数据整理到一�?
heatmap_data <- rbind(GOBP_data,GOCC_data,GOMF_data,HP_data)
total_genes <- heatmap_data$`total genes`
heatmap_data <- heatmap_data[,-1]
rownames(heatmap_data) <- heatmap_data[,1]
heatmap_data <- heatmap_data[,-1]
heatmap_data <- heatmap_data[,-5]

heatmap_data_1 <- -log10(heatmap_data)
heatmap_data_1$type <- c(rep("GO Biological Process",2),rep("GO Cellular Component",2),rep("GO Molecular Function",2),rep("Human Phenotype",2))
heatmap_data_1$type <- factor(heatmap_data_1$type)
# 假设你的数据框名为data，需要将�?4列转换为数值型
heatmap_data_1[, 1:4] <- lapply(heatmap_data_1[, 1:4], as.numeric)
#heatmap_data_1$type <- as.numeric(heatmap_data_1$type)
heatmap_data_2 <- as.matrix(heatmap_data_1[,1:4])


#首先，自定义颜色�?
col_fun = colorRamp2(c(0,15), c("white", "red"))

#继续尝试绘制热图，这里取消列聚类，添加cell边框�?
Heatmap(heatmap_data_2[,1:4],col = col_fun,
        name = "-log10(p_value)",
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        rect_gp = gpar(col= "white"),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        row_names_side = "left")


# �?5个柱子一个颜色，后面4个一个颜色，后面2个一个颜色，最�?3个一个颜�?
colors <- c(rep("#1f77b4", 2), rep("#ff7f0e", 2), rep("#9467bd", 2), rep("#2ca02c", 2))
#为type设置分组颜色                                                
colors_type <- c("#1f77b4", "#ff7f0e", "#9467bd", "#2ca02c")
names(colors_type) <- unique(heatmap_data_1$type)



bar2 = rowAnnotation(
  sum2 = anno_barplot(
    total_genes,
    bar_width = 0.9,
    #gp = gpar(col = "white", fill = "orange"),
    gp = gpar(col = "white", fill = colors),
    border = T,
    axis_param = list(at = c(0,392,786),
                      labels = c("0","393","786")),
    width = unit(1, "cm")), show_annotation_name = F)





###绘制热图
Heatmap(heatmap_data_2[,1:4],col = col_fun,
        width = unit(8, "cm"),
        height = unit(8, "cm"),
        name = "-log10(p_value)",
        cluster_columns = F,
        cluster_rows = F,
        show_row_dend = F,
        show_column_dend = F,
        rect_gp = gpar(col= "white",lwd = 1),
        #top_annotation = bar1,
        left_annotation =bar2,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        column_title = "Increased with age    Decreased with age",
        column_title_gp = gpar(fontsize = 12),
        cell_fun = function(j, i, x, y, width, height, fill, data=10^(-heatmap_data_2[,1:4])){
          value <- data[i, j]  # 使用行列名称来访问数�?
          formatted_value <- format(value, scientific = TRUE)
          # 将科学计数法字符串解析为数字
          scientific_notation_value <- as.numeric(formatted_value)
          # 四舍五入整数部分
          integer_part <- round(scientific_notation_value)
          # 提取小数部分
          decimal_part <- scientific_notation_value - integer_part
          # 将整数部分和小数部分重新格式化为字符�?
          formatted_value <- paste0(format(decimal_part, digits = 2))
          grid.text(formatted_value, x, y, gp = gpar(fontsize = 10))
        })+
  Heatmap(heatmap_data_1$type,name = "Ontology",width = unit(3, "mm"),col = colors_type)



#2000 CpG
# 读取normal_pos_data文件
normal_pos_data <- read.table("data/GREAT_heatmap/normal_positive_2000.tsv", header = F, sep = "\t")
normal_pos_data <- normal_pos_data[,c(1,2,3,4,5,20)]
colnames(normal_pos_data) <- c("Ontology","ID","Term Name","Rank","P-Value","Total Genes")

# 读取disease_pos_data文件
disease_pos_data <- read.table("data/GREAT_heatmap/disease_positive_2000.tsv", header = F, sep = "\t")
disease_pos_data <- disease_pos_data[,c(1,2,3,4,5,20)]
colnames(disease_pos_data) <- c("Ontology","ID","Term Name","Rank","P-Value","Total Genes")

# 读取normal_neg_data文件
normal_neg_data <- read.table("data/GREAT_heatmap/normal_negative_2000.tsv", header = F, sep = "\t")
normal_neg_data <- normal_neg_data[,c(1,2,3,4,5,20)]
colnames(normal_neg_data) <- c("Ontology","ID","Term Name","Rank","P-Value","Total Genes")

# 读取disease_neg_data文件
disease_neg_data <- read.table("data/GREAT_heatmap/disease_negative_2000.tsv", header = F, sep = "\t")
disease_neg_data <- disease_neg_data[,c(1,2,3,4,5,20)]
colnames(disease_neg_data) <- c("Ontology","ID","Term Name","Rank","P-Value","Total Genes")



########筛选GO Biological Process
normal_pos_GOBP_data <- normal_pos_data[normal_pos_data[,1] == "GO Biological Process", ]

disease_pos_GOBP_data <- disease_pos_data[disease_pos_data[,1] == "GO Biological Process", ]

normal_neg_GOBP_data <- normal_neg_data[normal_neg_data[,1] == "GO Biological Process",]

disease_neg_GOBP_data <- disease_neg_data[disease_neg_data[,1] == "GO Biological Process", ]


#四个数据框中第二列相同的数据
common_GOBP_data <- Reduce(intersect, list(normal_pos_GOBP_data[, 2], disease_pos_GOBP_data[, 2],normal_neg_GOBP_data[,2],disease_neg_GOBP_data[,2]))

#筛选交集数�? data[data[, 1] %in% c(12, 13), ]
normal_pos_GOBP_data <- normal_pos_data[normal_pos_data[,2] %in% common_GOBP_data,]

normal_pos_GOBP_data <- normal_pos_data[normal_pos_data[,2] %in% c("GO:0038202","GO:0032012"),]
normal_pos_GOBP_data <- normal_pos_GOBP_data[,c(1,3,5,6)]

disease_pos_GOBP_data <- disease_pos_data[disease_pos_data[,2] %in% c("GO:0038202","GO:0032012"), ]
#disease_pos_GOBP_data <- disease_pos_GOBP_data[c(2,1),]
disease_pos_GOBP_data <- disease_pos_GOBP_data[,c(5,6)]

normal_neg_GOBP_data <- normal_neg_data[normal_neg_data[,2] %in% c("GO:0038202","GO:0032012"),]
normal_neg_GOBP_data <- normal_neg_GOBP_data[,c(5,6)]

disease_neg_GOBP_data <- disease_neg_data[disease_neg_data[,2] %in% c("GO:0038202","GO:0032012"), ]
disease_neg_GOBP_data <- disease_neg_GOBP_data[c(2,1),]
disease_neg_GOBP_data <- disease_neg_GOBP_data[,c(5,6)]

#整理到一�?
GOBP_data <- cbind(normal_pos_GOBP_data,disease_pos_GOBP_data,normal_neg_GOBP_data,disease_neg_GOBP_data)
GOBP_data <- GOBP_data[,c(1,2,3,4,5,7,9)]
colnames(GOBP_data) <- c("ontology","term name","normal_pos_p_value","total genes","disease_pos_p_value","normal_neg_p_value","disease_neg_p_value")
#GOBP_data$normal_neg_p_value <- c(1)
#GOBP_data$disease_neg_p_value <- c(1,1)
GOBP_data <- GOBP_data[,c(1,2,3,5,6,7,4)]





########筛选GO Cellular Component
normal_pos_GOCC_data <- normal_pos_data[normal_pos_data[,1] == "GO Cellular Component", ]

disease_pos_GOCC_data <- disease_pos_data[disease_pos_data[,1] == "GO Cellular Component", ]

normal_neg_GOCC_data <- normal_neg_data[normal_neg_data[,1] == "GO Cellular Component",]

disease_neg_GOCC_data <- disease_neg_data[disease_neg_data[,1] == "GO Cellular Component", ]

#四个数据框中第二列相同的数据
common_GOCC_data <- Reduce(intersect, list(normal_pos_GOCC_data[, 2],disease_neg_GOCC_data[, 2]))


#筛选交集数�?              data[data[, 1] %in% c(12, 13), ]
normal_pos_GOCC_data <- normal_pos_data[normal_pos_data[,2] %in% common_GOCC_data,]

normal_pos_GOCC_data <- normal_pos_data[normal_pos_data[,2] %in% c("GO:0033291","GO:0031234"),]
normal_pos_GOCC_data <- normal_pos_GOCC_data[,c(1,3,5,6)]

#disease_pos_GOCC_data <- disease_pos_data[disease_pos_data[,2] %in% common_GOCC_data,]

#disease_pos_GOCC_data <- disease_pos_data[disease_pos_data[,2] %in% c("GO:0032040","GO:0005652","GO:0005793"), ]
#disease_pos_GOCC_data <- disease_pos_GOCC_data[c(2,1),]
#disease_pos_GOCC_data <- disease_pos_GOCC_data[,c(1,3,5,6)]

#normal_neg_GOCC_data <- normal_neg_data[normal_neg_data[,2] %in% c("GO:0032040","GO:0005652","GO:0005793"),]
#normal_neg_GOCC_data <- normal_neg_GOCC_data[c(2,3,1),]
#normal_neg_GOCC_data <- normal_neg_GOCC_data[,c(5,6)]

disease_neg_GOCC_data <- disease_neg_data[disease_neg_data[,2] %in% c("GO:0033291","GO:0031234","GO:0005793"), ]
#disease_neg_GOCC_data <- disease_neg_GOCC_data[c(3,2,1),]
disease_neg_GOCC_data <- disease_neg_GOCC_data[,c(5,6)]

#整理到一�?
GOCC_data <- cbind(normal_pos_GOCC_data,disease_neg_GOCC_data)
GOCC_data <- GOCC_data[,c(1,2,3,4,5)]
colnames(GOCC_data) <- c("ontology","term name","normal_pos_p_value","total genes","disease_neg_p_value")
GOCC_data$disease_pos_p_value <- c(1,1)
GOCC_data$normal_neg_p_value <- c(1,1)
GOCC_data <- GOCC_data[,c(1,2,3,6,7,5,4)]








########筛选GO Molecular Function
normal_pos_GOMF_data <- normal_pos_data[normal_pos_data[,1] == "GO Molecular Function", ]

disease_pos_GOMF_data <- disease_pos_data[disease_pos_data[,1] == "GO Molecular Function", ]

normal_neg_GOMF_data <- normal_neg_data[normal_neg_data[,1] == "GO Molecular Function",]

disease_neg_GOMF_data <- disease_neg_data[disease_neg_data[,1] == "GO Molecular Function", ]

#四个数据框中第二列相同的数据
common_GOMF_data <- Reduce(intersect, list(normal_pos_GOMF_data[, 2],normal_neg_GOMF_data[, 2],disease_neg_GOMF_data[, 2]))

#筛选交集数�?              data[data[, 1] %in% c(12, 13), ]
#normal_pos_GOMF_data <- normal_pos_GOMF_data[normal_pos_GOMF_data[,2] %in% common_GOMF_data,]

normal_pos_GOMF_data <- normal_pos_data[normal_pos_data[,2] %in% common_GOMF_data,]

normal_pos_GOMF_data <- normal_pos_data[normal_pos_data[,2] %in% c("GO:0000406","GO:0032181"),]
normal_pos_GOMF_data <- normal_pos_GOMF_data[,c(1,3,5,6)]

#disease_pos_GOMF_data <- disease_pos_data[disease_pos_data[,2] %in% c("GO:0047184","GO:0008170"), ]
#disease_pos_GOMF_data <- disease_pos_GOMF_data[c(2,1),]
#disease_pos_GOMF_data <- disease_pos_GOMF_data[,c(5,6)]

normal_neg_GOMF_data <- normal_neg_data[normal_neg_data[,2] %in% c("GO:0000406","GO:0032181"),]
normal_neg_GOMF_data <- normal_neg_GOMF_data[c(2,1),]
normal_neg_GOMF_data <- normal_neg_GOMF_data[,c(5,6)]

disease_neg_GOMF_data <- disease_neg_data[disease_neg_data[,2] %in% c("GO:0000406","GO:0032181"), ]
#disease_neg_GOMF_data <- disease_neg_GOMF_data[c(2,1),]
disease_neg_GOMF_data <- disease_neg_GOMF_data[,c(5,6)]

#整理到一�?
GOMF_data <- cbind(normal_pos_GOMF_data,normal_neg_GOMF_data,disease_neg_GOMF_data)
GOMF_data <- GOMF_data[,c(1,2,3,4,5,7)]
colnames(GOMF_data) <- c("ontology","term name","normal_pos_p_value","total genes","normal_neg_p_value","disease_neg_p_value")
GOMF_data$disease_pos_p_value <- c(1,1)
GOMF_data <- GOMF_data[,c(1,2,3,7,5,6,4)]








########筛选Human Phenotype
normal_pos_HP_data <- normal_pos_data[normal_pos_data[,1] == "Human Phenotype", ]

disease_pos_HP_data <- disease_pos_data[disease_pos_data[,1] == "Human Phenotype", ]

normal_neg_HP_data <- normal_neg_data[normal_neg_data[,1] == "Human Phenotype",]

disease_neg_HP_data <- disease_neg_data[disease_neg_data[,1] == "Human Phenotype", ]

#四个数据框中第二列相同的数据
#common_HP_data <- Reduce(intersect, list(normal_pos_HP_data[, 2], normal_neg_HP_data[, 2], disease_neg_HP_data[, 2]))


#筛选交集数�?              data[data[, 1] %in% c(12, 13), ]
#normal_pos_HP_data <- normal_pos_data[normal_pos_data[,2] %in% common_HP_data,]

normal_pos_HP_data <- normal_pos_data[normal_pos_data[,2] %in% c("HP:0010444"),]
normal_pos_HP_data <- normal_pos_HP_data[,c(1,3,5,6)]

#disease_pos_HP_data <- disease_pos_data[disease_pos_data[,2] %in% c("HP:0002205","HP:0001633"), ]
#disease_pos_HP_data <- disease_pos_HP_data[c(2,1),]
#disease_pos_HP_data <- disease_pos_HP_data[,c(5,6)]

#normal_neg_HP_data <- normal_neg_data[normal_neg_data[,2] %in% c("HP:0002877","HP:0002619"),]
#normal_neg_HP_data <- normal_neg_HP_data[c(2,1),]
#normal_neg_HP_data <- normal_neg_HP_data[,c(5,6)]

#disease_neg_HP_data <- disease_neg_data[disease_neg_data[,2] %in% c("HP:0002877","HP:0002619"), ]
#disease_neg_HP_data <- disease_neg_HP_data[,c(5,6)]


#整理到一�?
HP_data <- cbind(normal_pos_HP_data)
#HP_data <- HP_data[,c(1,2,3,4,5,7)]
colnames(HP_data) <- c("ontology","term name","normal_pos_p_value","total genes")
HP_data$disease_pos_p_value <- c(1)
HP_data$normal_neg_p_value <- c(1)
HP_data$disease_neg_p_value <- c(1)
HP_data <- HP_data[,c(1,2,3,5,6,7,4)]


#将上面的数据整理到一�?
heatmap_data <- rbind(GOBP_data,GOCC_data,GOMF_data,HP_data)
total_genes <- heatmap_data$`total genes`
heatmap_data <- heatmap_data[,-1]
rownames(heatmap_data) <- heatmap_data[,1]
heatmap_data <- heatmap_data[,-1]
heatmap_data <- heatmap_data[,-5]


heatmap_data_1 <- -log10(heatmap_data)
heatmap_data_1$type <- c(rep("GO Biological Process",2),rep("GO Cellular Component",2),rep("GO Molecular Function",2),rep("Human Phenotype",1))
heatmap_data_1$type <- factor(heatmap_data_1$type)
# 假设你的数据框名为data，需要将�?4列转换为数值型
heatmap_data_1[, 1:4] <- lapply(heatmap_data_1[, 1:4], as.numeric)
#heatmap_data_1$type <- as.numeric(heatmap_data_1$type)
heatmap_data_2 <- as.matrix(heatmap_data_1[,1:4])


#首先，自定义颜色�?
col_fun = colorRamp2(c(0,20), c("white", "red"))

#继续尝试绘制热图，这里取消列聚类，添加cell边框�?
Heatmap(heatmap_data_2[,1:4],col = col_fun,
        name = "-log10(p_value)",
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        rect_gp = gpar(col= "white"),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        row_names_side = "left")


# �?5个柱子一个颜色，后面4个一个颜色，后面2个一个颜色，最�?3个一个颜�?
colors <- c(rep("#1f77b4", 2), rep("#ff7f0e", 2), rep("#9467bd", 2), rep("#2ca02c", 1))
#为type设置分组颜色                                                
colors_type <- c("#1f77b4", "#ff7f0e", "#9467bd", "#2ca02c")
names(colors_type) <- unique(heatmap_data_1$type)



bar2 = rowAnnotation(
  sum2 = anno_barplot(
    total_genes,
    bar_width = 0.9,
    #gp = gpar(col = "white", fill = "orange"),
    gp = gpar(col = "white", fill = colors),
    border = T,
    axis_param = list(at = c(0,54,108),
                      labels = c("0","54","108")),
    width = unit(1, "cm")), show_annotation_name = F)





###绘制热图
Heatmap(heatmap_data_2[,1:4],col = col_fun,
        width = unit(8, "cm"),
        height = unit(8, "cm"),
        name = "-log10(p_value)",
        cluster_columns = F,
        cluster_rows = F,
        show_row_dend = F,
        show_column_dend = F,
        rect_gp = gpar(col= "white",lwd = 1),
        #top_annotation = bar1,
        left_annotation =bar2,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        column_title = "Increased with age    Decreased with age",
        column_title_gp = gpar(fontsize = 12),
        cell_fun = function(j, i, x, y, width, height, fill, data=10^(-heatmap_data_2[,1:4])){
          value <- data[i, j]  # 使用行列名称来访问数�?
          formatted_value <- format(value, scientific = TRUE)
          # 将科学计数法字符串解析为数字
          scientific_notation_value <- as.numeric(formatted_value)
          # 四舍五入整数部分
          integer_part <- round(scientific_notation_value)
          # 提取小数部分
          decimal_part <- scientific_notation_value - integer_part
          # 将整数部分和小数部分重新格式化为字符�?
          formatted_value <- paste0(format(decimal_part, digits = 2))
          grid.text(formatted_value, x, y, gp = gpar(fontsize = 10))
        })+
  Heatmap(heatmap_data_1$type,name = "Ontology",width = unit(3, "mm"),col = colors_type)


#3000 CpG
# 读取normal_pos_data文件
normal_pos_data <- read.table("data/GREAT_heatmap/normal_positive_3000.tsv", header = F, sep = "\t")
normal_pos_data <- normal_pos_data[,c(1,2,3,4,5,20)]
colnames(normal_pos_data) <- c("Ontology","ID","Term Name","Rank","P-Value","Total Genes")

# 读取disease_pos_data文件
disease_pos_data <- read.table("data/GREAT_heatmap/disease_positive_3000.tsv", header = F, sep = "\t")
disease_pos_data <- disease_pos_data[,c(1,2,3,4,5,20)]
colnames(disease_pos_data) <- c("Ontology","ID","Term Name","Rank","P-Value","Total Genes")

# 读取normal_neg_data文件
normal_neg_data <- read.table("data/GREAT_heatmap/normal_negative_3000.tsv", header = F, sep = "\t")
normal_neg_data <- normal_neg_data[,c(1,2,3,4,5,20)]
colnames(normal_neg_data) <- c("Ontology","ID","Term Name","Rank","P-Value","Total Genes")

# 读取disease_neg_data文件
disease_neg_data <- read.table("data/GREAT_heatmap/disease_negative_3000.tsv", header = F, sep = "\t")
disease_neg_data <- disease_neg_data[,c(1,2,3,4,5,20)]
colnames(disease_neg_data) <- c("Ontology","ID","Term Name","Rank","P-Value","Total Genes")



########筛选GO Biological Process
normal_pos_GOBP_data <- normal_pos_data[normal_pos_data[,1] == "GO Biological Process", ]

disease_pos_GOBP_data <- disease_pos_data[disease_pos_data[,1] == "GO Biological Process", ]

normal_neg_GOBP_data <- normal_neg_data[normal_neg_data[,1] == "GO Biological Process",]

disease_neg_GOBP_data <- disease_neg_data[disease_neg_data[,1] == "GO Biological Process", ]


#四个数据框中第二列相同的数据
common_GOBP_data <- Reduce(intersect, list(normal_pos_GOBP_data[, 2], disease_pos_GOBP_data[, 2],normal_neg_GOBP_data[,2],disease_neg_GOBP_data[,2]))

#筛选交集数�? data[data[, 1] %in% c(12, 13), ]
normal_pos_GOBP_data <- normal_pos_data[normal_pos_data[,2] %in% common_GOBP_data,]

normal_pos_GOBP_data <- normal_pos_data[normal_pos_data[,2] %in% c("GO:0038202","GO:0061052","GO:0032053"),]
normal_pos_GOBP_data <- normal_pos_GOBP_data[,c(1,3,5,6)]

disease_pos_GOBP_data <- disease_pos_data[disease_pos_data[,2] %in% c("GO:0038202","GO:0061052","GO:0032053"), ]
#disease_pos_GOBP_data <- disease_pos_GOBP_data[c(2,1),]
disease_pos_GOBP_data <- disease_pos_GOBP_data[,c(5,6)]

normal_neg_GOBP_data <- normal_neg_data[normal_neg_data[,2] %in% c("GO:0038202","GO:0061052","GO:0032053"),]
normal_neg_GOBP_data <- normal_neg_GOBP_data[c(1,3,2),]
normal_neg_GOBP_data <- normal_neg_GOBP_data[,c(5,6)]

disease_neg_GOBP_data <- disease_neg_data[disease_neg_data[,2] %in% c("GO:0038202","GO:0061052","GO:0032053"), ]
disease_neg_GOBP_data <- disease_neg_GOBP_data[c(2,3,1),]
disease_neg_GOBP_data <- disease_neg_GOBP_data[,c(5,6)]

#整理到一�?
GOBP_data <- cbind(normal_pos_GOBP_data,disease_pos_GOBP_data,normal_neg_GOBP_data,disease_neg_GOBP_data)
GOBP_data <- GOBP_data[,c(1,2,3,4,5,7,9)]
colnames(GOBP_data) <- c("ontology","term name","normal_pos_p_value","total genes","disease_pos_p_value","normal_neg_p_value","disease_neg_p_value")
#GOBP_data$disease_neg_p_value <- c(1)
#GOBP_data$disease_neg_p_value <- c(1,1)
GOBP_data <- GOBP_data[,c(1,2,3,5,6,7,4)]





########筛选GO Cellular Component
normal_pos_GOCC_data <- normal_pos_data[normal_pos_data[,1] == "GO Cellular Component", ]

disease_pos_GOCC_data <- disease_pos_data[disease_pos_data[,1] == "GO Cellular Component", ]

normal_neg_GOCC_data <- normal_neg_data[normal_neg_data[,1] == "GO Cellular Component",]

disease_neg_GOCC_data <- disease_neg_data[disease_neg_data[,1] == "GO Cellular Component", ]

#四个数据框中第二列相同的数据
#common_GOCC_data <- Reduce(intersect, list(disease_pos_GOCC_data[, 2], normal_neg_GOCC_data[,2],disease_neg_GOCC_data[, 2]))


#筛选交集数�?              data[data[, 1] %in% c(12, 13), ]
#normal_pos_GOCC_data <- normal_pos_data[normal_pos_data[,2] %in% common_GOCC_data,]

#normal_pos_GOCC_data <- normal_pos_data[normal_pos_data[,2] %in% c("GO:0005764","GO:0005765"),]
#normal_pos_GOCC_data <- normal_pos_GOCC_data[,c(1,3,5,6)]

#disease_pos_GOCC_data <- disease_pos_data[disease_pos_data[,2] %in% common_GOCC_data,]

#disease_pos_GOCC_data <- disease_pos_data[disease_pos_data[,2] %in% c("GO:0045323","GO:0097440","GO:0031519"), ]
#disease_pos_GOCC_data <- disease_pos_GOCC_data[c(2,1),]
#disease_pos_GOCC_data <- disease_pos_GOCC_data[,c(1,3,5,6)]

normal_neg_GOCC_data <- normal_neg_data[normal_neg_data[,2] %in% c("GO:0031931"),]
#normal_neg_GOCC_data <- normal_neg_GOCC_data[c(2,3,1),]
normal_neg_GOCC_data <- normal_neg_GOCC_data[,c(1,3,5,6)]

#disease_neg_GOCC_data <- disease_neg_data[disease_neg_data[,2] %in% c("GO:0045323","GO:0097440","GO:0031519"), ]
#disease_neg_GOCC_data <- disease_neg_GOCC_data[c(1,3,2),]
#disease_neg_GOCC_data <- disease_neg_GOCC_data[,c(5,6)]

#整理到一�?
GOCC_data <- cbind(normal_neg_GOCC_data)
#GOCC_data <- GOCC_data[,c(1,2,3,4,5,7)]
colnames(GOCC_data) <- c("ontology","term name","normal_neg_p_value","total genes")
GOCC_data$normal_pos_p_value <- c(1)
GOCC_data$disease_pos_p_value <- c(1)
GOCC_data$disease_neg_p_value <- c(1)
GOCC_data <- GOCC_data[,c(1,2,5,6,3,7,4)]








########筛选GO Molecular Function
normal_pos_GOMF_data <- normal_pos_data[normal_pos_data[,1] == "GO Molecular Function", ]

disease_pos_GOMF_data <- disease_pos_data[disease_pos_data[,1] == "GO Molecular Function", ]

normal_neg_GOMF_data <- normal_neg_data[normal_neg_data[,1] == "GO Molecular Function",]

disease_neg_GOMF_data <- disease_neg_data[disease_neg_data[,1] == "GO Molecular Function", ]

#四个数据框中第二列相同的数据
common_GOMF_data <- Reduce(intersect, list(normal_pos_GOMF_data[, 2],normal_neg_GOMF_data[,2],disease_neg_GOMF_data[,2]))

#筛选交集数�?              data[data[, 1] %in% c(12, 13), ]
#normal_pos_GOMF_data <- normal_pos_GOMF_data[normal_pos_GOMF_data[,2] %in% common_GOMF_data,]

normal_pos_GOMF_data <- normal_pos_data[normal_pos_data[,2] %in% common_GOMF_data,]

normal_pos_GOMF_data <- normal_pos_data[normal_pos_data[,2] %in% c("GO:0005179","GO:0016278"),]
normal_pos_GOMF_data <- normal_pos_GOMF_data[,c(1,3,5,6)]

#disease_pos_GOMF_data <- disease_pos_data[disease_pos_data[,2] %in% c("GO:0005179","GO:0016278"), ]
#disease_pos_GOMF_data <- disease_pos_GOMF_data[c(2,1),]
#disease_pos_GOMF_data <- disease_pos_GOMF_data[,c(5,6)]

normal_neg_GOMF_data <- normal_neg_data[normal_neg_data[,2] %in% c("GO:0005179","GO:0016278"),]
normal_neg_GOMF_data <- normal_neg_GOMF_data[c(2,1),]
normal_neg_GOMF_data <- normal_neg_GOMF_data[,c(5,6)]

disease_neg_GOMF_data <- disease_neg_data[disease_neg_data[,2] %in% c("GO:0005179","GO:0016278"), ]
disease_neg_GOMF_data <- disease_neg_GOMF_data[c(2,1),]
disease_neg_GOMF_data <- disease_neg_GOMF_data[,c(5,6)]

#整理到一�?
GOMF_data <- cbind(normal_pos_GOMF_data,normal_neg_GOMF_data,disease_neg_GOMF_data)
GOMF_data <- GOMF_data[,c(1,2,3,4,5,7)]
colnames(GOMF_data) <- c("ontology","term name","normal_pos_p_value","total genes","normal_neg_p_value","disease_neg_p_value")
GOMF_data$disease_pos_p_value <- c(1,1)
GOMF_data <- GOMF_data[,c(1,2,3,7,5,6,4)]








########筛选Human Phenotype
normal_pos_HP_data <- normal_pos_data[normal_pos_data[,1] == "Human Phenotype", ]

disease_pos_HP_data <- disease_pos_data[disease_pos_data[,1] == "Human Phenotype", ]

normal_neg_HP_data <- normal_neg_data[normal_neg_data[,1] == "Human Phenotype",]

disease_neg_HP_data <- disease_neg_data[disease_neg_data[,1] == "Human Phenotype", ]

#四个数据框中第二列相同的数据
#common_HP_data <- Reduce(intersect, list(normal_pos_HP_data[, 2],normal_neg_HP_data[, 2]))


#筛选交集数�?              data[data[, 1] %in% c(12, 13), ]
#normal_pos_HP_data <- normal_pos_data[normal_pos_data[,2] %in% common_HP_data,]

normal_pos_HP_data <- normal_pos_data[normal_pos_data[,2] %in% c("HP:0008214"),]
normal_pos_HP_data <- normal_pos_HP_data[,c(1,3,5,6)]

#disease_pos_HP_data <- disease_pos_data[disease_pos_data[,2] %in% c("HP:0002205","HP:0001633"), ]
#disease_pos_HP_data <- disease_pos_HP_data[c(2,1),]
#disease_pos_HP_data <- disease_pos_HP_data[,c(5,6)]

#normal_neg_HP_data <- normal_neg_data[normal_neg_data[,2] %in% c("HP:0045046","HP:0006705"),]
#normal_neg_HP_data <- normal_neg_HP_data[c(2,1),]
#normal_neg_HP_data <- normal_neg_HP_data[,c(5,6)]

#disease_neg_HP_data <- disease_neg_data[disease_neg_data[,2] %in% c("HP:0002877","HP:0002619"), ]
#disease_neg_HP_data <- disease_neg_HP_data[,c(5,6)]


#整理到一�?
HP_data <- cbind(normal_pos_HP_data)
#HP_data <- HP_data[,c(1,2,3,4,5)]
colnames(HP_data) <- c("ontology","term name","normal_pos_p_value","total genes")
HP_data$disease_pos_p_value <- c(1)
HP_data$normal_neg_p_value <- c(1)
HP_data$disease_neg_p_value <- c(1)
HP_data <- HP_data[,c(1,2,3,5,6,7,4)]

#将上面的数据整理到一�?
heatmap_data <- rbind(GOBP_data,GOCC_data,GOMF_data,HP_data)
total_genes <- heatmap_data$`total genes`
heatmap_data <- heatmap_data[,-1]
rownames(heatmap_data) <- heatmap_data[,1]
heatmap_data <- heatmap_data[,-1]
heatmap_data <- heatmap_data[,-5]

heatmap_data_1 <- -log10(heatmap_data)
heatmap_data_1$type <- c(rep("GO Biological Process",3),rep("GO Cellular Component",1),rep("GO Molecular Function",2),rep("Human Phenotype",1))
heatmap_data_1$type <- factor(heatmap_data_1$type)
# 假设你的数据框名为data，需要将�?4列转换为数值型
heatmap_data_1[, 1:4] <- lapply(heatmap_data_1[, 1:4], as.numeric)
#heatmap_data_1$type <- as.numeric(heatmap_data_1$type)
heatmap_data_2 <- as.matrix(heatmap_data_1[,1:4])


#首先，自定义颜色�?
col_fun = colorRamp2(c(0,30), c("white", "red"))

#继续尝试绘制热图，这里取消列聚类，添加cell边框�?
Heatmap(heatmap_data_2[,1:4],col = col_fun,
        name = "-log10(p_value)",
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        rect_gp = gpar(col= "white"),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        row_names_side = "left")


# �?5个柱子一个颜色，后面4个一个颜色，后面2个一个颜色，最�?3个一个颜�?
colors <- c(rep("#1f77b4", 3), rep("#ff7f0e", 1), rep("#9467bd", 2), rep("#2ca02c", 1))
#为type设置分组颜色                                                
colors_type <- c("#1f77b4", "#ff7f0e", "#9467bd", "#2ca02c")
names(colors_type) <- unique(heatmap_data_1$type)



bar2 = rowAnnotation(
  sum2 = anno_barplot(
    total_genes,
    bar_width = 0.9,
    #gp = gpar(col = "white", fill = "orange"),
    gp = gpar(col = "white", fill = colors),
    border = T,
    axis_param = list(at = c(0,59,117),
                      labels = c("0","59","117")),
    width = unit(1, "cm")), show_annotation_name = F)





###绘制热图
Heatmap(heatmap_data_2[,1:4],col = col_fun,
        width = unit(8, "cm"),
        height = unit(8, "cm"),
        name = "-log10(p_value)",
        cluster_columns = F,
        cluster_rows = F,
        show_row_dend = F,
        show_column_dend = F,
        rect_gp = gpar(col= "white",lwd = 1),
        #top_annotation = bar1,
        left_annotation =bar2,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        column_title = "Increased with age    Decreased with age",
        column_title_gp = gpar(fontsize = 12),
        cell_fun = function(j, i, x, y, width, height, fill, data=10^(-heatmap_data_2[,1:4])){
          value <- data[i, j]  # 使用行列名称来访问数�?
          formatted_value <- format(value, scientific = TRUE)
          # 将科学计数法字符串解析为数字
          scientific_notation_value <- as.numeric(formatted_value)
          # 四舍五入整数部分
          integer_part <- round(scientific_notation_value)
          # 提取小数部分
          decimal_part <- scientific_notation_value - integer_part
          # 将整数部分和小数部分重新格式化为字符�?
          formatted_value <- paste0(format(decimal_part, digits = 2))
          grid.text(formatted_value, x, y, gp = gpar(fontsize = 10))
        })+
  Heatmap(heatmap_data_1$type,name = "Ontology",width = unit(3, "mm"),col = colors_type)
