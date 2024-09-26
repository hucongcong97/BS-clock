# Manhattan plot
# import packages
#install.packages("CMplot")
library(CMplot)

# read blood data
data <- read.csv("Manhattan_data_Blood_new.csv")
data1 <- data[,c(2,3,5,7)]
data1 <- data1[,c(3,1,2,4)]
colnames(data1) <- c("SNP","CHR","BP","P")
data1 <- na.omit(data1)
data1$CHR <- gsub("chr", "", data1$CHR)
data1 <- data1[data1$CHR != "M", ]
#odata1$SNP <- paste("cpg_", data1$SNP, sep = "")

SNPs <- data1[data1[,4] < 1e-08, 1]
CMplot(data1,
       type = "p",
       plot.type="m",
       LOG10=TRUE,
       col= c("#3E0A52", "#423D77","#3F678B",
              "#468C8D", "#5FB47F", "#9FD55C","#F9E956"),
       highlight = SNPs,
       highlight.col = NULL,
       highlight.cex = 1,
       #highlight.pch = c(15:17), 
       highlight.text = SNPs,      
       highlight.text.col = "black",
       highlight.text.cex = 0.5,
       threshold = c(1e-06,1e-08),
       threshold.lt=c(2,1),
       amplify = FALSE,
       file = "pdf",
       dpi = 300,
       file.output = T,
       verbose = T,
       width = 13,height = 6,chr.labels.angle=45)


# read brain data
data <- read.csv("Manhattan_data_Brain_new.csv")
data1 <- data[,c(2,3,5,7)]
data1 <- data1[,c(3,1,2,4)]
colnames(data1) <- c("SNP","CHR","BP","P")
data1 <- na.omit(data1)
data1$CHR <- gsub("chr", "", data1$CHR)
data1 <- data1[data1$CHR != "M", ]
#odata1$SNP <- paste("cpg_", data1$SNP, sep = "")

SNPs <- data1[data1[,4] < 1e-10, 1]

chr_order <- data.frame(CHR = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"), 
                        order = 1:24)

data1$CHR <- factor(data1$CHR, levels = chr_order$CHR)
data1 <- data1[order(data1$CHR),]


CMplot(data1,
       type = "p",
       plot.type="m",
       LOG10=TRUE,
       col= c("#3E0A52", "#423D77","#3F678B",
              "#468C8D", "#5FB47F", "#9FD55C","#F9E956"),
       highlight = SNPs,
       highlight.col = NULL,
       highlight.cex = 1,
       #highlight.pch = c(15:17), 
       highlight.text = SNPs,      
       highlight.text.col = "black",
       highlight.text.cex = 0.5,
       #threshold = 5e-02,
       threshold=c(5e-02,1e-10),
       threshold.lt=c(2,1),
       amplify = FALSE,
       file = "pdf",
       dpi = 300,
       file.output = T,
       verbose = T,
       width = 13,height = 6,chr.labels.angle=45)

# read lung data
data <- read.csv("Manhattan_data_Lung_new.csv")
data1 <- data[,c(2,3,5,7)]
data1 <- data1[,c(3,1,2,4)]
colnames(data1) <- c("SNP","CHR","BP","P")
data1 <- na.omit(data1)
data1$CHR <- gsub("chr", "", data1$CHR)
data1 <- data1[data1$CHR != "M", ]
#data1$SNP <- paste("cpg_", data1$SNP, sep = "")

SNPs <- data1[data1[,4] < 1e-04, 1]

chr_order <- data.frame(CHR = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"), 
                        order = 1:24)

data1$CHR <- factor(data1$CHR, levels = chr_order$CHR)
data1 <- data1[order(data1$CHR),]

CMplot(data1,
       type = "p",
       plot.type="m",
       LOG10=TRUE,
       col= c("#3E0A52", "#423D77","#3F678B",
              "#468C8D", "#5FB47F", "#9FD55C","#F9E956"),
       highlight = SNPs,
       highlight.col = NULL,
       highlight.cex = 1,
       #highlight.pch = c(15:17), 
       highlight.text = SNPs,      
       highlight.text.col = "black",
       highlight.text.cex = 0.5 ,
       threshold = c(1e-02,1e-04),
       threshold.lt=c(2,1),
       amplify = FALSE,
       file = "pdf",
       dpi = 300,
       file.output = T,
       verbose = T,
       width = 13,height = 6,chr.labels.angle=45)

# read skin data
data <- read.csv("Manhattan_data_Skin_new.csv")
data1 <- data[,c(2,3,5,7)]
data1 <- data1[,c(3,1,2,4)]
colnames(data1) <- c("SNP","CHR","BP","P")
data1 <- na.omit(data1)
data1$CHR <- gsub("chr", "", data1$CHR)
data1 <- data1[data1$CHR != "M", ]
#data1$SNP <- paste("cpg_", data1$SNP, sep = "")

SNPs <- data1[data1[,4] < 1e-03, 1]

chr_order <- data.frame(CHR = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"), 
                        order = 1:24)

data1$CHR <- factor(data1$CHR, levels = chr_order$CHR)
data1 <- data1[order(data1$CHR),]


CMplot(data1,
       type = "p",
       plot.type="m",
       LOG10=TRUE,
       col= c("#3E0A52", "#423D77","#3F678B",
              "#468C8D", "#5FB47F", "#9FD55C","#F9E956"),
       highlight = SNPs,
       highlight.col = NULL,
       highlight.cex = 1,
       #highlight.pch = c(15:17), 
       highlight.text = SNPs,      
       highlight.text.col = "black",
       highlight.text.cex = 0.5,
       threshold = 1e-03,
       amplify = FALSE,
       file = "pdf",
       dpi = 300,
       file.output = T,
       verbose = T,
       width = 13,height = 6,chr.labels.angle=45)
