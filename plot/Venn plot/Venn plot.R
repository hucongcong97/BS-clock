# 加载ggvenn�?--绘图�?
library(ggvenn)
# 加载RColorBrewer--生成颜色�?
library(RColorBrewer)
data <- read.csv("Venn_data.csv")
set1 <- which(data$Blood == 1)
set2 <- which(data$Brain == 1)
set3 <- which(data$Lung == 1)
set4 <- which(data$Skin == 1)

# 查看四个集合的交集元�?
intersect(intersect(set1,set2),intersect(set3,set4))

# 把四个集合的向量合并成一个列�?
# 注意这里的PART1、PART2、PART3、PART4就是后续展现在图里的集合名称
#vennlist1 <- list(PART1=part1,PART2=part2,PART3=part3,PART4=part4)
vennlist1 <- list(Blood=set1,Brain=set2,Lung=set3,Skin=set4)
#绘制韦恩�?
p <- ggvenn(vennlist1,
            fill_color=c(brewer.pal(8, 'Set2')[2:5]), # 设置填充颜色
            stroke_size=0.5,  # 集合圆圈的线�?
            set_name_size=3.3, # 集合名称的文本大�?
            digits = 0, # 小数点后保留位数
            show_percentage = F # 是否展示每一部分所占的百分�?
)
p          
ggsave("venn.pdf",p,width = 4,height = 4)