# library packages
library(ggvenn)
library(RColorBrewer)

# read data
data <- read.csv("Venn_data.csv")
set1 <- which(data$Blood == 1)
set2 <- which(data$Brain == 1)
set3 <- which(data$Lung == 1)
set4 <- which(data$Skin == 1)

# Get the intersection of four data
intersect(intersect(set1,set2),intersect(set3,set4))

# Note that PART1, PART2, PART3, and PART4 are the names of the collections 
# shown in the figure later
# vennlist1 <- list(PART1=part1,PART2=part2,PART3=part3,PART4=part4)
vennlist1 <- list(Blood=set1,Brain=set2,Lung=set3,Skin=set4)
# plot venn
p <- ggvenn(vennlist1,
            fill_color=c(brewer.pal(8, 'Set2')[2:5]), 
            stroke_size=0.5,  
            set_name_size=3.3, 
            digits = 0, 
            show_percentage = F 
)
p          
ggsave("venn.pdf",p,width = 4,height = 4)