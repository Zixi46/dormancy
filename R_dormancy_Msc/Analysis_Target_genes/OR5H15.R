### This script analyses OR5H15 mutation data and links with quiescence and expressin levels in tumours.
##### 初步数据处理 aa mutaion status与Exp QS合并#####
### Read mutation data downloaded from cbioportal:
setwd("/Users/lizixi/Desktop/Final Project/5.18 task")
mutations_or5h15 <- read.delim("mutations_cbioportal_or5h15.txt")
head(mutations_or5h15)
## How many samples have a mutation?
table(mutations_or5h15$OR5H15)
tb_or5h15 <- table(mutations_or5h15$OR5H15)
sort(tb_or5h15) #descend升序排列PHIN1列的数字元素

# Add a new column that just denotes WT and MUT for OR5H15:
mutations_or5h15[which(mutations_or5h15$OR5H15 == "NS"),]$OR5H15 <- "WT" #****将OR5H15中的NS替换为WT
mutations_or5h15$OR5H15_status <- mutations_or5h15$OR5H15 #为mutations_or5h15表添加一列为OR5H15_status，同时将这一列定义为OR5H15（和OR5H15一样）
mutations_or5h15[which(mutations_or5h15$OR5H15 != "WT"),]$OR5H15_status <- "MUT" #选出mutation中OR5H15列中不等于“WT”的元素，并在OR5H15_status中将其全赋值为“MUT”
head(mutations_or5h15)
# How many samples with mutations_or5h15 in OR5H15?
table(mutations_or5h15$OR5H15_status)
# Read expression data for OR5H15:
expression_or5h15 <- read.delim("mRNA Expression_RSEM_cbioportal_or5h15.txt")
head(expression_or5h15)
# Plot histogram of expression values of OR5H15:
hist(expression_or5h15$OR5H15)
# It's good to take the log2 of expression values to avoid working with very large numbers:
expression_or5h15$OR5H15_log2expression <- log2(expression_or5h15$OR5H15+1)
head(expression_or5h15)
hist(expression_or5h15$OR5H15_log2expression) # more decent value

# Merge OR5H15 expression and mutation data:
df.merged_or5h15 <- merge(mutations_or5h15[,c("SAMPLE_ID","OR5H15_status","STUDY_ID")], #mutation中需要合并的纵列是sample id，gnaq_status
                          expression_or5h15[,c("SAMPLE_ID","OR5H15_log2expression")], #expresssion中需要合并的纵列是sample ID，OR5H15_log2expression
                          by.x="SAMPLE_ID", by.y="SAMPLE_ID", #x和y表合并的依据是相同的sample id
                          all.x=FALSE, all.y=FALSE) #x和y表中不需要所有的行都合并（只合并sample id一样的，只在其中一个表中才存在的sample id丢弃）
head(df.merged_or5h15)
# Merge TP53 mutation data+QS
df.merged_or5h15 <- merge(df.merged_or5h15[,c("SAMPLE_ID","OR5H15_status","STUDY_ID","OR5H15_log2expression")], #mutation中需要合并的纵列是sample id，or5h15_status
                          df.merged2[,c("SAMPLE_ID","TP53_status","QuiescenceScore","Group")], #expresssion中需要合并的纵列是sample ID，OR5H15_log2expression
                          by.x="SAMPLE_ID", by.y="SAMPLE_ID", #x和y表合并的依据是相同的sample id
                          all.x=FALSE, all.y=FALSE) #x和y表中不需要所有的行都合并（只合并sample id一样的，只在其中一个表中才存在的sample id丢弃）
head(df.merged_or5h15)

#####Expression and QS levels - OR5H15_aamutation#####
#QS
## by boxplot
ggboxplot(df.merged_or5h15, x = "OR5H15_status", y = "QuiescenceScore",
          color = "OR5H15_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")
#Exp
## by boxplot
ggboxplot(df.merged_or5h15, x = "OR5H15_status", y = "OR5H15_log2expression",
          color = "OR5H15_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")
