### This script analyses GTF2I mutation data and links with quiescence and expressin levels in tumours.
##### 初步数据处理 aa mutaion status与Exp QS合并#####
### Read mutation data downloaded from cbioportal:
setwd("/Users/lizixi/Desktop/Final Project/5.18 task")
mutations_gtf2i <- read.delim("mutations_cbioportal_gtf2i.txt")
head(mutations_gtf2i)
## How many samples have a mutation?
table(mutations_gtf2i$GTF2I)
tb_gtf2i <- table(mutations_gtf2i$GTF2I)
sort(tb_gtf2i) #descend升序排列PHIN1列的数字元素

# Add a new column that just denotes WT and MUT for GTF2I:
mutations_gtf2i[which(mutations_gtf2i$GTF2I == "NS"),]$GTF2I <- "WT" #****将GTF2I中的NS替换为WT
mutations_gtf2i$GTF2I_status <- mutations_gtf2i$GTF2I #为mutations_gtf2i表添加一列为GTF2I_status，同时将这一列定义为GTF2I（和GTF2I一样）
mutations_gtf2i[which(mutations_gtf2i$GTF2I != "WT"),]$GTF2I_status <- "MUT" #选出mutation中GTF2I列中不等于“WT”的元素，并在GTF2I_status中将其全赋值为“MUT”
head(mutations_gtf2i)
# How many samples with mutations_gtf2i in GTF2I?
table(mutations_gtf2i$GTF2I_status)
# Read expression data for GTF2I:
expression_gtf2i <- read.delim("mRNA Expression_RSEM_cbioportal_gtf2i.txt")
head(expression_gtf2i)
# Plot histogram of expression values of GTF2I:
hist(expression_gtf2i$GTF2I)
# It's good to take the log2 of expression values to avoid working with very large numbers:
expression_gtf2i$GTF2I_log2expression <- log2(expression_gtf2i$GTF2I+1)
head(expression_gtf2i)
hist(expression_gtf2i$GTF2I_log2expression) # more decent value

# Merge GTF2I expression and mutation data:
df.merged_gtf2i <- merge(mutations_gtf2i[,c("SAMPLE_ID","GTF2I_status","STUDY_ID")], #mutation中需要合并的纵列是sample id，gnaq_status
                          expression_gtf2i[,c("SAMPLE_ID","GTF2I_log2expression")], #expresssion中需要合并的纵列是sample ID，GTF2I_log2expression
                          by.x="SAMPLE_ID", by.y="SAMPLE_ID", #x和y表合并的依据是相同的sample id
                          all.x=FALSE, all.y=FALSE) #x和y表中不需要所有的行都合并（只合并sample id一样的，只在其中一个表中才存在的sample id丢弃）
head(df.merged_gtf2i)
# Merge TP53 mutation data+QS
df.merged_gtf2i <- merge(df.merged_gtf2i[,c("SAMPLE_ID","GTF2I_status","STUDY_ID","GTF2I_log2expression")], #mutation中需要合并的纵列是sample id，gtf2i_status
                          df.merged2[,c("SAMPLE_ID","TP53_status","QuiescenceScore","Group")], #expresssion中需要合并的纵列是sample ID，GTF2I_log2expression
                          by.x="SAMPLE_ID", by.y="SAMPLE_ID", #x和y表合并的依据是相同的sample id
                          all.x=FALSE, all.y=FALSE) #x和y表中不需要所有的行都合并（只合并sample id一样的，只在其中一个表中才存在的sample id丢弃）
head(df.merged_gtf2i)

#####Expression and QS levels - GTF2I_aamutation#####
#QS
## by boxplot
ggboxplot(df.merged_gtf2i, x = "GTF2I_status", y = "QuiescenceScore",
          color = "GTF2I_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")
#Exp
## by boxplot
ggboxplot(df.merged_gtf2i, x = "GTF2I_status", y = "GTF2I_log2expression",
          color = "GTF2I_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")
