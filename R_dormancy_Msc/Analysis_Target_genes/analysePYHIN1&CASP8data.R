### This script analyses PYHIN1 mutation data and links with quiescence levels in tumours.
### Read mutation data downloaded from cbioportal:
setwd("/Users/lizixi/Desktop/Final Project/5.18 task")
mutations <- read.delim("mutations_cbioportal.txt")
head(mutations)
## How many samples have a mutation?
table(mutations$PYHIN1)
tb <- table(mutations$PYHIN1)
sort(tb) #descend升序排列PHIN1列的数字元素
# Add a new column that just denotes WT and MUT for PYHIN1:
mutations[which(mutations$PYHIN1 == "NS"),]$PYHIN1 <- "WT" #****将PYHIN1中的NS替换为WT
mutations$PYHIN1_status <- mutations$PYHIN1 #为mutations表添加一列为PYHIN1_status，同时将这一列定义为PYHIN1（和PYHIN1一样）
mutations[which(mutations$PYHIN1 != "WT") ,]$PYHIN1_status <- "MUT" #选出mutation中PYHIN1列中不等于“WT”的元素，并在PYHIN1_status中将其全赋值为“MUT”
head(mutations)
table(mutations$PYHIN1_status)

## >>> Do the same for TP53:
## How many samples have a mutation?
table(mutations$TP53)
tb_TP53 <- table(mutations$TP53)
sort(tb_TP53) #descend升序排列TP53列的数字元素
getOption("max.print") #查看最大结果显示行，只有1000行，不能显示完整，
options(max.print = 5000) #改变最大的结果显示行
sort(tb_TP53) 
# Add a new column that just denotes WT and MUT for TP53:
mutations[which(mutations$TP53 == "NS"),]$TP53 <- "NS"#****将TP53中的NS替换为WT
mutations$TP53_status <- mutations$TP53 #为mutations表添加一列为TP53_status，同时将这一列定义为TP53（和TP53一样）
mutations[which(mutations$TP53 != "WT"),]$TP53_status <- "MUT" #选出mutation中TP53列中不等于“WT”的元素，并在TP53_status中将其全赋值为“MUT”
head(mutations)
# How many samples with mutations in PYHIN1?
table(mutations$TP53_status)



##### Do mutations in PYHIN1 increase or decrease PYHIN1 expression compared to samples that have a WT PYHIN1?#####
# Read expression data for PYHIN1:
expression <- read.delim("mRNAExpression_RSEM_cbioportal.txt")
head(expression)
# Plot histogram of expression values of PYHIN1:
hist(expression$PYHIN1)
# It's good to take the log2 of expression values to avoid working with very large numbers:
expression$PYHIN1_log2expression <- log2(expression$PYHIN1+1)
head(expression)
hist(expression$PYHIN1_log2expression) # more decent value
# Merge expression and mutation data:
df.merged <- merge(mutations[,c("SAMPLE_ID","TP53_status","PYHIN1_status")], #mutation中需要合并的纵列是sample id，tp53_status和pyhin1_status
                   expression[,c("SAMPLE_ID","PYHIN1_log2expression")], #expresssion中需要合并的纵列是sample ID，PYHIN1_log2expression
                   by.x="SAMPLE_ID", by.y="SAMPLE_ID", #x和y表合并的依据是相同的sample id
                   all.x=FALSE, all.y=FALSE) #x和y表中不需要所有的行都合并（只合并sample id一样的，只在其中一个表中才存在的sample id丢弃）
head(df.merged)


# Now compare the expression levels of PYHIN1 by the mutation status by boxplot:
library(ggplot2)
library(ggpubr)
my_comparisons <- list(c("WT","MUT"))
pdf("PYHIN1expression.byMutationStatus.pdf")
table(df.merged$PYHIN1_status)#看看MUT 和 WT的数量
ggboxplot(df.merged, x = "PYHIN1_status", y = "PYHIN1_log2expression",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+ #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)
dev.off()
### >>>Try and create a nicer violin plot for this.
ggviolin(df.merged, x = "PYHIN1_status", y = "PYHIN1_log2expression", fill = "PYHIN1_status",
         palette = c("#00AFBB", "#E7B800"),
         add = "boxplot", add.params = list(fill = "white"),order = c("WT", "MUT"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 12) # Add global the p-value
dev.off()


## Next, read quiescence data, merge it with the existing table and compare quiescence levels 
# between groups with/without PYHIN1 mutation and with/without TP53 mutation.

# Read quiescence data:
quiescence <- read.delim("PancancerQuiescenceGroups.txt")
head(quiescence)
# Plot histogram of quiescence values:
hist(quiescence$QuiescenceScore) # decent value
# Merge quiescence, expression and mutation data:
df.merged2 <- merge(df.merged[,c("SAMPLE_ID","TP53_status","PYHIN1_status","PYHIN1_log2expression")], #合并之前的表达合并表和静止分数表
                   quiescence[,c("Barcode","QuiescenceScore","Group")], 
                   by.x="SAMPLE_ID", by.y="Barcode", #x和y表合并的依据是相同的sample id ,在y表中列名为barcode
                   all.x=FALSE, all.y=FALSE) #x和y表中不需要所有的行都合并（只合并sample id一样的，只在其中一个表中才存在的sample id丢弃）
head(df.merged2) # 合并不的列名的表时，合并表的列名被前一表覆盖，都统一为SAMPLE_ID了
##### Compare quiescence levels between WT and MUT PYHIN1 cases. #####
# by boxplot
ggboxplot(df.merged2, x = "PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  labs(x=NULL)+
  stat_compare_means(comparisons = my_comparisons)
dev.off()
# by violin plot
ggviolin(df.merged2, x = "PYHIN1_status", y = "QuiescenceScore", fill = "PYHIN1_status",
         palette = c("#00AFBB", "#E7B800"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 40)
dev.off()


###### Repeat all of the above and what is otherwise included in the script for P53 #####
# Compare quiescence levels between WT and MUT TP53 cases.
# by boxplot
table(df.merged2$TP53_status)
ggboxplot(df.merged2, x = "TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter", order = c("WT", "MUT"))+
  stat_compare_means(comparisons = my_comparisons)
# by violin plot
ggviolin(df.merged2, x = "TP53_status", y = "QuiescenceScore", fill = "TP53_status",
         palette = c("#00AFBB", "#E7B800"),
         add = "boxplot", add.params = list(fill = "white"), order = c("WT", "MUT"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 40)

#Quiescence group的组间比例分布图 for PYHIN1 and TP53(因为评分不能反映静止或是快速增殖)
# for PYHIN1
table(df.merged2$PYHIN1_status) #看看wt和mut的数量，差距比较大
ggplot(df.merged2, mapping = aes(
  x = PYHIN1_status, fill = Group), order = c("WT", "MUT"))+
  geom_bar(stat="count",width=0.5,position='fill')+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  geom_text(stat='count',aes(label=scales::percent(..count../sum(..count..))), color="white", size=3.5,position=position_fill(0.5))+
  theme_minimal()
# for TP53
table(df.merged2$TP53_status) #数量差距正常
ggplot(df.merged2, mapping = aes(
  x = TP53_status, fill = Group), order = c("WT", "MUT"))+
  geom_bar(stat="count",width=0.5,position='fill')+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  geom_text(stat='count',aes(label=scales::percent(..count../sum(..count..))), color="white", size=3.5,position=position_fill(0.5))+
  theme_minimal()


##### for CASP8 from cbioportal TCGA Pan-cancer Atlas ##### 
#Repeat the same analysis you did for PYHIN1 for casp8.
### Read mutation data downloaded from cbioportal:
mutations2 <- read.delim("mutations_cbioportal_casp8.txt")
head(mutations2)
## How many samples have a mutation?
table(mutations2$CASP8)
tb_CASP8 <- table(mutations2$CASP8)
sort(tb_CASP8) #descend升序排列CASP8列的数字元素
# Add a new column that just denotes WT and MUT for CASP8:
mutations2[which(mutations2$CASP8 == "NS"),]$CASP8 <- "WT"#****将PYHIN1中的NS替换为WT
mutations2$CASP8_status <- mutations2$CASP8 #为mutations表添加一列为CASP8_status，同时将这一列定义为CASP8（和CASP8一样）
mutations2[which(mutations2$CASP8 != "WT"),]$CASP8_status <- "MUT" #选出mutation中CASP8列中不等于“WT”的元素，并在CASP8_status中将其全赋值为“MUT”
head(mutations2)
# How many samples with mutations in CASP8?
table(mutations2$CASP8_status)

## Do mutations in CASP8 increase or decrease CASP8 expression compared to samples that have a WT CASP8?
# Read expression data for PYHIN1:
expression2 <- read.delim("mRNA Expression_RSEM_cbioportal_casp8.txt")
head(expression2)
# Plot histogram of expression values of PYHIN1:
hist(expression2$CASP8)
# It's good to take the log2 of expression values to avoid working with very large numbers:
expression2$CASP8_log2expression <- log2(expression2$CASP8+1)
head(expression2)
hist(expression2$CASP8_log2expression) # more decent value
# Merge expression and mutation data of casp8:
df.merged3 <- merge(mutations2[,c("SAMPLE_ID","CASP8_status")], #mutation中需要合并的纵列是sample id，tp53_status和CASP8_status
                   expression2[,c("SAMPLE_ID","CASP8_log2expression")], #expresssion中需要合并的纵列是sample ID，CASP8_log2expression
                   by.x="SAMPLE_ID", by.y="SAMPLE_ID", #x和y表合并的依据是相同的sample id
                   all.x=FALSE, all.y=FALSE) #x和y表中不需要所有的行都合并（只合并sample id一样的，只在其中一个表中才存在的sample id丢弃）
head(df.merged3)
# Merge expression/mutation data of casp8+mutation data of TP53+mutation data of PYHIN1+quiencense level
df.merged4 <- merge(df.merged3[,c("SAMPLE_ID","CASP8_status","CASP8_log2expression")], 
                   df.merged2[,c("SAMPLE_ID","TP53_status","PYHIN1_status","PYHIN1_log2expression","QuiescenceScore","Group")],
                   by.x="SAMPLE_ID", by.y="SAMPLE_ID", #x和y表合并的依据是相同的sample id
                   all.x=FALSE, all.y=FALSE) 
head(df.merged4)

# Now compare the expression levels of CASP8 by the mutation status by boxplot:
table(df.merged4$CASP8_status) #数量差距正常
my_comparisons <- list(c("WT","MUT"))
pdf("CASP8expression.byMutationStatus.pdf")
ggboxplot(df.merged4, x = "PYHIN1_status", y = "PYHIN1_log2expression",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")
dev.off()
### >>>Try and create a nicer violin plot for this.
ggviolin(df.merged4, x = "CASP8_status", y = "CASP8_log2expression", fill = "CASP8_status",
         palette = c("#00AFBB", "#E7B800"),
         add = "boxplot", add.params = list(fill = "white"),order = c("WT", "MUT"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 15) # Add global the p-value
dev.off()
# Compare quiescence levels between WT and MUT CASP8 cases. 
## by boxplot
ggboxplot(df.merged4, x = "CASP8_status", y = "QuiescenceScore",
          color = "CASP8_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)
dev.off()
## by violin plot
ggviolin(df.merged4, x = "CASP8_status", y = "QuiescenceScore", fill = "CASP8_status",
         palette = c("#00AFBB", "#E7B800"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 50)
dev.off()
#Quiescence group的组间比例分布图 for CASP8(因为评分不能反映静止或是快速增殖)
table(df.merged4$CASP8_status) #看看wt和mut的数量，差距比较大
ggplot(df.merged4, mapping = aes(
  x = CASP8_status, fill = Group), order = c("WT", "MUT"))+
  geom_bar(stat="count",width=0.5,position='fill')+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  geom_text(stat='count',aes(label=scales::percent(..count../sum(..count..))), color="white", size=3.5,position=position_fill(0.5))+
  theme_minimal()


#####Combine PYHIN1+TP53 status and CASP8+TP53 status to create 4 groups that you then compare the quiescence levels for them #####
#Combine PYHIN1+TP53 to create 4 groups 。用tidyr::unite合并两列（PYHIN1_status，TP53_status）的元素为新的一列（PYHIN1_TP53_status），默认以_连接（可以修改）
library(tidyr)
df.merged5 <- tidyr::unite(df.merged4, "PYHIN1_TP53_status", PYHIN1_status, TP53_status, remove = FALSE)
head(df.merged5)
table(df.merged5$PYHIN1_TP53_status)
#Combine CASP8+TP53 to create 4 groups 用tidyr::unite合并两列的元素为新的一列，默认以_连接（可以修改）
df.merged6 <- tidyr::unite(df.merged5, "CASP8_TP53_status", CASP8_status, TP53_status, remove = FALSE)
head(df.merged6)
table(df.merged6$CASP8_TP53_status)
#Reorder the column 重新给列排序
df.merged6<-df.merged6[,c("SAMPLE_ID","TP53_status","PYHIN1_status","CASP8_status","PYHIN1_TP53_status","CASP8_TP53_status","PYHIN1_log2expression","CASP8_log2expression","QuiescenceScore","Group")]  # column_name 设置为期望的排列顺序
head(df.merged6)

# Compare quiescence levels between 4 groups of PYHIN1-TP53 cases
## by boxplot
table(df.merged6$PYHIN1_TP53_status)
my_comparisons2 <- list(c("MUT_WT","MUT_MUT"),c("WT_MUT","MUT_MUT"),c("WT_WT","WT_MUT"),c("WT_MUT","MUT_WT"),c("WT_WT","MUT_MUT"),c("WT_WT","MUT_WT")) #显示4组之间的p
ggboxplot(df.merged6, x = "PYHIN1_TP53_status", y = "QuiescenceScore",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons2)
dev.off()
## by violin plot
ggviolin(df.merged6, x = "PYHIN1_TP53_status", y = "QuiescenceScore", fill = "PYHIN1_TP53_status",
         palette = c("#00AFBB", "#E7B800","#D2691E","#008000"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons2, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 50)
dev.off()
##Quiescence group的组间比例分布图 (因为评分不能反映静止或是快速增殖)
table(df.merged6$PYHIN1_TP53_status) #看看wt和mut的数量，差距比较大
ggplot(df.merged6, mapping = aes(
  x = PYHIN1_TP53_status, fill = Group), order = c("WT_WT","MUT_WT","WT_MUT","MUT_MUT"))+
  geom_bar(stat="count",width=0.5,position='fill')+
  scale_fill_manual(values=c('#999999','#E69F00','#D2691E','#008000'))+
  geom_text(stat='count',aes(label=scales::percent(..count../sum(..count..))), color="white", size=3.5,position=position_fill(0.5))+
  theme_minimal()

# Compare quiescence levels between 4 groups of CASP8-TP53 cases
## by boxplot
table(df.merged6$CASP8_TP53_status)
ggboxplot(df.merged6, x = "CASP8_TP53_status", y = "QuiescenceScore",
          color = "CASP8_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons2)
dev.off()
## by violin plot
ggviolin(df.merged6, x = "CASP8_TP53_status", y = "QuiescenceScore", fill = "CASP8_TP53_status",
         palette = c("#00AFBB", "#E7B800","#D2691E","#008000"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons2, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 50)
dev.off()
##Quiescence group的组间比例分布图(因为评分不能反映静止或是快速增殖)
table(df.merged6$CASP8_TP53_status) #看看wt和mut的数量，差距比较大
ggplot(df.merged6, mapping = aes(
  x = CASP8_TP53_status, fill = Group), order = c("WT_WT","MUT_WT","WT_MUT","MUT_MUT"))+
  geom_bar(stat="count",width=0.5,position='fill')+
  scale_fill_manual(values=c('#999999','#E69F00','#D2691E','#008000'))+
  geom_text(stat='count',aes(label=scales::percent(..count../sum(..count..))), color="white", size=3.5,position=position_fill(0.5))+
  theme_minimal()


##Q:NS要去除吗？因为常常成对出现，被判定为MUT且数量多，在三个表里数量也一致为524（所以基本可以判断为NON SIGNIFICANT)
sort(tb)
sort(tb_TP53)
sort(tb_CASP8)

#**Compare quiescence levels between 4 groups of PYHIN1-CASP8 cases
##Combine PYHIN1+CASP8 to create 4 groups 用tidyr::unite合并两列的元素为新的一列，默认以_连接（可以修改）
df.merged7 <- tidyr::unite(df.merged5, "PYHIN1_CASP8_status", PYHIN1_status, CASP8_status, remove = FALSE)
head(df.merged7)
table(df.merged7$PYHIN1_CASP8_status)
##Compare quiescence levels between 4 groups of PYHIN1-CASP8 cases
### by boxplot
ggboxplot(df.merged7, x = "PYHIN1_CASP8_status", y = "QuiescenceScore",
          color = "PYHIN1_CASP8_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons2)
dev.off()
## by violin plot
ggviolin(df.merged7, x = "PYHIN1_CASP8_status", y = "QuiescenceScore", fill = "PYHIN1_CASP8_status",
         palette = c("#00AFBB", "#E7B800","#D2691E","#008000"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons2, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 50)
dev.off()

###### 增加PYHIN1突变类型分析列 “AMP”“DEL” #####
#读取cna文件
cna_PYHIN1 <- read.delim("PYHIN1_cna_fromcBioportal.txt")
head(cna_PYHIN1)
table(cna_PYHIN1$PYHIN1_CNA)
#将cna文件中的ID 和 CNA与总表合并
head(df.merged7)
head(cna_PYHIN1)
df.merged8 <- merge(df.merged7[,c("SAMPLE_ID","PYHIN1_CASP8_status","CASP8_status",
                                  "CASP8_log2expression","PYHIN1_TP53_status",
                                  "TP53_status","PYHIN1_status","PYHIN1_log2expression",
                                  "QuiescenceScore","Group")], 
                    cna_PYHIN1[,c("STUDY_ID","SAMPLE_ID","PYHIN1_CNA")],
                    by.x="SAMPLE_ID", by.y="SAMPLE_ID", #x和y表合并的依据是相同的sample id
                    all.x=FALSE, all.y=FALSE) 
head(df.merged8)
table(df.merged8$PYHIN1_CNA)
#数据整理，完整案例分析
##exclude cns-NP(not profiled) samples from the analysis 剔除包含NP的样本
##把NP都换为NA
df.merged8[which(df.merged8$PYHIN1_CNA == "NP"),]$PYHIN1_CNA <- "NA" #****将PYHIN1_CAN中的NP替换为NA
table(df.merged8$PYHIN1_CNA)
##把-2换为-1（数量太少了），把2换为1
  #"-2" is a deep loss, possibly a homozygous deletion
  #"-1" is a single-copy loss (heterozygous deletion)
  #"0" is diploid
  #"1" indicates a low-level gain
  #"2" is a high-level amplification.
##将1转化为AMP,-1转化为DEL，O转化为NCN
df.merged8[which(df.merged8$PYHIN1_CNA == "-2"),]$PYHIN1_CNA <- "-1" 
df.merged8[which(df.merged8$PYHIN1_CNA == "2"),]$PYHIN1_CNA <- "1" 
df.merged8[which(df.merged8$PYHIN1_CNA == "-1"),]$PYHIN1_CNA <- "DEL" 
df.merged8[which(df.merged8$PYHIN1_CNA == "1"),]$PYHIN1_CNA <- "AMP" 
df.merged8[which(df.merged8$PYHIN1_CNA == "0"),]$PYHIN1_CNA <- "NCN" 
table(df.merged8$PYHIN1_CNA)
##去掉NA样本，生成df.merged9
df.merged9 <- df.merged8[!grepl("NA", df.merged8$PYHIN1_CNA),]
table(df.merged9$PYHIN1_CNA)
##将pyhin1突变信息和CNA信息结合为新的一列Mutation_CNA,生成6个不同组别
df.merged9 <- tidyr::unite(df.merged9, "Pyhin1_Pyhin1CNA", PYHIN1_status, PYHIN1_CNA, remove = FALSE)
head(df.merged9)
table(df.merged9$Pyhin1_Pyhin1CNA)


#####*****以下为错误代码*****#####


#结合CNA重看PYhin1表达
###设置x轴分组排序
df.merged9$Pyhin1_Pyhin1CNA <- factor(df.merged9$Pyhin1_Pyhin1CNA, levels=c('WT_AMP','WT_NCN','WT_DEL','MUT_AMP','MUT_NCN','MUT_DEL'))
##boxplot
my_comparisons3 <- list(c("MUT_AMP","MUT_NCN","MUT_DEL","WT_AMP","WT_NCN","WT_DEL"))
ggboxplot(df.merged9, x = "Pyhin1_Pyhin1CNA", y = "PYHIN1_log2expression",
          color = "Pyhin1_Pyhin1CNA", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE"),
          add = "jitter")+ #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT_AMP","WT_NCN")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("WT_AMP","WT_DEL")), label = "p.signif")+
  stat_compare_means(comparisons = list(c("WT_NCN","WT_DEL")), label = "p.signif")+
  stat_compare_means(comparisons = list(c("MUT_AMP","MUT_NCN")), label = "p.signif")+
  stat_compare_means(comparisons = list(c("MUT_AMP","MUT_DEL")), label = "p.signif")+
  stat_compare_means(comparisons = list(c("MUT_NCN","MUT_DEL")), label = "p.signif")
##violin plot
ggviolin(df.merged9, x = "Pyhin1_Pyhin1CNA", y = "PYHIN1_log2expression", fill = "Pyhin1_Pyhin1CNA",
        palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE"),
        add = "boxplot", add.params = list(fill = "white"))+ #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 15)
#结合CNA重看PYHIN1的quiescence level
##boxplot
my_comparisons3 <- list(c("MUT_AMP","MUT_NCN","MUT_DEL","WT_AMP","WT_NCN","WT_DEL"))
ggboxplot(df.merged9, x = "Pyhin1_Pyhin1CNA", y = "QuiescenceScore",
          color = "Pyhin1_Pyhin1CNA", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE"),
          add = "jitter")+ #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT_AMP","WT_NCN")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("WT_AMP","WT_DEL")), label = "p.signif")+
  stat_compare_means(comparisons = list(c("WT_NCN","WT_DEL")), label = "p.signif")+
  stat_compare_means(comparisons = list(c("MUT_AMP","MUT_NCN")), label = "p.signif")+
  stat_compare_means(comparisons = list(c("MUT_AMP","MUT_DEL")), label = "p.signif")+
  stat_compare_means(comparisons = list(c("MUT_NCN","MUT_DEL")), label = "p.signif")
##violin plot
ggviolin(df.merged9, x = "Pyhin1_Pyhin1CNA", y = "QuiescenceScore", fill = "Pyhin1_Pyhin1CNA",
         palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE"),
         add = "boxplot", add.params = list(fill = "white"))+ #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 50)
#结合CNA重看PYHIN1的quiescence level
##将Mutation_CNA与TP53_status结合为新列,生成df.merged10
df.merged10 <- tidyr::unite(df.merged9, "Pyhin1_Pyhin1CNA_TP53", Pyhin1_Pyhin1CNA, TP53_status, remove = FALSE)
table(df.merged10$Pyhin1_Pyhin1CNA_TP53)
###设置x轴分组排序
df.merged10$Pyhin1_Pyhin1CNA_TP53 <- factor(df.merged10$Pyhin1_Pyhin1CNA_TP53, levels=c("WT_AMP_WT","WT_AMP_MUT","WT_NCN_WT","WT_NCN_MUT","WT_DEL_WT","WT_DEL_MUT",
                                                                                        "MUT_AMP_WT","MUT_AMP_MUT","MUT_NCN_WT","MUT_NCN_MUT","MUT_DEL_WT","MUT_DEL_MUT"))
##boxplot
my_comparisons4 <- list(c("WT_AMP_WT","WT_AMP_MUT","WT_NCN_WT","WT_NCN_MUT","WT_DEL_WT","WT_DEL_MUT",
                          "MUT_AMP_WT","MUT_AMP_MUT","MUT_NCN_WT","MUT_NCN_MUT","MUT_DEL_WT","MUT_DEL_MUT"))
plot1 <- ggboxplot(df.merged10, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
          color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                      "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
          add = "jitter") #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons4)
###调整横坐标的角度为45°
plot1 <- plot1+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
##violinplot
plot2 <- ggviolin(df.merged10, x = "Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore", fill = "Pyhin1_Pyhin1CNA_TP53",
         palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
         add = "boxplot", add.params = list(fill = "white"))+ #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons4, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 50)
###调整横坐标的角度为45°
plot2 <-plot2+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
#绘制p值的矩阵热图

#分癌症类型，重复以上步骤，比较quiescence score
##看看有哪些癌症种类
table(df.merged10$STUDY_ID) #29种

#####Quiescence levels - PYHIN1_PYHIN1CNA_TP53 - Cancer type#####
#把每种癌症都分为一个df
#再分类绘制quiescencescore plot（boxplot 感觉violinplot没有意义）
##1 acc
acc<- df.merged10[df.merged10$STUDY_ID=="acc_tcga_pan_can_atlas_2018",]
table(acc$Pyhin1_Pyhin1CNA_TP53)
my_comparisons4 <- list(c("WT_AMP_WT","WT_AMP_MUT","WT_NCN_WT","WT_NCN_MUT","WT_DEL_WT","WT_DEL_MUT",
                          "MUT_AMP_WT","MUT_AMP_MUT","MUT_NCN_WT","MUT_NCN_MUT","MUT_DEL_WT","MUT_DEL_MUT"))
###boxplot
plot1_acc <- ggboxplot(acc, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                   color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                               "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                   add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_acc <- plot1_acc+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_acc

##2 blac
blca<- df.merged10[df.merged10$STUDY_ID=="blca_tcga_pan_can_atlas_2018",]
table(blca$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_blca <- ggboxplot(blca, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                       color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                   "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                       add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_blca <- plot1_blca+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_blca

##3 brca
brca<- df.merged10[df.merged10$STUDY_ID=="brca_tcga_pan_can_atlas_2018",]
table(brca$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_brca <- ggboxplot(brca, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_brca <- plot1_brca+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_brca

##4 cesc
cesc<- df.merged10[df.merged10$STUDY_ID=="cesc_tcga_pan_can_atlas_2018",]
table(cesc$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_cesc <- ggboxplot(cesc, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_cesc <- plot1_cesc+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_cesc

##5 chol
chol<- df.merged10[df.merged10$STUDY_ID=="chol_tcga_pan_can_atlas_2018",]
table(chol$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_chol <- ggboxplot(chol, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_chol <- plot1_chol+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_chol

##6 coadread
coadread<- df.merged10[df.merged10$STUDY_ID=="coadread_tcga_pan_can_atlas_2018",]
table(coadread$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_coadread <- ggboxplot(coadread, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_coadread <- plot1_coadread+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_coadread

##7 esca
esca<- df.merged10[df.merged10$STUDY_ID=="esca_tcga_pan_can_atlas_2018",]
table(esca$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_esca <- ggboxplot(esca, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                            color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                        "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                            add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_esca <- plot1_esca+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_esca

##8gbm
gbm<- df.merged10[df.merged10$STUDY_ID=="gbm_tcga_pan_can_atlas_2018",]
table(gbm$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_gbm <- ggboxplot(gbm, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_gbm <- plot1_gbm+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_gbm

##9hnsc
hnsc<- df.merged10[df.merged10$STUDY_ID=="hnsc_tcga_pan_can_atlas_2018",]
table(hnsc$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_hnsc <- ggboxplot(hnsc, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                       color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                   "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                       add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_hnsc <- plot1_hnsc+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_hnsc

##10kich
kich<- df.merged10[df.merged10$STUDY_ID=="kich_tcga_pan_can_atlas_2018",]
table(kich$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_kich <- ggboxplot(kich, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_kich <- plot1_kich+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_kich

##11kirc
kirc<- df.merged10[df.merged10$STUDY_ID=="kirc_tcga_pan_can_atlas_2018",]
table(kirc$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_kirc <- ggboxplot(kirc, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_kirc <- plot1_kirc+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_kirc

##12kirp
kirp<- df.merged10[df.merged10$STUDY_ID=="kirp_tcga_pan_can_atlas_2018",]
table(kirp$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_kirp <- ggboxplot(kirp, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_kirp <- plot1_kirp+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_kirp

##13lgg
lgg<- df.merged10[df.merged10$STUDY_ID=="lgg_tcga_pan_can_atlas_2018",]
table(lgg$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_lgg <- ggboxplot(lgg, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_lgg <- plot1_lgg+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_lgg

##14lihc
lihc<- df.merged10[df.merged10$STUDY_ID=="lihc_tcga_pan_can_atlas_2018",]
table(lihc$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_lihc <- ggboxplot(lihc, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                       color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                   "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                       add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_lihc <- plot1_lihc+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_lihc

##15luad
luad<- df.merged10[df.merged10$STUDY_ID=="luad_tcga_pan_can_atlas_2018",]
table(luad$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_luad <- ggboxplot(luad, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_luad <- plot1_luad+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_luad

##16lusc
lusc<- df.merged10[df.merged10$STUDY_ID=="lusc_tcga_pan_can_atlas_2018",]
table(lusc$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_lusc <- ggboxplot(lusc, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_lusc <- plot1_lusc+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_lusc

##17meso
meso<- df.merged10[df.merged10$STUDY_ID=="meso_tcga_pan_can_atlas_2018",]
table(meso$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_meso <- ggboxplot(meso, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_meso <- plot1_meso+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_meso

##18ov
ov<- df.merged10[df.merged10$STUDY_ID=="ov_tcga_pan_can_atlas_2018",]
table(ov$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_ov <- ggboxplot(ov, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_ov <- plot1_ov+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_ov

##19paad
paad<- df.merged10[df.merged10$STUDY_ID=="paad_tcga_pan_can_atlas_2018",]
table(paad$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_paad <- ggboxplot(paad, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                      color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                  "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                      add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_paad <- plot1_paad+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_paad

##20pcpg
pcpg<- df.merged10[df.merged10$STUDY_ID=="pcpg_tcga_pan_can_atlas_2018",]
table(pcpg$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_pcpg <- ggboxplot(pcpg, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_pcpg <- plot1_pcpg+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_pcpg

##21prad
prad<- df.merged10[df.merged10$STUDY_ID=="prad_tcga_pan_can_atlas_2018",]
table(prad$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_prad <- ggboxplot(prad, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_prad <- plot1_prad+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_prad

##22sarc
sarc<- df.merged10[df.merged10$STUDY_ID=="sarc_tcga_pan_can_atlas_2018",]
table(sarc$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_sarc <- ggboxplot(sarc, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_sarc <- plot1_sarc+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_sarc

##23stad
stad<- df.merged10[df.merged10$STUDY_ID=="stad_tcga_pan_can_atlas_2018",]
table(stad$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_stad <- ggboxplot(stad, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_stad <- plot1_stad+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_stad

##24tgct
tgct<- df.merged10[df.merged10$STUDY_ID=="tgct_tcga_pan_can_atlas_2018",]
table(tgct$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_tgct <- ggboxplot(tgct, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_tgct <- plot1_tgct+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_tgct

##25thca
thca<- df.merged10[df.merged10$STUDY_ID=="thca_tcga_pan_can_atlas_2018",]
table(thca$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_thca <- ggboxplot(thca, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_thca <- plot1_thca+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_thca

##26thym
thym<- df.merged10[df.merged10$STUDY_ID=="thym_tcga_pan_can_atlas_2018",]
table(thym$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_thym <- ggboxplot(thym, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_thym <- plot1_thym+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_thym

##27ucec
ucec<- df.merged10[df.merged10$STUDY_ID=="ucec_tcga_pan_can_atlas_2018",]
table(ucec$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_ucec <- ggboxplot(ucec, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_ucec <- plot1_ucec+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_ucec

##28ucs
ucs<- df.merged10[df.merged10$STUDY_ID=="ucs_tcga_pan_can_atlas_2018",]
table(ucs$Pyhin1_Pyhin1CNA_TP53)
###boxplot
plot1_ucs <- ggboxplot(ucs, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                        color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                        add = "jitter") #添加jitter后，会自动隐藏异常值。
plot1_ucs <- plot1_ucs+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_ucs

##29uvm
uvm<- df.merged10[df.merged10$STUDY_ID=="uvm_tcga_pan_can_atlas_2018",]
table(uvm$Pyhin1_Pyhin1CNA_TP53)
###boxplot
my_comparisons5 <- list(c("WT_AMP_WT","WT_NCN_WT"))
plot1_uvm <- ggboxplot(uvm, x ="Pyhin1_Pyhin1CNA_TP53", y = "QuiescenceScore",
                       color = "Pyhin1_Pyhin1CNA_TP53", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                   "#D5D2C1","#00AFBB","#D2691E","#370D00","#D387D8","#A13E97"),
                       add = "jitter")+ #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)
plot1_uvm <- plot1_uvm+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#调整横坐标角度
plot1_uvm

#####Quiescence levels - PYHIN1CNA + Cancer type#####
table(df.merged10$PYHIN1_CNA)
#检查下expression
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
ggboxplot(df.merged10, x = "PYHIN1_CNA", y = "PYHIN1_log2expression",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter")+ #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)
#总三类boxplot
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
ggboxplot(df.merged10, x = "PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter")+ #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

#分癌症类型
##1 acc
table(acc$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(acc, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##2 blca
table(blca$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(blca, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##3 brca
table(brca$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(brca, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##4 cesc
table(cesc$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(cesc, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##5 chol
table(chol$PYHIN1_CNA)
###boxplot
ggboxplot(chol, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("AMP","NCN")))

##6 coadread
table(coadread$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(coadread, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##7 esca
table(esca$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(esca, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##8 gbm
table(gbm$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(gbm, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##9 hnsc
table(hnsc$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(hnsc, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##10 kich
table(kich$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(kich, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##11 kirc
table(kirc$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(kirc, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##12 kirp
table(kirp$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(kirp, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##13 lgg
table(lgg$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(lgg, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##14 lihc
table(lihc$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(lihc, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##15 luad
table(luad$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(luad, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##16 lusc
table(lusc$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(lusc, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##17 meso
table(meso$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(meso, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##18 ov
table(ov$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(ov, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##19 paad
table(paad$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(paad, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##20 pcpg
table(pcpg$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(pcpg, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##21 prad
table(prad$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(prad, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##22 sarc
table(sarc$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(sarc, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)


##23 stad
table(stad$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(stad, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##24 tgct
table(tgct$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(tgct, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##25 thca
table(thca$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(thca, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##26 thym
table(thym$PYHIN1_CNA)
###boxplot
ggboxplot(thym, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("AMP","NCN")))

##27 ucec
table(ucec$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(ucec, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##28 ucs
table(ucs$PYHIN1_CNA)
my_comparisons5 <- list(c("AMP","NCN"),c("AMP","DEL"),c("NCN","DEL"))
###boxplot
ggboxplot(ucs, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons5)

##29 uvm
table(uvm$PYHIN1_CNA)
###boxplot
ggboxplot(uvm, x ="PYHIN1_CNA", y = "QuiescenceScore",
          color = "PYHIN1_CNA", palette =c("#00AFBB", "#E7B800","#D2691E"),
          add = "jitter",order = c("AMP", "NCN","DEL"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("AMP","NCN")))

#####Quiescence levels - PYHIN1_status + Cancer type#####
##1 acc
table(acc$PYHIN1_status)
###boxplot
ggboxplot(acc, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##2 blca
table(blca$PYHIN1_status)
###boxplot
ggboxplot(blca, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##3 brca
table(brca$PYHIN1_status)
###boxplot
ggboxplot(brca, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##4 cesc
table(cesc$PYHIN1_status)
###boxplot
ggboxplot(cesc, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##5 chol
table(chol$PYHIN1_status)
###boxplot
ggboxplot(chol, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##6 coadread
table(coadread$PYHIN1_status)
###boxplot
ggboxplot(coadread, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##7 esca
table(esca$PYHIN1_status)
###boxplot
ggboxplot(esca, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)


##8 gbm
table(gbm$PYHIN1_status)
###boxplot
ggboxplot(gbm, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##9 hnsc
table(hnsc$PYHIN1_status)
###boxplot
ggboxplot(hnsc, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##10 kich
table(kich$PYHIN1_status)
###boxplot
ggboxplot(kich, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##11 kirc
table(kirc$PYHIN1_status)
###boxplot
ggboxplot(kirc, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##12 kirp
table(kirp$PYHIN1_status)
###boxplot
ggboxplot(kirp, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##13 lgg
table(lgg$PYHIN1_status)
###boxplot
ggboxplot(lgg, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##14 lihc
table(lihc$PYHIN1_status)
###boxplot
ggboxplot(lihc, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##15 luad
table(luad$PYHIN1_status)
###boxplot
ggboxplot(luad, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##16 lusc
table(lusc$PYHIN1_status)
###boxplot
ggboxplot(lusc, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##17 meso
table(meso$PYHIN1_status)
###boxplot
ggboxplot(meso, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##18 ov
table(ov$PYHIN1_status)
###boxplot
ggboxplot(ov, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##19 paad
table(paad$PYHIN1_status)
###boxplot
ggboxplot(paad, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##20 pcpg
table(pcpg$PYHIN1_status)
###boxplot
ggboxplot(pcpg, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##21 prad
table(prad$PYHIN1_status)
###boxplot
ggboxplot(prad, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##22 sarc
table(sarc$PYHIN1_status)
###boxplot
ggboxplot(sarc, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##23 stad
table(stad$PYHIN1_status)
###boxplot
ggboxplot(stad, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##24 tgct
table(tgct$PYHIN1_status)
###boxplot
ggboxplot(tgct, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##25 thca
table(thca$PYHIN1_status)
###boxplot
ggboxplot(thca, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##26 thym
table(thym$PYHIN1_status)
###boxplot
ggboxplot(thym, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##27 ucec
table(ucec$PYHIN1_status)
###boxplot
ggboxplot(ucec, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##28 ucs
table(ucs$PYHIN1_status)
###boxplot
ggboxplot(ucs, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##29 uvm
table(uvm$PYHIN1_status)
###boxplot
ggboxplot(uvm, x ="PYHIN1_status", y = "QuiescenceScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

#####Quiescence levels - TP53_status + Cancer type#####
##1 acc
table(acc$TP53_status)
###boxplot
ggboxplot(acc, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##2 blca
table(blca$TP53_status)
###boxplot
ggboxplot(blca, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##3 brca
table(brca$TP53_status)
###boxplot
ggboxplot(brca, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##4 cesc
table(cesc$TP53_status)
###boxplot
ggboxplot(cesc, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##5 chol
table(chol$TP53_status)
###boxplot
ggboxplot(chol, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##6 coadread
table(coadread$TP53_status)
###boxplot
ggboxplot(coadread, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##7 esca
table(esca$TP53_status)
###boxplot
ggboxplot(esca, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##8 gbm
table(gbm$TP53_status)
###boxplot
ggboxplot(gbm, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##9 hnsc
table(hnsc$TP53_status)
###boxplot
ggboxplot(hnsc, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##10 kich
table(kich$TP53_status)
###boxplot
ggboxplot(kich, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##11 kirc
table(kirc$TP53_status)
###boxplot
ggboxplot(kich, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##12 kirp
table(kirp$TP53_status)
###boxplot
ggboxplot(kirp, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##13 lgg
table(lgg$TP53_status)
###boxplot
ggboxplot(lgg, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##14 lihc
table(lihc$TP53_status)
###boxplot
ggboxplot(lihc, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##15 luad
table(luad$TP53_status)
###boxplot
ggboxplot(luad, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##16 lusc
table(lusc$TP53_status)
###boxplot
ggboxplot(lusc, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##17 meso
table(meso$TP53_status)
###boxplot
ggboxplot(meso, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##18 ov
table(ov$TP53_status)
###boxplot
ggboxplot(ov, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##19 paad
table(paad$TP53_status)
###boxplot
ggboxplot(paad, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##20 pcpg
table(pcpg$TP53_status)
###boxplot
ggboxplot(pcpg, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##21 prad
table(prad$TP53_status)
###boxplot
ggboxplot(prad, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##22 sarc
table(sarc$TP53_status)
###boxplot
ggboxplot(sarc, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##23 stad
table(stad$TP53_status)
###boxplot
ggboxplot(stad, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##24 tgct
table(tgct$TP53_status)
###boxplot
ggboxplot(tgct, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##25 thca
table(thca$TP53_status)
###boxplot
ggboxplot(thca, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##26 thym
table(thym$TP53_status)
###boxplot
ggboxplot(thym, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##27 ucec
table(ucec$TP53_status)
###boxplot
ggboxplot(ucec, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##28 ucs
table(ucs$TP53_status)
###boxplot
ggboxplot(ucs, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

##29 uvm
table(uvm$TP53_status)
###boxplot
ggboxplot(uvm, x ="TP53_status", y = "QuiescenceScore",
          color = "TP53_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT"))+#添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)


#####Quiescence levels - PYHIN1_TP53_status + Cancer type#####
#封装 table：f_Table_PYHIN1_TP53_status(df)
f_Table_PYHIN1_TP53_status <- function(a) {
  table(a$PYHIN1_TP53_status)}
#封装 绘图：f_QS_PYHIN1_TP53_status(df)
f_QS_PYHIN1_TP53_status <- function(a) {
  table(a$PYHIN1_TP53_status)
  plot_QS_PYHIN1_TP53_status<-ggboxplot(a, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+stat_compare_means(comparisons = my_comparisons2)
  return(plot_QS_PYHIN1_TP53_status)
}
#分Cancer type
##1 acc
f_Table_PYHIN1_TP53_status(acc)
###p比较和f中的不对应
ggboxplot(acc, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
           color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
           add = "jitter",order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+stat_compare_means(comparisons = list(c("WT_WT","WT_MUT")))
##2 blca
f_Table_PYHIN1_TP53_status(blca)
f_QS_PYHIN1_TP53_status(blca)
##3 brca
f_Table_PYHIN1_TP53_status(brca)
f_QS_PYHIN1_TP53_status(brca)
##4 cesc
f_Table_PYHIN1_TP53_status(brca)
f_QS_PYHIN1_TP53_status(brca)
##5 chol
f_Table_PYHIN1_TP53_status(chol)
ggboxplot(chol, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = list(c("WT_WT","WT_MUT")))
##6 coadread
f_Table_PYHIN1_TP53_status(coadread)
f_QS_PYHIN1_TP53_status(coadread)
##7 esca
f_Table_PYHIN1_TP53_status(esca)
ggboxplot(esca, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = list(c("WT_WT","MUT_MUT"),c("WT_WT","WT_MUT"),c("WT_MUT","MUT_MUT")))
##8 gbm
f_Table_PYHIN1_TP53_status(gbm)
ggboxplot(esca, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = list(c("WT_WT","WT_MUT"),c("WT_WT","MUT_MUT"),c("WT_MUT","MUT_MUT")))
##9 hnsc
f_Table_PYHIN1_TP53_status(hnsc)
ggboxplot(hnsc, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = list(c("WT_WT","WT_MUT"),c("WT_WT","MUT_MUT"),c("WT_MUT","MUT_MUT")))
##10 kich
f_Table_PYHIN1_TP53_status(kich)
ggboxplot(kich, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = list(c("WT_MUT","WT_WT")))
##11 kirc
f_Table_PYHIN1_TP53_status(kirc)
ggboxplot(kirc, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = list(c("WT_MUT","WT_WT")))
##12 kirp
f_Table_PYHIN1_TP53_status(kirp)
ggboxplot(kirp, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = list(c("MUT_WT","WT_MUT"),c("MUT_WT","WT_WT"),c("WT_MUT","WT_WT")))
##13 lgg
f_Table_PYHIN1_TP53_status(lgg)
ggboxplot(lgg, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
         color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
         add = "jitter",order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = list(c("WT_MUT","WT_WT")))
##14 lihc
f_Table_PYHIN1_TP53_status(lihc)
f_QS_PYHIN1_TP53_status(lihc)
##15 luad
f_Table_PYHIN1_TP53_status(luad)
f_QS_PYHIN1_TP53_status(luad)
##16 lusc
f_Table_PYHIN1_TP53_status(lusc)
f_QS_PYHIN1_TP53_status(lusc)
##17 meso
f_Table_PYHIN1_TP53_status(meso)
ggboxplot(meso, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = list(c("WT_MUT","WT_WT")))
##18 ov
f_Table_PYHIN1_TP53_status(ov)
ggboxplot(ov, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = list(c("WT_WT","WT_MUT"),c("WT_WT","MUT_MUT"),c("WT_MUT","MUT_MUT")))
##19 paad
f_Table_PYHIN1_TP53_status(paad)
ggboxplot(paad, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = list(c("WT_WT","WT_MUT"),c("WT_WT","MUT_WT"),c("WT_MUT","MUT_WT")))
##20 pcpg
f_Table_PYHIN1_TP53_status(pcpg)
ggboxplot(pcpg, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = list(c("WT_MUT","WT_WT")))
##21 prad
f_Table_PYHIN1_TP53_status(prad)
f_QS_PYHIN1_TP53_status(prad)
##22 sarc
f_Table_PYHIN1_TP53_status(sarc)
ggboxplot(sarc, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = list(c("WT_WT","WT_MUT"),c("WT_WT","MUT_MUT"),c("WT_MUT","MUT_MUT")))
##23 stad
f_Table_PYHIN1_TP53_status(stad)
f_QS_PYHIN1_TP53_status(stad)
##24 tgct
f_Table_PYHIN1_TP53_status(tgct)
ggboxplot(tgct, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = list(c("WT_MUT","WT_WT")))
##25 thca
f_Table_PYHIN1_TP53_status(thca)
ggboxplot(thca, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = list(c("WT_MUT","WT_WT")))
##26 thym
f_Table_PYHIN1_TP53_status(thym)
ggboxplot(thym, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = list(c("WT_MUT","WT_WT")))
##27 ucec
f_Table_PYHIN1_TP53_status(ucec)
f_QS_PYHIN1_TP53_status(ucec)
##28 ucs
f_Table_PYHIN1_TP53_status(ucs)
ggboxplot(ucs, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = list(c("WT_MUT","WT_WT")))
##28 uvm
f_Table_PYHIN1_TP53_status(uvm)
f_QS_PYHIN1_TP53_status(uvm)

#####Quiescence levels - PYHIN1_PYHIN1CNA + Cancer type#####
my_comparisons3 <- list(c("WT_AMP","WT_NCN"),c("WT_AMP","WT"))
#封装table 和 boxplot
f_Table_PYHIN1_PYHINCNA <- function(a){
  table(a$Pyhin1_Pyhin1CNA)
}
f_QS_PYHIN1_PYHINCNA <- function(a){
  ggboxplot(a, x = "Pyhin1_Pyhin1CNA", y = "QuiescenceScore",
            color = "Pyhin1_Pyhin1CNA", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE"),
            add = "jitter",order = c("WT_AMP","WT_NCN","WT_DEL","MUT_AMP","MUT_NCN","MUT_DEL"))+ #添加jitter后，会自动隐藏异常值。
    stat_compare_means(comparisons = list(c("WT_AMP","WT_NCN")))+
    stat_compare_means(comparisons = list(c("WT_AMP","WT_DEL")))+
    stat_compare_means(comparisons = list(c("WT_NCN","WT_DEL")))+
    stat_compare_means(comparisons = list(c("MUT_AMP","MUT_NCN")))+
    stat_compare_means(comparisons = list(c("MUT_AMP","MUT_DEL")))+
    stat_compare_means(comparisons = list(c("MUT_NCN","MUT_DEL")))
}
##1 acc
f_Table_PYHIN1_PYHINCNA(acc)
f_QS_PYHIN1_PYHINCNA(acc)
##2 blca
f_Table_PYHIN1_PYHINCNA(blca)
f_QS_PYHIN1_PYHINCNA(blca)
##3 brca
f_Table_PYHIN1_PYHINCNA(brca)
f_QS_PYHIN1_PYHINCNA(brca)
##4 cesc
f_Table_PYHIN1_PYHINCNA(cesc)
f_QS_PYHIN1_PYHINCNA(cesc)
##5 chol
f_Table_PYHIN1_PYHINCNA(chol)
f_QS_PYHIN1_PYHINCNA(chol)
##6 coadread
f_Table_PYHIN1_PYHINCNA(coadread)
f_QS_PYHIN1_PYHINCNA(coadread)
##7 esca
f_Table_PYHIN1_PYHINCNA(esca)
f_QS_PYHIN1_PYHINCNA(esca)
##8 gbm
f_Table_PYHIN1_PYHINCNA(gbm)
f_QS_PYHIN1_PYHINCNA(gbm)
##9 hnsc
f_Table_PYHIN1_PYHINCNA(hnsc)
f_QS_PYHIN1_PYHINCNA(hnsc)
##10 kich
f_Table_PYHIN1_PYHINCNA(kich)
f_QS_PYHIN1_PYHINCNA(kich)
##11 kirc
f_Table_PYHIN1_PYHINCNA(kirc)
f_QS_PYHIN1_PYHINCNA(kirc)
##12 kirp
f_Table_PYHIN1_PYHINCNA(kirp)
f_QS_PYHIN1_PYHINCNA(kirp)
##13 lgg
f_Table_PYHIN1_PYHINCNA(kirp)
f_QS_PYHIN1_PYHINCNA(kirp)
##14 lihc
f_Table_PYHIN1_PYHINCNA(lihc)
f_QS_PYHIN1_PYHINCNA(lihc)
##15 luad
f_Table_PYHIN1_PYHINCNA(luad)
f_QS_PYHIN1_PYHINCNA(luad)
##16 lusc
f_Table_PYHIN1_PYHINCNA(lusc)
f_QS_PYHIN1_PYHINCNA(lusc)
##17 meso
f_Table_PYHIN1_PYHINCNA(meso)
f_QS_PYHIN1_PYHINCNA(meso)
##18 ov
f_Table_PYHIN1_PYHINCNA(ov)
f_QS_PYHIN1_PYHINCNA(ov)
##19 paad
f_Table_PYHIN1_PYHINCNA(paad)
f_QS_PYHIN1_PYHINCNA(paad)
##20 pcpg
f_Table_PYHIN1_PYHINCNA(pcpg)
f_QS_PYHIN1_PYHINCNA(pcpg)
##21 prad
f_Table_PYHIN1_PYHINCNA(prad)
f_QS_PYHIN1_PYHINCNA(prad)
##22 sarc
f_Table_PYHIN1_PYHINCNA(sarc)
f_QS_PYHIN1_PYHINCNA(sarc)
##23 stad
f_Table_PYHIN1_PYHINCNA(stad)
f_QS_PYHIN1_PYHINCNA(stad)
##24 tgct
f_Table_PYHIN1_PYHINCNA(tgct)
f_QS_PYHIN1_PYHINCNA(tgct)
##25 thca
f_Table_PYHIN1_PYHINCNA(thca)
f_QS_PYHIN1_PYHINCNA(thca)
##26 thca
f_Table_PYHIN1_PYHINCNA(thym)
f_QS_PYHIN1_PYHINCNA(thym)
##27 ucec
f_Table_PYHIN1_PYHINCNA(ucec)
f_QS_PYHIN1_PYHINCNA(ucec)
##28 ucs
f_Table_PYHIN1_PYHINCNA(ucs)
f_QS_PYHIN1_PYHINCNA(ucs)
##28 uvm
f_Table_PYHIN1_PYHINCNA(uvm)
f_QS_PYHIN1_PYHINCNA(uvm)

#####*****以上为错误代码*****#####



#####修正：加入“AMP”和“DEL”后重新分类进行作图，之前的错了#####
#分为4组的判定
## WT:WT(snv)-NCN。这是完全的野生型，基因功能和水平都正常。
## MUT:MUT(snv)-AMP+NCN。这是蛋白突变的（基因功能发生了改变）。
## DEL:MUT(snv)-DEL。如果一个MUT的基因是DEL（完全缺失）那么这个MUT是没有意义的，所以可以直接划分为DEL。
## AMP:WT(snv)-AMP。这是基因功能没变，但水平可能上升的。
## DEL:WT(snv)-DEL。这是基因功能没变，但水平可能下降或缺失的。

# Add a new column 将以上SNV和CNA信息整理合并为PYHIN1_aamutation_cna列
##一个要避免的报错，会导致replace后的值变成NA
### 由于在cbind()一步中将两列原始数据变成了factor类型,0不属于自动判定的level,直接替换“WT_NCN”为"WT"会报错。
### 解决方法：把数据保存为csv格式，再重新在R中用read.csv()时将参数stringAsFactors设置为FALSE。或者在建立dataframe的时候设置同样的参数！
write.table (df.merged9, file ="AAmutation+CNA.csv", sep =",", row.names =FALSE)
df.merged9 <- read.csv("AAmutation+CNA.csv", header = TRUE, stringsAsFactors = FALSE)
## NV和CNA信息整理合并为新列
df.merged9$PYHIN1_aamutation_cna <- df.merged9$Pyhin1_Pyhin1CNA #添加新列
df.merged9[which(df.merged9$PYHIN1_aamutation_cna == "WT_NCN"),]$PYHIN1_aamutation_cna <- "WT"# WT-NCN → WT
df.merged9[which(df.merged9$PYHIN1_aamutation_cna == "MUT_AMP"),]$PYHIN1_aamutation_cna <- "MUT"# MUT-AMP → MUT
df.merged9[which(df.merged9$PYHIN1_aamutation_cna == "MUT_NCN"),]$PYHIN1_aamutation_cna <- "MUT"# WUT-NCN → MUT
df.merged9[which(df.merged9$PYHIN1_aamutation_cna == "MUT_DEL"),]$PYHIN1_aamutation_cna <- "DEL"# MUT-DEL → DEL
df.merged9[which(df.merged9$PYHIN1_aamutation_cna == "WT_AMP"),]$PYHIN1_aamutation_cna <- "AMP"# WT-AMP → AMP
df.merged9[which(df.merged9$PYHIN1_aamutation_cna == "WT_DEL"),]$PYHIN1_aamutation_cna <- "DEL"# WT-DEL → DEL
table(df.merged9$PYHIN1_aamutation_cna)

#####Expression and QS levels - PYHIN1_aamutation_cna#####
##表达plot
my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged9, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged9, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3)
##QS plot
my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
aggregate(df.merged9$QuiescenceScore, by=list(type=df.merged9$PYHIN1_aamutation_cna),mean)

ggboxplot(df.merged9, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged9, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3)

#####Expression and QS levels - PYHIN1_aamutation_cna + TP53_AAmutation#####
df.merged10 <- tidyr::unite(df.merged9, "PYHIN1_aamutation_cna_TP53_aamutation", PYHIN1_aamutation_cna, TP53_status, remove = FALSE)
head(df.merged10)
table(df.merged10$PYHIN1_aamutation_cna_TP53_aamutation)
##表达plot
my_comparisons4 <- list(c("AMP_MUT","AMP_WT"),c("DEL_MUT","DEL_WT"),c("MUT_MUT","MUT_WT"),c("WT_MUT","WT_WT"))
ggboxplot(df.merged10, x = "PYHIN1_aamutation_cna_TP53_aamutation", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna_TP53_aamutation", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                      "#D5D2C1","#00AFBB"),
          add = "jitter",order = c("WT_WT","WT_MUT","MUT_WT","MUT_MUT","AMP_WT", "AMP_MUT","DEL_WT","DEL_MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons4,label = "p.signif")
##QS plot
my_comparisons4 <- list(c("AMP_MUT","AMP_WT"),c("DEL_MUT","DEL_WT"),c("MUT_MUT","MUT_WT"),c("WT_MUT","WT_WT"))
ggboxplot(df.merged10, x = "PYHIN1_aamutation_cna_TP53_aamutation", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna_TP53_aamutation", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                      "#D5D2C1","#00AFBB"),
          add = "jitter",order = c("WT_WT","WT_MUT","MUT_WT","MUT_MUT","AMP_WT", "AMP_MUT","DEL_WT","DEL_MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons4,label = "p.signif")
my_comparisons4 <- list(c("AMP_MUT","AMP_WT"),c("DEL_MUT","DEL_WT"),c("MUT_MUT","MUT_WT"),c("WT_MUT","WT_WT"))
ggboxplot(df.merged10, x = "PYHIN1_aamutation_cna_TP53_aamutation", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna_TP53_aamutation", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                      "#D5D2C1","#00AFBB"),
          add = "jitter",order = c("WT_WT","WT_MUT","MUT_WT","MUT_MUT","AMP_WT", "AMP_MUT","DEL_WT","DEL_MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons4)

#####Expression and QS levels - PYHIN1_aamutation_cna - Cancer type#####
#把每种癌症都分为一个df
#再分类绘制QS & Expression plot（boxplot 感觉violinplot没有意义）

##1 acc
acc<- df.merged10[df.merged10$STUDY_ID=="acc_tcga_pan_can_atlas_2018",]
table(acc$PYHIN1_aamutation_cna)
aggregate(acc$QuiescenceScore, by=list(type=acc$PYHIN1_aamutation_cna),mean)
ggboxplot(acc, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP"),c("WT","DEL"),c("AMP","DEL")),label = "p.signif")

ggboxplot(acc, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP"),c("WT","DEL"),c("AMP","DEL")))

##2 blca
blca<- df.merged10[df.merged10$STUDY_ID=="blca_tcga_pan_can_atlas_2018",]
table(blca$PYHIN1_aamutation_cna)
aggregate(blca$QuiescenceScore, by=list(type=blca$PYHIN1_aamutation_cna),mean)
ggboxplot(blca, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
ggboxplot(blca, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")

##3 brca
brca<- df.merged10[df.merged10$STUDY_ID=="brca_tcga_pan_can_atlas_2018",]
table(brca$PYHIN1_aamutation_cna)
aggregate(brca$QuiescenceScore, by=list(type=brca$PYHIN1_aamutation_cna),mean)
ggboxplot(brca, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
ggboxplot(brca, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
##4 cesc
cesc<- df.merged10[df.merged10$STUDY_ID=="cesc_tcga_pan_can_atlas_2018",]
table(cesc$PYHIN1_aamutation_cna)
aggregate(cesc$QuiescenceScore, by=list(type=cesc$PYHIN1_aamutation_cna),mean)
ggboxplot(cesc, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
ggboxplot(cesc, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
##5 chol
chol<- df.merged10[df.merged10$STUDY_ID=="chol_tcga_pan_can_atlas_2018",]
table(chol$PYHIN1_aamutation_cna)
aggregate(chol$QuiescenceScore, by=list(type=chol$PYHIN1_aamutation_cna),mean)
ggboxplot(chol, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP")),label = "p.signif")
ggboxplot(chol, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP")),label = "p.signif")
##6 coadread
coadread<- df.merged10[df.merged10$STUDY_ID=="coadread_tcga_pan_can_atlas_2018",]
table(coadread$PYHIN1_aamutation_cna)
aggregate(coadread$QuiescenceScore, by=list(type=coadread$PYHIN1_aamutation_cna),mean)
ggboxplot(coadread, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
ggboxplot(coadread, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
##7 esca
esca<- df.merged10[df.merged10$STUDY_ID=="esca_tcga_pan_can_atlas_2018",]
table(esca$PYHIN1_aamutation_cna)
aggregate(esca$QuiescenceScore, by=list(type=esca$PYHIN1_aamutation_cna),mean)
ggboxplot(esca, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
ggboxplot(esca, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
##8 gbm
gbm<- df.merged10[df.merged10$STUDY_ID=="gbm_tcga_pan_can_atlas_2018",]
table(gbm$PYHIN1_aamutation_cna)
aggregate(gbm$QuiescenceScore, by=list(type=gbm$PYHIN1_aamutation_cna),mean)
ggboxplot(gbm, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
ggboxplot(gbm, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
##9 hnsc
hnsc<- df.merged10[df.merged10$STUDY_ID=="hnsc_tcga_pan_can_atlas_2018",]
table(hnsc$PYHIN1_aamutation_cna)
aggregate(hnsc$QuiescenceScore, by=list(type=hnsc$PYHIN1_aamutation_cna),mean)
ggboxplot(hnsc, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
ggboxplot(hnsc, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
##10 kich
kich<- df.merged10[df.merged10$STUDY_ID=="kich_tcga_pan_can_atlas_2018",]
table(kich$PYHIN1_aamutation_cna)
aggregate(kich$QuiescenceScore, by=list(type=kich$PYHIN1_aamutation_cna),mean)
ggboxplot(kich, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP"),c("WT","DEL"),c("AMP","DEL")),label = "p.signif")
ggboxplot(kich, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP"),c("WT","DEL"),c("AMP","DEL")),label = "p.signif")
##11 kirc
kirc<- df.merged10[df.merged10$STUDY_ID=="kirc_tcga_pan_can_atlas_2018",]
table(kirc$PYHIN1_aamutation_cna)
aggregate(kirc$QuiescenceScore, by=list(type=kirc$PYHIN1_aamutation_cna),mean)
ggboxplot(kirc, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP"),c("WT","DEL"),c("AMP","DEL")),label = "p.signif")
ggboxplot(kirc, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP"),c("WT","DEL"),c("AMP","DEL")),label = "p.signif")
##12 kirp
kirp<- df.merged10[df.merged10$STUDY_ID=="kirp_tcga_pan_can_atlas_2018",]
table(kirp$PYHIN1_aamutation_cna)
aggregate(kirp$QuiescenceScore, by=list(type=kirp$PYHIN1_aamutation_cna),mean)
ggboxplot(kirp, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
ggboxplot(kirp, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
##13 lgg
lgg<- df.merged10[df.merged10$STUDY_ID=="lgg_tcga_pan_can_atlas_2018",]
table(lgg$PYHIN1_aamutation_cna)
aggregate(lgg$QuiescenceScore, by=list(type=lgg$PYHIN1_aamutation_cna),mean)
ggboxplot(lgg, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP"),c("WT","DEL"),c("AMP","DEL")),label = "p.signif")
ggboxplot(lgg, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP"),c("WT","DEL"),c("AMP","DEL")),label = "p.signif")
##14 lihc
lihc<- df.merged10[df.merged10$STUDY_ID=="lihc_tcga_pan_can_atlas_2018",]
table(lihc$PYHIN1_aamutation_cna)
aggregate(lihc$QuiescenceScore, by=list(type=lihc$PYHIN1_aamutation_cna),mean)
ggboxplot(lihc, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
ggboxplot(lihc, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
##15 luad
luad<- df.merged10[df.merged10$STUDY_ID=="luad_tcga_pan_can_atlas_2018",]
table(luad$PYHIN1_aamutation_cna)
aggregate(luad$QuiescenceScore, by=list(type=luad$PYHIN1_aamutation_cna),mean)
ggboxplot(luad, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
ggboxplot(luad, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
##16 lusc
lusc<- df.merged10[df.merged10$STUDY_ID=="lusc_tcga_pan_can_atlas_2018",]
table(lusc$PYHIN1_aamutation_cna)
aggregate(lusc$QuiescenceScore, by=list(type=lusc$PYHIN1_aamutation_cna),mean)
ggboxplot(lusc, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
ggboxplot(lusc, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
##17 meso
meso<- df.merged10[df.merged10$STUDY_ID=="meso_tcga_pan_can_atlas_2018",]
table(meso$PYHIN1_aamutation_cna)
aggregate(meso$QuiescenceScore, by=list(type=meso$PYHIN1_aamutation_cna),mean)
ggboxplot(meso, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP"),c("WT","DEL"),c("AMP","DEL")),label = "p.signif")
ggboxplot(meso, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP"),c("WT","DEL"),c("AMP","DEL")),label = "p.signif")
##18 ov
ov<- df.merged10[df.merged10$STUDY_ID=="ov_tcga_pan_can_atlas_2018",]
table(ov$PYHIN1_aamutation_cna)
aggregate(ov$QuiescenceScore, by=list(type=ov$PYHIN1_aamutation_cna),mean)
ggboxplot(ov, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
ggboxplot(ov, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
##19 paad
paad<- df.merged10[df.merged10$STUDY_ID=="paad_tcga_pan_can_atlas_2018",]
table(paad$PYHIN1_aamutation_cna)
aggregate(paad$QuiescenceScore, by=list(type=paad$PYHIN1_aamutation_cna),mean)
ggboxplot(paad, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
ggboxplot(paad, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
##20 pcpg
pcpg<- df.merged10[df.merged10$STUDY_ID=="pcpg_tcga_pan_can_atlas_2018",]
table(pcpg$PYHIN1_aamutation_cna)
aggregate(pcpg$QuiescenceScore, by=list(type=pcpg$PYHIN1_aamutation_cna),mean)
ggboxplot(pcpg, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP"),c("WT","DEL"),c("AMP","DEL")),label = "p.signif")
ggboxplot(pcpg, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP"),c("WT","DEL"),c("AMP","DEL")),label = "p.signif")
##21 prad
prad<- df.merged10[df.merged10$STUDY_ID=="prad_tcga_pan_can_atlas_2018",]
table(prad$PYHIN1_aamutation_cna)
aggregate(prad$QuiescenceScore, by=list(type=prad$PYHIN1_aamutation_cna),mean)
ggboxplot(prad, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
ggboxplot(prad, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
##22 sarc
sarc<- df.merged10[df.merged10$STUDY_ID=="sarc_tcga_pan_can_atlas_2018",]
table(sarc$PYHIN1_aamutation_cna)
aggregate(sarc$QuiescenceScore, by=list(type=sarc$PYHIN1_aamutation_cna),mean)
ggboxplot(sarc, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
ggboxplot(sarc, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
##23 stad
stad<- df.merged10[df.merged10$STUDY_ID=="stad_tcga_pan_can_atlas_2018",]
table(stad$PYHIN1_aamutation_cna)
aggregate(stad$QuiescenceScore, by=list(type=stad$PYHIN1_aamutation_cna),mean)
ggboxplot(stad, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
ggboxplot(stad, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
##24 tgct
tgct<- df.merged10[df.merged10$STUDY_ID=="tgct_tcga_pan_can_atlas_2018",]
table(tgct$PYHIN1_aamutation_cna)
aggregate(tgct$QuiescenceScore, by=list(type=tgct$PYHIN1_aamutation_cna),mean)
ggboxplot(tgct, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP"),c("WT","DEL"),c("AMP","DEL")),label = "p.signif")
ggboxplot(tgct, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP"),c("WT","DEL"),c("AMP","DEL")),label = "p.signif")
##25 thca
thca<- df.merged10[df.merged10$STUDY_ID=="thca_tcga_pan_can_atlas_2018",]
table(thca$PYHIN1_aamutation_cna)
aggregate(thca$QuiescenceScore, by=list(type=thca$PYHIN1_aamutation_cna),mean)
ggboxplot(thca, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP"),c("WT","DEL"),c("AMP","DEL")),label = "p.signif")
ggboxplot(thca, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP"),c("WT","DEL"),c("AMP","DEL")),label = "p.signif")
##26 thym
thym<- df.merged10[df.merged10$STUDY_ID=="thym_tcga_pan_can_atlas_2018",]
table(thym$PYHIN1_aamutation_cna)
aggregate(thym$QuiescenceScore, by=list(type=thym$PYHIN1_aamutation_cna),mean)
ggboxplot(thym, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP")),label = "p.signif")
ggboxplot(thym, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP")),label = "p.signif")
##27 ucec
ucec<- df.merged10[df.merged10$STUDY_ID=="ucec_tcga_pan_can_atlas_2018",]
table(ucec$PYHIN1_aamutation_cna)
aggregate(ucec$QuiescenceScore, by=list(type=ucec$PYHIN1_aamutation_cna),mean)
ggboxplot(ucec, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
ggboxplot(ucec, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
##28 ucs
ucs<- df.merged10[df.merged10$STUDY_ID=="ucs_tcga_pan_can_atlas_2018",]
table(ucs$PYHIN1_aamutation_cna)
aggregate(ucs$QuiescenceScore, by=list(type=ucs$PYHIN1_aamutation_cna),mean)
ggboxplot(ucs, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP"),c("WT","DEL"),c("AMP","DEL")),label = "p.signif")
ggboxplot(ucs, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP"),c("WT","DEL"),c("AMP","DEL")),label = "p.signif")
##29 uvm
uvm<- df.merged10[df.merged10$STUDY_ID=="uvm_tcga_pan_can_atlas_2018",]
table(uvm$PYHIN1_aamutation_cna)
aggregate(uvm$QuiescenceScore, by=list(type=uvm$PYHIN1_aamutation_cna),mean)
ggboxplot(uvm, x = "PYHIN1_aamutation_cna", y = "PYHIN1_log2expression",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP")),label = "p.signif")
ggboxplot(uvm, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","AMP")),label = "p.signif")

#####Expression and QS levels - PYHIN1_aamutation_cna_TP53_aamutation - Cancer type#####
##expression plot function
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation <- function(a){
  ggboxplot(a, x = "PYHIN1_aamutation_cna_TP53_aamutation", y = "PYHIN1_log2expression",
            color = "PYHIN1_aamutation_cna_TP53_aamutation", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                        "#D5D2C1","#00AFBB"),
            add = "jitter",order = c("WT_WT","WT_MUT","MUT_WT","MUT_MUT","AMP_WT", "AMP_MUT","DEL_WT","DEL_MUT")) + #添加jitter后，会自动隐藏异常值。
    stat_compare_means(comparisons = list(c("AMP_WT","AMP_MUT")),label = "p.signif")+
    stat_compare_means(comparisons = list(c("DEL_WT","DEL_MUT")),label = "p.signif")+
    stat_compare_means(comparisons = list(c("MUT_WT","MUT_MUT")),label = "p.signif")+
    stat_compare_means(comparisons = list(c("WT_WT","WT_MUT")),label = "p.signif")
}

f_QS_PYHIN1_aamutation_cna_TP53_aamutation <- function(a){
  ggboxplot(a, x = "PYHIN1_aamutation_cna_TP53_aamutation", y = "QuiescenceScore",
            color = "PYHIN1_aamutation_cna_TP53_aamutation", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                        "#D5D2C1","#00AFBB"),
            add = "jitter",order = c("WT_WT","WT_MUT","MUT_WT","MUT_MUT","AMP_WT", "AMP_MUT","DEL_WT","DEL_MUT")) + #添加jitter后，会自动隐藏异常值。
    stat_compare_means(comparisons = list(c("AMP_WT","AMP_MUT")),label = "p.signif")+
    stat_compare_means(comparisons = list(c("DEL_WT","DEL_MUT")),label = "p.signif")+
    stat_compare_means(comparisons = list(c("MUT_WT","MUT_MUT")),label = "p.signif")+
    stat_compare_means(comparisons = list(c("WT_WT","WT_MUT")),label = "p.signif")
}
##1 acc
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(acc)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(acc)
##2 blca
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(blca)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(blca)
##3 brca
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(brca)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(brca)
##4 cesc
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(cesc)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(cesc)
##5 chol
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(chol)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(chol)
##6 coadread
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(coadread)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(coadread)
##7 esca
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(esca)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(esca)
##8 gbm
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(gbm)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(gbm)
##9 hnsc
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(hnsc)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(hnsc)
##10 kich
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(kich)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(kich)
##11 kirc
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(kirc)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(kirc)
##12 kirp
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(kirp)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(kirp)
##13 lgg
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(lgg)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(lgg)
##14 lihc
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(lihc)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(lihc)
##15 luad
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(luad)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(luad)
##16 lusc
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(lusc)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(lusc)
##17 meso
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(meso)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(meso)
##18 ov
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(ov)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(ov)
##19 paad
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(paad)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(paad)
##20 pcpg
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(pcpg)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(pcpg)
##21 prad
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(prad)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(prad)
##22 sarc
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(sarc)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(sarc)
##23 stad
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(stad)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(stad)
##24 tgct
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(tgct)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(tgct)
##25 thca
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(thca)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(thca)
##26 thym
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(thym)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(thym)
##27 ucec
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(ucec)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(ucec)
##28 ucs
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(ucs)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(ucs)
##29 uvm
f_QS_PYHIN1_aamutation_cna_TP53_aamutation(uvm)
f_Ex_PYHIN1_aamutation_cna_TP53_aamutation(uvm)

#####导入整倍体数据AS#####
aneuploidy <- read.csv("TCGA_aneuploidy.csv")
head(aneuploidy)
hist(aneuploidy$AneuploidyScore)

df.merged11 <- merge(df.merged10[,c("SAMPLE_ID","PYHIN1_CASP8_status","CASP8_status",
                                  "CASP8_log2expression","PYHIN1_TP53_status",
                                  "TP53_status","PYHIN1_status","PYHIN1_log2expression",
                                  "QuiescenceScore","Group",
                                  "Pyhin1_Pyhin1CNA",
                                  "STUDY_ID",
                                  "PYHIN1_CNA",
                                  "PYHIN1_aamutation_cna",
                                  "PYHIN1_aamutation_cna_TP53_aamutation")], 
                    aneuploidy[,c("Sample","Type","AneuploidyScore",
                                  "AS_del","AS_amp",
                                  "Genome_doublings","Leuk","Purity","Stroma",
                                  "Stroma_notLeukocyte","Stroma_notLeukocyte_Floor",
                                  "SilentMutationspeMb",
                                  "Non.silentMutationsperMb")],
                    by.x="SAMPLE_ID", by.y="Sample", #x和y表合并的依据是相同的sample id
                    all.x=FALSE, all.y=FALSE) 
head(df.merged11)

##### AS - PYHIN1_aamutation #####
table(df.merged11$PYHIN1_status)
ggboxplot(df.merged11, x = "PYHIN1_status", y = "AneuploidyScore",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### AS_amp - PYHIN1_aamutation #####
table(df.merged11$PYHIN1_status)
ggboxplot(df.merged11, x = "PYHIN1_status", y = "AS_amp",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### AS_del - PYHIN1_aamutation #####
table(df.merged11$PYHIN1_status)
ggboxplot(df.merged11, x = "PYHIN1_status", y = "AS_del",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Gd - PYHIN1_aamutation #####
table(df.merged11$PYHIN1_status)
ggboxplot(df.merged11, x = "PYHIN1_status", y = "Genome_doublings",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Leuk - PYHIN1_aamutation #####
table(df.merged11$PYHIN1_status)
df.merged11$Leuk[which(df.merged11$Leuk =='#N/A')] <- 'NA'
df.merged12 <- df.merged11[!grepl("NA", df.merged11$Leuk),]
class(df.merged12$Leuk)
df.merged12$Leuk <- as.numeric(df.merged12$Leuk)
class(df.merged12$Leuk)
hist(df.merged12$Leuk)

ggboxplot(df.merged12, x = "PYHIN1_status", y = "Leuk",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Purity - PYHIN1_aamutation #####
table(df.merged11$PYHIN1_status)
df.merged11$Purity[which(df.merged11$Purity =='#N/A')] <- 'NA'
class(df.merged11$Purity)
df.merged11$Purity <- as.numeric(df.merged11$Purity)
class(df.merged11$Purity)
hist(df.merged11$Purity)

ggboxplot(df.merged11, x = "PYHIN1_status", y = "Purity",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Stroma - PYHIN1_aamutation #####
df.merged11$Stroma[which(df.merged11$Stroma =='#N/A')] <- 'NA'
class(df.merged11$Stroma)
df.merged11$Stroma <- as.numeric(df.merged11$Stroma)
class(df.merged11$Stroma)
hist(df.merged11$Stroma)

ggboxplot(df.merged11, x = "PYHIN1_status", y = "Stroma",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Stroma_notLeukocyte - PYHIN1_aamutation #####
df.merged11$Stroma_notLeukocyte[which(df.merged11$Stroma_notLeukocyte =='#N/A')] <- 'NA'
class(df.merged11$Stroma_notLeukocyte)
df.merged11$Stroma_notLeukocyte <- as.numeric(df.merged11$Stroma_notLeukocyte)
class(df.merged11$Stroma_notLeukocyte)
hist(df.merged11$Stroma_notLeukocyte)

ggboxplot(df.merged11, x = "PYHIN1_status", y = "Stroma_notLeukocyte",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### SilentMutationspeMb - PYHIN1_aamutation #####
df.merged11$SilentMutationspeMb[which(df.merged11$SilentMutationspeMb =='#N/A')] <- 'NA'
class(df.merged11$SilentMutationspeMb)
df.merged11$SilentMutationspeMb <- as.numeric(df.merged11$SilentMutationspeMb)
class(df.merged11$SilentMutationspeMb)
hist(df.merged11$SilentMutationspeMb)

ggboxplot(df.merged11, x = "PYHIN1_status", y = "SilentMutationspeMb",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Non.silentMutationsperMb - PYHIN1_aamutation #####
df.merged11$Non.silentMutationsperMb[which(df.merged11$Non.silentMutationsperMb =='#N/A')] <- 'NA'
class(df.merged11$Non.silentMutationsperMb)
df.merged11$Non.silentMutationsperMb <- as.numeric(df.merged11$Non.silentMutationsperMb)
class(df.merged11$Non.silentMutationsperMb)
hist(df.merged11$Non.silentMutationsperMb)

ggboxplot(df.merged11, x = "PYHIN1_status", y = "Non.silentMutationsperMb",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")


##### QS - AS #####

# 导入所需包和数据
options(scipen=999)  # 关闭科学记数法，如 1e+48
library(ggplot2)
theme_set(theme_bw())  # 设置 theme_bw()为默认主题.
data("midwest", package = "ggplot2") # 导入数据
# midwest <- read.csv("http://goo.gl/G1K41K")  # 数据源

# Scatterplot
gg1 <- ggplot(df.merged11, aes(x=AneuploidyScore, y=QuiescenceScore)) +
  geom_point(aes(col=Type)) +
  geom_smooth(method="loess", se=F) + 
  xlim(c(0, 40)) + 
  ylim(c(-30, 30)) + 
  labs(subtitle="QS Vs AS", 
       title="Scatterplot")+
  stat_cor(method = 'pearson', aes(x = AneuploidyScore, y = QuiescenceScore))
plot(gg1)

gg2 <- ggplot(df.merged11, aes(x=AneuploidyScore, y=QuiescenceScore)) +
  geom_point(aes(col=Type)) +
  geom_smooth(method="lm", se=F) + 
  xlim(c(0, 40)) + 
  ylim(c(-30, 30)) + 
  labs(subtitle="QS Vs AS", 
       title="Scatterplot")+
  stat_cor(method = 'pearson', aes(x = AneuploidyScore, y = QuiescenceScore))
plot(gg2)

gg3 <- ggplot(df.merged11, aes(x=AS_del, y=QuiescenceScore)) +
  geom_point(aes(col=Type)) +
  geom_smooth(method="lm", se=F) + 
  xlim(c(0, 40)) + 
  ylim(c(-30, 30)) + 
  labs(subtitle="QS Vs AS_del", 
       title="Scatterplot")+
  stat_cor(method = 'pearson', aes(x = AS_del, y = QuiescenceScore))
plot(gg3)

gg4 <- ggplot(df.merged11, aes(x=AS_amp, y=QuiescenceScore)) +
  geom_point(aes(col=Type)) +
  geom_smooth(method="lm", se=F) + 
  xlim(c(0, 40)) + 
  ylim(c(-30, 30)) + 
  labs(subtitle="QS Vs AS_amp", 
       title="Scatterplot")+
  stat_cor(method = 'pearson', aes(x = AS_amp, y = QuiescenceScore))
plot(gg4)

ggplot(data = df.merged11, aes(x = AneuploidyScore, y = QuiescenceScore,  color = PYHIN1_aamutation_cna)) + 
  geom_point() + geom_smooth(method = lm) +
  scale_color_manual(values = c('#FF7400', '#009999', '#3914AF',"#000000"))  +
  labs(title = 'QS Vs AS - PYHIN1_aamutation_cna') + 
  stat_cor(method = 'pearson', aes(x = AneuploidyScore, y = QuiescenceScore,  color = PYHIN1_aamutation_cna))

#####AS- PYHIN1_aamutation_cna#####
#mean of groups
aggregate(df.merged11$AneuploidyScore, by=list(type=df.merged11$PYHIN1_aamutation_cna),mean)

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged11, x = "PYHIN1_aamutation_cna", y = "AneuploidyScore",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3)
#####AS_amp- PYHIN1_aamutation_cna#####
#mean of groups
aggregate(df.merged11$AS_amp, by=list(type=df.merged11$PYHIN1_aamutation_cna),mean)

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged11, x = "PYHIN1_aamutation_cna", y = "AS_amp",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3)
#####AS_del- PYHIN1_aamutation_cna#####
#mean of groups
aggregate(df.merged11$AS_del, by=list(type=df.merged11$PYHIN1_aamutation_cna),mean)

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged11, x = "PYHIN1_aamutation_cna", y = "AS_del",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####Gd- PYHIN1_aamutation_cna#####
#mean of groups
aggregate(df.merged11$Genome_doublings, by=list(type=df.merged11$PYHIN1_aamutation_cna),mean)

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged11, x = "PYHIN1_aamutation_cna", y = "Genome_doublings",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####Leuk- PYHIN1_aamutation_cna#####
#mean of groups
aggregate(df.merged12$Leuk, by=list(type=df.merged12$PYHIN1_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged12, x = "PYHIN1_aamutation_cna", y = "Leuk",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####Purity- PYHIN1_aamutation_cna#####
#mean of groups
aggregate(df.merged11$Purity, by=list(type=df.merged11$PYHIN1_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged11, x = "PYHIN1_aamutation_cna", y = "Purity",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####Stroma- PYHIN1_aamutation_cna#####
#mean of groups
aggregate(df.merged11$Stroma, by=list(type=df.merged11$PYHIN1_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged11, x = "PYHIN1_aamutation_cna", y = "Stroma",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####Stroma_notLeukocyte- PYHIN1_aamutation_cna#####
#mean of groups
aggregate(df.merged11$Stroma_notLeukocyte, by=list(type=df.merged11$PYHIN1_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged11, x = "PYHIN1_aamutation_cna", y = "Stroma",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####SilentMutationspeMb- PYHIN1_aamutation_cna#####
my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged11, x = "PYHIN1_aamutation_cna", y = "SilentMutationspeMb",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####Non.silentMutationsperMb- PYHIN1_aamutation_cna#####
my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged11, x = "PYHIN1_aamutation_cna", y = "Non.silentMutationsperMb",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")

#####导入突变总数数据NM#####
pgen <- read.csv("pgen.csv")
head(pgen)
hist(pgen$number_of_mutations)

df.merged13 <- merge(df.merged11[,c("SAMPLE_ID","PYHIN1_CASP8_status","CASP8_status",
                                    "CASP8_log2expression","PYHIN1_TP53_status",
                                    "TP53_status","PYHIN1_status","PYHIN1_log2expression",
                                    "QuiescenceScore","Group",
                                    "Pyhin1_Pyhin1CNA",
                                    "STUDY_ID",
                                    "PYHIN1_CNA",
                                    "PYHIN1_aamutation_cna",
                                    "PYHIN1_aamutation_cna_TP53_aamutation",
                                    "Type","AneuploidyScore",
                                    "AS_del","AS_amp",
                                    "Genome_doublings","Leuk","Purity","Stroma",
                                    "Stroma_notLeukocyte","Stroma_notLeukocyte_Floor",
                                    "SilentMutationspeMb",
                                    "Non.silentMutationsperMb")], 
                     pgen[,c("sample_name","tumor_type",
                             "number_of_mutations")],
                     by.x="SAMPLE_ID", by.y="sample_name", #x和y表合并的依据是相同的sample id
                     all.x=FALSE, all.y=FALSE) 
head(df.merged13)
#####Number of Mutations- PYHIN1_aamutation_cna#####
#mean of groups
aggregate(df.merged13$number_of_mutations, by=list(type=df.merged13$PYHIN1_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged13, x = "PYHIN1_aamutation_cna", y = "number_of_mutations",
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3)
#####Signature 导入癌症的标签数据#####
#Combine 合并所有癌症类型的标签数据#
sigs_ACC <- sigs.defaultnounknown
sigs_BLCA <- sigs.defaultnounknown
sigs_BRCA <- sigs.defaultnounknown
sigs_CESC <- sigs.defaultnounknown
sigs_CHOL <- sigs.defaultnounknown
sigs_COAD <- sigs.defaultnounknown
sigs_ESCA <- sigs.defaultnounknown
sigs_GBM <- sigs.defaultnounknown
sigs_HNSC <- sigs.defaultnounknown
sigs_KICH <- sigs.defaultnounknown
sigs_KIRC <- sigs.defaultnounknown
sigs_KIRP <- sigs.defaultnounknown
sigs_LGG <- sigs.defaultnounknown
sigs_LIHC <- sigs.defaultnounknown
sigs_LUAD <- sigs.defaultnounknown
sigs_LUSC <- sigs.defaultnounknown
sigs_MESO <- sigs.defaultnounknown
sigs_OV <- sigs.defaultnounknown
sigs_PAAD <- sigs.defaultnounknown
sigs_PRAD <- sigs.defaultnounknown
sigs_READ <- sigs.defaultnounknown
sigs_SARC <- sigs.defaultnounknown
sigs_SKCM <- sigs.defaultnounknown
sigs_STAD <- sigs.defaultnounknown
sigs_TGCT <- sigs.defaultnounknown
sigs_THCA <- sigs.defaultnounknown
sigs_THYM <- sigs.defaultnounknown
sigs_UCEC <- sigs.defaultnounknown
sigs_UCS <- sigs.defaultnounknown

#new column，sample id 将行名id转化为一列#
library(data.table)
setDT(sigs_ACC, keep.rownames = "sig_sample_id")
setDT(sigs_BLCA, keep.rownames = "sig_sample_id")
setDT(sigs_BRCA, keep.rownames = "sig_sample_id")
setDT(sigs_CESC, keep.rownames = "sig_sample_id")
setDT(sigs_CHOL, keep.rownames = "sig_sample_id")
setDT(sigs_COAD, keep.rownames = "sig_sample_id")
setDT(sigs_ESCA, keep.rownames = "sig_sample_id")
setDT(sigs_GBM, keep.rownames = "sig_sample_id")
setDT(sigs_HNSC, keep.rownames = "sig_sample_id")
setDT(sigs_KICH, keep.rownames = "sig_sample_id")
setDT(sigs_KIRC, keep.rownames = "sig_sample_id")
setDT(sigs_KIRP, keep.rownames = "sig_sample_id")
setDT(sigs_LGG, keep.rownames = "sig_sample_id")
setDT(sigs_LIHC, keep.rownames = "sig_sample_id")
setDT(sigs_LUAD, keep.rownames = "sig_sample_id")
setDT(sigs_LUSC, keep.rownames = "sig_sample_id")
setDT(sigs_MESO, keep.rownames = "sig_sample_id")
setDT(sigs_OV, keep.rownames = "sig_sample_id")
setDT(sigs_PAAD, keep.rownames = "sig_sample_id")
setDT(sigs_PRAD, keep.rownames = "sig_sample_id")
setDT(sigs_READ, keep.rownames = "sig_sample_id")
setDT(sigs_SARC, keep.rownames = "sig_sample_id")
setDT(sigs_SKCM, keep.rownames = "sig_sample_id")
setDT(sigs_STAD, keep.rownames = "sig_sample_id")
setDT(sigs_TGCT, keep.rownames = "sig_sample_id")
setDT(sigs_THCA, keep.rownames = "sig_sample_id")
setDT(sigs_THYM, keep.rownames = "sig_sample_id")
setDT(sigs_UCEC, keep.rownames = "sig_sample_id")
setDT(sigs_UCS, keep.rownames = "sig_sample_id")

#Complete signature合并为完整的标签数据集#
library(plyr)
sigs_all <- rbind.fill(sigs_ACC,sigs_BLCA,sigs_BRCA,sigs_CESC,sigs_CHOL,sigs_COAD,
                       sigs_ESCA,sigs_GBM,sigs_HNSC,sigs_KICH,sigs_KIRC,sigs_KIRP,
                       sigs_LGG,sigs_LIHC,sigs_LUAD,sigs_LUSC,sigs_MESO,sigs_OV,
                       sigs_PAAD,sigs_PRAD,sigs_READ,sigs_SARC,sigs_SKCM,sigs_STAD,
                       sigs_TGCT,sigs_THCA,sigs_THYM,sigs_UCEC,sigs_UCS)
#将标签数据集中sample id多余部分删除，统一id格式.只保留前15个字符#
library("stringr")
rownames(sigs_all)
sigs_all_new<-sigs_all
sigs_all_new$sig_sample_id<-str_sub(sigs_all$sig_sample_id,end=15)
write.csv(sigs_all_new,file="sigs_all_new")
#将整理好的标签数据 与样本之前的数据 合并#
df.merged14 <- merge(df.merged13[,c("SAMPLE_ID","PYHIN1_CASP8_status","CASP8_status",
                                    "CASP8_log2expression","PYHIN1_TP53_status",
                                    "TP53_status","PYHIN1_status","PYHIN1_log2expression",
                                    "QuiescenceScore","Group",
                                    "Pyhin1_Pyhin1CNA",
                                    "STUDY_ID",
                                    "PYHIN1_CNA",
                                    "PYHIN1_aamutation_cna",
                                    "PYHIN1_aamutation_cna_TP53_aamutation",
                                    "Type","AneuploidyScore",
                                    "AS_del","AS_amp",
                                    "Genome_doublings","Leuk","Purity","Stroma",
                                    "Stroma_notLeukocyte","Stroma_notLeukocyte_Floor",
                                    "SilentMutationspeMb",
                                    "Non.silentMutationsperMb",
                                    "tumor_type",
                                    "number_of_mutations")], 
                     sigs_all_new[,c("sig_sample_id","SBS1","SBS2","SBS3",
                                "SBS4","SBS5","SBS6","SBS7a",
                                "SBS7b","SBS7c","SBS7d","SBS8",
                                "SBS9","SBS10a","SBS10b","SBS11",
                                "SBS12","SBS13","SBS14","SBS15",
                                "SBS16","SBS17a","SBS17b","SBS18",
                                "SBS19","SBS20","SBS21","SBS22",
                                "SBS23","SBS24","SBS25","SBS26",
                                "SBS27","SBS28","SBS29","SBS30",
                                "SBS31","SBS32","SBS33","SBS34",
                                "SBS35","SBS36","SBS37","SBS38",
                                "SBS39","SBS40","SBS41","SBS42",
                                "SBS43","SBS44","SBS45","SBS46",
                                "SBS47","SBS48","SBS49","SBS50",
                                "SBS51","SBS52","SBS53","SBS54",
                                "SBS55","SBS56","SBS57","SBS58",
                                "SBS59","SBS60")],
                     by.x="SAMPLE_ID", by.y="sig_sample_id", #x和y表合并的依据是相同的sample id
                     all.x=FALSE, all.y=FALSE) 


#####Signatures- PYHIN1_aamutation_cna 完整案例分析#####
#mean of groups
SBS_y = list(df.merged14$SBS1,df.merged14$SBS2,df.merged14$SBS3,
             df.merged14$SBS4,df.merged14$SBS5,df.merged14$SBS6,df.merged14$SBS7a,
             df.merged14$SBS7b,df.merged14$SBS7c,df.merged14$SBS7d,df.merged14$SBS8,
             df.merged14$SBS9,df.merged14$SBS10a,df.merged14$SBS10b,df.merged14$SBS11,
             df.merged14$SBS12,df.merged14$SBS13,df.merged14$SBS14,df.merged14$SBS15,
             df.merged14$SBS16,df.merged14$SBS17a,df.merged14$SBS17b,df.merged14$SBS18,
             df.merged14$SBS19,df.merged14$SBS20,df.merged14$SBS21,df.merged14$SBS22,
             df.merged14$SBS23,df.merged14$SBS24,df.merged14$SBS25,df.merged14$SBS26,
             df.merged14$SBS27,df.merged14$SBS28,df.merged14$SBS29,df.merged14$SBS30,
             df.merged14$SBS31,df.merged14$SBS32,df.merged14$SBS33,df.merged14$SBS34,
             df.merged14$SBS35,df.merged14$SBS36,df.merged14$SBS37,df.merged14$SBS38,
             df.merged14$SBS39,df.merged14$SBS40,df.merged14$SBS41,df.merged14$SBS42,
             df.merged14$SBS43,df.merged14$SBS44,df.merged14$SBS45,df.merged14$SBS46,
             df.merged14$SBS47,df.merged14$SBS48,df.merged14$SBS49,df.merged14$SBS50,
             df.merged14$SBS51,df.merged14$SBS52,df.merged14$SBS53,df.merged14$SBS54,
             df.merged14$SBS55,df.merged14$SBS56,df.merged14$SBS57,df.merged14$SBS58,
             df.merged14$SBS59,df.merged14$SBS60)
SBS_y_name<-list("SBS1","SBS2","SBS3",
"SBS4","SBS5","SBS6","SBS7a",
"SBS7b","SBS7c","SBS7d","SBS8",
"SBS9","SBS10a","SBS10b","SBS11",
"SBS12","SBS13","SBS14","SBS15",
"SBS16","SBS17a","SBS17b","SBS18",
"SBS19","SBS20","SBS21","SBS22",
"SBS23","SBS24","SBS25","SBS26",
"SBS27","SBS28","SBS29","SBS30",
"SBS31","SBS32","SBS33","SBS34",
"SBS35","SBS36","SBS37","SBS38",
"SBS39","SBS40","SBS41","SBS42",
"SBS43","SBS44","SBS45","SBS46",
"SBS47","SBS48","SBS49","SBS50",
"SBS51","SBS52","SBS53","SBS54",
"SBS55","SBS56","SBS57","SBS58",
"SBS59","SBS60"
)

aggregate(df.merged14$SBS1,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS1",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS2,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS2",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS3,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS3",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS4,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS4",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS5,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS5",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS6,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS6",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS7a,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS7a",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS7b,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS7b",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS7c,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS7c",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS7d,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS7d",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS8,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS8",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS9,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS9",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS10a,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS10a",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS10b,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS10b",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS11,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS11",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS12,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS12",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS13,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS13",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS14,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS14",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS15,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS15",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS16,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS16",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS17a,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS17a",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS17b,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS17b",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS18,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS18",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS19,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS19",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS20,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS20",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS21,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS21",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS22,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS22",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS23,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS23",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS24,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS24",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS25,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS25",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS26,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS26",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS27,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS27",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS28,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS28",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS29,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS29",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS30,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS30",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS31,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS31",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS32,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS32",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS33,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS33",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS34,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS34",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS35,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS35",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS36,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS36",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS37,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS37",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS38,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS38",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS39,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS39",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS40,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS40",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS41,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS41",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS42,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS42",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS43,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS43",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS44,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS44",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS45,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS45",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS46,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS46",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS47,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS47",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS48,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS48",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS49,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS49",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS50,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS50",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS51,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS51",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS52,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS52",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS53,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS53",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS54,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS54",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS55,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS55",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS56,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS56",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS57,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS57",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS58,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS58",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS59,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS59",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14$SBS60,by=list(type=df.merged14$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14, x = "PYHIN1_status", y = "SBS60",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)
#####Ballon plot-Signature_aa#####
install.packages("reshape2") 
library(reshape2)
library(scales)

sig_ballon_aa <- read.csv("sig_ballon_data.csv",row.names = 1)
  
ggballoonplot(sig_ballon_aa, 
                fill = "value", #气球填充颜色
                ggtheme = theme_bw(),#画板主题
                size = "value",#气球大小
                size.range = c(1,7),
                color = "grey",#气球边框颜色
                shape = 23,#shape可以改变显示形状
                show.label = F)+#是否显示标签
    scale_fill_viridis_c(option = "C")+
    guides(size = FALSE)+#气球图例是否显示
    #设置颜色
    scale_fill_gradientn(
      colors=c("blue","white","red"),
      values=rescale(c(-0.04,0,0.13)),
      limits=c(-0.04,0.13))+
    
    labs(fill="Difference
(MUT-WT)")
  
