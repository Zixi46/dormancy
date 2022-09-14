### This script analyses FGFR2 mutation data and links with quiescence and expressin levels in tumours.
##### 初步数据处理 aa mutaion status与Exp QS合并#####
### Read mutation data downloaded from cbioportal:
setwd("/Users/lizixi/Desktop/Final Project/5.18 task")
mutations_fgfr2 <- read.delim("mutations_cbioportal_fgfr2.txt")
head(mutations_fgfr2)
## How many samples have a mutation?
table(mutations_fgfr2$FGFR2)
tb_fgfr2 <- table(mutations_fgfr2$FGFR2)
sort(tb_fgfr2) #descend升序排列PHIN1列的数字元素

# Add a new column that just denotes WT and MUT for FGFR2:
mutations_fgfr2[which(mutations_fgfr2$FGFR2 == "NS"),]$FGFR2 <- "WT" #****将FGFR2中的NS替换为WT
mutations_fgfr2$FGFR2_status <- mutations_fgfr2$FGFR2 #为mutations_fgfr2表添加一列为FGFR2_status，同时将这一列定义为FGFR2（和FGFR2一样）
mutations_fgfr2[which(mutations_fgfr2$FGFR2 != "WT"),]$FGFR2_status <- "MUT" #选出mutation中FGFR2列中不等于“WT”的元素，并在FGFR2_status中将其全赋值为“MUT”
head(mutations_fgfr2)
# How many samples with mutations_fgfr2 in FGFR2?
table(mutations_fgfr2$FGFR2_status)
# Read expression data for FGFR2:
expression_fgfr2 <- read.delim("mRNA Expression_RSEM_cbioportal_fgfr2.txt")
head(expression_fgfr2)
# Plot histogram of expression values of FGFR2:
hist(expression_fgfr2$FGFR2)
# It's good to take the log2 of expression values to avoid working with very large numbers:
expression_fgfr2$FGFR2_log2expression <- log2(expression_fgfr2$FGFR2+1)
head(expression_fgfr2)
hist(expression_fgfr2$FGFR2_log2expression) # more decent value

# Merge FGFR2 expression and mutation data:
df.merged_fgfr2 <- merge(mutations_fgfr2[,c("SAMPLE_ID","FGFR2_status","STUDY_ID")], #mutation中需要合并的纵列是sample id，fgfr2__status
                        expression_fgfr2[,c("SAMPLE_ID","FGFR2_log2expression")], #expresssion中需要合并的纵列是sample ID，FGFR2_log2expression
                        by.x="SAMPLE_ID", by.y="SAMPLE_ID", #x和y表合并的依据是相同的sample id
                        all.x=FALSE, all.y=FALSE) #x和y表中不需要所有的行都合并（只合并sample id一样的，只在其中一个表中才存在的sample id丢弃）
head(df.merged_fgfr2)
# Merge TP53 mutation data+QS
df.merged_fgfr2 <- merge(df.merged_fgfr2[,c("SAMPLE_ID","FGFR2_status","STUDY_ID","FGFR2_log2expression")], #mutation中需要合并的纵列是sample id，fgfr2_status
                        df.merged2[,c("SAMPLE_ID","TP53_status","QuiescenceScore","Group")], #expresssion中需要合并的纵列是sample ID，FGFR2_log2expression
                        by.x="SAMPLE_ID", by.y="SAMPLE_ID", #x和y表合并的依据是相同的sample id
                        all.x=FALSE, all.y=FALSE) #x和y表中不需要所有的行都合并（只合并sample id一样的，只在其中一个表中才存在的sample id丢弃）
head(df.merged_fgfr2)
###### 增加FGFR2突变类型分析列 “AMP”“DEL” #####
#读取cna文件
cna_FGFR2 <- read.delim("cna_fromcBioportal_fgfr2.txt")
names(cna_FGFR2)[3] <- "FGFR2_CNA"#给第三列重命名
head(cna_FGFR2)
table(cna_FGFR2$FGFR2_CNA)
#将cna文件中的ID 和 CNA与总表合并
head(df.merged_fgfr2)
df.merged_fgfr2_2 <- merge(df.merged_fgfr2[,c("SAMPLE_ID","STUDY_ID","FGFR2_status",
                                           "TP53_status","FGFR2_log2expression",
                                           "QuiescenceScore","Group")], 
                         cna_FGFR2[,c("SAMPLE_ID","FGFR2_CNA")],
                         by.x="SAMPLE_ID", by.y="SAMPLE_ID", #x和y表合并的依据是相同的sample id
                         all.x=FALSE, all.y=FALSE) 
head(df.merged_fgfr2_2)
table(df.merged_fgfr2_2$FGFR2_CNA)
#数据整理，完整案例分析
##exclude cns-NP(not profiled) samples from the analysis 剔除包含NP的样本
##把NP都换为NA
df.merged_fgfr2_2[which(df.merged_fgfr2_2$FGFR2_CNA == "NP"),]$FGFR2_CNA <- "NA" #****将FGFR2_CAN中的NP替换为NA
table(df.merged_fgfr2_2$FGFR2_CNA)
##把-2换为-1（数量太少了），把2换为1
#"-2" is a deep loss, possibly a homozygous deletion
#"-1" is a single-copy loss (heterozygous deletion)
#"0" is diploid
#"1" indicates a low-level gain
#"2" is a high-level amplification.
##将1转化为AMP,-1转化为DEL，O转化为NCN
df.merged_fgfr2_2[which(df.merged_fgfr2_2$FGFR2_CNA == "-2"),]$FGFR2_CNA <- "-1" 
df.merged_fgfr2_2[which(df.merged_fgfr2_2$FGFR2_CNA == "2"),]$FGFR2_CNA <- "1" 
df.merged_fgfr2_2[which(df.merged_fgfr2_2$FGFR2_CNA == "-1"),]$FGFR2_CNA <- "DEL" 
df.merged_fgfr2_2[which(df.merged_fgfr2_2$FGFR2_CNA == "1"),]$FGFR2_CNA <- "AMP" 
df.merged_fgfr2_2[which(df.merged_fgfr2_2$FGFR2_CNA == "0"),]$FGFR2_CNA <- "NCN" 
table(df.merged_fgfr2_2$FGFR2_CNA)
##去掉NA样本
df.merged_fgfr2_2 <- df.merged_fgfr2_2[!grepl("NA", df.merged_fgfr2_2$FGFR2_CNA),]
table(df.merged_fgfr2_2$FGFR2_CNA)
##将FGFR2 AAmutation信息和CNA信息结合为新的一列Mutation_CNA,生成6个不同组别
df.merged_fgfr2_2 <- tidyr::unite(df.merged_fgfr2_2, "FGFR2_FGFR2CNA", FGFR2_status, FGFR2_CNA, remove = FALSE)
head(df.merged_fgfr2_2)
table(df.merged_fgfr2_2$FGFR2_FGFR2CNA)
#重排顺序
df.merged_fgfr2_2 <- df.merged_fgfr2_2[,c("SAMPLE_ID","STUDY_ID","FGFR2_status","FGFR2_CNA","FGFR2_FGFR2CNA","TP53_status","FGFR2_log2expression","QuiescenceScore","Group")]
table(df.merged_fgfr2$FGFR2_status)
table(df.merged_fgfr2$TP53_status)
#根据aamutaiton+CNA进行重新分组
df.merged_fgfr2_2$FGFR2_aamutation_cna <- df.merged_fgfr2_2$FGFR2_FGFR2CNA #添加新列
df.merged_fgfr2_2[which(df.merged_fgfr2_2$FGFR2_aamutation_cna == "WT_NCN"),]$FGFR2_aamutation_cna <- "WT"# WT-NCN → WT
df.merged_fgfr2_2[which(df.merged_fgfr2_2$FGFR2_aamutation_cna == "MUT_AMP"),]$FGFR2_aamutation_cna <- "MUT"# MUT-AMP → MUT
df.merged_fgfr2_2[which(df.merged_fgfr2_2$FGFR2_aamutation_cna == "MUT_NCN"),]$FGFR2_aamutation_cna <- "MUT"# WUT-NCN → MUT
df.merged_fgfr2_2[which(df.merged_fgfr2_2$FGFR2_aamutation_cna == "MUT_DEL"),]$FGFR2_aamutation_cna <- "DEL"# MUT-DEL → DEL
df.merged_fgfr2_2[which(df.merged_fgfr2_2$FGFR2_aamutation_cna == "WT_AMP"),]$FGFR2_aamutation_cna <- "AMP"# WT-AMP → AMP
df.merged_fgfr2_2[which(df.merged_fgfr2_2$FGFR2_aamutation_cna == "WT_DEL"),]$FGFR2_aamutation_cna <- "DEL"# WT-DEL → DEL
table(df.merged_fgfr2_2$FGFR2_aamutation_cna)
#####合并 FGFR2+TP53 status and FGFR2+TP53 status #####
#Combine FGFR2+TP53 to create 4 groups 。用tidyr::unite合并两列（FGFR2_status，TP53_status）的元素为新的一列（FGFR2_TP53_status），默认以_连接（可以修改）
library(tidyr)
df.merged_fgfr2_2 <- tidyr::unite(df.merged_fgfr2_2, "FGFR2_TP53_status", FGFR2_status, TP53_status, remove = FALSE)
head(df.merged_fgfr2_2)
table(df.merged_fgfr2_2$FGFR2_TP53_status)
#####合并 FGFR2_aamutation_cna+TP53_status#####
df.merged_fgfr2_2 <- tidyr::unite(df.merged_fgfr2_2, "FGFR2_aamutation_cna_TP53_aamutation", FGFR2_aamutation_cna, TP53_status, remove = FALSE)
table(df.merged_fgfr2_2$FGFR2_aamutation_cna_TP53_aamutation)

#####Expression and QS levels - FGFR2_aamutation#####
#QS
## by boxplot
ggboxplot(df.merged_fgfr2, x = "FGFR2_status", y = "QuiescenceScore",
          color = "FGFR2_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")
ggboxplot(df.merged_fgfr2, x = "FGFR2_status", y = "QuiescenceScore",
          color = "FGFR2_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)
#Exp
## by boxplot
ggboxplot(df.merged_fgfr2, x = "FGFR2_status", y = "FGFR2_log2expression",
          color = "FGFR2_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")
ggboxplot(df.merged_fgfr2, x = "FGFR2_status", y = "FGFR2_log2expression",
          color = "FGFR2_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)

#####Expression and QS levels - FGFR2_aamutation_cna#####
#QS
## by boxplot
ggboxplot(df.merged_fgfr2_2, x = "FGFR2_aamutation_cna", y = "QuiescenceScore",
          color = "FGFR2_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons3,label="p.signif")
ggboxplot(df.merged_fgfr2_2, x = "FGFR2_aamutation_cna", y = "QuiescenceScore",
          color = "FGFR2_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons3)
#Exp
## by boxplot
ggboxplot(df.merged_fgfr2_2, x = "FGFR2_aamutation_cna", y = "FGFR2_log2expression",
          color = "FGFR2_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons3,label="p.signif")
ggboxplot(df.merged_fgfr2_2, x = "FGFR2_aamutation_cna", y = "FGFR2_log2expression",
          color = "FGFR2_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons3)

#####Quiescence group的组间比例分布图 for FGFR2 (因为评分不能反映静止或是快速增殖)#####
table(df.merged_fgfr2$FGFR2_status) #看看wt和mut的数量，差距比较大
ggplot(df.merged_fgfr2, mapping = aes(
  x = FGFR2_status, fill = Group), order = c("WT", "MUT"))+
  geom_bar(stat="count",width=0.5,position='fill')+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  geom_text(stat='count',aes(label=scales::percent(..count../sum(..count..))),
            color="white", size=3.5,position=position_fill(0.5))+
  theme_minimal()

##### QS - FGFR2_aamutation_TP53_aamutation #####
## by boxplot
my_comparisons2 <- list(c("MUT_WT","MUT_MUT"),c("WT_MUT","MUT_MUT"),c("WT_WT","WT_MUT"),c("WT_MUT","MUT_WT"),c("WT_WT","MUT_MUT"),c("WT_WT","MUT_WT")) #显示4组之间的p
ggboxplot(df.merged_fgfr2_2, x = "FGFR2_TP53_status", y = "QuiescenceScore",
          color = "FGFR2_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons2)

#合并QSgroup和aamutaiton
df.merged_fgfr2_2 <- tidyr::unite(df.merged_fgfr2_2, "FGFR2_status_group",TP53_status, Group,remove = FALSE)
head(df.merged_fgfr2_2)
table(df.merged_fgfr2_2$FGFR2_status_group)
##### QS - FGFR2_status_group #####
## by boxplot

my_comparisons_group <- list(c("MUT_Fast Cycling","MUT_Highly Quiescent"),c("MUT_Fast Cycling","WT_Fast Cycling"),c("MUT_Fast Cycling","WT_Highly Quiescent"),
                        c("MUT_Highly Quiescent","WT_Fast Cycling"),c("MUT_Highly Quiescent","WT_Highly Quiescent"),
                        c("WT_Fast Cycling","WT_Highly Quiescent")) #显示4组之间的p
ggboxplot(df.merged_fgfr2_2, x = "FGFR2_status_group", y = "QuiescenceScore",
          color = "FGFR2_status_group", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons_group)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))#横坐标倾斜45度

aggregate(df.merged_fgfr2_2$QuiescenceScore, by=list(type=df.merged_fgfr2_2$FGFR2_aamutation_cna),mean)
#####QS - FGFR2_aamutation_cna + TP53_AAmutation #####
ggboxplot(df.merged_fgfr2_2, x = "FGFR2_aamutation_cna_TP53_aamutation", y = "QuiescenceScore",
          color = "FGFR2_aamutation_cna_TP53_aamutation", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB"),
          add = "jitter",order = c("WT_WT","WT_MUT","MUT_WT","MUT_MUT","AMP_WT", "AMP_MUT","DEL_WT","DEL_MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("AMP_WT","AMP_MUT")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("DEL_WT","DEL_MUT")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("MUT_WT","MUT_MUT")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("WT_WT","WT_MUT")),label = "p.signif")
ggboxplot(df.merged_fgfr2_2, x = "FGFR2_aamutation_cna_TP53_aamutation", y = "QuiescenceScore",
          color = "FGFR2_aamutation_cna_TP53_aamutation", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                     "#D5D2C1","#00AFBB"),
          add = "jitter",order = c("WT_WT","WT_MUT","MUT_WT","MUT_MUT","AMP_WT", "AMP_MUT","DEL_WT","DEL_MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("AMP_WT","AMP_MUT")))+
  stat_compare_means(comparisons = list(c("DEL_WT","DEL_MUT")))+
  stat_compare_means(comparisons = list(c("MUT_WT","MUT_MUT")))+
  stat_compare_means(comparisons = list(c("WT_WT","WT_MUT")))

#####Expression and QS levels - FGFR2_aamutation_cna - Cancer type#####
#把每种癌症都分为一个df
#再分类绘制QS & Expression plot（boxplot 感觉violinplot没有意义）

## QS - Cancer type 一种种替换 eg.“lihc”全部替换为“ov”
nhsc_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="lihc_tcga_pan_can_atlas_2018",]
table(lihc_fgfr2$FGFR2_aamutation_cna)
ggboxplot(lihc_fgfr2, x = "FGFR2_aamutation_cna", y = "QuiescenceScore",
          color = "FGFR2_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("AMP","DEL")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("AMP","MUT")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("DEL","MUT")),label = "p.signif")

#ACC
acc_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="acc_tcga_pan_can_atlas_2018",]
table(acc_fgfr2$FGFR2_aamutation_cna)
aggregate(acc_fgfr2$QuiescenceScore, by=list(type=acc_fgfr2$FGFR2_aamutation_cna),mean)
#BLCA
blca_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="blca_tcga_pan_can_atlas_2018",]
table(blca_fgfr2$FGFR2_aamutation_cna)
aggregate(blca_fgfr2$QuiescenceScore, by=list(type=blca_fgfr2$FGFR2_aamutation_cna),mean)
#BRCA
brca_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="brca_tcga_pan_can_atlas_2018",]
table(brca_fgfr2$FGFR2_aamutation_cna)
aggregate(brca_fgfr2$QuiescenceScore, by=list(type=brca_fgfr2$FGFR2_aamutation_cna),mean)
#CESC
cesc_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="cesc_tcga_pan_can_atlas_2018",]
table(cesc_fgfr2$FGFR2_aamutation_cna)
aggregate(cesc_fgfr2$QuiescenceScore, by=list(type=cesc_fgfr2$FGFR2_aamutation_cna),mean)
#CHOL
chol_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="chol_tcga_pan_can_atlas_2018",]
table(chol_fgfr2$FGFR2_aamutation_cna)
aggregate(chol_fgfr2$QuiescenceScore, by=list(type=chol_fgfr2$FGFR2_aamutation_cna),mean)
#COAD
coadread_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="coadread_tcga_pan_can_atlas_2018",]
table(coadread_fgfr2$FGFR2_aamutation_cna)
aggregate(coadread_fgfr2$QuiescenceScore, by=list(type=coadread_fgfr2$FGFR2_aamutation_cna),mean)
#ESCA
esca_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="esca_tcga_pan_can_atlas_2018",]
table(esca_fgfr2$FGFR2_aamutation_cna)
aggregate(esca_fgfr2$QuiescenceScore, by=list(type=esca_fgfr2$FGFR2_aamutation_cna),mean)
#GBM
gbm_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="gbm_tcga_pan_can_atlas_2018",]
table(gbm_fgfr2$FGFR2_aamutation_cna)
aggregate(gbm_fgfr2$QuiescenceScore, by=list(type=gbm_fgfr2$FGFR2_aamutation_cna),mean)
#HNSC
hnsc_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="hnsc_tcga_pan_can_atlas_2018",]
table(hnsc_fgfr2$FGFR2_aamutation_cna)
aggregate(hnsc_fgfr2$QuiescenceScore, by=list(type=hnsc_fgfr2$FGFR2_aamutation_cna),mean)
#KICH
kich_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="kich_tcga_pan_can_atlas_2018",]
table(kich_fgfr2$FGFR2_aamutation_cna)
aggregate(kich_fgfr2$QuiescenceScore, by=list(type=kich_fgfr2$FGFR2_aamutation_cna),mean)
#KIRC
kirc_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="kirc_tcga_pan_can_atlas_2018",]
table(kirc_fgfr2$FGFR2_aamutation_cna)
aggregate(kirc_fgfr2$QuiescenceScore, by=list(type=kirc_fgfr2$FGFR2_aamutation_cna),mean)
#KIRP
kirp_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="kirp_tcga_pan_can_atlas_2018",]
table(kirp_fgfr2$FGFR2_aamutation_cna)
aggregate(kirp_fgfr2$QuiescenceScore, by=list(type=kirp_fgfr2$FGFR2_aamutation_cna),mean)
#LGG
lgg_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="lgg_tcga_pan_can_atlas_2018",]
table(lgg_fgfr2$FGFR2_aamutation_cna)
aggregate(lgg_fgfr2$QuiescenceScore, by=list(type=lgg_fgfr2$FGFR2_aamutation_cna),mean)
#LIHC
lihc_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="lihc_tcga_pan_can_atlas_2018",]
table(lihc_fgfr2$FGFR2_aamutation_cna)
aggregate(lihc_fgfr2$QuiescenceScore, by=list(type=lihc_fgfr2$FGFR2_aamutation_cna),mean)
#LUAD
luad_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="luad_tcga_pan_can_atlas_2018",]
table(luad_fgfr2$FGFR2_aamutation_cna)
aggregate(luad_fgfr2$QuiescenceScore, by=list(type=luad_fgfr2$FGFR2_aamutation_cna),mean)
#LUSC
lusc_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="lusc_tcga_pan_can_atlas_2018",]
table(lusc_fgfr2$FGFR2_aamutation_cna)
aggregate(lusc_fgfr2$QuiescenceScore, by=list(type=lusc_fgfr2$FGFR2_aamutation_cna),mean)
#MESO
meso_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="meso_tcga_pan_can_atlas_2018",]
table(meso_fgfr2$FGFR2_aamutation_cna)
aggregate(meso_fgfr2$QuiescenceScore, by=list(type=meso_fgfr2$FGFR2_aamutation_cna),mean)
#OV
ov_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="ov_tcga_pan_can_atlas_2018",]
table(ov_fgfr2$FGFR2_aamutation_cna)
aggregate(ov_fgfr2$QuiescenceScore, by=list(type=ov_fgfr2$FGFR2_aamutation_cna),mean)
#PAAD
paad_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="paad_tcga_pan_can_atlas_2018",]
table(paad_fgfr2$FGFR2_aamutation_cna)
aggregate(paad_fgfr2$QuiescenceScore, by=list(type=paad_fgfr2$FGFR2_aamutation_cna),mean)
#PCPG
pcpg_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="pcpg_tcga_pan_can_atlas_2018",]
table(pcpg_fgfr2$FGFR2_aamutation_cna)
aggregate(pcpg_fgfr2$QuiescenceScore, by=list(type=pcpg_fgfr2$FGFR2_aamutation_cna),mean)
#PRAD
prad_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="prad_tcga_pan_can_atlas_2018",]
table(prad_fgfr2$FGFR2_aamutation_cna)
aggregate(prad_fgfr2$QuiescenceScore, by=list(type=prad_fgfr2$FGFR2_aamutation_cna),mean)
#SARC
sarc_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="sarc_tcga_pan_can_atlas_2018",]
table(sarc_fgfr2$FGFR2_aamutation_cna)
aggregate(sarc_fgfr2$QuiescenceScore, by=list(type=sarc_fgfr2$FGFR2_aamutation_cna),mean)
#STAD
stad_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="stad_tcga_pan_can_atlas_2018",]
table(stad_fgfr2$FGFR2_aamutation_cna)
aggregate(stad_fgfr2$QuiescenceScore, by=list(type=stad_fgfr2$FGFR2_aamutation_cna),mean)
#TGCT
tgct_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="tgct_tcga_pan_can_atlas_2018",]
table(tgct_fgfr2$FGFR2_aamutation_cna)
aggregate(tgct_fgfr2$QuiescenceScore, by=list(type=tgct_fgfr2$FGFR2_aamutation_cna),mean)
#THCA
thca_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="thca_tcga_pan_can_atlas_2018",]
table(thca_fgfr2$FGFR2_aamutation_cna)
aggregate(thca_fgfr2$QuiescenceScore, by=list(type=thca_fgfr2$FGFR2_aamutation_cna),mean)
#THYM
thym_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="thym_tcga_pan_can_atlas_2018",]
table(thym_fgfr2$FGFR2_aamutation_cna)
aggregate(thym_fgfr2$QuiescenceScore, by=list(type=thym_fgfr2$FGFR2_aamutation_cna),mean)
#UCEC
ucec_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="ucec_tcga_pan_can_atlas_2018",]
table(ucec_fgfr2$FGFR2_aamutation_cna)
aggregate(ucec_fgfr2$QuiescenceScore, by=list(type=ucec_fgfr2$FGFR2_aamutation_cna),mean)
#UCS
ucs_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="ucs_tcga_pan_can_atlas_2018",]
table(ucs_fgfr2$FGFR2_aamutation_cna)
aggregate(ucs_fgfr2$QuiescenceScore, by=list(type=ucs_fgfr2$FGFR2_aamutation_cna),mean)
#UVM
uvm_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="uvm_tcga_pan_can_atlas_2018",]
table(uvm_fgfr2$FGFR2_aamutation_cna)
aggregate(uvm_fgfr2$QuiescenceScore, by=list(type=uvm_fgfr2$FGFR2_aamutation_cna),mean)
## Exp - Cancer type 一种种替换 eg.“ov”全部替换为“ov”
ggboxplot(lihc_fgfr2, x = "FGFR2_aamutation_cna", y = "FGFR2_log2expression",
          color = "FGFR2_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("AMP","DEL")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("AMP","MUT")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("DEL","MUT")),label = "p.signif")
#####QS - FGFR2_aamutation_cna + TP53_AAmutation - Cancer type#####
lihc_fgfr2<- df.merged_fgfr2_2[df.merged_fgfr2_2$STUDY_ID=="lihc_tcga_pan_can_atlas_2018",]
ggboxplot(lihc_fgfr2, x = "FGFR2_aamutation_cna_TP53_aamutation", y = "QuiescenceScore",
          color = "FGFR2_aamutation_cna_TP53_aamutation", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB"),
          add = "jitter",order = c("WT_WT","WT_MUT","MUT_WT","MUT_MUT","AMP_WT", "AMP_MUT","DEL_WT","DEL_MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("AMP_WT","AMP_MUT")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("DEL_WT","DEL_MUT")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("MUT_WT","MUT_MUT")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("WT_WT","WT_MUT")),label = "p.signif")


#####导入整倍体数据AS---已下为整倍体所有代码 可直接复制后分别替换基因大小写#####
aneuploidy <- read.csv("TCGA_aneuploidy.csv")
head(aneuploidy)
hist(aneuploidy$AneuploidyScore)

df.merged_fgfr2_3 <- merge(df.merged_fgfr2_2[,c("SAMPLE_ID",
                                             "FGFR2_TP53_status",
                                             "TP53_status","FGFR2_status","FGFR2_log2expression",
                                             "QuiescenceScore","Group",
                                             "FGFR2_FGFR2CNA",
                                             "STUDY_ID",
                                             "FGFR2_CNA",
                                             "FGFR2_aamutation_cna",
                                             "FGFR2_aamutation_cna_TP53_aamutation")], 
                          aneuploidy[,c("Sample","Type","AneuploidyScore",
                                        "AS_del","AS_amp",
                                        "Genome_doublings","Leuk","Purity","Stroma",
                                        "Stroma_notLeukocyte","Stroma_notLeukocyte_Floor",
                                        "SilentMutationspeMb",
                                        "Non.silentMutationsperMb")],
                          by.x="SAMPLE_ID", by.y="Sample", #x和y表合并的依据是相同的sample id
                          all.x=FALSE, all.y=FALSE) 
head(df.merged_fgfr2_3)

##### AS - FGFR2_aamutation #####
table(df.merged_fgfr2_3$FGFR2_status)
ggboxplot(df.merged_fgfr2_3, x = "FGFR2_status", y = "AneuploidyScore",
          color = "FGFR2_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### AS_amp - FGFR2_aamutation #####
table(df.merged_fgfr2_3$FGFR2_status)
ggboxplot(df.merged_fgfr2_3, x = "FGFR2_status", y = "AS_amp",
          color = "FGFR2_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### AS_del - FGFR2_aamutation #####
table(df.merged_fgfr2_3$FGFR2_status)
ggboxplot(df.merged_fgfr2_3, x = "FGFR2_status", y = "AS_del",
          color = "FGFR2_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Gb - FGFR2_aamutation #####
table(df.merged_fgfr2_3$FGFR2_status)
ggboxplot(df.merged_fgfr2_3, x = "FGFR2_status", y = "Genome_doublings",
          color = "FGFR2_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Leuk - FGFR2_aamutation #####
df.merged_fgfr2_3$Leuk[which(df.merged_fgfr2_3$Leuk =='#N/A')] <- 'NA'
class(df.merged_fgfr2_3$Leuk) #检查Leuk的数据类型，原数据框中为字符
df.merged_fgfr2_3$Leuk <- as.numeric(df.merged_fgfr2_3$Leuk)#使用as.numeric将其数据类型换为数值型
class(df.merged_fgfr2_3$Leuk)#再次检查
hist(df.merged_fgfr2_3$Leuk)

ggboxplot(df.merged_fgfr2_3, x = "FGFR2_status", y = "Leuk",
          color = "FGFR2_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Purity - FGFR2_aamutation #####
df.merged_fgfr2_3$Purity[which(df.merged_fgfr2_3$Purity =='#N/A')] <- 'NA'
class(df.merged_fgfr2_3$Purity) #检查Leuk的数据类型，原数据框中为字符
df.merged_fgfr2_3$Purity <- as.numeric(df.merged_fgfr2_3$Purity)#使用as.numeric将其数据类型换为数值型
class(df.merged_fgfr2_3$Purity)#再次检查
hist(df.merged_fgfr2_3$Purity)

ggboxplot(df.merged_fgfr2_3, x = "FGFR2_status", y = "Purity",
          color = "FGFR2_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Stroma - FGFR2_aamutation #####
df.merged_fgfr2_3$Stroma[which(df.merged_fgfr2_3$Stroma =='#N/A')] <- 'NA'
class(df.merged_fgfr2_3$Stroma) #检查Leuk的数据类型，原数据框中为字符
df.merged_fgfr2_3$Stroma <- as.numeric(df.merged_fgfr2_3$Stroma)#使用as.numeric将其数据类型换为数值型
class(df.merged_fgfr2_3$Stroma)#再次检查
hist(df.merged_fgfr2_3$Stroma)

ggboxplot(df.merged_fgfr2_3, x = "FGFR2_status", y = "Stroma",
          color = "FGFR2_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Stroma_notLeukocyte - FGFR2_aamutation #####
df.merged_fgfr2_3$Stroma_notLeukocyte[which(df.merged_fgfr2_3$Stroma_notLeukocyte =='#N/A')] <- 'NA'
class(df.merged_fgfr2_3$Stroma_notLeukocyte) #检查Leuk的数据类型，原数据框中为字符
df.merged_fgfr2_3$Stroma_notLeukocyte <- as.numeric(df.merged_fgfr2_3$Stroma_notLeukocyte)#使用as.numeric将其数据类型换为数值型
class(df.merged_fgfr2_3$Stroma_notLeukocyte)#再次检查
hist(df.merged_fgfr2_3$Stroma_notLeukocyte)

ggboxplot(df.merged_fgfr2_3, x = "FGFR2_status", y = "Stroma_notLeukocyte",
          color = "FGFR2_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### SilentMutationspeMb - FGFR2_aamutation #####
df.merged_fgfr2_3$SilentMutationspeMb[which(df.merged_fgfr2_3$SilentMutationspeMb =='#N/A')] <- 'NA'
class(df.merged_fgfr2_3$SilentMutationspeMb) #检查Leuk的数据类型，原数据框中为字符
df.merged_fgfr2_3$SilentMutationspeMb <- as.numeric(df.merged_fgfr2_3$SilentMutationspeMb)#使用as.numeric将其数据类型换为数值型
class(df.merged_fgfr2_3$SilentMutationspeMb)#再次检查
hist(df.merged_fgfr2_3$SilentMutationspeMb)

ggboxplot(df.merged_fgfr2_3, x = "FGFR2_status", y = "SilentMutationspeMb",
          color = "FGFR2_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Non.silentMutationsperMb - FGFR2_aamutation #####
df.merged_fgfr2_3$Non.silentMutationsperMb[which(df.merged_fgfr2_3$Non.silentMutationsperMb =='#N/A')] <- 'NA'
class(df.merged_fgfr2_3$Non.silentMutationsperMb) #检查Leuk的数据类型，原数据框中为字符
df.merged_fgfr2_3$Non.silentMutationsperMb <- as.numeric(df.merged_fgfr2_3$Non.silentMutationsperMb)#使用as.numeric将其数据类型换为数值型
class(df.merged_fgfr2_3$Non.silentMutationsperMb)#再次检查
hist(df.merged_fgfr2_3$Non.silentMutationsperMb)

ggboxplot(df.merged_fgfr2_3, x = "FGFR2_status", y = "Non.silentMutationsperMb",
          color = "FGFR2_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
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
ggplot(df.merged_fgfr2_3, aes(x=AneuploidyScore, y=QuiescenceScore)) +
  geom_point(aes(col=Type)) +
  geom_smooth(method="loess", se=F) + 
  xlim(c(0, 40)) + 
  ylim(c(-30, 30)) + 
  labs(subtitle="QS Vs AS", 
       title="Scatterplot_FGFR2")+
  stat_cor(method = 'pearson', aes(x = AneuploidyScore, y = QuiescenceScore))

ggplot(df.merged_fgfr2_3, aes(x=AneuploidyScore, y=QuiescenceScore)) +
  geom_point(aes(col=Type)) +
  geom_smooth(method="lm", se=F) + 
  xlim(c(0, 40)) + 
  ylim(c(-30, 30)) + 
  labs(subtitle="QS Vs AS", 
       title="Scatterplot_FGFR2")+
  stat_cor(method = 'pearson', aes(x = AneuploidyScore, y = QuiescenceScore))

ggplot(df.merged_fgfr2_3, aes(x=AS_del, y=QuiescenceScore)) +
  geom_point(aes(col=Type)) +
  geom_smooth(method="lm", se=F) + 
  xlim(c(0, 40)) + 
  ylim(c(-30, 30)) + 
  labs(subtitle="QS Vs AS_del", 
       title="Scatterplot_FGFR2")+
  stat_cor(method = 'pearson', aes(x = AS_del, y = QuiescenceScore))

ggplot(df.merged_fgfr2_3, aes(x=AS_amp, y=QuiescenceScore)) +
  geom_point(aes(col=Type)) +
  geom_smooth(method="lm", se=F) + 
  xlim(c(0, 40)) + 
  ylim(c(-30, 30)) + 
  labs(subtitle="QS Vs AS_amp", 
       title="Scatterplot_FGFR2")+
  stat_cor(method = 'pearson', aes(x = AS_amp, y = QuiescenceScore))

ggplot(data = df.merged_fgfr2_3, aes(x = AneuploidyScore, y = QuiescenceScore,  color = FGFR2_aamutation_cna)) + 
  geom_point() + geom_smooth(method = lm) +
  scale_color_manual(values = c('#FF7400', '#009999', '#3914AF',"#000000"))  +
  labs(title = 'QS Vs AS - FGFR2_aamutation_cna') + 
  stat_cor(method = 'pearson', aes(x = AneuploidyScore, y = QuiescenceScore,  color = FGFR2_aamutation_cna))
#####AS- FGFR2_aamutation_cna#####
aggregate(df.merged_fgfr2_3$AneuploidyScore, by=list(type=df.merged_fgfr2_3$FGFR2_aamutation_cna),mean)

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_fgfr2_3, x = "FGFR2_aamutation_cna", y = "AneuploidyScore",
          color = "FGFR2_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####AS_amp- FGFR2_aamutation_cna#####
#mean of groups
aggregate(df.merged_fgfr2_3$AS_amp, by=list(type=df.merged_fgfr2_3$FGFR2_aamutation_cna),mean)

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_fgfr2_3, x = "FGFR2_aamutation_cna", y = "AS_amp",
          color = "FGFR2_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####AS_del- FGFR2_aamutation_cna#####
#mean of groups
aggregate(df.merged_fgfr2_3$AS_del, by=list(type=df.merged_fgfr2_3$FGFR2_aamutation_cna),mean)

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_fgfr2_3, x = "FGFR2_aamutation_cna", y = "AS_del",
          color = "FGFR2_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####Gb- FGFR2_aamutation_cna#####
#mean of groups
aggregate(df.merged_fgfr2_3$Genome_doublings, by=list(type=df.merged_fgfr2_3$FGFR2_aamutation_cna),mean)

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_fgfr2_3, x = "FGFR2_aamutation_cna", y = "Genome_doublings",
          color = "FGFR2_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####Leuk- FGFR2_aamutation_cna#####
#mean of groups
aggregate(df.merged_fgfr2_3$Leuk, by=list(type=df.merged_fgfr2_3$FGFR2_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_fgfr2_3, x = "FGFR2_aamutation_cna", y = "Leuk",
          color = "FGFR2_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####Purity- FGFR2_aamutation_cna#####
#mean of groups
aggregate(df.merged_fgfr2_3$Purity, by=list(type=df.merged_fgfr2_3$FGFR2_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_fgfr2_3, x = "FGFR2_aamutation_cna", y = "Purity",
          color = "FGFR2_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####Stroma- FGFR2_aamutation_cna#####
#mean of groups
aggregate(df.merged_fgfr2_3$Stroma, by=list(type=df.merged_fgfr2_3$FGFR2_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_fgfr2_3, x = "FGFR2_aamutation_cna", y = "Stroma",
          color = "FGFR2_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####Stroma_notLeukocyte- FGFR2_aamutation_cna#####
#mean of groups
aggregate(df.merged_fgfr2_3$Stroma_notLeukocyte, by=list(type=df.merged_fgfr2_3$FGFR2_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_fgfr2_3, x = "FGFR2_aamutation_cna", y = "Stroma_notLeukocyte",
          color = "FGFR2_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####SilentMutationspeMb- FGFR2_aamutation_cna#####
my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_fgfr2_3, x = "FGFR2_aamutation_cna", y = "SilentMutationspeMb",
          color = "FGFR2_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####Non.silentMutationsperMb- FGFR2_aamutation_cna#####
my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_fgfr2_3, x = "FGFR2_aamutation_cna", y = "Non.silentMutationsperMb",
          color = "FGFR2_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####导入样本总突变数据NM---已下为突变数体所有代码 可直接复制后分别替换基因大小写#####
pgen <- read.csv("pgen.csv")
head(pgen)
hist(pgen$number_of_mutations)

df.merged_fgfr2_4 <- merge(df.merged_fgfr2_3[,c("SAMPLE_ID",
                                                "FGFR2_TP53_status",
                                                "TP53_status","FGFR2_status","FGFR2_log2expression",
                                                "QuiescenceScore","Group",
                                                "FGFR2_FGFR2CNA",
                                                "STUDY_ID",
                                                "FGFR2_CNA",
                                                "FGFR2_aamutation_cna",
                                                "FGFR2_aamutation_cna_TP53_aamutation",
                                                "Type",
                                                "AneuploidyScore",
                                                "AS_del","AS_amp",
                                                "Genome_doublings","Leuk","Purity","Stroma",
                                                "Stroma_notLeukocyte","Stroma_notLeukocyte_Floor",
                                                "SilentMutationspeMb",
                                                "Non.silentMutationsperMb")], 
                           pgen[,c("sample_name","tumor_type",
                                   "number_of_mutations")],
                           by.x="SAMPLE_ID", by.y="sample_name", #x和y表合并的依据是相同的sample id
                           all.x=FALSE, all.y=FALSE) 
head(df.merged_fgfr2_4)
#####Number of Mutations- FGFR2_aamutation_cna#####
#mean of groups
aggregate(df.merged_fgfr2_4$number_of_mutations, by=list(type=df.merged_fgfr2_4$FGFR2_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_fgfr2_4, x = "FGFR2_aamutation_cna", y = "number_of_mutations",
          color = "FGFR2_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3)
