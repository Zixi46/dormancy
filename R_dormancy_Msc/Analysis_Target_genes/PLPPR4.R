### This script analyses PLPPR4 mutation data and links with quiescence and expressin levels in tumours.
##### 初步数据处理 aa mutaion status与Exp QS合并#####
### Read mutation data downloaded from cbioportal:
setwd("/Users/lizixi/Desktop/Final Project/5.18 task")
mutations_plppr4 <- read.delim("mutations_cbioportal_plppr4.txt")
head(mutations_plppr4)
## How many samples have a mutation?
table(mutations_plppr4$PLPPR4)
tb_plppr4 <- table(mutations_plppr4$PLPPR4)
sort(tb_plppr4) #descend升序排列PHIN1列的数字元素

# Add a new column that just denotes WT and MUT for PLPPR4:
mutations_plppr4[which(mutations_plppr4$PLPPR4 == "NS"),]$PLPPR4 <- "WT" #****将PLPPR4中的NS替换为WT
mutations_plppr4$PLPPR4_status <- mutations_plppr4$PLPPR4 #为mutations_plppr4表添加一列为PLPPR4_status，同时将这一列定义为PLPPR4（和PLPPR4一样）
mutations_plppr4[which(mutations_plppr4$PLPPR4 != "WT"),]$PLPPR4_status <- "MUT" #选出mutation中PLPPR4列中不等于“WT”的元素，并在PLPPR4_status中将其全赋值为“MUT”
head(mutations_plppr4)
# How many samples with mutations_plppr4 in PLPPR4?
table(mutations_plppr4$PLPPR4_status)
# Read expression data for PLPPR4:
expression_plppr4 <- read.delim("mRNA Expression_RSEM_cbioportal_plppr4.txt")
head(expression_plppr4)
# Plot histogram of expression values of PLPPR4:
hist(expression_plppr4$PLPPR4)
# It's good to take the log2 of expression values to avoid working with very large numbers:
expression_plppr4$PLPPR4_log2expression <- log2(expression_plppr4$PLPPR4+1)
head(expression_plppr4)
hist(expression_plppr4$PLPPR4_log2expression) # more decent value

# Merge PLPPR4 expression and mutation data:
df.merged_plppr4 <- merge(mutations_plppr4[,c("SAMPLE_ID","PLPPR4_status","STUDY_ID")], #mutation中需要合并的纵列是sample id，plppr4_status
                         expression_plppr4[,c("SAMPLE_ID","PLPPR4_log2expression")], #expresssion中需要合并的纵列是sample ID，PLPPR4_log2expression
                         by.x="SAMPLE_ID", by.y="SAMPLE_ID", #x和y表合并的依据是相同的sample id
                         all.x=FALSE, all.y=FALSE) #x和y表中不需要所有的行都合并（只合并sample id一样的，只在其中一个表中才存在的sample id丢弃）
head(df.merged_plppr4)
# Merge TP53 mutation data+QS
df.merged_plppr4 <- merge(df.merged_plppr4[,c("SAMPLE_ID","PLPPR4_status","STUDY_ID","PLPPR4_log2expression")], #mutation中需要合并的纵列是sample id，plppr4_status
                         df.merged2[,c("SAMPLE_ID","TP53_status","QuiescenceScore","Group")], #expresssion中需要合并的纵列是sample ID，PLPPR4_log2expression
                         by.x="SAMPLE_ID", by.y="SAMPLE_ID", #x和y表合并的依据是相同的sample id
                         all.x=FALSE, all.y=FALSE) #x和y表中不需要所有的行都合并（只合并sample id一样的，只在其中一个表中才存在的sample id丢弃）
head(df.merged_plppr4)
###### 增加PLPPR4突变类型分析列 “AMP”“DEL” #####
#读取cna文件
cna_PLPPR4 <- read.delim("cna_fromcBioportal_plppr4.txt")
names(cna_PLPPR4)[3] <- "PLPPR4_CNA"#给第三列重命名
head(cna_PLPPR4)
table(cna_PLPPR4$PLPPR4_CNA)
#将cna文件中的ID 和 CNA与总表合并
head(df.merged_plppr4)
df.merged_plppr4_2 <- merge(df.merged_plppr4[,c("SAMPLE_ID","STUDY_ID","PLPPR4_status",
                                              "TP53_status","PLPPR4_log2expression",
                                              "QuiescenceScore","Group")], 
                           cna_PLPPR4[,c("SAMPLE_ID","PLPPR4_CNA")],
                           by.x="SAMPLE_ID", by.y="SAMPLE_ID", #x和y表合并的依据是相同的sample id
                           all.x=FALSE, all.y=FALSE) 
head(df.merged_plppr4_2)
table(df.merged_plppr4_2$PLPPR4_CNA)
#数据整理，完整案例分析
##exclude cns-NP(not profiled) samples from the analysis 剔除包含NP的样本
##把NP都换为NA
df.merged_plppr4_2[which(df.merged_plppr4_2$PLPPR4_CNA == "NP"),]$PLPPR4_CNA <- "NA" #****将PLPPR4_CAN中的NP替换为NA
table(df.merged_plppr4_2$PLPPR4_CNA)
##把-2换为-1（数量太少了），把2换为1
#"-2" is a deep loss, possibly a homozygous deletion
#"-1" is a single-copy loss (heterozygous deletion)
#"0" is diploid
#"1" indicates a low-level gain
#"2" is a high-level amplification.
##将1转化为AMP,-1转化为DEL，O转化为NCN
df.merged_plppr4_2[which(df.merged_plppr4_2$PLPPR4_CNA == "-2"),]$PLPPR4_CNA <- "-1" 
df.merged_plppr4_2[which(df.merged_plppr4_2$PLPPR4_CNA == "2"),]$PLPPR4_CNA <- "1" 
df.merged_plppr4_2[which(df.merged_plppr4_2$PLPPR4_CNA == "-1"),]$PLPPR4_CNA <- "DEL" 
df.merged_plppr4_2[which(df.merged_plppr4_2$PLPPR4_CNA == "1"),]$PLPPR4_CNA <- "AMP" 
df.merged_plppr4_2[which(df.merged_plppr4_2$PLPPR4_CNA == "0"),]$PLPPR4_CNA <- "NCN" 
table(df.merged_plppr4_2$PLPPR4_CNA)
##去掉NA样本
df.merged_plppr4_2 <- df.merged_plppr4_2[!grepl("NA", df.merged_plppr4_2$PLPPR4_CNA),]
table(df.merged_plppr4_2$PLPPR4_CNA)
##将PLPPR4 AAmutation信息和CNA信息结合为新的一列Mutation_CNA,生成6个不同组别
df.merged_plppr4_2 <- tidyr::unite(df.merged_plppr4_2, "PLPPR4_PLPPR4CNA", PLPPR4_status, PLPPR4_CNA, remove = FALSE)
head(df.merged_plppr4_2)
table(df.merged_plppr4_2$PLPPR4_PLPPR4CNA)
#重排顺序
df.merged_plppr4_2 <- df.merged_plppr4_2[,c("SAMPLE_ID","STUDY_ID","PLPPR4_status","PLPPR4_CNA","PLPPR4_PLPPR4CNA","TP53_status","PLPPR4_log2expression","QuiescenceScore","Group")]
table(df.merged_plppr4$PLPPR4_status)
table(df.merged_plppr4$TP53_status)
#根据aamutaiton+CNA进行重新分组
df.merged_plppr4_2$PLPPR4_aamutation_cna <- df.merged_plppr4_2$PLPPR4_PLPPR4CNA #添加新列
df.merged_plppr4_2[which(df.merged_plppr4_2$PLPPR4_aamutation_cna == "WT_NCN"),]$PLPPR4_aamutation_cna <- "WT"# WT-NCN → WT
df.merged_plppr4_2[which(df.merged_plppr4_2$PLPPR4_aamutation_cna == "MUT_AMP"),]$PLPPR4_aamutation_cna <- "MUT"# MUT-AMP → MUT
df.merged_plppr4_2[which(df.merged_plppr4_2$PLPPR4_aamutation_cna == "MUT_NCN"),]$PLPPR4_aamutation_cna <- "MUT"# WUT-NCN → MUT
df.merged_plppr4_2[which(df.merged_plppr4_2$PLPPR4_aamutation_cna == "MUT_DEL"),]$PLPPR4_aamutation_cna <- "DEL"# MUT-DEL → DEL
df.merged_plppr4_2[which(df.merged_plppr4_2$PLPPR4_aamutation_cna == "WT_AMP"),]$PLPPR4_aamutation_cna <- "AMP"# WT-AMP → AMP
df.merged_plppr4_2[which(df.merged_plppr4_2$PLPPR4_aamutation_cna == "WT_DEL"),]$PLPPR4_aamutation_cna <- "DEL"# WT-DEL → DEL
table(df.merged_plppr4_2$PLPPR4_aamutation_cna)
#####合并 PLPPR4+TP53 status and PLPPR4+TP53 status #####
#Combine PLPPR4+TP53 to create 4 groups 。用tidyr::unite合并两列（PLPPR4_status，TP53_status）的元素为新的一列（PLPPR4_TP53_status），默认以_连接（可以修改）
library(tidyr)
df.merged_plppr4_2 <- tidyr::unite(df.merged_plppr4_2, "PLPPR4_TP53_status", PLPPR4_status, TP53_status, remove = FALSE)
head(df.merged_plppr4_2)
table(df.merged_plppr4_2$PLPPR4_TP53_status)
#####合并 PLPPR4_aamutation_cna+TP53_status#####
df.merged_plppr4_2 <- tidyr::unite(df.merged_plppr4_2, "PLPPR4_aamutation_cna_TP53_aamutation", PLPPR4_aamutation_cna, TP53_status, remove = FALSE)
table(df.merged_plppr4_2$PLPPR4_aamutation_cna_TP53_aamutation)

#####Expression and QS levels - PLPPR4_aamutation#####
#QS
## by boxplot
ggboxplot(df.merged_plppr4, x = "PLPPR4_status", y = "QuiescenceScore",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")
ggboxplot(df.merged_plppr4, x = "PLPPR4_status", y = "QuiescenceScore",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)
#Exp
## by boxplot
ggboxplot(df.merged_plppr4, x = "PLPPR4_status", y = "PLPPR4_log2expression",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")
ggboxplot(df.merged_plppr4, x = "PLPPR4_status", y = "PLPPR4_log2expression",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)

#####Expression and QS levels - PLPPR4_aamutation_cna#####
#QS
## by boxplot
aggregate(df.merged_plppr4_2$QuiescenceScore, by=list(type=df.merged_plppr4_2$PLPPR4_aamutation_cna),mean)
ggboxplot(df.merged_plppr4_2, x = "PLPPR4_aamutation_cna", y = "QuiescenceScore",
          color = "PLPPR4_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons3,label="p.signif")
ggboxplot(df.merged_plppr4_2, x = "PLPPR4_aamutation_cna", y = "QuiescenceScore",
          color = "PLPPR4_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons3)
#Exp
## by boxplot
ggboxplot(df.merged_plppr4_2, x = "PLPPR4_aamutation_cna", y = "PLPPR4_log2expression",
          color = "PLPPR4_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons3,label="p.signif")
ggboxplot(df.merged_plppr4_2, x = "PLPPR4_aamutation_cna", y = "PLPPR4_log2expression",
          color = "PLPPR4_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons3)

#####Quiescence group的组间比例分布图 for PLPPR4 (因为评分不能反映静止或是快速增殖)#####
table(df.merged_plppr4$PLPPR4_status) #看看wt和mut的数量，差距比较大
ggplot(df.merged_plppr4, mapping = aes(
  x = PLPPR4_status, fill = Group), order = c("WT", "MUT"))+
  geom_bar(stat="count",width=0.5,position='fill')+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  geom_text(stat='count',aes(label=scales::percent(..count../sum(..count..))),
            color="white", size=3.5,position=position_fill(0.5))+
  theme_minimal()

##### QS - PLPPR4_aamutation_TP53_aamutation #####
## by boxplot
my_comparisons2 <- list(c("MUT_WT","MUT_MUT"),c("WT_MUT","MUT_MUT"),c("WT_WT","WT_MUT"),c("WT_MUT","MUT_WT"),c("WT_WT","MUT_MUT"),c("WT_WT","MUT_WT")) #显示4组之间的p
ggboxplot(df.merged_plppr4_2, x = "PLPPR4_TP53_status", y = "QuiescenceScore",
          color = "PLPPR4_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons2)
#####QS - PLPPR4_aamutation_cna + TP53_AAmutation #####
ggboxplot(df.merged_plppr4_2, x = "PLPPR4_aamutation_cna_TP53_aamutation", y = "QuiescenceScore",
          color = "PLPPR4_aamutation_cna_TP53_aamutation", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                     "#D5D2C1","#00AFBB"),
          add = "jitter",order = c("WT_WT","WT_MUT","MUT_WT","MUT_MUT","AMP_WT", "AMP_MUT","DEL_WT","DEL_MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("AMP_WT","AMP_MUT")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("DEL_WT","DEL_MUT")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("MUT_WT","MUT_MUT")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("WT_WT","WT_MUT")),label = "p.signif")
ggboxplot(df.merged_plppr4_2, x = "PLPPR4_aamutation_cna_TP53_aamutation", y = "QuiescenceScore",
          color = "PLPPR4_aamutation_cna_TP53_aamutation", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                      "#D5D2C1","#00AFBB"),
          add = "jitter",order = c("WT_WT","WT_MUT","MUT_WT","MUT_MUT","AMP_WT", "AMP_MUT","DEL_WT","DEL_MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("AMP_WT","AMP_MUT")))+
  stat_compare_means(comparisons = list(c("DEL_WT","DEL_MUT")))+
  stat_compare_means(comparisons = list(c("MUT_WT","MUT_MUT")))+
  stat_compare_means(comparisons = list(c("WT_WT","WT_MUT")))


#####Expression and QS levels - PLPPR4_aamutation_cna - Cancer type热气球图#####
#把每种癌症都分为一个df
#再分类绘制QS & Expression plot（boxplot 感觉violinplot没有意义）

## QS - Cancer type 一种种替换 eg.“ov”全部替换为“ov”
uvm_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="uvm_tcga_pan_can_atlas_2018",]
table(uvm_plppr4$PLPPR4_aamutation_cna)
ggboxplot(uvm_plppr4, x = "PLPPR4_aamutation_cna", y = "QuiescenceScore",
          color = "PLPPR4_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("AMP","MUT")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("AMP","DEL")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("DEL","MUT")),label = "p.signif")

#ACC
acc_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="acc_tcga_pan_can_atlas_2018",]
table(acc_plppr4$PLPPR4_aamutation_cna)
aggregate(acc_plppr4$QuiescenceScore, by=list(type=acc_plppr4$PLPPR4_aamutation_cna),mean)
#BLCA
blca_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="blca_tcga_pan_can_atlas_2018",]
table(blca_plppr4$PLPPR4_aamutation_cna)
aggregate(blca_plppr4$QuiescenceScore, by=list(type=blca_plppr4$PLPPR4_aamutation_cna),mean)
#BRCA
brca_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="brca_tcga_pan_can_atlas_2018",]
table(brca_plppr4$PLPPR4_aamutation_cna)
aggregate(brca_plppr4$QuiescenceScore, by=list(type=brca_plppr4$PLPPR4_aamutation_cna),mean)
#CESC
cesc_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="cesc_tcga_pan_can_atlas_2018",]
table(cesc_plppr4$PLPPR4_aamutation_cna)
aggregate(cesc_plppr4$QuiescenceScore, by=list(type=cesc_plppr4$PLPPR4_aamutation_cna),mean)
#CHOL
chol_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="chol_tcga_pan_can_atlas_2018",]
table(chol_plppr4$PLPPR4_aamutation_cna)
aggregate(chol_plppr4$QuiescenceScore, by=list(type=chol_plppr4$PLPPR4_aamutation_cna),mean)
#COAD
coadread_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="coadread_tcga_pan_can_atlas_2018",]
table(coadread_plppr4$PLPPR4_aamutation_cna)
aggregate(coadread_plppr4$QuiescenceScore, by=list(type=coadread_plppr4$PLPPR4_aamutation_cna),mean)
#ESCA
esca_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="esca_tcga_pan_can_atlas_2018",]
table(esca_plppr4$PLPPR4_aamutation_cna)
aggregate(esca_plppr4$QuiescenceScore, by=list(type=esca_plppr4$PLPPR4_aamutation_cna),mean)
#GBM
gbm_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="gbm_tcga_pan_can_atlas_2018",]
table(gbm_plppr4$PLPPR4_aamutation_cna)
aggregate(gbm_plppr4$QuiescenceScore, by=list(type=gbm_plppr4$PLPPR4_aamutation_cna),mean)
#HNSC
hnsc_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="hnsc_tcga_pan_can_atlas_2018",]
table(hnsc_plppr4$PLPPR4_aamutation_cna)
aggregate(hnsc_plppr4$QuiescenceScore, by=list(type=hnsc_plppr4$PLPPR4_aamutation_cna),mean)
#KICH
kich_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="kich_tcga_pan_can_atlas_2018",]
table(kich_plppr4$PLPPR4_aamutation_cna)
aggregate(kich_plppr4$QuiescenceScore, by=list(type=kich_plppr4$PLPPR4_aamutation_cna),mean)
#KIRC
kirc_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="kirc_tcga_pan_can_atlas_2018",]
table(kirc_plppr4$PLPPR4_aamutation_cna)
aggregate(kirc_plppr4$QuiescenceScore, by=list(type=kirc_plppr4$PLPPR4_aamutation_cna),mean)
#KIRP
kirp_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="kirp_tcga_pan_can_atlas_2018",]
table(kirp_plppr4$PLPPR4_aamutation_cna)
aggregate(kirp_plppr4$QuiescenceScore, by=list(type=kirp_plppr4$PLPPR4_aamutation_cna),mean)
#LGG
lgg_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="lgg_tcga_pan_can_atlas_2018",]
table(lgg_plppr4$PLPPR4_aamutation_cna)
aggregate(lgg_plppr4$QuiescenceScore, by=list(type=lgg_plppr4$PLPPR4_aamutation_cna),mean)
#LIHC
lihc_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="lihc_tcga_pan_can_atlas_2018",]
table(lihc_plppr4$PLPPR4_aamutation_cna)
aggregate(lihc_plppr4$QuiescenceScore, by=list(type=lihc_plppr4$PLPPR4_aamutation_cna),mean)
#LUAD
luad_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="luad_tcga_pan_can_atlas_2018",]
table(luad_plppr4$PLPPR4_aamutation_cna)
aggregate(luad_plppr4$QuiescenceScore, by=list(type=luad_plppr4$PLPPR4_aamutation_cna),mean)
#LUSC
lusc_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="lusc_tcga_pan_can_atlas_2018",]
table(lusc_plppr4$PLPPR4_aamutation_cna)
aggregate(lusc_plppr4$QuiescenceScore, by=list(type=lusc_plppr4$PLPPR4_aamutation_cna),mean)
#MESO
meso_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="meso_tcga_pan_can_atlas_2018",]
table(meso_plppr4$PLPPR4_aamutation_cna)
aggregate(meso_plppr4$QuiescenceScore, by=list(type=meso_plppr4$PLPPR4_aamutation_cna),mean)
#OV
ov_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="ov_tcga_pan_can_atlas_2018",]
table(ov_plppr4$PLPPR4_aamutation_cna)
aggregate(ov_plppr4$QuiescenceScore, by=list(type=ov_plppr4$PLPPR4_aamutation_cna),mean)
#PAAD
paad_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="paad_tcga_pan_can_atlas_2018",]
table(paad_plppr4$PLPPR4_aamutation_cna)
aggregate(paad_plppr4$QuiescenceScore, by=list(type=paad_plppr4$PLPPR4_aamutation_cna),mean)
#PCPG
pcpg_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="pcpg_tcga_pan_can_atlas_2018",]
table(pcpg_plppr4$PLPPR4_aamutation_cna)
aggregate(pcpg_plppr4$QuiescenceScore, by=list(type=pcpg_plppr4$PLPPR4_aamutation_cna),mean)
#PRAD
prad_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="prad_tcga_pan_can_atlas_2018",]
table(prad_plppr4$PLPPR4_aamutation_cna)
aggregate(prad_plppr4$QuiescenceScore, by=list(type=prad_plppr4$PLPPR4_aamutation_cna),mean)
#SARC
sarc_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="sarc_tcga_pan_can_atlas_2018",]
table(sarc_plppr4$PLPPR4_aamutation_cna)
aggregate(sarc_plppr4$QuiescenceScore, by=list(type=sarc_plppr4$PLPPR4_aamutation_cna),mean)
#STAD
stad_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="stad_tcga_pan_can_atlas_2018",]
table(stad_plppr4$PLPPR4_aamutation_cna)
aggregate(stad_plppr4$QuiescenceScore, by=list(type=stad_plppr4$PLPPR4_aamutation_cna),mean)
#TGCT
tgct_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="tgct_tcga_pan_can_atlas_2018",]
table(tgct_plppr4$PLPPR4_aamutation_cna)
aggregate(tgct_plppr4$QuiescenceScore, by=list(type=tgct_plppr4$PLPPR4_aamutation_cna),mean)
#THCA
thca_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="thca_tcga_pan_can_atlas_2018",]
table(thca_plppr4$PLPPR4_aamutation_cna)
aggregate(thca_plppr4$QuiescenceScore, by=list(type=thca_plppr4$PLPPR4_aamutation_cna),mean)
#THYM
thym_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="thym_tcga_pan_can_atlas_2018",]
table(thym_plppr4$PLPPR4_aamutation_cna)
aggregate(thym_plppr4$QuiescenceScore, by=list(type=thym_plppr4$PLPPR4_aamutation_cna),mean)
#UCEC
ucec_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="ucec_tcga_pan_can_atlas_2018",]
table(ucec_plppr4$PLPPR4_aamutation_cna)
aggregate(ucec_plppr4$QuiescenceScore, by=list(type=ucec_plppr4$PLPPR4_aamutation_cna),mean)
#UCS
ucs_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="ucs_tcga_pan_can_atlas_2018",]
table(ucs_plppr4$PLPPR4_aamutation_cna)
aggregate(ucs_plppr4$QuiescenceScore, by=list(type=ucs_plppr4$PLPPR4_aamutation_cna),mean)
#UVM
uvm_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="uvm_tcga_pan_can_atlas_2018",]
table(uvm_plppr4$PLPPR4_aamutation_cna)
aggregate(uvm_plppr4$QuiescenceScore, by=list(type=uvm_plppr4$PLPPR4_aamutation_cna),mean)

#####Ballon plot:QS-GENE_aamutation_cna + Cancer type#####
install.packages("reshape2") 
library(reshape2)
library(scales)

QS_GENE_ballon_aa_cna <- read.csv("QS_GENE_ballon_data.CSV",row.names = 1)
QS_GENE_ballon_aa_cna
ggballoonplot(QS_GENE_ballon_aa_cna)

ggballoonplot(QS_GENE_ballon_aa_cna, fill = "value")+
  scale_fill_gradientn(colors = my_cols)

ggballoonplot(QS_GENE_ballon_aa_cna, 
              fill = "value", #气球填充颜色
              ggtheme = theme_bw(),#画板主题
              size = "value",#气球大小
              size.range = c(1,5),
              color = "grey",#气球边框颜色
              shape = 23,#shape可以改变显示形状
              show.label = F)+#是否显示标签
  scale_fill_viridis_c(option = "C")+
  guides(size = FALSE)+#气球图例是否显示
  #设置颜色
  scale_fill_gradientn(
    colors=c("blue","white","red"),
    values=rescale(c(-15.5,0,5.5)),
    limits=c(-15.5,7.51))+
  labs(x="Gene status (AA mutation + CNV)",y="Cancer type ( p < 0.05)",fill="Difference of Quiescence Score 
(MUT/AMP/DEL - WT)")



## Exp - Cancer type 一种种替换 eg.“ov”全部替换为“ov”
ggboxplot(uvm_plppr4, x = "PLPPR4_aamutation_cna", y = "PLPPR4_log2expression",
          color = "PLPPR4_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("DEL","AMP")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("AMP","MUT")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("DEL","MUT")),label = "p.signif")

#####QS - PLPPR4_aamutation_cna + TP53_AAmutation - Cancer type#####
uvm_plppr4<- df.merged_plppr4_2[df.merged_plppr4_2$STUDY_ID=="uvm_tcga_pan_can_atlas_2018",]
ggboxplot(uvm_plppr4, x = "PLPPR4_aamutation_cna_TP53_aamutation", y = "QuiescenceScore",
          color = "PLPPR4_aamutation_cna_TP53_aamutation", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
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

df.merged_plppr4_3 <- merge(df.merged_plppr4_2[,c("SAMPLE_ID",
                                                "PLPPR4_TP53_status",
                                                "TP53_status","PLPPR4_status","PLPPR4_log2expression",
                                                "QuiescenceScore","Group",
                                                "PLPPR4_PLPPR4CNA",
                                                "STUDY_ID",
                                                "PLPPR4_CNA",
                                                "PLPPR4_aamutation_cna",
                                                "PLPPR4_aamutation_cna_TP53_aamutation")], 
                           aneuploidy[,c("Sample","Type","AneuploidyScore",
                                         "AS_del","AS_amp",
                                         "Genome_doublings","Leuk","Purity","Stroma",
                                         "Stroma_notLeukocyte","Stroma_notLeukocyte_Floor",
                                         "SilentMutationspeMb",
                                         "Non.silentMutationsperMb")],
                           by.x="SAMPLE_ID", by.y="Sample", #x和y表合并的依据是相同的sample id
                           all.x=FALSE, all.y=FALSE) 
head(df.merged_plppr4_3)

##### AS - PLPPR4_aamutation #####
table(df.merged_plppr4_3$PLPPR4_status)
ggboxplot(df.merged_plppr4_3, x = "PLPPR4_status", y = "AneuploidyScore",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### AS_amp - PLPPR4_aamutation #####
table(df.merged_plppr4_3$PLPPR4_status)
ggboxplot(df.merged_plppr4_3, x = "PLPPR4_status", y = "AS_amp",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### AS_del - PLPPR4_aamutation #####
table(df.merged_plppr4_3$PLPPR4_status)
ggboxplot(df.merged_plppr4_3, x = "PLPPR4_status", y = "AS_del",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Gb - PLPPR4_aamutation #####
table(df.merged_plppr4_3$PLPPR4_status)
ggboxplot(df.merged_plppr4_3, x = "PLPPR4_status", y = "Genome_doublings",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Leuk - PLPPR4_aamutation #####
df.merged_plppr4_3$Leuk[which(df.merged_plppr4_3$Leuk =='#N/A')] <- 'NA'
class(df.merged_plppr4_3$Leuk) #检查Leuk的数据类型，原数据框中为字符
df.merged_plppr4_3$Leuk <- as.numeric(df.merged_plppr4_3$Leuk)#使用as.numeric将其数据类型换为数值型
class(df.merged_plppr4_3$Leuk)#再次检查
hist(df.merged_plppr4_3$Leuk)

ggboxplot(df.merged_plppr4_3, x = "PLPPR4_status", y = "Leuk",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Purity - PLPPR4_aamutation #####
df.merged_plppr4_3$Purity[which(df.merged_plppr4_3$Purity =='#N/A')] <- 'NA'
class(df.merged_plppr4_3$Purity) #检查Leuk的数据类型，原数据框中为字符
df.merged_plppr4_3$Purity <- as.numeric(df.merged_plppr4_3$Purity)#使用as.numeric将其数据类型换为数值型
class(df.merged_plppr4_3$Purity)#再次检查
hist(df.merged_plppr4_3$Purity)

ggboxplot(df.merged_plppr4_3, x = "PLPPR4_status", y = "Purity",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Stroma - PLPPR4_aamutation #####
df.merged_plppr4_3$Stroma[which(df.merged_plppr4_3$Stroma =='#N/A')] <- 'NA'
class(df.merged_plppr4_3$Stroma) #检查Leuk的数据类型，原数据框中为字符
df.merged_plppr4_3$Stroma <- as.numeric(df.merged_plppr4_3$Stroma)#使用as.numeric将其数据类型换为数值型
class(df.merged_plppr4_3$Stroma)#再次检查
hist(df.merged_plppr4_3$Stroma)

ggboxplot(df.merged_plppr4_3, x = "PLPPR4_status", y = "Stroma",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Stroma_notLeukocyte - PLPPR4_aamutation #####
df.merged_plppr4_3$Stroma_notLeukocyte[which(df.merged_plppr4_3$Stroma_notLeukocyte =='#N/A')] <- 'NA'
class(df.merged_plppr4_3$Stroma_notLeukocyte) #检查Leuk的数据类型，原数据框中为字符
df.merged_plppr4_3$Stroma_notLeukocyte <- as.numeric(df.merged_plppr4_3$Stroma_notLeukocyte)#使用as.numeric将其数据类型换为数值型
class(df.merged_plppr4_3$Stroma_notLeukocyte)#再次检查
hist(df.merged_plppr4_3$Stroma_notLeukocyte)

ggboxplot(df.merged_plppr4_3, x = "PLPPR4_status", y = "Stroma_notLeukocyte",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### SilentMutationspeMb - PLPPR4_aamutation #####
df.merged_plppr4_3$SilentMutationspeMb[which(df.merged_plppr4_3$SilentMutationspeMb =='#N/A')] <- 'NA'
class(df.merged_plppr4_3$SilentMutationspeMb) #检查Leuk的数据类型，原数据框中为字符
df.merged_plppr4_3$SilentMutationspeMb <- as.numeric(df.merged_plppr4_3$SilentMutationspeMb)#使用as.numeric将其数据类型换为数值型
class(df.merged_plppr4_3$SilentMutationspeMb)#再次检查
hist(df.merged_plppr4_3$SilentMutationspeMb)

ggboxplot(df.merged_plppr4_3, x = "PLPPR4_status", y = "SilentMutationspeMb",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Non.silentMutationsperMb - PLPPR4_aamutation #####
df.merged_plppr4_3$Non.silentMutationsperMb[which(df.merged_plppr4_3$Non.silentMutationsperMb =='#N/A')] <- 'NA'
class(df.merged_plppr4_3$Non.silentMutationsperMb) #检查Leuk的数据类型，原数据框中为字符
df.merged_plppr4_3$Non.silentMutationsperMb <- as.numeric(df.merged_plppr4_3$Non.silentMutationsperMb)#使用as.numeric将其数据类型换为数值型
class(df.merged_plppr4_3$Non.silentMutationsperMb)#再次检查
hist(df.merged_plppr4_3$Non.silentMutationsperMb)

ggboxplot(df.merged_plppr4_3, x = "PLPPR4_status", y = "Non.silentMutationsperMb",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
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
ggplot(df.merged_plppr4_3, aes(x=AneuploidyScore, y=QuiescenceScore)) +
  geom_point(aes(col=Type)) +
  geom_smooth(method="loess", se=F) + 
  xlim(c(0, 40)) + 
  ylim(c(-30, 30)) + 
  labs(subtitle="QS Vs AS", 
       title="Scatterplot_PLPPR4")+
  stat_cor(method = 'pearson', aes(x = AneuploidyScore, y = QuiescenceScore))

ggplot(df.merged_plppr4_3, aes(x=AneuploidyScore, y=QuiescenceScore)) +
  geom_point(aes(col=Type)) +
  geom_smooth(method="lm", se=F) + 
  xlim(c(0, 40)) + 
  ylim(c(-30, 30)) + 
  labs(subtitle="QS Vs AS", 
       title="Scatterplot_PLPPR4")+
  stat_cor(method = 'pearson', aes(x = AneuploidyScore, y = QuiescenceScore))

ggplot(df.merged_plppr4_3, aes(x=AS_del, y=QuiescenceScore)) +
  geom_point(aes(col=Type)) +
  geom_smooth(method="lm", se=F) + 
  xlim(c(0, 40)) + 
  ylim(c(-30, 30)) + 
  labs(subtitle="QS Vs AS_del", 
       title="Scatterplot_PLPPR4")+
  stat_cor(method = 'pearson', aes(x = AS_del, y = QuiescenceScore))

ggplot(df.merged_plppr4_3, aes(x=AS_amp, y=QuiescenceScore)) +
  geom_point(aes(col=Type)) +
  geom_smooth(method="lm", se=F) + 
  xlim(c(0, 40)) + 
  ylim(c(-30, 30)) + 
  labs(subtitle="QS Vs AS_amp", 
       title="Scatterplot_PLPPR4")+
  stat_cor(method = 'pearson', aes(x = AS_amp, y = QuiescenceScore))

ggplot(data = df.merged_plppr4_3, aes(x = AneuploidyScore, y = QuiescenceScore,  color = PLPPR4_aamutation_cna)) + 
  geom_point() + geom_smooth(method = lm) +
  scale_color_manual(values = c('#FF7400', '#009999', '#3914AF',"#000000"))  +
  labs(title = 'QS Vs AS - PLPPR4_aamutation_cna') + 
  stat_cor(method = 'pearson', aes(x = AneuploidyScore, y = QuiescenceScore,  color = PLPPR4_aamutation_cna))
#####AS- PLPPR4_aamutation_cna#####
aggregate(df.merged_plppr4_3$AneuploidyScore, by=list(type=df.merged_plppr4_3$PLPPR4_aamutation_cna),mean)

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_plppr4_3, x = "PLPPR4_aamutation_cna", y = "AneuploidyScore",
          color = "PLPPR4_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#Ballon plot
my_cols <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
             "#E16462FF", "#FCA636FF", "#F0F921FF")
AS_ballon_aa_cna <- read.csv("AS_ballon_data.csv",row.names = 1)
AS_ballon_aa_cna
ggballoonplot(AS_ballon_aa_cna, fill = "value")+
     scale_fill_gradientn(colors = my_cols)

#####AS_amp- PLPPR4_aamutation_cna + Ballon plot#####
#mean of groups
aggregate(df.merged_plppr4_3$AS_amp, by=list(type=df.merged_plppr4_3$PLPPR4_aamutation_cna),mean)

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_plppr4_3, x = "PLPPR4_aamutation_cna", y = "AS_amp",
          color = "PLPPR4_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#Ballon plot
AS_amp_ballon_aa_cna <- read.csv("AS_amp_ballon_data.csv",row.names = 1)
AS_amp_ballon_aa_cna
ggballoonplot(AS_amp_ballon_aa_cna, fill = "value")+
  scale_fill_gradientn(colors = my_cols)
#####AS_del- PLPPR4_aamutation_cna#####
#mean of groups
aggregate(df.merged_plppr4_3$AS_del, by=list(type=df.merged_plppr4_3$PLPPR4_aamutation_cna),mean)

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_plppr4_3, x = "PLPPR4_aamutation_cna", y = "AS_del",
          color = "PLPPR4_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#Ballon plot
AS_del_ballon_aa_cna <- read.csv("AS_del_ballon_data.csv",row.names = 1)
AS_del_ballon_aa_cna
ggballoonplot(AS_del_ballon_aa_cna, fill = "value")+
  scale_fill_gradientn(colors = my_cols)
#####Gb- PLPPR4_aamutation_cna#####
#mean of groups
aggregate(df.merged_plppr4_3$Genome_doublings, by=list(type=df.merged_plppr4_3$PLPPR4_aamutation_cna),mean)

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_plppr4_3, x = "PLPPR4_aamutation_cna", y = "Genome_doublings",
          color = "PLPPR4_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#Ballon plot
GB_ballon_aa_cna <- read.csv("GB_ballon_data.csv",row.names = 1)
GB_ballon_aa_cna
ggballoonplot(GB_ballon_aa_cna, fill = "value")+
  scale_fill_gradientn(colors = my_cols)
#####Leuk- PLPPR4_aamutation_cna#####
#mean of groups
aggregate(df.merged_plppr4_3$Leuk, by=list(type=df.merged_plppr4_3$PLPPR4_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_plppr4_3, x = "PLPPR4_aamutation_cna", y = "Leuk",
          color = "PLPPR4_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#Ballon plot
Leuk_ballon_aa_cna <- read.csv("Leuk_ballon_data.csv",row.names = 1)
Leuk_ballon_aa_cna
ggballoonplot(Leuk_ballon_aa_cna, fill = "value")+
  scale_fill_gradientn(colors = my_cols)
#####Purity- PLPPR4_aamutation_cna#####
#mean of groups
aggregate(df.merged_plppr4_3$Purity, by=list(type=df.merged_plppr4_3$PLPPR4_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_plppr4_3, x = "PLPPR4_aamutation_cna", y = "Purity",
          color = "PLPPR4_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#Ballon plot
Purity_ballon_aa_cna <- read.csv("Purity_ballon_data.csv",row.names = 1)
Purity_ballon_aa_cna
ggballoonplot(Purity_ballon_aa_cna, fill = "value")+
  scale_fill_gradientn(colors = my_cols)
#####Stroma- PLPPR4_aamutation_cna#####
#mean of groups
aggregate(df.merged_plppr4_3$Stroma, by=list(type=df.merged_plppr4_3$PLPPR4_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_plppr4_3, x = "PLPPR4_aamutation_cna", y = "Stroma",
          color = "PLPPR4_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#Ballon plot
Stroma_ballon_aa_cna <- read.csv("Stroma_ballon_data.csv",row.names = 1)
Stroma_ballon_aa_cna
ggballoonplot(Stroma_ballon_aa_cna, fill = "value")+
  scale_fill_gradientn(colors = my_cols)
#####Stroma_notLeukocyte- PLPPR4_aamutation_cna#####
#mean of groups
aggregate(df.merged_plppr4_3$Stroma_notLeukocyte, by=list(type=df.merged_plppr4_3$PLPPR4_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_plppr4_3, x = "PLPPR4_aamutation_cna", y = "Stroma_notLeukocyte",
          color = "PLPPR4_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#Ballon plot
Stroma_notLeukocyte_ballon_aa_cna <- read.csv("Stroma_notLeukocyte_ballon_data.csv",row.names = 1)
Stroma_notLeukocyte_ballon_aa_cna
ggballoonplot(Stroma_notLeukocyte_ballon_aa_cna, fill = "value")+
  scale_fill_gradientn(colors = my_cols)
#####SilentMutationspeMb- PLPPR4_aamutation_cna#####
my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_plppr4_3, x = "PLPPR4_aamutation_cna", y = "SilentMutationspeMb",
          color = "PLPPR4_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####Non.silentMutationsperMb- PLPPR4_aamutation_cna#####
my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_plppr4_3, x = "PLPPR4_aamutation_cna", y = "Non.silentMutationsperMb",
          color = "PLPPR4_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")

#####导入样本总突变数据NM---已下为突变数所有代码 可直接复制后分别替换基因大小写#####
pgen <- read.csv("pgen.csv")
head(pgen)
hist(pgen$number_of_mutations)

df.merged_plppr4_4 <- merge(df.merged_plppr4_3[,c("SAMPLE_ID",
                                                  "PLPPR4_TP53_status",
                                                  "TP53_status","PLPPR4_status","PLPPR4_log2expression",
                                                  "QuiescenceScore","Group",
                                                  "PLPPR4_PLPPR4CNA",
                                                  "STUDY_ID",
                                                  "PLPPR4_CNA",
                                                  "PLPPR4_aamutation_cna",
                                                  "PLPPR4_aamutation_cna_TP53_aamutation",
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
head(df.merged_plppr4_4)
#####Number of Mutations- PLPPR4_aamutation_cna#####
#mean of groups
aggregate(df.merged_plppr4_4$number_of_mutations, by=list(type=df.merged_plppr4_4$PLPPR4_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_plppr4_4, x = "PLPPR4_aamutation_cna", y = "number_of_mutations",
          color = "PLPPR4_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3)
#####Ballon plot-NM_aa_cna#####
NM_ballon_aa_cna <- read.csv("NM_ballon_data.csv",row.names = 1)
NM_ballon_aa_cna
ggballoonplot(NM_ballon_aa_cna, fill = "value")+
  scale_fill_gradientn(colors = my_cols)
