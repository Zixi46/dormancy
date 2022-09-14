### This script analyses GNAQ mutation data and links with quiescence and expressin levels in tumours.
##### 初步数据处理 aa mutaion status与Exp QS合并#####
### Read mutation data downloaded from cbioportal:
setwd("/Users/lizixi/Desktop/Final Project/5.18 task")
mutations_gnaq <- read.delim("mutations_cbioportal_gnaq.txt")
head(mutations_gnaq)
## How many samples have a mutation?
table(mutations_gnaq_gnaq$GNAQ)
tb_gnaq <- table(mutations_gnaq$GNAQ)
sort(tb_gnaq) #descend升序排列PHIN1列的数字元素

# Add a new column that just denotes WT and MUT for GNAQ:
mutations_gnaq[which(mutations_gnaq$GNAQ == "NS"),]$GNAQ <- "WT" #****将GNAQ中的NS替换为WT
mutations_gnaq$GNAQ_status <- mutations_gnaq$GNAQ #为mutations_gnaq表添加一列为GNAQ_status，同时将这一列定义为GNAQ（和GNAQ一样）
mutations_gnaq[which(mutations_gnaq$GNAQ != "WT"),]$GNAQ_status <- "MUT" #选出mutation中GNAQ列中不等于“WT”的元素，并在GNAQ_status中将其全赋值为“MUT”
head(mutations_gnaq)
# How many samples with mutations_gnaq in GNAQ?
table(mutations_gnaq$GNAQ_status)
# Read expression data for GNAQ:
expression_gnaq <- read.delim("mRNA Expression_RSEM_cbioportal_gnaq.txt")
head(expression_gnaq)
# Plot histogram of expression values of GNAQ:
hist(expression_gnaq$GNAQ)
# take the log2 of expression values to avoid working with very large numbers:
expression_gnaq$GNAQ_log2expression <- log2(expression_gnaq$GNAQ+1)
head(expression_gnaq)
hist(expression_gnaq$GNAQ_log2expression) # more decent value

# Merge GNAQ expression data with its mutation data:
df.merged_gnaq <- merge(mutations_gnaq[,c("SAMPLE_ID","GNAQ_status","STUDY_ID")], #mutation中需要合并的纵列是sample id，gnaq_status
                   expression_gnaq[,c("SAMPLE_ID","GNAQ_log2expression")], #expresssion中需要合并的纵列是sample ID，GNAQ_log2expression
                   by.x="SAMPLE_ID", by.y="SAMPLE_ID", #x和y表合并的依据是相同的sample id
                   all.x=FALSE, all.y=FALSE) #x和y表中不需要所有的行都合并（只合并sample id一样的，只在其中一个表中才存在的sample id丢弃）
head(df.merged_gnaq)
# Merge TP53 mutation data with QS data
df.merged_gnaq <- merge(df.merged_gnaq[,c("SAMPLE_ID","GNAQ_status","STUDY_ID",
                                          "GNAQ_log2expression")], #mutation中需要合并的纵列是sample id，gnaq_status
                        df.merged2[,c("SAMPLE_ID","TP53_status","QuiescenceScore",
                                      "Group")], #expresssion中需要合并的纵列是sample ID，GNAQ_log2expression
                        by.x="SAMPLE_ID", by.y="SAMPLE_ID", #x和y表合并的依据是相同的sample id
                        all.x=FALSE, all.y=FALSE) #x和y表中不需要所有的行都合并（只合并sample id一样的，只在其中一个表中才存在的sample id丢弃）
head(df.merged_gnaq)


###### 增加GNAQ突变类型分析列 “AMP”“DEL” #####
#read cna file
cna_GNAQ <- read.delim("cna_fromcBioportal_gnaq.txt")
names(cna_GNAQ)[3] <- "GNAQ_CNA"#给第三列重命名
head(cna_GNAQ)
table(cna_GNAQ$GNAQ_CNA)
#将cna文件中的ID 和 CNA与总表合并
head(df.merged_gnaq)
head(cna_GNAQ)
df.merged_gnaq2 <- merge(df.merged_gnaq[,c("SAMPLE_ID","STUDY_ID","GNAQ_status",
                                           "TP53_status","GNAQ_log2expression",
                                           "QuiescenceScore","Group")], 
                         cna_GNAQ[,c("SAMPLE_ID","GNAQ_CNA")],
                         by.x="SAMPLE_ID", by.y="SAMPLE_ID", #x和y表合并的依据是相同的sample id
                         all.x=FALSE, all.y=FALSE) 
head(df.merged_gnaq2)
table(df.merged_gnaq2$GNAQ_CNA)
#数据整理，完整案例分析
##exclude cns-NP(not profiled) samples from the analysis 剔除包含NP的样本
##把NP都换为NA
df.merged_gnaq2[which(df.merged_gnaq2$GNAQ_CNA == "NP"),]$GNAQ_CNA <- "NA" #****将GNAQ_CAN中的NP替换为NA
table(df.merged_gnaq2$GNAQ_CNA)
##把-2换为-1（数量太少了），把2换为1
#"-2" is a deep loss, possibly a homozygous deletion
#"-1" is a single-copy loss (heterozygous deletion)
#"0" is diploid
#"1" indicates a low-level gain
#"2" is a high-level amplification.
##将1转化为AMP,-1转化为DEL，O转化为NCN
df.merged_gnaq2[which(df.merged_gnaq2$GNAQ_CNA == "-2"),]$GNAQ_CNA <- "-1" 
df.merged_gnaq2[which(df.merged_gnaq2$GNAQ_CNA == "2"),]$GNAQ_CNA <- "1" 
df.merged_gnaq2[which(df.merged_gnaq2$GNAQ_CNA == "-1"),]$GNAQ_CNA <- "DEL" 
df.merged_gnaq2[which(df.merged_gnaq2$GNAQ_CNA == "1"),]$GNAQ_CNA <- "AMP" 
df.merged_gnaq2[which(df.merged_gnaq2$GNAQ_CNA == "0"),]$GNAQ_CNA <- "NCN" 
table(df.merged_gnaq2$GNAQ_CNA)
##去掉NA样本
df.merged_gnaq2 <- df.merged_gnaq2[!grepl("NA", df.merged_gnaq2$GNAQ_CNA),]
table(df.merged_gnaq2$GNAQ_CNA)
##将GNAQ AAmutation信息和CNA信息结合为新的一列Mutation_CNA,生成6个不同组别
df.merged_gnaq2 <- tidyr::unite(df.merged_gnaq2, "GNAQ_GNAQCNA", GNAQ_status, GNAQ_CNA, remove = FALSE)
head(df.merged_gnaq2)
table(df.merged_gnaq2$GNAQ_GNAQCNA)
#重排顺序
df.merged_gnaq2 <- df.merged_gnaq2[,c("SAMPLE_ID","STUDY_ID","GNAQ_status","GNAQ_CNA","GNAQ_GNAQCNA","TP53_status","GNAQ_log2expression","QuiescenceScore","Group")]
table(df.merged_gnaq$GNAQ_status)
table(df.merged_gnaq$TP53_status)
#根据aamutaiton+CNA进行重新分组
df.merged_gnaq2$GNAQ_aamutation_cna <- df.merged_gnaq2$GNAQ_GNAQCNA #添加新列
df.merged_gnaq2[which(df.merged_gnaq2$GNAQ_aamutation_cna == "WT_NCN"),]$GNAQ_aamutation_cna <- "WT"# WT-NCN → WT
df.merged_gnaq2[which(df.merged_gnaq2$GNAQ_aamutation_cna == "MUT_AMP"),]$GNAQ_aamutation_cna <- "MUT"# MUT-AMP → MUT
df.merged_gnaq2[which(df.merged_gnaq2$GNAQ_aamutation_cna == "MUT_NCN"),]$GNAQ_aamutation_cna <- "MUT"# WUT-NCN → MUT
df.merged_gnaq2[which(df.merged_gnaq2$GNAQ_aamutation_cna == "MUT_DEL"),]$GNAQ_aamutation_cna <- "DEL"# MUT-DEL → DEL
df.merged_gnaq2[which(df.merged_gnaq2$GNAQ_aamutation_cna == "WT_AMP"),]$GNAQ_aamutation_cna <- "AMP"# WT-AMP → AMP
df.merged_gnaq2[which(df.merged_gnaq2$GNAQ_aamutation_cna == "WT_DEL"),]$GNAQ_aamutation_cna <- "DEL"# WT-DEL → DEL
table(df.merged_gnaq2$GNAQ_aamutation_cna)

#####合并 GNAQ+TP53 status and GNAQ+TP53 status #####
#Combine GNAQ+TP53 to create 4 groups 。用tidyr::unite合并两列（GNAQ_status，TP53_status）的元素为新的一列（GNAQ_TP53_status），默认以_连接（可以修改）
library(tidyr)
df.merged_gnaq2 <- tidyr::unite(df.merged_gnaq2, "GNAQ_TP53_status", GNAQ_status, TP53_status, remove = FALSE)
head(df.merged_gnaq2)
table(df.merged_gnaq2$GNAQ_TP53_status)

#####合并 GNAQ_aamutation_cna+TP53_status#####
df.merged_gnaq2 <- tidyr::unite(df.merged_gnaq2, "GNAQ_aamutation_cna_TP53_aamutation", GNAQ_aamutation_cna, TP53_status, remove = FALSE)
table(df.merged_gnaq2$GNAQ_aamutation_cna_TP53_aamutation)

#####Expression and QS levels - GNAQ_aamutation#####
#QS
## by boxplot
ggboxplot(df.merged_gnaq, x = "GNAQ_status", y = "QuiescenceScore",
          color = "GNAQ_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")
ggboxplot(df.merged_gnaq, x = "GNAQ_status", y = "QuiescenceScore",
          color = "GNAQ_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)
#Exp
## by boxplot
ggboxplot(df.merged_gnaq, x = "GNAQ_status", y = "GNAQ_log2expression",
          color = "GNAQ_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")
ggboxplot(df.merged_gnaq, x = "GNAQ_status", y = "GNAQ_log2expression",
          color = "GNAQ_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)
#####Expression and QS levels - GNAQ_aamutation_cna#####
#QS
## by boxplot
ggboxplot(df.merged_gnaq2, x = "GNAQ_aamutation_cna", y = "QuiescenceScore",
          color = "GNAQ_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL"))+
  stat_compare_means(comparisons = my_comparisons3,label="p.signif")
ggboxplot(df.merged_gnaq2, x = "GNAQ_aamutation_cna", y = "QuiescenceScore",
          color = "GNAQ_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL"))+
  stat_compare_means(comparisons = my_comparisons3)
#Exp
## by boxplot
ggboxplot(df.merged_gnaq2, x = "GNAQ_aamutation_cna", y = "GNAQ_log2expression",
          color = "GNAQ_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL"))+
  stat_compare_means(comparisons = my_comparisons3,label="p.signif")
ggboxplot(df.merged_gnaq2, x = "GNAQ_aamutation_cna", y = "GNAQ_log2expression",
          color = "GNAQ_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL"))+
  stat_compare_means(comparisons = my_comparisons3)

#####Quiescence group 组间比例分布图 by GNAQ_aamutation #####
table(df.merged_gnaq$GNAQ_status) #看看wt和mut的数量，差距比较大
ggplot(df.merged_gnaq, mapping = aes(
  x = GNAQ_status, fill = Group), order = c("WT", "MUT"))+
  geom_bar(stat="count",width=0.5,position='fill')+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  geom_text(stat='count',aes(label=scales::percent(..count../sum(..count..))),
            color="white", size=3.5,position=position_fill(0.5))+
  theme_minimal()

##### QS - GNAQ_aamutation_TP53_aamutation #####
## by boxplot
my_comparisons2 <- list(c("MUT_WT","MUT_MUT"),c("WT_MUT","MUT_MUT"),c("WT_WT","WT_MUT"),c("WT_MUT","MUT_WT"),c("WT_WT","MUT_MUT"),c("WT_WT","MUT_WT")) #显示4组之间的p
ggboxplot(df.merged_gnaq2, x = "GNAQ_TP53_status", y = "QuiescenceScore",
          color = "GNAQ_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons2,label = "p.signif")
#####QS - GNAQ_aamutation_cna + TP53_AAmutation #####
aggregate(df.merged_gnaq2$QuiescenceScore, by=list(type=df.merged_gnaq2$GNAQ_aamutation_cna),mean)
ggboxplot(df.merged_gnaq2, x = "GNAQ_aamutation_cna_TP53_aamutation", y = "QuiescenceScore",
          color = "GNAQ_aamutation_cna_TP53_aamutation", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB"),
          add = "jitter",order = c("WT_WT","WT_MUT","MUT_WT","MUT_MUT","AMP_WT", "AMP_MUT","DEL_WT","DEL_MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("AMP_WT","AMP_MUT")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("DEL_WT","DEL_MUT")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("MUT_WT","MUT_MUT")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("WT_WT","WT_MUT")),label = "p.signif")
ggboxplot(df.merged_gnaq2, x = "GNAQ_aamutation_cna_TP53_aamutation", y = "QuiescenceScore",
          color = "GNAQ_aamutation_cna_TP53_aamutation", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB"),
          add = "jitter",order = c("WT_WT","WT_MUT","MUT_WT","MUT_MUT","AMP_WT", "AMP_MUT","DEL_WT","DEL_MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("AMP_WT","AMP_MUT")))+
  stat_compare_means(comparisons = list(c("DEL_WT","DEL_MUT")))+
  stat_compare_means(comparisons = list(c("MUT_WT","MUT_MUT")))+
  stat_compare_means(comparisons = list(c("WT_WT","WT_MUT")))

#####Expression and QS levels - GNAQ_aamutation_cna - Cancer type#####
#把每种癌症都分为一个df
#再分类绘制QS & Expression plot（boxplot 感觉violinplot没有意义）

## QS - Cancer type 一种种替换 eg.“ov”全部替换为“ov”
#ACC
acc_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="acc_tcga_pan_can_atlas_2018",]
table(acc_gnaq$GNAQ_aamutation_cna)
aggregate(acc_gnaq$QuiescenceScore, by=list(type=acc_gnaq$GNAQ_aamutation_cna),mean)
#BLCA
blca_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="blca_tcga_pan_can_atlas_2018",]
table(blca_gnaq$GNAQ_aamutation_cna)
aggregate(blca_gnaq$QuiescenceScore, by=list(type=blca_gnaq$GNAQ_aamutation_cna),mean)
#BRCA
brca_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="brca_tcga_pan_can_atlas_2018",]
table(brca_gnaq$GNAQ_aamutation_cna)
aggregate(brca_gnaq$QuiescenceScore, by=list(type=brca_gnaq$GNAQ_aamutation_cna),mean)
#CESC
cesc_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="cesc_tcga_pan_can_atlas_2018",]
table(cesc_gnaq$GNAQ_aamutation_cna)
aggregate(cesc_gnaq$QuiescenceScore, by=list(type=cesc_gnaq$GNAQ_aamutation_cna),mean)
#CHOL
chol_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="chol_tcga_pan_can_atlas_2018",]
table(chol_gnaq$GNAQ_aamutation_cna)
aggregate(chol_gnaq$QuiescenceScore, by=list(type=chol_gnaq$GNAQ_aamutation_cna),mean)
#COAD
coadread_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="coadread_tcga_pan_can_atlas_2018",]
table(coadread_gnaq$GNAQ_aamutation_cna)
aggregate(coadread_gnaq$QuiescenceScore, by=list(type=coadread_gnaq$GNAQ_aamutation_cna),mean)
#ESCA
esca_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="esca_tcga_pan_can_atlas_2018",]
table(esca_gnaq$GNAQ_aamutation_cna)
aggregate(esca_gnaq$QuiescenceScore, by=list(type=esca_gnaq$GNAQ_aamutation_cna),mean)
#GBM
gbm_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="gbm_tcga_pan_can_atlas_2018",]
table(gbm_gnaq$GNAQ_aamutation_cna)
aggregate(gbm_gnaq$QuiescenceScore, by=list(type=gbm_gnaq$GNAQ_aamutation_cna),mean)
#HNSC
hnsc_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="hnsc_tcga_pan_can_atlas_2018",]
table(hnsc_gnaq$GNAQ_aamutation_cna)
aggregate(hnsc_gnaq$QuiescenceScore, by=list(type=hnsc_gnaq$GNAQ_aamutation_cna),mean)
#KICH
kich_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="kich_tcga_pan_can_atlas_2018",]
table(kich_gnaq$GNAQ_aamutation_cna)
aggregate(kich_gnaq$QuiescenceScore, by=list(type=kich_gnaq$GNAQ_aamutation_cna),mean)
#KIRC
kirc_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="kirc_tcga_pan_can_atlas_2018",]
table(kirc_gnaq$GNAQ_aamutation_cna)
aggregate(kirc_gnaq$QuiescenceScore, by=list(type=kirc_gnaq$GNAQ_aamutation_cna),mean)
#KIRP
kirp_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="kirp_tcga_pan_can_atlas_2018",]
table(kirp_gnaq$GNAQ_aamutation_cna)
aggregate(kirp_gnaq$QuiescenceScore, by=list(type=kirp_gnaq$GNAQ_aamutation_cna),mean)
#LGG
lgg_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="lgg_tcga_pan_can_atlas_2018",]
table(lgg_gnaq$GNAQ_aamutation_cna)
aggregate(lgg_gnaq$QuiescenceScore, by=list(type=lgg_gnaq$GNAQ_aamutation_cna),mean)
#LIHC
lihc_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="lihc_tcga_pan_can_atlas_2018",]
table(lihc_gnaq$GNAQ_aamutation_cna)
aggregate(lihc_gnaq$QuiescenceScore, by=list(type=lihc_gnaq$GNAQ_aamutation_cna),mean)
#LUAD
luad_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="luad_tcga_pan_can_atlas_2018",]
table(luad_gnaq$GNAQ_aamutation_cna)
aggregate(luad_gnaq$QuiescenceScore, by=list(type=luad_gnaq$GNAQ_aamutation_cna),mean)
#LUSC
lusc_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="lusc_tcga_pan_can_atlas_2018",]
table(lusc_gnaq$GNAQ_aamutation_cna)
aggregate(lusc_gnaq$QuiescenceScore, by=list(type=lusc_gnaq$GNAQ_aamutation_cna),mean)
#MESO
meso_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="meso_tcga_pan_can_atlas_2018",]
table(meso_gnaq$GNAQ_aamutation_cna)
aggregate(meso_gnaq$QuiescenceScore, by=list(type=meso_gnaq$GNAQ_aamutation_cna),mean)
#OV
ov_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="ov_tcga_pan_can_atlas_2018",]
table(ov_gnaq$GNAQ_aamutation_cna)
aggregate(ov_gnaq$QuiescenceScore, by=list(type=ov_gnaq$GNAQ_aamutation_cna),mean)
#PAAD
paad_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="paad_tcga_pan_can_atlas_2018",]
table(paad_gnaq$GNAQ_aamutation_cna)
aggregate(paad_gnaq$QuiescenceScore, by=list(type=paad_gnaq$GNAQ_aamutation_cna),mean)
#PCPG
pcpg_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="pcpg_tcga_pan_can_atlas_2018",]
table(pcpg_gnaq$GNAQ_aamutation_cna)
aggregate(pcpg_gnaq$QuiescenceScore, by=list(type=pcpg_gnaq$GNAQ_aamutation_cna),mean)
#PRAD
prad_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="prad_tcga_pan_can_atlas_2018",]
table(prad_gnaq$GNAQ_aamutation_cna)
aggregate(prad_gnaq$QuiescenceScore, by=list(type=prad_gnaq$GNAQ_aamutation_cna),mean)
#SARC
sarc_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="sarc_tcga_pan_can_atlas_2018",]
table(sarc_gnaq$GNAQ_aamutation_cna)
aggregate(sarc_gnaq$QuiescenceScore, by=list(type=sarc_gnaq$GNAQ_aamutation_cna),mean)
#STAD
stad_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="stad_tcga_pan_can_atlas_2018",]
table(stad_gnaq$GNAQ_aamutation_cna)
aggregate(stad_gnaq$QuiescenceScore, by=list(type=stad_gnaq$GNAQ_aamutation_cna),mean)
#TGCT
tgct_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="tgct_tcga_pan_can_atlas_2018",]
table(tgct_gnaq$GNAQ_aamutation_cna)
aggregate(tgct_gnaq$QuiescenceScore, by=list(type=tgct_gnaq$GNAQ_aamutation_cna),mean)
#THCA
thca_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="thca_tcga_pan_can_atlas_2018",]
table(thca_gnaq$GNAQ_aamutation_cna)
aggregate(thca_gnaq$QuiescenceScore, by=list(type=thca_gnaq$GNAQ_aamutation_cna),mean)
#THYM
thym_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="thym_tcga_pan_can_atlas_2018",]
table(thym_gnaq$GNAQ_aamutation_cna)
aggregate(thym_gnaq$QuiescenceScore, by=list(type=thym_gnaq$GNAQ_aamutation_cna),mean)
#UCEC
ucec_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="ucec_tcga_pan_can_atlas_2018",]
table(ucec_gnaq$GNAQ_aamutation_cna)
aggregate(ucec_gnaq$QuiescenceScore, by=list(type=ucec_gnaq$GNAQ_aamutation_cna),mean)
#UCS
ucs_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="ucs_tcga_pan_can_atlas_2018",]
table(ucs_gnaq$GNAQ_aamutation_cna)
aggregate(ucs_gnaq$QuiescenceScore, by=list(type=ucs_gnaq$GNAQ_aamutation_cna),mean)
#UVM
uvm_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="uvm_tcga_pan_can_atlas_2018",]
table(uvm_gnaq$GNAQ_aamutation_cna)
aggregate(uvm_gnaq$QuiescenceScore, by=list(type=uvm_gnaq$GNAQ_aamutation_cna),mean)


ggboxplot(ov_gnaq, x = "GNAQ_aamutation_cna", y = "QuiescenceScore",
          color = "GNAQ_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("MUT","AMP")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("MUT","DEL")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("AMP","DEL")),label = "p.signif")
## Exp - Cancer type 一种种替换 eg.“ov”全部替换为“ov”
head(uvm_gnaq)
ggboxplot(uvm_gnaq, x = "GNAQ_aamutation_cna", y = "GNAQ_log2expression",
          color = "GNAQ_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("MUT","AMP")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("MUT","DEL")),label = "p.signif")+
  stat_compare_means(comparisons = list(c("AMP","DEL")),label = "p.signif")

#####QS - GNAQ_aamutation_cna + TP53_AAmutation - Cancer type#####

uvm_gnaq<- df.merged_gnaq2[df.merged_gnaq2$STUDY_ID=="uvm_tcga_pan_can_atlas_2018",]
ggboxplot(uvm_gnaq, x = "GNAQ_aamutation_cna_TP53_aamutation", y = "QuiescenceScore",
          color = "GNAQ_aamutation_cna_TP53_aamutation", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
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

df.merged_gnaq_3 <- merge(df.merged_gnaq2[,c("SAMPLE_ID",
                                    "GNAQ_TP53_status",
                                    "TP53_status","GNAQ_status","GNAQ_log2expression",
                                    "QuiescenceScore","Group",
                                    "GNAQ_GNAQCNA",
                                    "STUDY_ID",
                                    "GNAQ_CNA",
                                    "GNAQ_aamutation_cna",
                                    "GNAQ_aamutation_cna_TP53_aamutation")], 
                     aneuploidy[,c("Sample","Type","AneuploidyScore",
                                   "AS_del","AS_amp",
                                   "Genome_doublings","Leuk","Purity","Stroma",
                                   "Stroma_notLeukocyte","Stroma_notLeukocyte_Floor",
                                   "SilentMutationspeMb",
                                   "Non.silentMutationsperMb")],
                     by.x="SAMPLE_ID", by.y="Sample", #x和y表合并的依据是相同的sample id
                     all.x=FALSE, all.y=FALSE) 
head(df.merged_gnaq_3)

##### AS - GNAQ_aamutation #####
table(df.merged_gnaq_3$GNAQ_status)
ggboxplot(df.merged_gnaq_3, x = "GNAQ_status", y = "AneuploidyScore",
          color = "GNAQ_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### AS_amp - GNAQ_aamutation #####
table(df.merged_gnaq_3$GNAQ_status)
ggboxplot(df.merged_gnaq_3, x = "GNAQ_status", y = "AS_amp",
          color = "GNAQ_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### AS_del - GNAQ_aamutation #####
table(df.merged_gnaq_3$GNAQ_status)
ggboxplot(df.merged_gnaq_3, x = "GNAQ_status", y = "AS_del",
          color = "GNAQ_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Gb - GNAQ_aamutation #####
table(df.merged_gnaq_3$GNAQ_status)
ggboxplot(df.merged_gnaq_3, x = "GNAQ_status", y = "Genome_doublings",
          color = "GNAQ_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Leuk - GNAQ_aamutation #####
df.merged_gnaq_3$Leuk[which(df.merged_gnaq_3$Leuk =='#N/A')] <- 'NA'
class(df.merged_gnaq_3$Leuk) #检查Leuk的数据类型，原数据框中为字符
df.merged_gnaq_3$Leuk <- as.numeric(df.merged_gnaq_3$Leuk)#使用as.numeric将其数据类型换为数值型
class(df.merged_gnaq_3$Leuk)#再次检查
hist(df.merged_gnaq_3$Leuk)

ggboxplot(df.merged_gnaq_3, x = "GNAQ_status", y = "Leuk",
          color = "GNAQ_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Purity - GNAQ_aamutation #####
df.merged_gnaq_3$Purity [which(df.merged_gnaq_3$Purity  =='#N/A')] <- 'NA'
class(df.merged_gnaq_3$Purity) #检查数据类型，原数据框中为字符
df.merged_gnaq_3$Purity  <- as.numeric(df.merged_gnaq_3$Purity )#使用as.numeric将其数据类型换为数值型
class(df.merged_gnaq_3$Purity )#再次检查
hist(df.merged_gnaq_3$Purity )

ggboxplot(df.merged_gnaq_3, x = "GNAQ_status", y = "Purity",
          color = "GNAQ_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Stroma - GNAQ_aamutation #####
df.merged_gnaq_3$Stroma [which(df.merged_gnaq_3$Stroma  =='#N/A')] <- 'NA'
class(df.merged_gnaq_3$Stroma) #检查数据类型，原数据框中为字符
df.merged_gnaq_3$Stroma  <- as.numeric(df.merged_gnaq_3$Stroma)#使用as.numeric将其数据类型换为数值型
class(df.merged_gnaq_3$Stroma)#再次检查
hist(df.merged_gnaq_3$Stroma)

ggboxplot(df.merged_gnaq_3, x = "GNAQ_status", y = "Stroma",
          color = "GNAQ_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Stroma_notLeukocyte - GNAQ_aamutation #####
df.merged_gnaq_3$Stroma_notLeukocyte [which(df.merged_gnaq_3$Stroma_notLeukocyte  =='#N/A')] <- 'NA'
class(df.merged_gnaq_3$Stroma_notLeukocyte) #检查数据类型，原数据框中为字符
df.merged_gnaq_3$Stroma_notLeukocyte  <- as.numeric(df.merged_gnaq_3$Stroma_notLeukocyte)#使用as.numeric将其数据类型换为数值型
class(df.merged_gnaq_3$Stroma_notLeukocyte)#再次检查
hist(df.merged_gnaq_3$Stroma_notLeukocyte)

ggboxplot(df.merged_gnaq_3, x = "GNAQ_status", y = "Stroma_notLeukocyte",
          color = "GNAQ_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### SilentMutationspeMb - GNAQ_aamutation #####
df.merged_gnaq_3$SilentMutationspeMb [which(df.merged_gnaq_3$SilentMutationspeMb  =='#N/A')] <- 'NA'
class(df.merged_gnaq_3$SilentMutationspeMb) #检查数据类型，原数据框中为字符
df.merged_gnaq_3$SilentMutationspeMb  <- as.numeric(df.merged_gnaq_3$SilentMutationspeMb)#使用as.numeric将其数据类型换为数值型
class(df.merged_gnaq_3$SilentMutationspeMb)#再次检查
hist(df.merged_gnaq_3$SilentMutationspeMb)

ggboxplot(df.merged_gnaq_3, x = "GNAQ_status", y = "SilentMutationspeMb",
          color = "GNAQ_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = list(c("WT","MUT")),label = "p.signif")
##### Non.silentMutationsperMb - GNAQ_aamutation #####
df.merged_gnaq_3$Non.silentMutationsperMb [which(df.merged_gnaq_3$Non.silentMutationsperMb  =='#N/A')] <- 'NA'
class(df.merged_gnaq_3$Non.silentMutationsperMb) #检查数据类型，原数据框中为字符
df.merged_gnaq_3$Non.silentMutationsperMb  <- as.numeric(df.merged_gnaq_3$Non.silentMutationsperMb)#使用as.numeric将其数据类型换为数值型
class(df.merged_gnaq_3$Non.silentMutationsperMb)#再次检查
hist(df.merged_gnaq_3$Non.silentMutationsperMb)

ggboxplot(df.merged_gnaq_3, x = "GNAQ_status", y = "Non.silentMutationsperMb",
          color = "GNAQ_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
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
ggplot(df.merged_gnaq_3, aes(x=AneuploidyScore, y=QuiescenceScore)) +
  geom_point(aes(col=Type)) +
  geom_smooth(method="loess", se=F) + 
  xlim(c(0, 40)) + 
  ylim(c(-30, 30)) + 
  labs(subtitle="QS Vs AS", 
       title="Scatterplot_GNAQ")+
  stat_cor(method = 'pearson', aes(x = AneuploidyScore, y = QuiescenceScore))

ggplot(df.merged_gnaq_3, aes(x=AneuploidyScore, y=QuiescenceScore)) +
  geom_point(aes(col=Type)) +
  geom_smooth(method="lm", se=F) + 
  xlim(c(0, 40)) + 
  ylim(c(-30, 30)) + 
  labs(subtitle="QS Vs AS", 
       title="Scatterplot_GNAQ")+
  stat_cor(method = 'pearson', aes(x = AneuploidyScore, y = QuiescenceScore))

ggplot(df.merged_gnaq_3, aes(x=AS_del, y=QuiescenceScore)) +
  geom_point(aes(col=Type)) +
  geom_smooth(method="lm", se=F) + 
  xlim(c(0, 40)) + 
  ylim(c(-30, 30)) + 
  labs(subtitle="QS Vs AS_del", 
       title="Scatterplot_GNAQ")+
  stat_cor(method = 'pearson', aes(x = AS_del, y = QuiescenceScore))

ggplot(df.merged_gnaq_3, aes(x=AS_amp, y=QuiescenceScore)) +
  geom_point(aes(col=Type)) +
  geom_smooth(method="lm", se=F) + 
  xlim(c(0, 40)) + 
  ylim(c(-30, 30)) + 
  labs(subtitle="QS Vs AS_amp", 
       title="Scatterplot_GNAQ")+
  stat_cor(method = 'pearson', aes(x = AS_amp, y = QuiescenceScore))

ggplot(data = df.merged_gnaq_3, aes(x = AneuploidyScore, y = QuiescenceScore,  color = GNAQ_aamutation_cna)) + 
  geom_point() + geom_smooth(method = lm) +
  scale_color_manual(values = c('#FF7400', '#009999', '#3914AF',"#000000"))  +
  labs(title = 'QS Vs AS - GNAQ_aamutation_cna') + 
  stat_cor(method = 'pearson', aes(x = AneuploidyScore, y = QuiescenceScore,  color = GNAQ_aamutation_cna))
#####AS- GNAQ_aamutation_cna#####
#mean of groups
aggregate(df.merged_gnaq_3$AneuploidyScore, by=list(type=df.merged_gnaq_3$GNAQ_aamutation_cna),mean)

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_gnaq_3, x = "GNAQ_aamutation_cna", y = "AneuploidyScore",
          color = "GNAQ_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####AS_amp- GNAQ_aamutation_cna#####
#mean of groups
aggregate(df.merged_gnaq_3$AS_amp, by=list(type=df.merged_gnaq_3$GNAQ_aamutation_cna),mean)

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_gnaq_3, x = "GNAQ_aamutation_cna", y = "AS_amp",
          color = "GNAQ_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####AS_del- GNAQ_aamutation_cna#####
#mean of groups
aggregate(df.merged_gnaq_3$AS_del, by=list(type=df.merged_gnaq_3$GNAQ_aamutation_cna),mean)

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_gnaq_3, x = "GNAQ_aamutation_cna", y = "AS_del",
          color = "GNAQ_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####Gb- GNAQ_aamutation_cna#####
#mean of groups
aggregate(df.merged_gnaq_3$Genome_doublings, by=list(type=df.merged_gnaq_3$GNAQ_aamutation_cna),mean)

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_gnaq_3, x = "GNAQ_aamutation_cna", y = "Genome_doublings",
          color = "GNAQ_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####Leuk- GNAQ_aamutation_cna#####
#mean of groups
aggregate(df.merged_gnaq_3$Leuk, by=list(type=df.merged_gnaq_3$GNAQ_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_gnaq_3, x = "GNAQ_aamutation_cna", y = "Leuk",
          color = "GNAQ_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####Purity- GNAQ_aamutation_cna#####
#mean of groups
aggregate(df.merged_gnaq_3$Purity, by=list(type=df.merged_gnaq_3$GNAQ_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_gnaq_3, x = "GNAQ_aamutation_cna", y = "Purity",
          color = "GNAQ_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3)
#####Stroma - GNAQ_aamutation_cna#####
#mean of groups
aggregate(df.merged_gnaq_3$Stroma, by=list(type=df.merged_gnaq_3$GNAQ_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_gnaq_3, x = "GNAQ_aamutation_cna", y = "Stroma",
          color = "GNAQ_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3)
#####Stroma_notLeukocyte - GNAQ_aamutation_cna#####
#mean of groups
aggregate(df.merged_gnaq_3$Stroma_notLeukocyte, by=list(type=df.merged_gnaq_3$GNAQ_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_gnaq_3, x = "GNAQ_aamutation_cna", y = "Stroma_notLeukocyte",
          color = "GNAQ_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####SilentMutationspeMb - GNAQ_aamutation_cna#####
my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_gnaq_3, x = "GNAQ_aamutation_cna", y = "SilentMutationspeMb",
          color = "GNAQ_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####Non.silentMutationsperMb - GNAQ_aamutation_cna#####
my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_gnaq_3, x = "GNAQ_aamutation_cna", y = "Non.silentMutationsperMb",
          color = "GNAQ_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")
#####导入样本总突变数据NM---已下为突变数体所有代码 可直接复制后分别替换基因大小写#####
pgen <- read.csv("pgen.csv")
head(pgen)
hist(pgen$number_of_mutations)

df.merged_gnaq_4 <- merge(df.merged_gnaq_3[,c("SAMPLE_ID",
                                                "GNAQ_TP53_status",
                                                "TP53_status","GNAQ_status","GNAQ_log2expression",
                                                "QuiescenceScore","Group",
                                                "GNAQ_GNAQCNA",
                                                "STUDY_ID",
                                                "GNAQ_CNA",
                                                "GNAQ_aamutation_cna",
                                                "GNAQ_aamutation_cna_TP53_aamutation",
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
head(df.merged_gnaq_4)
#####Number of Mutations- GNAQ_aamutation_cna#####
#mean of groups
aggregate(df.merged_gnaq_4$number_of_mutations, by=list(type=df.merged_gnaq_4$GNAQ_aamutation_cna),mean,na.rm = TRUE)#exclude NA

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged_gnaq_4, x = "GNAQ_aamutation_cna", y = "number_of_mutations",
          color = "GNAQ_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          add = "jitter",order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3)

