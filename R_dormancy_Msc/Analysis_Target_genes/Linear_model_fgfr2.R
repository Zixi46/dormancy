#linear model_FGFR2
#####QS ——— AA mutaiton and CNV#####
######QS --- FGFR2_AA mutation######
library(tidyverse)
sample_n(df.merged_fgfr2_2, 3)
model_qs_fgfr2_aa<- lm(QuiescenceScore ~ FGFR2_status,data=df.merged_fgfr2_2)#建立一个简单的线性回归模型来解释
summary(model_qs_fgfr2_aa)$coef#关联，继续构建
df.merged_fgfr2_2_model <- df.merged_fgfr2_2 %>% transmute_all(as.factor)#所有变量转化为因素
contrasts(df.merged_fgfr2_2_model$FGFR2_status)#使用contrasts()函数可以返回R创建的哑变量编码：
df.merged_fgfr2_2_model <- df.merged_fgfr2_2_model %>%
  mutate(FGFR2_status = relevel(FGFR2_status, ref = "WT"))#基线类别设置为WT, WT为0
df.merged_fgfr2_2_model$QuiescenceScore <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_fgfr2_2_model$QuiescenceScore)) #QuiescenceScore向量有因子类。这种数据类型不适合进行数值计算(^平方根)。
model_qs_fgfr2_aa<- lm(QuiescenceScore ~ FGFR2_status,data=df.merged_fgfr2_2_model)
summary(model_qs_fgfr2_aa)$coef

######QS --- FGFR2_TP53_AA mutation######
sample_n(df.merged_fgfr2_2, 3)
model_qs_fgfr2_tp53_aa<- lm(QuiescenceScore ~ FGFR2_TP53_status,data=df.merged_fgfr2_2)#建立一个简单的线性回归模型来解释
summary(model_qs_fgfr2_tp53_aa)$coef#存在强关联，可以构建线性模型
df.merged_fgfr2_2_model <- df.merged_fgfr2_2 %>% transmute_all(as.factor)#所有变量转化为因素
contrasts(df.merged_fgfr2_2_model$FGFR2_TP53_status)#使用contrasts()函数可以返回R创建的哑变量编码：
df.merged_fgfr2_2_model <- df.merged_fgfr2_2_model %>%
  mutate(FGFR2_TP53_status = relevel(FGFR2_TP53_status, ref = "WT_WT"))#基线类别设置为WT, WT为0
df.merged_fgfr2_2_model$QuiescenceScore <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_fgfr2_2_model$QuiescenceScore)) #QuiescenceScore向量有因子类。这种数据类型不适合进行数值计算(^平方根)。
model_qs_fgfr2_tp53_aa<- lm(QuiescenceScore ~ FGFR2_TP53_status,data=df.merged_fgfr2_2_model)
summary(model_qs_fgfr2_tp53_aa)$coef
######EXP --- FGFR2_TP53_AA mutation######
sample_n(df.merged_fgfr2_2, 3)
model_exp_fgfr2_tp53_aa<- lm(FGFR2_log2expression ~ FGFR2_TP53_status,data=df.merged_fgfr2_2)#建立一个简单的线性回归模型来解释
summary(model_exp_fgfr2_tp53_aa)$coef#存在强关联，可以构建线性模型
df.merged_fgfr2_2_model <- df.merged_fgfr2_2 %>% transmute_all(as.factor)#所有变量转化为因素
contrasts(df.merged_fgfr2_2_model$FGFR2_TP53_status)#使用contrasts()函数可以返回R创建的哑变量编码：
df.merged_fgfr2_2_model <- df.merged_fgfr2_2_model %>%
  mutate(FGFR2_TP53_status = relevel(FGFR2_TP53_status, ref = "WT_WT"))#基线类别设置为WT, WT为0
df.merged_fgfr2_2_model$FGFR2_log2expression <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_fgfr2_2_model$FGFR2_log2expression)) #FGFR2_log2expression向量有因子类。这种数据类型不适合进行数值计算(^平方根)。
model_exp_fgfr2_tp53_aa<- lm(FGFR2_log2expression ~ FGFR2_TP53_status,data=df.merged_fgfr2_2_model)
summary(model_exp_fgfr2_tp53_aa)$coef

######QS --- FGFR2_AA mutation + CNA######
head(df.merged_fgfr2_2)
model_qs_fgfr2_aa_cna<- lm(QuiescenceScore ~ FGFR2_aamutation_cna,data=df.merged_fgfr2_2)
summary(model_qs_fgfr2_aa_cna)$coef#有关联，可以继续
df.merged_fgfr2_2_model <- df.merged_fgfr2_2 %>% transmute_all(as.factor)
contrasts(df.merged_fgfr2_2_model$FGFR2_aamutation_cna)
df.merged_fgfr2_2_model <- df.merged_fgfr2_2_model %>%
  mutate(FGFR2_aamutation_cna = relevel(FGFR2_aamutation_cna, ref = "WT"))
df.merged_fgfr2_2_model$QuiescenceScore <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_fgfr2_2_model$QuiescenceScore))
model_qs_fgfr2_aa_cna<- lm(QuiescenceScore ~ FGFR2_aamutation_cna,data=df.merged_fgfr2_2_model)
summary(model_qs_fgfr2_aa_cna)$coef
######EXP --- FGFR2_AA mutation######
hist(df.merged_fgfr2_2$FGFR2_log2expression)#正态分布
model_exp_fgfr2_aa<- lm(FGFR2_log2expression ~ FGFR2_status,data=df.merged_fgfr2_2)
summary(model_exp_fgfr2_aa)$coef
contrasts(df.merged_fgfr2_2_model$FGFR2_status)
df.merged_fgfr2_2_model <- df.merged_fgfr2_2_model %>%
  mutate(FGFR2_status = relevel(FGFR2_status, ref = "WT"))
df.merged_fgfr2_2_model$FGFR2_log2expression <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_fgfr2_2_model$FGFR2_log2expression))
model_exp_fgfr2_aa<- lm(FGFR2_log2expression ~ FGFR2_status,data=df.merged_fgfr2_2_model)
summary(model_exp_fgfr2_aa)$coef
######EXP --- FGFR2_AA mutation + CNA######
hist(df.merged_fgfr2_2$FGFR2_log2expression)#正态分布
model_exp_fgfr2_aa_cna<- lm(FGFR2_log2expression ~ FGFR2_aamutation_cna,data=df.merged_fgfr2_2)
summary(model_exp_fgfr2_aa_cna)$coef
contrasts(df.merged_fgfr2_2_model$FGFR2_aamutation_cna)
df.merged_fgfr2_2_model <- df.merged_fgfr2_2_model %>%
  mutate(FGFR2_aamutation_cna = relevel(FGFR2_aamutation_cna, ref = "WT"))
df.merged_fgfr2_2_model$FGFR2_log2expression <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_fgfr2_2_model$FGFR2_log2expression))
model_exp_fgfr2_aa_cna<- lm(FGFR2_log2expression ~ FGFR2_aamutation_cna,data=df.merged_fgfr2_2_model)
summary(model_exp_fgfr2_aa_cna)$coef
