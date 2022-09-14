#linear model_ERBB3
#####QS ——— AA mutaiton and CNV#####
######QS --- ERBB3_AA mutation######
library(tidyverse)
sample_n(df.merged_erbb3_2, 3)
model_qs_erbb3_aa<- lm(QuiescenceScore ~ ERBB3_status,data=df.merged_erbb3_2)#建立一个简单的线性回归模型来解释
summary(model_qs_erbb3_aa)$coef#关联，继续构建
df.merged_erbb3_2_model <- df.merged_erbb3_2 %>% transmute_all(as.factor)#所有变量转化为因素
contrasts(df.merged_erbb3_2_model$ERBB3_status)#使用contrasts()函数可以返回R创建的哑变量编码：
df.merged_erbb3_2_model <- df.merged_erbb3_2_model %>%
  mutate(ERBB3_status = relevel(ERBB3_status, ref = "WT"))#基线类别设置为WT, WT为0
df.merged_erbb3_2_model$QuiescenceScore <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_erbb3_2_model$QuiescenceScore)) #QuiescenceScore向量有因子类。这种数据类型不适合进行数值计算(^平方根)。
model_qs_erbb3_aa<- lm(QuiescenceScore ~ ERBB3_status,data=df.merged_erbb3_2_model)
summary(model_qs_erbb3_aa)$coef

######QS --- ERBB3_TP53_AA mutation######
sample_n(df.merged_erbb3_2, 3)
model_qs_erbb3_tp53_aa<- lm(QuiescenceScore ~ ERBB3_TP53_status,data=df.merged_erbb3_2)#建立一个简单的线性回归模型来解释
summary(model_qs_erbb3_tp53_aa)$coef#存在强关联，可以构建线性模型
df.merged_erbb3_2_model <- df.merged_erbb3_2 %>% transmute_all(as.factor)#所有变量转化为因素
contrasts(df.merged_erbb3_2_model$ERBB3_TP53_status)#使用contrasts()函数可以返回R创建的哑变量编码：
df.merged_erbb3_2_model <- df.merged_erbb3_2_model %>%
  mutate(ERBB3_TP53_status = relevel(ERBB3_TP53_status, ref = "WT_WT"))#基线类别设置为WT, WT为0
df.merged_erbb3_2_model$QuiescenceScore <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_erbb3_2_model$QuiescenceScore)) #QuiescenceScore向量有因子类。这种数据类型不适合进行数值计算(^平方根)。
model_qs_erbb3_tp53_aa<- lm(QuiescenceScore ~ ERBB3_TP53_status,data=df.merged_erbb3_2_model)
summary(model_qs_erbb3_tp53_aa)$coef
######EXP --- ERBB3_TP53_AA mutation######
sample_n(df.merged_erbb3_2, 3)
model_exp_erbb3_tp53_aa<- lm(ERBB3_log2expression ~ ERBB3_TP53_status,data=df.merged_erbb3_2)#建立一个简单的线性回归模型来解释
summary(model_exp_erbb3_tp53_aa)$coef#存在强关联，可以构建线性模型
df.merged_erbb3_2_model <- df.merged_erbb3_2 %>% transmute_all(as.factor)#所有变量转化为因素
contrasts(df.merged_erbb3_2_model$ERBB3_TP53_status)#使用contrasts()函数可以返回R创建的哑变量编码：
df.merged_erbb3_2_model <- df.merged_erbb3_2_model %>%
  mutate(ERBB3_TP53_status = relevel(ERBB3_TP53_status, ref = "WT_WT"))#基线类别设置为WT, WT为0
df.merged_erbb3_2_model$ERBB3_log2expression <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_erbb3_2_model$ERBB3_log2expression)) #ERBB3_log2expression向量有因子类。这种数据类型不适合进行数值计算(^平方根)。
model_exp_erbb3_tp53_aa<- lm(ERBB3_log2expression ~ ERBB3_TP53_status,data=df.merged_erbb3_2_model)
summary(model_exp_erbb3_tp53_aa)$coef

######QS --- ERBB3_AA mutation + CNA######
head(df.merged_erbb3_2)
model_qs_erbb3_aa_cna<- lm(QuiescenceScore ~ ERBB3_aamutation_cna,data=df.merged_erbb3_2)
summary(model_qs_erbb3_aa_cna)$coef#有关联，可以继续
df.merged_erbb3_2_model <- df.merged_erbb3_2 %>% transmute_all(as.factor)
contrasts(df.merged_erbb3_2_model$ERBB3_aamutation_cna)
df.merged_erbb3_2_model <- df.merged_erbb3_2_model %>%
  mutate(ERBB3_aamutation_cna = relevel(ERBB3_aamutation_cna, ref = "WT"))
df.merged_erbb3_2_model$QuiescenceScore <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_erbb3_2_model$QuiescenceScore))
model_qs_erbb3_aa_cna<- lm(QuiescenceScore ~ ERBB3_aamutation_cna,data=df.merged_erbb3_2_model)
summary(model_qs_erbb3_aa_cna)$coef
######EXP --- ERBB3_AA mutation######
hist(df.merged_erbb3_2$ERBB3_log2expression)#正态分布
model_exp_erbb3_aa<- lm(ERBB3_log2expression ~ ERBB3_status,data=df.merged_erbb3_2)
summary(model_exp_erbb3_aa)$coef
contrasts(df.merged_erbb3_2_model$ERBB3_status)
df.merged_erbb3_2_model <- df.merged_erbb3_2_model %>%
  mutate(ERBB3_status = relevel(ERBB3_status, ref = "WT"))
df.merged_erbb3_2_model$ERBB3_log2expression <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_erbb3_2_model$ERBB3_log2expression))
model_exp_erbb3_aa<- lm(ERBB3_log2expression ~ ERBB3_status,data=df.merged_erbb3_2_model)
summary(model_exp_erbb3_aa)$coef
######EXP --- ERBB3_AA mutation + CNA######
hist(df.merged_erbb3_2$ERBB3_log2expression)#正态分布
model_exp_erbb3_aa_cna<- lm(ERBB3_log2expression ~ ERBB3_aamutation_cna,data=df.merged_erbb3_2)
summary(model_exp_erbb3_aa_cna)$coef
contrasts(df.merged_erbb3_2_model$ERBB3_aamutation_cna)
df.merged_erbb3_2_model <- df.merged_erbb3_2_model %>%
  mutate(ERBB3_aamutation_cna = relevel(ERBB3_aamutation_cna, ref = "WT"))
df.merged_erbb3_2_model$ERBB3_log2expression <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_erbb3_2_model$ERBB3_log2expression))
model_exp_erbb3_aa_cna<- lm(ERBB3_log2expression ~ ERBB3_aamutation_cna,data=df.merged_erbb3_2_model)
summary(model_exp_erbb3_aa_cna)$coef