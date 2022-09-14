#linear model_PLPPR4
#####QS ——— AA mutaiton and CNV#####
######QS --- PLPPR4_AA mutation######
library(tidyverse)
sample_n(df.merged_plppr4_2, 3)
model_qs_plppr4_aa<- lm(QuiescenceScore ~ PLPPR4_status,data=df.merged_plppr4_2)#建立一个简单的线性回归模型来解释
summary(model_qs_plppr4_aa)$coef#关联，继续构建
df.merged_plppr4_2_model <- df.merged_plppr4_2 %>% transmute_all(as.factor)#所有变量转化为因素
contrasts(df.merged_plppr4_2_model$PLPPR4_status)#使用contrasts()函数可以返回R创建的哑变量编码：
df.merged_plppr4_2_model <- df.merged_plppr4_2_model %>%
  mutate(PLPPR4_status = relevel(PLPPR4_status, ref = "WT"))#基线类别设置为WT, WT为0
df.merged_plppr4_2_model$QuiescenceScore <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_plppr4_2_model$QuiescenceScore)) #QuiescenceScore向量有因子类。这种数据类型不适合进行数值计算(^平方根)。
model_qs_plppr4_aa<- lm(QuiescenceScore ~ PLPPR4_status,data=df.merged_plppr4_2_model)
summary(model_qs_plppr4_aa)$coef
######QS --- PLPPR4_TP53_AA mutation######
sample_n(df.merged_plppr4_2, 3)
model_qs_plppr4_tp53_aa<- lm(QuiescenceScore ~ PLPPR4_TP53_status,data=df.merged_plppr4_2)#建立一个简单的线性回归模型来解释
summary(model_qs_plppr4_tp53_aa)$coef#存在强关联，可以构建线性模型
df.merged_plppr4_2_model <- df.merged_plppr4_2 %>% transmute_all(as.factor)#所有变量转化为因素
contrasts(df.merged_plppr4_2_model$PLPPR4_TP53_status)#使用contrasts()函数可以返回R创建的哑变量编码：
df.merged_plppr4_2_model <- df.merged_plppr4_2_model %>%
  mutate(PLPPR4_TP53_status = relevel(PLPPR4_TP53_status, ref = "WT_WT"))#基线类别设置为WT, WT为0
df.merged_plppr4_2_model$QuiescenceScore <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_plppr4_2_model$QuiescenceScore)) #QuiescenceScore向量有因子类。这种数据类型不适合进行数值计算(^平方根)。
model_qs_plppr4_tp53_aa<- lm(QuiescenceScore ~ PLPPR4_TP53_status,data=df.merged_plppr4_2_model)
summary(model_qs_plppr4_tp53_aa)$coef
######EXP --- PLPPR4_TP53_AA mutation######
sample_n(df.merged_plppr4_2, 3)
model_exp_plppr4_tp53_aa<- lm(PLPPR4_log2expression ~ PLPPR4_TP53_status,data=df.merged_plppr4_2)#建立一个简单的线性回归模型来解释
summary(model_exp_plppr4_tp53_aa)$coef#存在强关联，可以构建线性模型
df.merged_plppr4_2_model <- df.merged_plppr4_2 %>% transmute_all(as.factor)#所有变量转化为因素
contrasts(df.merged_plppr4_2_model$PLPPR4_TP53_status)#使用contrasts()函数可以返回R创建的哑变量编码：
df.merged_plppr4_2_model <- df.merged_plppr4_2_model %>%
  mutate(PLPPR4_TP53_status = relevel(PLPPR4_TP53_status, ref = "WT_WT"))#基线类别设置为WT, WT为0
df.merged_plppr4_2_model$PLPPR4_log2expression <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_plppr4_2_model$PLPPR4_log2expression)) #PLPPR4_log2expression向量有因子类。这种数据类型不适合进行数值计算(^平方根)。
model_exp_plppr4_tp53_aa<- lm(PLPPR4_log2expression ~ PLPPR4_TP53_status,data=df.merged_plppr4_2_model)
summary(model_exp_plppr4_tp53_aa)$coef

######QS --- PLPPR4_AA mutation + CNA######
head(df.merged_plppr4_2)
model_qs_plppr4_aa_cna<- lm(QuiescenceScore ~ PLPPR4_aamutation_cna,data=df.merged_plppr4_2)
summary(model_qs_plppr4_aa_cna)$coef#有关联，可以继续
df.merged_plppr4_2_model <- df.merged_plppr4_2 %>% transmute_all(as.factor)
contrasts(df.merged_plppr4_2_model$PLPPR4_aamutation_cna)
df.merged_plppr4_2_model <- df.merged_plppr4_2_model %>%
  mutate(PLPPR4_aamutation_cna = relevel(PLPPR4_aamutation_cna, ref = "WT"))
df.merged_plppr4_2_model$QuiescenceScore <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_plppr4_2_model$QuiescenceScore))
model_qs_plppr4_aa_cna<- lm(QuiescenceScore ~ PLPPR4_aamutation_cna,data=df.merged_plppr4_2_model)
summary(model_qs_plppr4_aa_cna)$coef
######EXP --- PLPPR4_AA mutation######
hist(df.merged_plppr4_2$PLPPR4_log2expression)#正态分布
model_exp_plppr4_aa<- lm(PLPPR4_log2expression ~ PLPPR4_status,data=df.merged_plppr4_2)
summary(model_exp_plppr4_aa)$coef
contrasts(df.merged_plppr4_2_model$PLPPR4_status)
df.merged_plppr4_2_model <- df.merged_plppr4_2_model %>%
  mutate(PLPPR4_status = relevel(PLPPR4_status, ref = "WT"))
df.merged_plppr4_2_model$PLPPR4_log2expression <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_plppr4_2_model$PLPPR4_log2expression))
model_exp_plppr4_aa<- lm(PLPPR4_log2expression ~ PLPPR4_status,data=df.merged_plppr4_2_model)
summary(model_exp_plppr4_aa)$coef
######EXP --- PLPPR4_AA mutation + CNA######
hist(df.merged_plppr4_2$PLPPR4_log2expression)#正态分布
model_exp_plppr4_aa_cna<- lm(PLPPR4_log2expression ~ PLPPR4_aamutation_cna,data=df.merged_plppr4_2)
summary(model_exp_plppr4_aa_cna)$coef
contrasts(df.merged_plppr4_2_model$PLPPR4_aamutation_cna)
df.merged_plppr4_2_model <- df.merged_plppr4_2_model %>%
  mutate(PLPPR4_aamutation_cna = relevel(PLPPR4_aamutation_cna, ref = "WT"))
df.merged_plppr4_2_model$PLPPR4_log2expression <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_plppr4_2_model$PLPPR4_log2expression))
model_exp_plppr4_aa_cna<- lm(PLPPR4_log2expression ~ PLPPR4_aamutation_cna,data=df.merged_plppr4_2_model)
summary(model_exp_plppr4_aa_cna)$coef