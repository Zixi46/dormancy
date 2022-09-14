#linear model_GNAQ
#####QS ——— AA mutaiton and CNV#####
######QS --- GNAQ_AA mutation######
library(tidyverse)
sample_n(df.merged_gnaq2, 3)
model_qs_gnaq_aa<- lm(QuiescenceScore ~ GNAQ_status,data=df.merged_gnaq2)#建立一个简单的线性回归模型来解释
summary(model_qs_gnaq_aa)$coef#无关联，不继续构建
df.merged_gnaq2_model <- df.merged_gnaq2 %>% transmute_all(as.factor)#所有变量转化为因素
contrasts(df.merged_gnaq2_model$GNAQ_status)#使用contrasts()函数可以返回R创建的哑变量编码：
df.merged_gnaq2_model <- df.merged_gnaq2_model %>%
  mutate(GNAQ_status = relevel(GNAQ_status, ref = "WT"))#基线类别设置为WT, WT为0
df.merged_gnaq2_model$QuiescenceScore <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_gnaq2_model$QuiescenceScore)) #QuiescenceScore向量有因子类。这种数据类型不适合进行数值计算(^平方根)。
model_qs_gnaq_aa<- lm(QuiescenceScore ~ GNAQ_status,data=df.merged_gnaq2_model)
summary(model_qs_gnaq_aa)$coef

######QS --- GNAQ_TP53_AA mutation######
sample_n(df.merged_gnaq2, 3)
model_qs_gnaq_tp53_aa<- lm(QuiescenceScore ~ GNAQ_TP53_status,data=df.merged_gnaq2)#建立一个简单的线性回归模型来解释
summary(model_qs_gnaq_tp53_aa)$coef#存在强关联，可以构建线性模型
df.merged_gnaq2_model <- df.merged_gnaq2 %>% transmute_all(as.factor)#所有变量转化为因素
contrasts(df.merged_gnaq2_model$GNAQ_TP53_status)#使用contrasts()函数可以返回R创建的哑变量编码：
df.merged_gnaq2_model <- df.merged_gnaq2_model %>%
  mutate(GNAQ_TP53_status = relevel(GNAQ_TP53_status, ref = "WT_WT"))#基线类别设置为WT, WT为0
df.merged_gnaq2_model$QuiescenceScore <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_gnaq2_model$QuiescenceScore)) #QuiescenceScore向量有因子类。这种数据类型不适合进行数值计算(^平方根)。
model_qs_gnaq_tp53_aa<- lm(QuiescenceScore ~ GNAQ_TP53_status,data=df.merged_gnaq2_model)
summary(model_qs_gnaq_tp53_aa)$coef
######EXP --- GNAQ_TP53_AA mutation######
sample_n(df.merged_gnaq2, 3)
model_exp_gnaq_tp53_aa<- lm(GNAQ_log2expression ~ GNAQ_TP53_status,data=df.merged_gnaq2)#建立一个简单的线性回归模型来解释
summary(model_exp_gnaq_tp53_aa)$coef#存在强关联，可以构建线性模型
df.merged_gnaq2_model <- df.merged_gnaq2 %>% transmute_all(as.factor)#所有变量转化为因素
contrasts(df.merged_gnaq2_model$GNAQ_TP53_status)#使用contrasts()函数可以返回R创建的哑变量编码：
df.merged_gnaq2_model <- df.merged_gnaq2_model %>%
  mutate(GNAQ_TP53_status = relevel(GNAQ_TP53_status, ref = "WT_WT"))#基线类别设置为WT, WT为0
df.merged_gnaq2_model$GNAQ_log2expression <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_gnaq2_model$GNAQ_log2expression)) #GNAQ_log2expression向量有因子类。这种数据类型不适合进行数值计算(^平方根)。
model_exp_gnaq_tp53_aa<- lm(GNAQ_log2expression ~ GNAQ_TP53_status,data=df.merged_gnaq2_model)
summary(model_exp_gnaq_tp53_aa)$coef

######QS --- GNAQ_AA mutation + CNA######
head(df.merged_gnaq2)
model_qs_gnaq_aa_cna<- lm(QuiescenceScore ~ GNAQ_aamutation_cna,data=df.merged_gnaq2)
summary(model_qs_gnaq_aa_cna)$coef#有关联，可以继续
df.merged_gnaq2_model <- df.merged_gnaq2 %>% transmute_all(as.factor)
contrasts(df.merged_gnaq2_model$GNAQ_aamutation_cna)
df.merged_gnaq2_model <- df.merged_gnaq2_model %>%
  mutate(GNAQ_aamutation_cna = relevel(GNAQ_aamutation_cna, ref = "WT"))
df.merged_gnaq2_model$QuiescenceScore <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_gnaq2_model$QuiescenceScore))
model_qs_gnaq_aa_cna<- lm(QuiescenceScore ~ GNAQ_aamutation_cna,data=df.merged_gnaq2_model)
summary(model_qs_gnaq_aa_cna)$coef
######EXP --- GNAQ_AA mutation######
hist(df.merged_gnaq2$GNAQ_log2expression)#正态分布
model_exp_gnaq_aa<- lm(GNAQ_log2expression ~ GNAQ_status,data=df.merged_gnaq2)
summary(model_exp_gnaq_aa)$coef
contrasts(df.merged_gnaq2_model$GNAQ_status)
df.merged_gnaq2_model <- df.merged_gnaq2_model %>%
  mutate(GNAQ_status = relevel(GNAQ_status, ref = "WT"))
df.merged_gnaq2_model$GNAQ_log2expression <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_gnaq2_model$GNAQ_log2expression))
model_exp_gnaq_aa<- lm(GNAQ_log2expression ~ GNAQ_status,data=df.merged_gnaq2_model)
summary(model_exp_gnaq_aa)$coef
######EXP --- GNAQ_AA mutation + CNA######
hist(df.merged_gnaq2$GNAQ_log2expression)#正态分布
model_exp_gnaq_aa_cna<- lm(GNAQ_log2expression ~ GNAQ_aamutation_cna,data=df.merged_gnaq2)
summary(model_exp_gnaq_aa_cna)$coef
contrasts(df.merged_gnaq2_model$GNAQ_aamutation_cna)
df.merged_gnaq2_model <- df.merged_gnaq2_model %>%
  mutate(GNAQ_aamutation_cna = relevel(GNAQ_aamutation_cna, ref = "WT"))
df.merged_gnaq2_model$GNAQ_log2expression <- as.numeric(#平方转化为可计算数值
  as.character(df.merged_gnaq2_model$GNAQ_log2expression))
model_exp_gnaq_aa_cna<- lm(GNAQ_log2expression ~ GNAQ_aamutation_cna,data=df.merged_gnaq2_model)
summary(model_exp_gnaq_aa_cna)$coef
