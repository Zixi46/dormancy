#linear model_PYHIN1
#####QS ——— AA mutaiton and CNV#####
######QS --- PYHIN1_AA mutation######
library(tidyverse)
sample_n(df.merged2, 3)
model<- lm(QuiescenceScore ~ PYHIN1_status,data=df.merged2)#建立一个简单的线性回归模型来解释
summary(model)$coef#存在强关联，可以构建线性模型
df.merged2_model <- df.merged2 %>% transmute_all(as.factor)#所有变量转化为因素
contrasts(df.merged2_model$PYHIN1_status)#使用contrasts()函数可以返回R创建的哑变量编码：
df.merged2_model <- df.merged2_model %>%
  mutate(PYHIN1_status = relevel(PYHIN1_status, ref = "WT"))#基线类别设置为WT, WT为0
df.merged2_model$QuiescenceScore <- as.numeric(#平方转化为可计算数值
  as.character(df.merged2_model$QuiescenceScore)) #QuiescenceScore向量有因子类。这种数据类型不适合进行数值计算(^平方根)。
model<- lm(QuiescenceScore ~ PYHIN1_status,data=df.merged2_model)
summary(model)$coef

######QS --- PYHIN1_TP53_AA mutation######
library(tidyverse)
sample_n(df.merged9, 3)
model_qs_pyhin1_tp53_aa<- lm(QuiescenceScore ~ PYHIN1_TP53_status,data=df.merged9)#建立一个简单的线性回归模型来解释
summary(model_qs_pyhin1_tp53_aa)$coef#存在强关联，可以构建线性模型
df.merged9_model <- df.merged9 %>% transmute_all(as.factor)#所有变量转化为因素
contrasts(df.merged9_model$PYHIN1_TP53_status)#使用contrasts()函数可以返回R创建的哑变量编码：
df.merged9_model <- df.merged9_model %>%
  mutate(PYHIN1_TP53_status = relevel(PYHIN1_TP53_status, ref = "WT_WT"))#基线类别设置为WT, WT为0
df.merged9_model$QuiescenceScore <- as.numeric(#平方转化为可计算数值
  as.character(df.merged9_model$QuiescenceScore)) #QuiescenceScore向量有因子类。这种数据类型不适合进行数值计算(^平方根)。
model_qs_pyhin1_tp53_aa<- lm(QuiescenceScore ~ PYHIN1_TP53_status,data=df.merged9_model)
summary(model_qs_pyhin1_tp53_aa)$coef
######EXP --- PYHIN1_TP53_AA mutation######
sample_n(df.merged9, 3)
model_exp_pyhin1_tp53_aa<- lm(PYHIN1_log2expression ~ PYHIN1_TP53_status,data=df.merged9)#建立一个简单的线性回归模型来解释
summary(model_exp_pyhin1_tp53_aa)$coef#存在强关联，可以构建线性模型
df.merged9_model <- df.merged9 %>% transmute_all(as.factor)#所有变量转化为因素
contrasts(df.merged9_model$PYHIN1_TP53_status)#使用contrasts()函数可以返回R创建的哑变量编码：
df.merged9_model <- df.merged9_model %>%
  mutate(PYHIN1_TP53_status = relevel(PYHIN1_TP53_status, ref = "WT_WT"))#基线类别设置为WT, WT为0
df.merged9_model$PYHIN1_log2expression <- as.numeric(#平方转化为可计算数值
  as.character(df.merged9_model$PYHIN1_log2expression)) #PYHIN1_log2expression向量有因子类。这种数据类型不适合进行数值计算(^平方根)。
model_exp_pyhin1_tp53_aa<- lm(PYHIN1_log2expression ~ PYHIN1_TP53_status,data=df.merged9_model)
summary(model_exp_pyhin1_tp53_aa)$coef

######QS --- PYHIN1_AA mutation + CNA######
head(df.merged9)
model_qs_pyhin1_aa_cna<- lm(QuiescenceScore ~ PYHIN1_aamutation_cna,data=df.merged9)
summary(model_qs_pyhin1_aa_cna)$coef
df.merged9_model <- df.merged9 %>% transmute_all(as.factor)
contrasts(df.merged9_model$PYHIN1_aamutation_cna)
df.merged9_model <- df.merged9_model %>%
  mutate(PYHIN1_aamutation_cna = relevel(PYHIN1_aamutation_cna, ref = "WT"))
df.merged9_model$QuiescenceScore <- as.numeric(#平方转化为可计算数值
  as.character(df.merged9_model$QuiescenceScore))
model_qs_pyhin1_aa_cna<- lm(QuiescenceScore ~ PYHIN1_aamutation_cna,data=df.merged9_model)
summary(model_qs_pyhin1_aa_cna)$coef



######EXP --- PYHIN1_AA mutation######
hist(df.merged9$PYHIN1_log2expression)#正态分布
model_exp_pyhin1_aa<- lm(PYHIN1_log2expression ~ PYHIN1_status,data=df.merged9)
summary(model_exp_pyhin1_aa)$coef
contrasts(df.merged9_model$PYHIN1_status)
df.merged9_model <- df.merged9_model %>%
  mutate(PYHIN1_status = relevel(PYHIN1_status, ref = "WT"))
df.merged9_model$PYHIN1_log2expression <- as.numeric(#平方转化为可计算数值
  as.character(df.merged9_model$PYHIN1_log2expression))
model_exp_pyhin1_aa<- lm(PYHIN1_log2expression ~ PYHIN1_status,data=df.merged9_model)
summary(model_exp_pyhin1_aa)$coef
######EXP --- PYHIN1_AA mutation + CNA######
hist(df.merged9$PYHIN1_log2expression)#正态分布
model_exp_pyhin1_aa_cna<- lm(PYHIN1_log2expression ~ PYHIN1_aamutation_cna,data=df.merged9)
summary(model_exp_pyhin1_aa_cna)$coef
contrasts(df.merged9_model$PYHIN1_aamutation_cna)
df.merged9_model <- df.merged9_model %>%
  mutate(PYHIN1_aamutation_cna = relevel(PYHIN1_aamutation_cna, ref = "WT"))
df.merged9_model$PYHIN1_log2expression <- as.numeric(#平方转化为可计算数值
  as.character(df.merged9_model$PYHIN1_log2expression))
model_exp_pyhin1_aa_cna<- lm(PYHIN1_log2expression ~ PYHIN1_aamutation_cna,data=df.merged9_model)
summary(model_exp_pyhin1_aa_cna)$coef

