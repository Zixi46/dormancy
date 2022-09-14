#4a
my_cols <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
             "#E16462FF", "#FCA636FF", "#F0F921FF")
AS_ballon_aa_cna <- read.csv("AS_ballon_data.csv",row.names = 1)
ggballoonplot(AS_ballon_aa_cna, fill = "value")+
  scale_fill_gradientn(colors = my_cols)

ggballoonplot(AS_ballon_aa_cna, 
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
    colors=c("white","red"),
    values=rescale(c(0,5)),
    limits=c(0,5))+
  labs(x="Gene status (AA mutation + CNV)",y="Cancer type ( p < 0.05)",fill="Difference of 
Aneuploidy Score
(ref: WT)")


AS_amp_ballon_aa_cna <- read.csv("AS_amp_ballon_data.csv",row.names = 1)
ggballoonplot(AS_amp_ballon_aa_cna, 
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
    values=rescale(c(-2,0,10)),
    limits=c(-2,11))+
  labs(x="Gene status (AA mutation + CNV)",y="Cancer type ( p < 0.05)",fill="Difference of 
number of AS_amp
(ref: WT)")


AS_del_ballon_aa_cna <- read.csv("AS_del_ballon_data.csv",row.names = 1)
ggballoonplot(AS_del_ballon_aa_cna, 
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
    values=rescale(c(-2,0,7)),
    limits=c(-2,7))+
  labs(x="Gene status (AA mutation + CNV)",y="Cancer type ( p < 0.05)",fill="Difference of 
number of AS_del
(ref: WT)")


GB_ballon_aa_cna <- read.csv("GB_ballon_data.csv",row.names = 1)
ggballoonplot(GB_ballon_aa_cna, 
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
    values=rescale(c(-0.5,0,0.7)),
    limits=c(-0.5,0.7))+
  labs(x="Gene status (AA mutation + CNV)",y="Cancer type ( p < 0.05)",fill="Difference of 
Genome doublings
(ref: WT)")


Purity_ballon_aa_cna <- read.csv("Purity_ballon_data.csv",row.names = 1)
ggballoonplot(Purity_ballon_aa_cna, 
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
    values=rescale(c(-0.05,0,0.02)),
    limits=c(-0.05,0.12))+
  labs(x="Gene status (AA mutation + CNV)",y="Cancer type ( p < 0.05)",fill="Difference of 
Purity
(ref: WT)")


Stroma_ballon_aa_cna <- read.csv("Stroma_ballon_data.csv",row.names = 1)
ggballoonplot(Stroma_ballon_aa_cna, 
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
    values=rescale(c(-0.2,0,0.05)),
    limits=c(-0.2,0.05))+
  labs(x="Gene status (AA mutation + CNV)",y="Cancer type ( p < 0.05)",fill="Difference of 
Stromal content
(ref: WT)")


Stroma_notLeukocyte_ballon_aa_cna <- read.csv("Stroma_notLeukocyte_ballon_data.csv",row.names = 1)
ggballoonplot(Stroma_notLeukocyte_ballon_aa_cna, 
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
    values=rescale(c(-0.07,0,0.03)),
    limits=c(-0.07,0.03))+
  labs(x="Gene status (AA mutation + CNV)",y="Cancer type ( p < 0.05)",fill="Difference of 
non-Stromal content
(ref: WT)")


NM_ballon_aa_cna <- read.csv("NM_ballon_data.csv",row.names = 1)
ggballoonplot(NM_ballon_aa_cna, 
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
    colors=c("blue","red"),
    values=rescale(c(0,1100)))+
  labs(x="Gene status (AA mutation + CNV)",y="Cancer type ( p < 0.05)",fill="Difference of 
total number of mutations
(ref: WT)")

