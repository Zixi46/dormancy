my_comparisons4 <- list(c("WT_WT","WT_MUT"),c("WT_WT","MUT_MUT"),c("WT_WT","MUT_WT"))
ggboxplot(df.merged9, x ="PYHIN1_TP53_status", y = "PYHIN1_log2expression",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          width=0.3,
          order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = my_comparisons4,label = "p.signif")+
  labs(x="PYHIN1_TP53",y="mRNA Expression (log2)",fill="新图例名称")+
  scale_y_continuous(limits=c(0,20))        ##修改y轴标签的名称

ggboxplot(df.merged_gnaq2, x ="GNAQ_TP53_status", y = "GNAQ_log2expression",
          color = "GNAQ_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          width=0.3,
          order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = my_comparisons4,label = "p.signif")+
  labs(x="GNAQ_TP53",y="mRNA Expression (log2)",fill="新图例名称")+
  scale_y_continuous(limits=c(0,20))  

ggboxplot(df.merged_fgfr2_2, x ="FGFR2_TP53_status", y = "FGFR2_log2expression",
          color = "FGFR2_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          width=0.3,
          order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = my_comparisons4,label = "p.signif")+
  labs(x="FGFR2_TP53",y="mRNA Expression (log2)",fill="新图例名称")+
  scale_y_continuous(limits=c(0,20))  

ggboxplot(df.merged_erbb3_2, x ="ERBB3_TP53_status", y = "ERBB3_log2expression",
          color = "ERBB3_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          width=0.3,
          order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = my_comparisons4,label = "p.signif")+
  labs(x="ERBB3_TP53",y="mRNA Expression (log2)",fill="新图例名称")+
  scale_y_continuous(limits=c(0,20))  

ggboxplot(df.merged_plppr4_2, x ="PLPPR4_TP53_status", y = "PLPPR4_log2expression",
          color = "PLPPR4_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          width=0.3,
          order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = my_comparisons4,label = "p.signif")+
  labs(x="PLPPR4_TP53",y="mRNA Expression (log2)",fill="GENEstatus_TP53status")+
  guides(colour=guide_legend(title="GENEstatus_TP53status"))+
  scale_y_continuous(limits=c(0,20))  

table(df.merged9$PYHIN1_status)
table(df.merged9$TP53_status)
table(df.merged_gnaq2$GNAQ_status)
table(df.merged_fgfr2_2$FGFR2_status)
table(df.merged_erbb3_2$ERBB3_status)
table(df.merged_fgfr2_2$FGFR2_status)
table(df.merged_plppr4_2$PLPPR4_status)

table(df.merged9$PYHIN1_TP53_status)

