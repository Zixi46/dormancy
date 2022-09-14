my_comparisons4 <- list(c("WT_WT","WT_MUT"),c("WT_WT","MUT_MUT"),c("WT_WT","MUT_WT"))
ggboxplot(df.merged9, x ="PYHIN1_TP53_status", y = "QuiescenceScore",
          color = "PYHIN1_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          width=0.4,
          order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = my_comparisons4,label = "p.signif")+
  labs(x="PYHIN1_TP53",y="Quiescence Score",fill="新图例名称")

ggboxplot(df.merged_gnaq2, x ="GNAQ_TP53_status", y = "QuiescenceScore",
          color = "GNAQ_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          width=0.4,
          order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = my_comparisons4,label = "p.signif")+
  labs(x="GNAQ_TP53",y="Quiescence Score",fill="新图例名称")

ggboxplot(df.merged_fgfr2_2, x ="FGFR2_TP53_status", y = "QuiescenceScore",
          color = "FGFR2_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          width=0.4,
          order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = my_comparisons4,label = "p.signif")+
  labs(x="FGFR2_TP53",y="Quiescence Score",fill="新图例名称")

ggboxplot(df.merged_erbb3_2, x ="ERBB3_TP53_status", y = "QuiescenceScore",
          color = "ERBB3_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          width=0.4,
          order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = my_comparisons4,label = "p.signif")+
  labs(x="ERBB3_TP53",y="Quiescence Score",fill="新图例名称")

ggboxplot(df.merged_plppr4_2, x ="PLPPR4_TP53_status", y = "QuiescenceScore",
          color = "PLPPR4_TP53_status", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          width=0.4,
          order=c("WT_WT","WT_MUT","MUT_MUT","MUT_WT"))+
  stat_compare_means(comparisons = my_comparisons4,label = "p.signif")+
  labs(x="PLPPR4_TP53",y="Quiescence Score",fill="GENEstatus_TP53status")+
  guides(colour=guide_legend(title="GENEstatus_TP53status"))
