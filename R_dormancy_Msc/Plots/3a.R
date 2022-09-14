#3a

my_comparisons3 <- list(c("WT","MUT"),c("WT","AMP"),c("WT","DEL"),c("MUT","AMP"),c("MUT","DEL"),c("AMP","DEL"))
ggboxplot(df.merged9, x = "PYHIN1_aamutation_cna", y = "QuiescenceScore",
          width=0.3,
          color = "PYHIN1_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3)+
  labs(x="PYHIN1 status",y="Quiescence Score",fill="GENEstatus_TP53status")+
  guides(colour=guide_legend(title="GENE status (AA mutation and CNV)"))

ggboxplot(df.merged_gnaq2, x = "GNAQ_aamutation_cna", y = "QuiescenceScore",
          width=0.3,
          color = "GNAQ_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")+
  labs(x="GNAQ status",y="Quiescence Score",fill="GENEstatus_TP53status")+
  guides(colour=guide_legend(title="GENE status (AA mutation and CNV)"))

ggboxplot(df.merged_fgfr2_2, x = "FGFR2_aamutation_cna", y = "QuiescenceScore",
          width=0.3,
          color = "FGFR2_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")+
  labs(x="FGFR2 status",y="Quiescence Score",fill="GENEstatus_TP53status")+
  guides(colour=guide_legend(title="GENE status (AA mutation and CNV)"))

ggboxplot(df.merged_erbb3_2, x = "ERBB3_aamutation_cna", y = "QuiescenceScore",
          width=0.3,
          color = "ERBB3_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")+
  labs(x="ERBB3 status",y="Quiescence Score",fill="GENEstatus_TP53status")+
  guides(colour=guide_legend(title="GENE status (AA mutation and CNV)"))

ggboxplot(df.merged_plppr4_2, x = "PLPPR4_aamutation_cna", y = "QuiescenceScore",
          width=0.3,
          color = "PLPPR4_aamutation_cna", palette =c("#00AFBB", "#E7B800","#D2691E","#008000"),
          order = c("WT", "MUT","AMP","DEL")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons3,label = "p.signif")+
  labs(x="PLPPR4 status",y="Quiescence Score",fill="GENEstatus_TP53status")+
  guides(colour=guide_legend(title="GENE status (AA mutation and CNV)"))
