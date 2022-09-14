#3d

ggboxplot(df.merged10, x = "PYHIN1_aamutation_cna_TP53_aamutation", y = "QuiescenceScore",
          color = "PYHIN1_aamutation_cna_TP53_aamutation", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                      "#D5D2C1","#00AFBB"),
          width=0.3,order = c("WT_WT","WT_MUT","MUT_WT","MUT_MUT","AMP_WT", "AMP_MUT","DEL_WT","DEL_MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons4)+
  labs(x="PYHIN1_TP53",y="Quiescence Score",fill="GENEstatus_TP53status")+
  guides(colour=guide_legend(title=NULL))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

ggboxplot(df.merged_gnaq2, x = "GNAQ_aamutation_cna_TP53_aamutation", y = "QuiescenceScore",
          color = "GNAQ_aamutation_cna_TP53_aamutation", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                      "#D5D2C1","#00AFBB"),
          width=0.3,order = c("WT_WT","WT_MUT","MUT_WT","MUT_MUT","AMP_WT", "AMP_MUT","DEL_WT","DEL_MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons4)+
  labs(x="GNAQ_TP53",y="Quiescence Score",fill="GENEstatus_TP53status")+
  guides(colour=guide_legend(title=NULL))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

ggboxplot(df.merged_fgfr2_2, x = "FGFR2_aamutation_cna_TP53_aamutation", y = "QuiescenceScore",
          color = "FGFR2_aamutation_cna_TP53_aamutation", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                    "#D5D2C1","#00AFBB"),
          width=0.3,order = c("WT_WT","WT_MUT","MUT_WT","MUT_MUT","AMP_WT", "AMP_MUT","DEL_WT","DEL_MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons4)+
  labs(x="FGFR2_TP53",y="Quiescence Score",fill="GENEstatus_TP53status")+
  guides(colour=guide_legend(title=NULL))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

ggboxplot(df.merged_erbb3_2, x = "ERBB3_aamutation_cna_TP53_aamutation", y = "QuiescenceScore",
          color = "ERBB3_aamutation_cna_TP53_aamutation", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                     "#D5D2C1","#00AFBB"),
          width=0.3,order = c("WT_WT","WT_MUT","MUT_WT","MUT_MUT","AMP_WT", "AMP_MUT","DEL_WT","DEL_MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons4)+
  labs(x="ERBB3_TP53",y="Quiescence Score",fill="GENEstatus_TP53status")+
  guides(colour=guide_legend(title=NULL))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

ggboxplot(df.merged_plppr4_2, x = "PLPPR4_aamutation_cna_TP53_aamutation", y = "QuiescenceScore",
          color = "PLPPR4_aamutation_cna_TP53_aamutation", palette =c("#00AFBB", "#E7B800","#D2691E","#008000","#6F6F6F","#8E8BFE",
                                                                     "#D5D2C1","#00AFBB"),
          width=0.3,order = c("WT_WT","WT_MUT","MUT_WT","MUT_MUT","AMP_WT", "AMP_MUT","DEL_WT","DEL_MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons4)+
  labs(x="PLPPR4_TP53",y="Quiescence Score",fill="GENEstatus_TP53status")+
  guides(colour=guide_legend(title=NULL))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
