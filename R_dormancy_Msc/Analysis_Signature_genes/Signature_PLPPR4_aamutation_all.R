######Signature合并——最多样本——PLPPR4#####
df.merged_plppr4_sig <- merge(mutations_plppr4[,c("SAMPLE_ID","PLPPR4_status")], 
                             sigs_all_new[,c("sig_sample_id","SBS1","SBS2","SBS3",
                                             "SBS4","SBS5","SBS6","SBS7a",
                                             "SBS7b","SBS7c","SBS7d","SBS8",
                                             "SBS9","SBS10a","SBS10b","SBS11",
                                             "SBS12","SBS13","SBS14","SBS15",
                                             "SBS16","SBS17a","SBS17b","SBS18",
                                             "SBS19","SBS20","SBS21","SBS22",
                                             "SBS23","SBS24","SBS25","SBS26",
                                             "SBS27","SBS28","SBS29","SBS30",
                                             "SBS31","SBS32","SBS33","SBS34",
                                             "SBS35","SBS36","SBS37","SBS38",
                                             "SBS39","SBS40","SBS41","SBS42",
                                             "SBS43","SBS44","SBS45","SBS46",
                                             "SBS47","SBS48","SBS49","SBS50",
                                             "SBS51","SBS52","SBS53","SBS54",
                                             "SBS55","SBS56","SBS57","SBS58",
                                             "SBS59","SBS60")],
                             by.x="SAMPLE_ID", by.y="sig_sample_id", #x和y表合并的依据是相同的sample id
                             all.x=FALSE, all.y=FALSE) 
View(df.merged_plppr4_sig)

aggregate(df.merged_plppr4_sig$SBS1,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS1",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS2,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS2",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS3,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS3",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS4,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS4",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS5,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS5",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS6,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS6",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS7a,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS7a",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS7b,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS7b",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS7c,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS7c",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS7d,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS7d",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS8,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS8",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS9,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS9",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS10a,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS10a",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS10b,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS10b",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS11,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS11",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS12,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS12",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS13,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS13",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS14,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS14",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS15,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS15",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS16,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS16",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS17a,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS17a",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS17b,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS17b",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS18,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS18",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS19,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS19",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS20,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS20",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS21,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS21",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS22,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS22",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS23,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS23",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS24,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS24",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS25,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS25",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS26,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS26",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS27,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS27",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS28,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS28",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS29,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS29",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS30,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS30",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS31,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS31",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS32,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS32",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS33,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS33",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS34,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS34",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS35,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS35",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS36,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS36",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS37,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS37",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS38,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS38",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS39,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS39",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS40,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS40",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS41,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS41",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS42,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS42",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS43,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS43",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS44,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS44",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS45,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS45",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS46,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS46",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS47,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS47",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS48,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS48",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS49,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS49",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS50,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS50",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS51,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS51",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS52,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS52",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS53,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS53",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS54,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS54",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS55,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS55",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS56,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS56",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS57,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS57",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS58,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS58",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS59,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS59",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged_plppr4_sig$SBS60,by=list(type=df.merged_plppr4_sig$PLPPR4_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged_plppr4_sig, x = "PLPPR4_status", y = "SBS60",
          color = "PLPPR4_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)