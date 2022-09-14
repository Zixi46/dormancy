
######Signature合并——最多样本——PYHIN1#####
df.merged14_new <- merge(df.merged2[,c("SAMPLE_ID","PYHIN1_status")], 
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

aggregate(df.merged14_new$SBS1,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS1",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS2,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS2",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS3,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS3",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS4,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS4",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS5,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS5",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS6,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS6",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS7a,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS7a",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS7b,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS7b",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS7c,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS7c",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS7d,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS7d",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS8,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS8",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS9,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS9",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS10a,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS10a",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS10b,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS10b",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS11,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS11",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS12,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS12",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS13,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS13",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS14,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS14",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS15,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS15",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS16,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS16",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS17a,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS17a",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS17b,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS17b",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS18,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS18",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS19,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS19",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS20,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS20",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS21,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS21",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS22,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS22",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS23,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS23",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS24,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS24",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS25,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS25",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS26,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS26",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS27,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS27",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS28,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS28",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS29,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS29",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS30,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS30",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS31,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS31",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS32,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS32",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS33,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS33",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS34,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS34",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS35,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS35",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS36,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS36",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS37,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS37",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS38,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS38",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS39,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS39",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS40,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS40",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS41,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS41",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS42,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS42",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS43,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS43",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS44,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS44",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS45,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS45",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS46,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS46",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS47,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS47",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS48,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS48",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS49,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS49",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS50,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS50",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS51,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS51",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS52,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS52",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS53,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS53",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS54,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS54",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS55,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS55",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS56,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS56",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS57,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS57",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS58,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS58",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS59,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS59",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)

aggregate(df.merged14_new$SBS60,by=list(type=df.merged14_new$PYHIN1_status),mean,na.rm = TRUE)#exclude NA
ggboxplot(df.merged14_new, x = "PYHIN1_status", y = "SBS60",
          color = "PYHIN1_status", palette =c("#00AFBB", "#E7B800"),
          add = "jitter",order = c("WT", "MUT")) + #添加jitter后，会自动隐藏异常值。
  stat_compare_means(comparisons = my_comparisons)


