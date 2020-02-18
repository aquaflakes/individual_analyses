fjComm::clear_()
# well filter
#(total is 1487 TFs
batch_4_5_guide_200= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.5CAP2r_lig200_curate_v1.txt") %>% dplyr::filter(use==1 & order==1) %>% select(TF,pos,class,family)
batch_4_4_guide_200= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.4CAP2_lig200_curate_v1.txt") %>% dplyr::filter(use==1 & order==1)%>% select(TF,pos,class,family)
batch_4_5_guide_147= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.5CAP2r_lig147_curate_v1.txt") %>% dplyr::filter(use==1 & order==1)%>% select(TF,pos,class,family)
batch_4_4_guide_147= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.4CAP2_lig147_curate_v1.txt") %>% dplyr::filter(use==1 & order==1)%>% select(TF,pos,class,family)

all_TFs= rbind(batch_4_5_guide_200,
               batch_4_4_guide_200,
               batch_4_5_guide_147,
               batch_4_4_guide_147)
all_TFs$TF %<>% str_replace("ooo.*$","")
all_TFs %<>% dplyr::distinct(TF,.keep_all = T)

  # # the "tried" amount
  #   batch_4_5_guide_200_t= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.5CAP2r_lig200_curate_v1.txt") %>% select(TF,pos,class,family)
  #   batch_4_4_guide_200_t= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.4CAP2_lig200_curate_v1.txt") %>% select(TF,pos,class,family)
  #   batch_4_5_guide_147_t= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.5CAP2r_lig147_curate_v1.txt") %>% select(TF,pos,class,family)
  #   batch_4_4_guide_147_t= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.4CAP2_lig147_curate_v1.txt") %>% select(TF,pos,class,family)
  #
  #   all_TFs_t= rbind(batch_4_5_guide_200_t,
  #                  batch_4_4_guide_200_t,
  #                  batch_4_5_guide_147_t,
  #                  batch_4_4_guide_147_t)
  #   all_TFs_t$TF %<>% str_replace("ooo.*$","")
  #   all_TFs_t %<>% dplyr::distinct(TF,.keep_all = T)


all_TFs$newFam= TF_classify(all_TFs$TF)
family_stat=all_TFs$newFam %>% table() %>% as.data.frame() %>% set_colnames(c("family","exp_cnt"))

total= rio::import("~/Nut_zhuData/DefaultFilesforProc/TF_classification/TF_family_total.xlsx")
total= merge(total, family_stat,by="family") %>% arrange(desc(total)) %>% mutate(family= factor(family,levels = family))

gg_axis_x_labels<- function(x_labels=c("tick1","tick2","tick3"),...){scale_x_discrete(labels=x_labels,...)}

p=ggplot()+geom_bar(data = total, aes(family,total,fill="grey"),stat = "identity",width = 0.5) + geom_bar(data = total, aes(family,exp_cnt,fill="orange"),stat = "identity",width = 0.5,alpha=1)+
  gg_axis_y_noExp(breaks=gg_breaks(limits = c(0,70)))+
  coord_cartesian(ylim = c(0,77)) +gg_theme_Publication()+gg_axis_x_label_angle(45)+ ylab("Number of TF") + theme(axis.title.x = element_blank())+
  scale_fill_manual(values = c("grey","orange"),labels=c("Total","Successful"),name=NULL)+ theme(legend.position = c(0.8,0.8),legend.key.size = unit(0.3,"cm") )
print(p)

gg_save_pdf(p,11,6.5,filename = "TF_family_stat")




# total Number
batch_4_5_guide_200= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.5CAP2r_lig200_curate_v1.txt") %>% select(TF,pos,class,family)
batch_4_4_guide_200= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.4CAP2_lig200_curate_v1.txt") %>% select(TF,pos,class,family)
batch_4_5_guide_147= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.5CAP2r_lig147_curate_v1.txt") %>% select(TF,pos,class,family)
batch_4_4_guide_147= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.4CAP2_lig147_curate_v1.txt") %>% select(TF,pos,class,family)

all_TFs= rbind(batch_4_5_guide_200,
               batch_4_4_guide_200,
               batch_4_5_guide_147,
               batch_4_4_guide_147)
all_TFs$TF %<>% str_replace("ooo.*$","")
all_TFs %<>% dplyr::distinct(TF,.keep_all = T)
print(nrow(all_TFs))

