fjComm::clear_()
# well filter
#(total is 1487 TFs

lig147=fjComm::nusel_get_147() %>% select(TF,pos,class,family,stdFamily,use,order)
lig200=fjComm::nusel_get_200() %>% select(TF,pos,class,family,stdFamily,use,order)

all_TFs= rbind(lig147 %>% dplyr::filter(use==1 & order==1), lig200 %>% dplyr::filter(use==1 & order==1))
all_TFs$TF %<>% str_replace("ooo.*$","")
all_TFs %<>% dplyr::distinct(TF,.keep_all = T) %>% mutate(stdFamily=TF_classify(.$TF)) %>% dplyr::filter(!is.na(stdFamily))

  # the "tried" amount
  all_TFs_t= rbind(lig147, lig200)
  all_TFs_t$TF %<>% str_replace("ooo.*$","")
  all_TFs_t %<>% dplyr::distinct(TF,.keep_all = T) %>% mutate(stdFamily=TF_classify(.$TF)) %>% dplyr::filter(!is.na(stdFamily))


family_stat=all_TFs$stdFamily %>% table() %>% as.data.frame() %>% set_colnames(c("family","exp_cnt"))
family_stat_t=all_TFs_t$stdFamily %>% table() %>% as.data.frame() %>% set_colnames(c("family","tried"))

total= rio::import("~/Nut_zhuData/DefaultFilesforProc/TF_classification/TF_family_total.xlsx")
total= merge(total, family_stat,by="family") %>% merge(family_stat_t, by="family") %>% arrange(desc(total)) %>% mutate(family= factor(family,levels = family))

gg_axis_x_labels<- function(x_labels=c("tick1","tick2","tick3"),...){scale_x_discrete(labels=x_labels,...)}

p=ggplot()+geom_bar(data = total, aes(family,total,fill="grey"),stat = "identity",width = 0.5) +
   geom_bar(data = total, aes(family,exp_cnt,fill="orange"),stat = "identity",width = 0.5,alpha=1)+
  gg_axis_y_noExp(breaks=gg_breaks(limits = c(0,70)))+
  coord_cartesian(ylim = c(0,77)) +gg_theme_Publication()+gg_axis_x_label_angle(45)+ ylab("Number of TF") + theme(axis.title.x = element_blank())+
  scale_fill_manual(values = c("grey","orange"),labels=c("Total","Successful"),name=NULL)+ theme(legend.position = c(0.8,0.8),legend.key.size = unit(0.3,"cm") )
print(p)

gg_save_pdf(p,11,6.5,filename = "TF_family_stat")

p1=ggplot()+geom_bar(data = total, aes(family,total,fill="grey"),stat = "identity",width = 0.5) +
  geom_bar(data = total, aes(family,tried,fill="yellow"),stat = "identity",width = 0.5,alpha=1)+
  gg_axis_y_noExp(breaks=c(0,30,60,90))+
  coord_cartesian(ylim = c(0,100)) +gg_theme_Publication()+gg_axis_x_label_angle(45)+ ylab("Number of TF") + theme(axis.title.x = element_blank())+
  scale_fill_manual(values = c("grey","yellow"),labels=c("Total","Tried"),name=NULL)+ theme(legend.position = c(0.8,0.8),legend.key.size = unit(0.3,"cm") )
print(p1)

gg_save_pdf(p1,11,8.5,filename = "TF_family_stat_tried")







# # total Number
# batch_4_5_guide_200= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.5CAP2r_lig200_curate_v1.txt") %>% select(TF,pos,class,family)
# batch_4_4_guide_200= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.4CAP2_lig200_curate_v1.txt") %>% select(TF,pos,class,family)
# batch_4_5_guide_147= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.5CAP2r_lig147_curate_v1.txt") %>% select(TF,pos,class,family)
# batch_4_4_guide_147= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.4CAP2_lig147_curate_v1.txt") %>% select(TF,pos,class,family)
#
# all_TFs= rbind(batch_4_5_guide_200,
#                batch_4_4_guide_200,
#                batch_4_5_guide_147,
#                batch_4_4_guide_147)
# all_TFs$TF %<>% str_replace("ooo.*$","")
# all_TFs %<>% dplyr::distinct(TF,.keep_all = T)
# print(nrow(all_TFs))

