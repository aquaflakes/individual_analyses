fjComm::clear_()

lig147=nusel_get_147() %>% dplyr::filter(use==1&order==1)
temp=lig147 %>% dplyr::filter(TF=="RFX5")

df1=read_csv(temp$c4_EMI_file) %>% mutate(TF="NCAP") %>% dplyr::filter(pos2-pos1==3)
df2=read_csv(temp$c4_EMI_TFctrl_file) %>% mutate(TF="HT") %>% dplyr::filter(pos2-pos1==3)

MIdf=rbind(df1,df2)
MI_plot_RFX5=ggplot(MIdf)+geom_line(aes(pos1,topMIsum,color=TF))+xlab("(bp)")+ylab("E-MI (bit)")+scale_color_manual(values = c("steelblue","orange"),name="E-MI")+scale_x_continuous(expand = c(0,0),breaks = c(1,96))+scale_y_continuous(breaks = c(0,0.06)) 
gg_save_pdf(MI_plot_RFX5,7,3)



TFname= "RFX5"
pfmfile="RFX5_motif_FJ4.4.txt"

file=temp$c4_file
motif_map=fjComm::motif_plot_moods_singlePWM_hit(readsfile = file,pfmfile = pfmfile,dup_rm = T,combine_2strand_hits = T)
 readsNum=file %>% read_csv(col_names=F) %>% rmdup() %>% nrow()
df1=motif_map$data %>% mutate(TF="NCAP",count=count/readsNum)

file=temp$TFctrl_file
motif_map=fjComm::motif_plot_moods_singlePWM_hit(readsfile = file,pfmfile = pfmfile,dup_rm = T,combine_2strand_hits = T)
  readsNum=file %>% read_csv(col_names=F) %>% rmdup() %>% nrow()
df2=motif_map$data %>% mutate(TF="HT",count=count/readsNum)
Motifdf=rbind(df1,df2)
Motif_plot_RFX5=ggplot(Motifdf)+geom_line(aes(pos+1,count,color=TF))+xlab("(bp)")+ylab("Motif Freq.")+scale_color_manual(values = c("steelblue","orange"),name="")+scale_x_continuous(expand = c(0,0),breaks = c(1,92))+scale_y_continuous(breaks = c(0,0.01)) 
gg_save_pdf(Motif_plot_RFX5,7,3)

