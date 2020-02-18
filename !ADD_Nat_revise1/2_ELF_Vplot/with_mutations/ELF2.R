fjComm::clear_()

for (i in c("1stGtoA","1stGtoC","1stGtoT"))
{
  df=read_tsv(paste0("1nt/LoVoFixed.fragment_length.closest_dist_to_c_within_500.ELF2_modified_",i,"_TGCAAG20NAAC_AL_AACCCGGAAGTRN_m2_c2_short_p0.0001_S5.zero_fragment_cov.c_point.no_bl.txt"),col_names = F)
  p=ggplot(df)+stat_bin_2d(aes(X11,X4),binwidth = 2)+scale_fill_gradientn(colours =c("white","#ead100","red","red"),name="Count")+scale_y_continuous(limits = c(0,200),expand = c(0,0),breaks = c(0,100,200))+scale_x_continuous(limits = c(-200,200),expand = c(0,0),breaks = c(-200,0,200))+
    xlab("Distance from motif (bp)")+ ylab("Fragment length (bp)")+gg_theme_Publication()+theme(legend.direction = "horizontal",legend.position = c(0.8,0.2))
  print(p)
  # ,limits=c(0,15),breaks=c(0,15),oob=scales::squish
  gg_save_pdf(p+ggtitle(i),6,5.5,filename = "ELF2_Vplot_" %>% paste0(i))
}


  df=read_tsv("1nt/LoVoFixed.fragment_length.closest_dist_to_c_within_500.ELF2_TGCAAG20NAAC_AL_AACCCGGAAGTRN_m2_c2_short_p0.0001_S5.no_overlap_LoVo_ChIP.zero_fragment_cov.c_point.no_bl.txt",col_names = F)
  p=ggplot(df)+stat_bin_2d(aes(X11,X4),binwidth = 2)+scale_fill_gradientn(colours =c("white","#ead100","red","red"),name="Count")+scale_y_continuous(limits = c(0,200),expand = c(0,0),breaks = c(0,100,200))+scale_x_continuous(limits = c(-200,200),expand = c(0,0),breaks = c(-200,0,200))+
    xlab("Distance from motif (bp)")+ ylab("Fragment length (bp)")+gg_theme_Publication()+theme(legend.direction = "horizontal",legend.position = c(0.8,0.2))
  print(p)
  # ,limits=c(0,15),breaks=c(0,15),oob=scales::squish
  gg_save_pdf(p+ggtitle("no mutation"),6,5.5,filename = "ELF2_Vplot_no_mutation")


