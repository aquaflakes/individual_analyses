fjComm::clear_()

# pallete=RColorBrewer::brewer.pal(n = 11, name = "Spectral") %>% rev()
df=read_tsv("LoVoFixed.fragment_length.closest_dist_to_center_within_500.ELF1_Control_TAAGCA40NGTC_KW_NACCCGGAAGTN_m1_c3_short_overlap_LoVo_ChIP.no_blacklist.txt",col_names = F)
p=ggplot(df)+stat_bin_2d(aes(X11,X4),binwidth = 2)+scale_fill_gradientn(colours =c("white","#ead100","red","red"),limits=c(0,40),name="Count",breaks=c(0,40))+scale_y_continuous(limits = c(0,200),expand = c(0,0),breaks = c(0,100,200))+scale_x_continuous(limits = c(-200,200),expand = c(0,0),breaks = c(-200,0,200))+
    xlab("Distance from motif (bp)")+ ylab("Fragment length (bp)")+gg_theme_Publication()+theme(legend.direction = "horizontal",legend.position = c(0.8,0.2))
print(p)
gg_save_pdf(p,6,5.5,filename = "ELF1_Vplot")

