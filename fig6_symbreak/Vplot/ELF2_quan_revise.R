fjComm::clear_()
p_load(IRanges)
# pallete=RColorBrewer::brewer.pal(n = 11, name = "Spectral") %>% rev()
df=read_tsv("LoVoFixed.fragment_length.closest_dist_to_center_within_500.ELF2_TGCAAG20NAAC_AL_AACCCGGAAGTRN_m2_c2_short_overlap_LoVo_ChIP.no_blacklist.txt",col_names = F)
df=data.frame(center=df$X11,lens=df$X4) %>% mutate(left=center-floor(lens/2),right=center+floor(lens/2))

df_cut=rbind(data.frame(ori="left",pos=df$left),data.frame(ori="right",pos=df$right))
cut_density=ggplot(df_cut)+stat_density(aes(pos,color=ori),bw=2,geom = "line",position = "identity")+ scale_x_continuous(limits = c(-200,200),expand = c(0,0),breaks = c(-200,0,200))+ scale_y_continuous(expand = c(0,0),breaks = c(0.002,0.003))+xlab("Distance from motif (bp)")+ylab("Density")+
  scale_color_manual(values = c("#F8766D", "#00BFC4"),name="End",labels=c("5'","3'"))+theme(legend.position = c(0.6,0.3))
gg_save_pdf(cut_density,7,3,filename = "ELF2_cut_density")

df1=df %>% dplyr::filter(lens<170&lens>140)
irng=IRanges(start = df1$left+201,end = df1$right+201)
cvg=irng %>% coverage()
cvg=cvg[1:401]

cvg_plot=ggplot()+geom_line(aes(-200:200,cvg))+ scale_x_continuous(limits = c(-200,200),expand = c(0,0),breaks = c(-200,0,200))+ scale_y_continuous(limits = c(10000,14000),breaks = c(10000,12000,14000))+xlab("Distance from motif (bp)")+ylab("Coverage")
gg_save_pdf(cvg_plot,6,5,filename = "ELF2_cut_coverage")
print(cvg_plot)
#
# df_curr=df %>% dplyr::filter(lens >=120 & lens <150)
# df_cut_left=data.frame(ori="left",pos=df_curr$left); df_cut_right=data.frame(ori="right",pos=df_curr$right)
# density_=density.default(df_cut_left$pos,bw = 0.5)
# ggplot()+geom_line(aes(density_$x,density_$y))+geom_line(aes(density$x,fjComm::math_sub_baseline(density_$y,baseline_hwm = 5)),color="yellow")
#
# # density_$y=density_$y %>% fjComm::math_sub_baseline(baseline_hwm = 5)
# peakVect=fjComm::math_sub_baseline(density_$y,baseline_hwm = 5) %>% set_names(density_$x)
#
# sub_bline=fjComm::math_sub_baseline(df_cut_left$pos)
# bline_=fjComm::math_get_baseline(df_cut_left$pos)
# ggplot()+stat_density(aes(df_cut_left$pos),bw=1,geom = "line")

# p=ggplot(df)+stat_bin_2d(aes(X11,X4),binwidth = 2)+scale_fill_gradientn(colours =c("white","#ead100","red","red"),limits=c(0,25),name="Count",breaks=c(0,25))+scale_y_continuous(limits = c(0,200),expand = c(0,0),breaks = c(0,100,200))+scale_x_continuous(limits = c(-200,200),expand = c(0,0),breaks = c(-200,0,200))+
#     xlab("Distance from motif (bp)")+ ylab("Fragment length (bp)")+gg_theme_Publication()+theme(legend.direction = "horizontal",legend.position = c(0.8,0.2))
# print(p)
# gg_save_pdf(p,6,5.5,filename = "ELF2_Vplot_quan")

